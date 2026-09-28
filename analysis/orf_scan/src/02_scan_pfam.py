#!/usr/bin/env python3
"""hmmscan each species' non-core ORF protein set against the full Pfam-A
library (already pressed/binary on this cluster, so this is fast), one
process invocation per species (not per ORF) to avoid reloading the
~1.7GB database thousands of times.

Adds pfam hit columns to each per-species orfs.tsv written by
01_find_orfs.py (in place). Only the best (lowest full-sequence E-value)
domain hit per ORF is kept in this table - see the raw .domtblout in
work/ if you need every domain hit.

This is the heaviest step in the whole analysis suite: hmmscan against
Pfam-A's ~20,000 profiles is measured at ~15-90s per species (mito runs
with more ORFs cost more than plastid runs with very few) even against a
pressed database, so a full ~1250-species run is a multi-hour job, not a
few minutes like everything else here. Use --jobs to run several hmmscan
processes concurrently (each still using --cpu threads internally) - e.g.
inside `bsub -n 16`, `--jobs 4 --cpu 4` uses the full allocation.

Output (augmented in place): analysis/orf_scan/results/per_species/<Species>.<organelle>.orfs.tsv
Adds columns: pfam_acc pfam_name evalue score description
"""
from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import os
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

# Defaults to the old same-tree layout (code and data siblings under one
# repo root); set PLANT_ORGANELLE_DATA_ROOT to point at a separate data
# checkout instead (e.g. when this code is installed from its own repo).
ROOT_DIR = Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT", str(Path(__file__).resolve().parents[3])))
ANALYSIS_DIR = ROOT_DIR / "analysis"
HMMSCAN = shutil.which("hmmscan") or "/software/team301/hmmer-3.4/src/hmmscan"
PFAM_DB = "/software/team301/pfam_hmm/Pfam-A.hmm"

FILENAME_RE = re.compile(r"^(.+)\.(mito|pltd)\.noncore_orfs\.faa$")
PFAM_COLS = ["pfam_acc", "pfam_name", "evalue", "score", "description"]


def parse_domtblout(path: Path) -> pd.DataFrame:
    rows = []
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split(None, 22)
        if len(fields) < 22:
            continue
        rows.append({
            "pfam_name": fields[0], "pfam_acc": fields[1], "orf_query": fields[3],
            "evalue": float(fields[6]), "score": float(fields[7]),
            "description": fields[22] if len(fields) > 22 else "",
        })
    if not rows:
        return pd.DataFrame(columns=["orf_query"] + PFAM_COLS)
    df = pd.DataFrame(rows)
    best = df.loc[df.groupby("orf_query")["evalue"].idxmin()].copy()
    # avoid float round-trip artifacts like 1.4999999999999998e-29 in the output TSV
    best["evalue"] = best["evalue"].apply(lambda x: f"{x:.3g}")
    best["score"] = best["score"].apply(lambda x: f"{x:.1f}")
    return best[["orf_query"] + PFAM_COLS]


def scan_one(species: str, organelle: str, faa_path: Path, work_dir: Path, per_species_dir: Path,
             cpu: str, force: bool) -> str:
    tsv_path = per_species_dir / f"{species}.{organelle}.orfs.tsv"
    domtblout_path = work_dir / f"{species}.{organelle}.domtblout"
    if not force and domtblout_path.exists() and domtblout_path.stat().st_mtime >= faa_path.stat().st_mtime:
        return "skipped"

    cmd = [HMMSCAN, "--cut_ga", "--domtblout", str(domtblout_path), "--cpu", cpu, PFAM_DB, str(faa_path)]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=1200)
    if proc.returncode != 0:
        print(f"[warn] hmmscan failed for {species}.{organelle}: {proc.stderr.strip()[:300]}", file=sys.stderr)
        return "failed"

    hits = parse_domtblout(domtblout_path)
    orfs = pd.read_csv(tsv_path, sep="\t")
    # drop any pfam/category columns from a previous run of this script (or of
    # 03_orf_scan_summary.py) before merging - otherwise pandas silently
    # appends _x/_y suffixed duplicate columns instead of overwriting them.
    orfs = orfs.drop(columns=[c for c in PFAM_COLS + ["category"] if c in orfs.columns])
    orfs["orf_query"] = f"{species}|{organelle}|" + orfs["orf_id"].astype(str)
    merged = orfs.merge(hits, on="orf_query", how="left").drop(columns="orf_query")
    merged.to_csv(tsv_path, sep="\t", index=False)
    return "scanned"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--species-list", default=None)
    ap.add_argument("--cpu", default="2", help="threads per hmmscan process (default 2)")
    ap.add_argument("--jobs", type=int, default=1, help="number of hmmscan processes to run concurrently (default 1)")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    species_filter = None
    if args.species_list:
        species_filter = {l.strip() for l in open(args.species_list) if l.strip() and not l.startswith("#")}

    work_dir = ANALYSIS_DIR / "orf_scan" / "work"
    per_species_dir = ANALYSIS_DIR / "orf_scan" / "results" / "per_species"

    tasks = []
    n_no_orfs = 0
    for faa_path in sorted(work_dir.glob("*.noncore_orfs.faa")):
        m = FILENAME_RE.match(faa_path.name)
        if not m:
            continue
        species, organelle = m.groups()
        if species_filter is not None and species not in species_filter:
            continue
        if faa_path.stat().st_size == 0:
            n_no_orfs += 1
            continue
        tasks.append((species, organelle, faa_path))

    counts = {"scanned": 0, "skipped": 0, "failed": 0}
    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        futures = {pool.submit(scan_one, sp, org, faa, work_dir, per_species_dir, args.cpu, args.force): sp
                   for sp, org, faa in tasks}
        for fut in as_completed(futures):
            counts[fut.result()] += 1

    print(f"[info] scan_pfam: scanned={counts['scanned']} skipped={counts['skipped']} "
          f"failed={counts['failed']} no_noncore_orfs={n_no_orfs}", file=sys.stderr)


if __name__ == "__main__":
    main()
