#!/usr/bin/env python3
"""Whole-contig minimap2 alignments, two modes.

MTPT mode (default): automates the manual example already in the top-level
README across the whole dataset - for every QC-eligible species with both
a mito and a plastid resolved assembly, align mito (ref) vs plastid (query)
with `minimap2 -x asm5` and save the PAF.

Pairwise mode: for a student-supplied --species-list (e.g. one genus) and a
single --organelle, align every pair within that list. Deliberately guarded
(small species-list only, or --allow-all-pairs) so an accidental
all-vs-all across ~1250 species can't happen - that would no longer be
"lightweight."

Output:
  MTPT mode:     analysis/synteny/results/mtpt/<Species>.mito_vs_pltd.paf
  Pairwise mode: analysis/synteny/results/paf/<A>__<B>.<organelle>.paf
"""
from __future__ import annotations

import os
import argparse
import itertools
import shutil
import subprocess
import sys
from pathlib import Path

# CODE_DIR is this checkout's analysis/ (shared modules, bundled reference
# data); data/ and every module's results/ and work/ live under the data
# root - PLANT_ORGANELLE_DATA_ROOT, or this checkout if unset.
CODE_DIR = Path(__file__).resolve().parents[2]
ROOT_DIR = Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT", str(CODE_DIR.parent)))
ANALYSIS_DIR = ROOT_DIR / "analysis"
sys.path.insert(0, str(CODE_DIR / "common"))

import species_discovery as sd  # noqa: E402

import pandas as pd  # noqa: E402

MINIMAP2 = shutil.which("minimap2") or "/software/team301/minimap2/minimap2"
MAX_PAIRWISE_SPECIES = 30


def qc_eligible(qc_status: list[str]) -> set[tuple[str, str]]:
    qc_path = ANALYSIS_DIR / "qc_basic_stats" / "results" / "qc_summary.tsv"
    qc = pd.read_csv(qc_path, sep="\t")
    sub = qc[qc["status"].isin(qc_status)]
    return set(zip(sub["species"], sub["organelle"]))


def run_minimap2(ref_fasta: str, query_fasta: str, out_paf: Path, force: bool) -> bool:
    if not force and out_paf.exists() and out_paf.stat().st_mtime >= max(
            Path(ref_fasta).stat().st_mtime, Path(query_fasta).stat().st_mtime):
        return False
    out_paf.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_paf.with_suffix(out_paf.suffix + ".tmp")
    with open(tmp, "w") as fh:
        proc = subprocess.run([MINIMAP2, "-x", "asm5", ref_fasta, query_fasta],
                               stdout=fh, stderr=subprocess.PIPE, text=True, timeout=300)
    if proc.returncode != 0:
        tmp.unlink(missing_ok=True)
        print(f"[warn] minimap2 failed for {ref_fasta} vs {query_fasta}: {proc.stderr.strip()}", file=sys.stderr)
        return False
    tmp.rename(out_paf)
    return True


def mtpt_mode(species_filter, qc_status: list[str], force: bool):
    eligible = qc_eligible(qc_status)
    mito_root = sd.repo_data_root(ROOT_DIR, "mito")
    pltd_root = sd.repo_data_root(ROOT_DIR, "pltd")
    species_list = species_filter or {p.name for p in mito_root.iterdir() if p.is_dir()}

    out_dir = ANALYSIS_DIR / "synteny" / "results" / "mtpt"
    n_computed = n_skipped = n_missing = 0
    for species in sorted(species_list):
        if (species, "mito") not in eligible or (species, "pltd") not in eligible:
            n_missing += 1
            continue
        mito_r = sd.resolve_species(mito_root / species, "mito")
        pltd_r = sd.resolve_species(pltd_root / species, "pltd")
        if mito_r.status != "ok" or pltd_r.status != "ok":
            n_missing += 1
            continue
        out_paf = out_dir / f"{species}.mito_vs_pltd.paf"
        if run_minimap2(mito_r.ctg_fasta, pltd_r.ctg_fasta, out_paf, force):
            n_computed += 1
        else:
            n_skipped += 1
    print(f"[info] pairwise_paf (mtpt): computed={n_computed} skipped={n_skipped} "
          f"not_eligible_or_missing={n_missing} -> {out_dir}", file=sys.stderr)


def pairwise_mode(species_filter, organelle: str, qc_status: list[str], force: bool, allow_all_pairs: bool):
    if organelle == "both":
        print("[err] pairwise mode requires a single --organelle (mito or pltd)", file=sys.stderr)
        sys.exit(1)
    if not species_filter:
        print("[err] pairwise mode requires --species-list", file=sys.stderr)
        sys.exit(1)
    if len(species_filter) > MAX_PAIRWISE_SPECIES and not allow_all_pairs:
        print(f"[err] --species-list has {len(species_filter)} species (max {MAX_PAIRWISE_SPECIES} "
              f"for pairwise mode without --allow-all-pairs) - an all-vs-all run at dataset scale "
              f"is not lightweight. Pass --allow-all-pairs to override.", file=sys.stderr)
        sys.exit(1)

    eligible = qc_eligible(qc_status)
    data_root = sd.repo_data_root(ROOT_DIR, organelle)
    resolved = {}
    for species in sorted(species_filter):
        if (species, organelle) not in eligible:
            continue
        r = sd.resolve_species(data_root / species, organelle)
        if r.status == "ok":
            resolved[species] = r.ctg_fasta

    out_dir = ANALYSIS_DIR / "synteny" / "results" / "paf"
    n_computed = n_skipped = 0
    for a, b in itertools.combinations(sorted(resolved), 2):
        out_paf = out_dir / f"{a}__{b}.{organelle}.paf"
        if run_minimap2(resolved[a], resolved[b], out_paf, force):
            n_computed += 1
        else:
            n_skipped += 1
    print(f"[info] pairwise_paf ({organelle}): {len(resolved)} species, "
          f"{len(resolved) * (len(resolved) - 1) // 2} pairs, computed={n_computed} skipped={n_skipped} -> {out_dir}",
          file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mode", default="mtpt", choices=["mtpt", "pairwise"])
    ap.add_argument("--organelle", default="both", choices=["mito", "pltd", "both"])
    ap.add_argument("--species-list", default=None)
    ap.add_argument("--qc-status", default="pass")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--allow-all-pairs", action="store_true")
    args = ap.parse_args()

    species_filter = sd.load_species_list(Path(args.species_list)) if args.species_list else None
    qc_status = args.qc_status.split(",")

    if args.mode == "mtpt":
        mtpt_mode(species_filter, qc_status, args.force)
    else:
        pairwise_mode(species_filter, args.organelle, qc_status, args.force, args.allow_all_pairs)


if __name__ == "__main__":
    main()
