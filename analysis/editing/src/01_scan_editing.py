#!/usr/bin/env python3
"""Run `orfedit` over every (species, gene) anchor that has a reference
profile, correcting gene boundaries for C->U RNA editing.

Builds one manifest TSV from `annotation/results/gene_calls.tsv` (already
QC-gated - see `annotation/src/01_gene_table.py`) joined against the
profiles `00_build_reference_profiles.py` built, then makes a single call
to `orfedit scan-batch`, which rayon-parallelizes across every row itself -
deliberately not the `--jobs`/`--cpu` multi-process pattern `orf_scan` uses
for `hmmscan`. That pattern exists because `hmmscan` isn't internally
parallel; `orfedit` is (via rayon), so one process using `--threads` is the
right shape here, and it sidesteps the LSF-array-of-many-small-jobs
scaling pain `orf_scan/README.md` already documents for its own heavy step.

Output: analysis/editing/results/edited_gene_calls.tsv
        analysis/editing/results/predicted_edits.tsv
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import os
import sys
from pathlib import Path

import pandas as pd

# Defaults to the old same-tree layout (code and data siblings under one
# repo root); set PLANT_ORGANELLE_DATA_ROOT to point at a separate data
# checkout instead (e.g. when this code is installed from its own repo).
ROOT_DIR = Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT", str(Path(__file__).resolve().parents[3])))
ANALYSIS_DIR = ROOT_DIR / "analysis"
sys.path.insert(0, str(ANALYSIS_DIR / "common"))

import species_discovery as sd  # noqa: E402

ORFEDIT = shutil.which("orfedit") or str(Path.home() / ".cargo" / "bin" / "orfedit")


def best_hits_per_species(gene_calls: pd.DataFrame, organelle: str, genes: list[str]) -> pd.DataFrame:
    sub = gene_calls[(gene_calls["organelle"] == organelle) & (gene_calls["gene"].isin(genes))]
    if sub.empty:
        return sub
    idx = sub.groupby(["species", "gene"])["score"].idxmax()
    return sub.loc[idx]


def build_manifest(gene_calls: pd.DataFrame, organelles: list[str], species_filter: set | None,
                    profiles_dir: Path) -> list[dict]:
    rows = []
    for organelle in organelles:
        organelle_profiles_dir = profiles_dir / organelle
        genes = sorted(p.stem for p in organelle_profiles_dir.glob("*.pssm"))
        if not genes:
            print(f"[warn] {organelle}: no profiles in {organelle_profiles_dir} - run "
                  f"00_build_reference_profiles.py first", file=sys.stderr)
            continue
        hits = best_hits_per_species(gene_calls, organelle, genes)
        if hits.empty:
            continue
        data_root = sd.repo_data_root(ROOT_DIR, organelle)
        # Resolve each species' ctg.fasta path once, reused across all its genes.
        fasta_cache: dict[str, str | None] = {}
        for row in hits.itertuples():
            if species_filter is not None and row.species not in species_filter:
                continue
            if row.species not in fasta_cache:
                r = sd.resolve_species(data_root / row.species, organelle)
                fasta_cache[row.species] = r.ctg_fasta if r.status == "ok" else None
            fasta_path = fasta_cache[row.species]
            if fasta_path is None:
                continue
            rows.append({
                "species": row.species, "organelle": organelle, "gene": row.gene,
                "fasta": fasta_path, "profile": str(organelle_profiles_dir / f"{row.gene}.pssm"),
                "contig": row.contig_id, "start": int(row.start), "end": int(row.end),
                "strand": row.strand,
            })
    return rows


def write_manifest(rows: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    cols = ["species", "organelle", "gene", "fasta", "profile", "contig", "start", "end", "strand"]
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organelle", default="both", choices=["mito", "pltd", "both"])
    ap.add_argument("--species-list", default=None)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--flank", type=int, default=300)
    ap.add_argument("--edit-penalty", type=float, default=2.5)
    ap.add_argument("--gene-calls", default=None, help="path to a gene_calls.tsv-shaped table to anchor "
                     "on (default: annotation/results/gene_calls.tsv, oatk's own calls) - e.g. point this "
                     "at denovo_annotation/results/gene_calls.tsv to use the oatk-independent "
                     "nhmmscan-based calls instead")
    ap.add_argument("--force", action="store_true", help="accepted for CLI consistency; orfedit "
                     "always does a full recompute over the manifest (no incremental skip yet)")
    args = ap.parse_args()

    if not Path(ORFEDIT).exists():
        print(f"[err] orfedit binary not found at {ORFEDIT} - see orfedit's own README "
              f"(`cargo install --path .`)", file=sys.stderr)
        sys.exit(1)

    gene_calls_path = Path(args.gene_calls) if args.gene_calls else ANALYSIS_DIR / "annotation" / "results" / "gene_calls.tsv"
    if not gene_calls_path.exists():
        print(f"[err] missing {gene_calls_path} - run annotation/src/01_gene_table.py first "
              f"(or pass --gene-calls)", file=sys.stderr)
        sys.exit(1)
    gene_calls = pd.read_csv(gene_calls_path, sep="\t")

    organelles = ["mito", "pltd"] if args.organelle == "both" else [args.organelle]
    species_filter = sd.load_species_list(Path(args.species_list)) if args.species_list else None
    profiles_dir = ANALYSIS_DIR / "editing" / "results" / "profiles"

    rows = build_manifest(gene_calls, organelles, species_filter, profiles_dir)
    if not rows:
        print("[err] empty manifest - nothing to scan", file=sys.stderr)
        sys.exit(1)

    work_dir = ANALYSIS_DIR / "editing" / "work"
    manifest_path = work_dir / "manifest.tsv"
    write_manifest(rows, manifest_path)
    print(f"[info] scan_editing: manifest has {len(rows)} rows -> {manifest_path}", file=sys.stderr)

    results_dir = ANALYSIS_DIR / "editing" / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        ORFEDIT, "scan-batch",
        "--manifest", str(manifest_path),
        "--gene-calls-out", str(results_dir / "edited_gene_calls.tsv"),
        "--edits-out", str(results_dir / "predicted_edits.tsv"),
        "--threads", str(args.threads),
        "--flank", str(args.flank),
        "--edit-penalty", str(args.edit_penalty),
    ]
    proc = subprocess.run(cmd)
    if proc.returncode != 0:
        print(f"[err] orfedit scan-batch failed (exit {proc.returncode})", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
