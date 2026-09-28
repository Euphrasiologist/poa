#!/usr/bin/env python3
"""Reshape oatk's own .ctg.bed gene calls into one tidy long table.

No new annotation is performed here - oatk already annotates every contig
(that's what .ctg.bed / .annot_*.txt already are); this just adds
species/organelle provenance columns and applies the QC filter so
downstream modules read one trusted table instead of re-deriving this
themselves.

Output: analysis/annotation/results/gene_calls.tsv
Columns: species organelle contig_id start end gene score strand
"""
from __future__ import annotations

import argparse
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

BED_COLS = ["contig_id", "start", "end", "gene", "score", "strand"]
COLUMNS = ["species", "organelle"] + BED_COLS


def qc_eligible_species(qc_status: list[str]) -> set[tuple[str, str]]:
    qc_path = ANALYSIS_DIR / "qc_basic_stats" / "results" / "qc_summary.tsv"
    if not qc_path.exists():
        print(f"[err] missing {qc_path} - run qc_basic_stats/src/05_qc_summary.py first", file=sys.stderr)
        sys.exit(1)
    qc = pd.read_csv(qc_path, sep="\t")
    sub = qc[qc["status"].isin(qc_status)]
    return set(zip(sub["species"], sub["organelle"]))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organelle", default="both", choices=["mito", "pltd", "both"])
    ap.add_argument("--species-list", default=None)
    ap.add_argument("--qc-status", default="pass", help="comma-separated subset of pass,flag,fail")
    ap.add_argument("--force", action="store_true", help="accepted for CLI consistency with run_all.sh; "
                     "this script always does a full recompute (cheap), so it's a no-op")
    args = ap.parse_args()

    qc_status = args.qc_status.split(",")
    eligible = qc_eligible_species(qc_status)

    organelles = ["mito", "pltd"] if args.organelle == "both" else [args.organelle]
    species_filter = sd.load_species_list(Path(args.species_list)) if args.species_list else None

    all_rows = []
    for organelle in organelles:
        data_root = sd.repo_data_root(ROOT_DIR, organelle)
        resolved = sd.discover_all(data_root, organelle, species_filter=species_filter)
        n_used = 0
        for r in resolved:
            if (r.species, organelle) not in eligible or not r.ctg_bed:
                continue
            df = pd.read_csv(r.ctg_bed, sep="\t", comment="#", names=BED_COLS)
            if df.empty:
                continue
            df["species"] = r.species
            df["organelle"] = organelle
            all_rows.append(df[COLUMNS])
            n_used += 1
        print(f"[info] gene_table: {organelle}: {n_used} QC-eligible species with gene calls", file=sys.stderr)

    out_path = ANALYSIS_DIR / "annotation" / "results" / "gene_calls.tsv"
    out_df = pd.concat(all_rows, ignore_index=True) if all_rows else pd.DataFrame(columns=COLUMNS)
    out_df.to_csv(out_path, sep="\t", index=False)
    print(f"[info] gene_table: {len(out_df)} rows -> {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
