#!/usr/bin/env python3
"""Compute per-gene prevalence across eligible (non-empty-assembly) species.

Two derived, data-driven thresholds are produced (never hardcoded gene
lists) and used downstream:
  - QC core-gene set:        pct_present   >= 0.90  (used for completeness scoring)
  - phylogeny marker-gene set: pct_present >= 0.90 AND pct_single_copy >= 0.90
    (this mechanically excludes multi-exon/trans-spliced genes like
    nad1/2/4/5/7, cox2, ccmFc/Fn, rps3 for mito or clpP/ycf3/rps12 for
    plastid, since they show up with multiple BED hits per species -
    no hardcoded blacklist needed).

"Eligible" = species whose gfa_stats.tsv parse_status is 'ok' for that
organelle (i.e. a non-empty assembly graph exists) - this is a cheap,
already-computed filter and intentionally does NOT require QC pass/fail
status, since this script's own output feeds into computing QC status.

Output: analysis/qc_basic_stats/results/core_genes_{mito,pltd}.tsv
Columns: gene pct_present pct_single_copy n_species_present n_species_single_copy n_eligible_species
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
COLUMNS = ["gene", "pct_present", "pct_single_copy", "n_species_present", "n_species_single_copy", "n_eligible_species"]


def eligible_species(gfa_stats: pd.DataFrame, organelle: str) -> set:
    sub = gfa_stats[(gfa_stats["organelle"] == organelle) & (gfa_stats["parse_status"] == "ok")]
    return set(sub["species"].unique())


def compute_core_genes(gene_matrix: pd.DataFrame, eligible: set) -> pd.DataFrame:
    n_eligible = len(eligible)
    sub = gene_matrix[gene_matrix["species"].isin(eligible)]
    if sub.empty or n_eligible == 0:
        return pd.DataFrame(columns=COLUMNS)

    sub = sub.copy()
    sub["is_single_copy"] = (sub["n_hits"] == 1) & (sub["n_contigs"] == 1)

    grouped = sub.groupby("gene").agg(
        n_species_present=("species", "nunique"),
        n_species_single_copy=("is_single_copy", "sum"),
    ).reset_index()
    grouped["n_eligible_species"] = n_eligible
    grouped["pct_present"] = grouped["n_species_present"] / n_eligible
    grouped["pct_single_copy"] = grouped["n_species_single_copy"] / n_eligible
    return grouped[COLUMNS].sort_values("pct_present", ascending=False)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organelle", default="both", choices=["mito", "pltd", "both"])
    args = ap.parse_args()

    results_dir = ANALYSIS_DIR / "qc_basic_stats" / "results"
    gfa_stats_path = results_dir / "gfa_stats.tsv"
    if not gfa_stats_path.exists():
        print(f"[err] missing {gfa_stats_path} - run 01_gfa_stats.py first", file=sys.stderr)
        sys.exit(1)
    gfa_stats = pd.read_csv(gfa_stats_path, sep="\t", dtype=str)

    organelles = ["mito", "pltd"] if args.organelle == "both" else [args.organelle]
    for organelle in organelles:
        gm_path = results_dir / f"gene_matrix_{organelle}.tsv"
        if not gm_path.exists():
            print(f"[err] missing {gm_path} - run 03_gene_matrix.py first", file=sys.stderr)
            continue
        gene_matrix = pd.read_csv(gm_path, sep="\t")
        eligible = eligible_species(gfa_stats, organelle)
        core = compute_core_genes(gene_matrix, eligible)

        out_path = results_dir / f"core_genes_{organelle}.tsv"
        core.to_csv(out_path, sep="\t", index=False)

        n_qc_core = (core["pct_present"] >= 0.90).sum()
        n_marker = ((core["pct_present"] >= 0.90) & (core["pct_single_copy"] >= 0.90)).sum()
        print(f"[info] core_genes_{organelle}: n_eligible_species={len(eligible)} "
              f"qc_core_genes(>=90% present)={n_qc_core} phylogeny_marker_genes(>=90% present & single-copy)={n_marker} "
              f"-> {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
