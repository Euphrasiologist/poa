#!/usr/bin/env python3
"""Finalize the marker-gene list actually used for the phylogeny, from
qc_basic_stats/results/core_genes_{organelle}.tsv (same threshold as
annotation/02_gene_fasta_extract.py: pct_present>=0.90 AND
pct_single_copy>=0.90). Also writes an explicit "excluded because
multi-exon/fragmented" note for teaching value, rather than silently
dropping those genes.

Output:
  analysis/phylogeny/results/marker_genes_{mito,pltd}.txt
  analysis/phylogeny/results/excluded_genes_{mito,pltd}.txt
"""
from __future__ import annotations

import os
import argparse
import sys
from pathlib import Path

import pandas as pd

# CODE_DIR is this checkout's analysis/ (shared modules, bundled reference
# data); data/ and every module's results/ and work/ live under the data
# root - PLANT_ORGANELLE_DATA_ROOT, or this checkout if unset.
CODE_DIR = Path(__file__).resolve().parents[2]
ROOT_DIR = Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT", str(CODE_DIR.parent)))
ANALYSIS_DIR = ROOT_DIR / "analysis"
PRESENT_THRESHOLD = 0.90
# See analysis/annotation/src/02_gene_fasta_extract.py for why this is 0.80, not 0.90:
# verified at full-dataset scale that mito genes top out around 89% single-copy
# (real biology - multipartite duplication - not a bug), so 90% excluded every
# mito gene entirely.
SINGLE_COPY_THRESHOLD = 0.80


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organelle", default="both", choices=["mito", "pltd", "both"])
    args = ap.parse_args()

    organelles = ["mito", "pltd"] if args.organelle == "both" else [args.organelle]
    results_dir = ANALYSIS_DIR / "phylogeny" / "results"
    results_dir.mkdir(parents=True, exist_ok=True)

    for organelle in organelles:
        core_path = ANALYSIS_DIR / "qc_basic_stats" / "results" / f"core_genes_{organelle}.tsv"
        if not core_path.exists():
            print(f"[err] missing {core_path} - run qc_basic_stats/src/04_core_gene_list.py first", file=sys.stderr)
            continue
        df = pd.read_csv(core_path, sep="\t")

        present = df[df["pct_present"] >= PRESENT_THRESHOLD]
        marker = present[present["pct_single_copy"] >= SINGLE_COPY_THRESHOLD].sort_values("gene")
        excluded = present[present["pct_single_copy"] < SINGLE_COPY_THRESHOLD].sort_values("gene")

        marker_path = results_dir / f"marker_genes_{organelle}.txt"
        marker_path.write_text("\n".join(marker["gene"]) + "\n" if len(marker) else "")

        excluded_path = results_dir / f"excluded_genes_{organelle}.txt"
        with open(excluded_path, "w") as fh:
            fh.write("# genes present in >=90% of species but excluded from the phylogeny marker set\n")
            fh.write("# because they are multi-copy/fragmented in >=10% of species where present\n")
            fh.write("# (typically multi-exon/trans-spliced genes - naive concatenation across\n")
            fh.write("# exon fragments would be wrong, so they're left out rather than mishandled)\n")
            fh.write("gene\tpct_present\tpct_single_copy\n")
            for _, row in excluded.iterrows():
                fh.write(f"{row['gene']}\t{row['pct_present']:.3f}\t{row['pct_single_copy']:.3f}\n")

        print(f"[info] select_core_genes: {organelle}: {len(marker)} marker genes, "
              f"{len(excluded)} excluded (multi-exon/fragmented) -> {marker_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
