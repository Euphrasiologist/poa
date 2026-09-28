#!/usr/bin/env python3
"""Summarize each species' mito-vs-plastid PAF into MTPT burden metrics.

Operationalizes the exercises already suggested in the top-level README
("quantify total aligned length... compare across species... test whether
MTPT burden correlates with genome size") as one reusable table instead of
a one-off snippet.

Output: analysis/synteny/results/mtpt_summary.tsv
Columns: species total_aligned_length_bp n_alignment_blocks longest_block_bp
         pct_plastid_genome_in_mito mito_genome_length pltd_genome_length
"""
from __future__ import annotations

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
import minimap_utils  # noqa: E402

MTPT_DIR = ANALYSIS_DIR / "synteny" / "results" / "mtpt"


def genome_lengths() -> pd.Series:
    """Length of the actual RESOLVED contig(s) that were aligned (not the raw,
    pre-path-resolution unitig graph total from gfa_stats.tsv, which can
    legitimately differ - path resolution may revisit/duplicate unitigs)."""
    path = ANALYSIS_DIR / "qc_basic_stats" / "results" / "contig_stats.tsv"
    df = pd.read_csv(path, sep="\t")
    resolved = df[df["contig_source"] == "ctg_fasta"]
    return resolved.groupby(["species", "organelle"])["length"].sum()


def main():
    lengths = genome_lengths()
    rows = []
    for paf_path in sorted(MTPT_DIR.glob("*.mito_vs_pltd.paf")):
        species = paf_path.name.removesuffix(".mito_vs_pltd.paf")
        records = []
        for line in paf_path.read_text().splitlines():
            fields = line.split("\t")
            if len(fields) < 12:
                continue
            rec = dict(zip(minimap_utils.PAF_COLS, fields[:12]))
            for k in ("query_start", "query_end", "target_start", "target_end"):
                rec[k] = int(rec[k])
            records.append(rec)

        total_aligned = minimap_utils.total_aligned_length(records, side="query")  # plastid side (query)
        n_blocks = len(records)
        longest = max((r["query_end"] - r["query_start"] for r in records), default=0)
        pltd_len = lengths.get((species, "pltd"))
        mito_len = lengths.get((species, "mito"))
        pct_in_mito = (total_aligned / pltd_len) if pltd_len else None

        rows.append({
            "species": species, "total_aligned_length_bp": total_aligned, "n_alignment_blocks": n_blocks,
            "longest_block_bp": longest, "pct_plastid_genome_in_mito": pct_in_mito,
            "mito_genome_length": mito_len, "pltd_genome_length": pltd_len,
        })

    out_path = ANALYSIS_DIR / "synteny" / "results" / "mtpt_summary.tsv"
    columns = ["species", "total_aligned_length_bp", "n_alignment_blocks", "longest_block_bp",
               "pct_plastid_genome_in_mito", "mito_genome_length", "pltd_genome_length"]
    pd.DataFrame(rows, columns=columns).to_csv(out_path, sep="\t", index=False)
    print(f"[info] mtpt_summary: {len(rows)} species -> {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
