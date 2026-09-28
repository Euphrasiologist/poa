#!/usr/bin/env python3
"""Aggregate per-species repeat tables into one summary row per species x organelle.

Meant to be joined against qc_basic_stats/results/gfa_stats.tsv (dead-end/
subgraph counts) for the natural "does repeat content predict assembly
graph complexity" mini-project.

Output: analysis/repeats/results/recomb_repeats_summary.tsv
Columns: species organelle n_repeat_families n_repeat_copies
         total_repeat_length pct_genome_repetitive n_putative_recomb_repeats
"""
from __future__ import annotations

import re
import os
import sys
from pathlib import Path

import pandas as pd

# Defaults to the old same-tree layout (code and data siblings under one
# repo root); set PLANT_ORGANELLE_DATA_ROOT to point at a separate data
# checkout instead (e.g. when this code is installed from its own repo).
ROOT_DIR = Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT", str(Path(__file__).resolve().parents[3])))
ANALYSIS_DIR = ROOT_DIR / "analysis"
PER_SPECIES_DIR = ANALYSIS_DIR / "repeats" / "results" / "per_species"
FILENAME_RE = re.compile(r"^(.+)\.(mito|pltd)\.repeats\.tsv$")


def genome_lengths() -> pd.DataFrame:
    path = ANALYSIS_DIR / "qc_basic_stats" / "results" / "gfa_stats.tsv"
    df = pd.read_csv(path, sep="\t")
    ok = df[df["parse_status"] == "ok"]
    return ok.groupby(["species", "organelle"])["total_length"].sum().rename("genome_length").reset_index()


def main():
    lengths = genome_lengths().set_index(["species", "organelle"])["genome_length"]

    rows = []
    for path in sorted(PER_SPECIES_DIR.glob("*.repeats.tsv")):
        m = FILENAME_RE.match(path.name)
        if not m:
            continue
        species, organelle = m.groups()
        df = pd.read_csv(path, sep="\t")

        n_families = df["repeat_family"].nunique() if not df.empty else 0
        n_copies = len(df)
        total_len = int(df["copy_length"].sum()) if not df.empty else 0
        genome_len = lengths.get((species, organelle))
        pct_repetitive = (total_len / genome_len) if genome_len else None
        n_recomb = int(df["putative_recomb_repeat"].sum()) if "putative_recomb_repeat" in df.columns else None

        rows.append({
            "species": species, "organelle": organelle, "n_repeat_families": n_families,
            "n_repeat_copies": n_copies, "total_repeat_length": total_len,
            "pct_genome_repetitive": pct_repetitive, "n_putative_recomb_repeats": n_recomb,
        })

    out_path = ANALYSIS_DIR / "repeats" / "results" / "recomb_repeats_summary.tsv"
    out_df = pd.DataFrame(rows, columns=["species", "organelle", "n_repeat_families", "n_repeat_copies",
                                          "total_repeat_length", "pct_genome_repetitive", "n_putative_recomb_repeats"])
    out_df.to_csv(out_path, sep="\t", index=False)
    print(f"[info] repeats_summary: {len(out_df)} species x organelle rows -> {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
