#!/usr/bin/env python3
"""Long-format per-(species, gene) hit table built from oatk's own .ctg.bed gene calls.

.ctg.bed columns (headerless except '#' comment line): seq_name align_from
align_to gene_name score strand. Isoacceptor tRNA names (e.g. trnP-UGG) are
kept distinct, never collapsed to a "trn" family, since isoacceptor identity
matters for a real completeness signal.

Output: analysis/qc_basic_stats/results/gene_matrix_{mito,pltd}.tsv
Columns: species organelle gene n_hits max_score total_aligned_length n_contigs
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

import io_utils  # noqa: E402
import species_discovery as sd  # noqa: E402

COLUMNS = ["species", "organelle", "gene", "n_hits", "max_score", "total_aligned_length", "n_contigs"]
KEY_COLS = ["species", "organelle"]
BED_COLS = ["seq_name", "align_from", "align_to", "gene", "score", "strand"]


def rows_for_species(species: str, organelle: str, ctg_bed: str) -> list[dict]:
    df = pd.read_csv(ctg_bed, sep="\t", comment="#", names=BED_COLS)
    if df.empty:
        return []
    df["aligned_length"] = (df["align_to"] - df["align_from"]).abs()
    grouped = df.groupby("gene").agg(
        n_hits=("gene", "size"),
        max_score=("score", "max"),
        total_aligned_length=("aligned_length", "sum"),
        n_contigs=("seq_name", "nunique"),
    ).reset_index()
    grouped["species"] = species
    grouped["organelle"] = organelle
    return grouped[COLUMNS].to_dict("records")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organelle", default="both", choices=["mito", "pltd", "both"])
    ap.add_argument("--species-list", default=None)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    organelles = ["mito", "pltd"] if args.organelle == "both" else [args.organelle]
    species_filter = sd.load_species_list(Path(args.species_list)) if args.species_list else None

    for organelle in organelles:
        out_path = ANALYSIS_DIR / "qc_basic_stats" / "results" / f"gene_matrix_{organelle}.tsv"
        existing_df, existing_keys = io_utils.load_existing(out_path, KEY_COLS)

        data_root = sd.repo_data_root(ROOT_DIR, organelle)
        resolved = sd.discover_all(data_root, organelle, species_filter=species_filter)

        new_rows = []
        stale_keys = []  # species that HAD a ctg_bed in a prior run but don't now
                          # (e.g. archived by a promotion that didn't re-annotate) -
                          # must be purged explicitly, since merge_rows only replaces
                          # keys present in new_rows and would otherwise leave this
                          # species' old gene calls in gene_matrix.tsv forever, silently
                          # feeding a stale core_gene_pct into qc_summary.
        n_computed = n_skipped = n_no_bed = 0
        for r in resolved:
            key = (r.species, organelle)
            if not r.ctg_bed:
                n_no_bed += 1
                if key in existing_keys:
                    stale_keys.append(key)
                continue
            if io_utils.needs_recompute(key, existing_keys, out_path, Path(r.ctg_bed), args.force):
                new_rows.extend(rows_for_species(r.species, organelle, r.ctg_bed))
                n_computed += 1
            else:
                n_skipped += 1

        if stale_keys and not existing_df.empty:
            existing_df = existing_df.set_index(KEY_COLS).drop(index=stale_keys, errors="ignore").reset_index()

        merged = io_utils.merge_rows(existing_df, new_rows, KEY_COLS, COLUMNS)
        io_utils.atomic_write_tsv(out_path, merged.to_dict("records"), COLUMNS)
        print(f"[info] gene_matrix_{organelle}: computed={n_computed} skipped={n_skipped} "
              f"no_bed={n_no_bed} purged_stale={len(stale_keys)} total_rows={len(merged)} -> {out_path}",
              file=sys.stderr)


if __name__ == "__main__":
    main()
