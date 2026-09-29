#!/usr/bin/env python3
"""Flag repeat copies near a contig's end as putative recombination-mediating
repeats - the classic mechanism for plant mitochondrial multipartite/
flip-flop structure. Rewrites each per-species repeats.tsv in place with
two extra columns.

A repeat near the end of a LINEAR contig is a much stronger candidate than
one near the end of a CIRCULAR contig, where "the end" is just an arbitrary
linearisation point rather than a real sequence boundary - both are flagged,
but circular ones get a note saying so, so students don't over-interpret it.

Output (in place): analysis/repeats/results/per_species/<Species>.<organelle>.repeats.tsv
Adds columns: distance_to_contig_end putative_recomb_repeat note
"""
from __future__ import annotations

import os
import argparse
import re
import sys
from pathlib import Path

import pandas as pd

# CODE_DIR is this checkout's analysis/ (shared modules, bundled reference
# data); data/ and every module's results/ and work/ live under the data
# root - PLANT_ORGANELLE_DATA_ROOT, or this checkout if unset.
CODE_DIR = Path(__file__).resolve().parents[2]
ROOT_DIR = Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT", str(CODE_DIR.parent)))
ANALYSIS_DIR = ROOT_DIR / "analysis"
PER_SPECIES_DIR = ANALYSIS_DIR / "repeats" / "results" / "per_species"

FILENAME_RE = re.compile(r"^(.+)\.(mito|pltd)\.repeats\.tsv$")


def load_contig_lengths() -> pd.DataFrame:
    path = ANALYSIS_DIR / "qc_basic_stats" / "results" / "contig_stats.tsv"
    df = pd.read_csv(path, sep="\t")
    return df[df["contig_source"] == "ctg_fasta"][["species", "organelle", "contig_id", "length", "circular"]]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--species-list", default=None)
    ap.add_argument("--end-distance", type=int, default=2000,
                     help="max distance (bp) from a contig end to flag as putative_recomb_repeat (default 2000)")
    args = ap.parse_args()

    species_filter = None
    if args.species_list:
        species_filter = {l.strip() for l in open(args.species_list) if l.strip() and not l.startswith("#")}

    contig_lengths = load_contig_lengths()
    n_processed = 0
    for path in sorted(PER_SPECIES_DIR.glob("*.repeats.tsv")):
        m = FILENAME_RE.match(path.name)
        if not m:
            continue
        species, organelle = m.groups()
        if species_filter is not None and species not in species_filter:
            continue

        df = pd.read_csv(path, sep="\t")
        if df.empty or "distance_to_contig_end" in df.columns:
            continue  # already flagged, or nothing to flag

        cl = contig_lengths[(contig_lengths["species"] == species) & (contig_lengths["organelle"] == organelle)]
        cl = cl.set_index("contig_id")

        def compute(row):
            if row.contig_id not in cl.index:
                return pd.Series({"distance_to_contig_end": None, "putative_recomb_repeat": False, "note": "contig length unknown"})
            length = cl.loc[row.contig_id, "length"]
            circular = str(cl.loc[row.contig_id, "circular"]).lower() == "true"
            dist = min(row.start, length - row.end)
            flagged = dist <= args.end_distance
            note = ""
            if flagged:
                note = ("circular contig - 'end' is an arbitrary linearisation point, weaker signal"
                        if circular else
                        "linear contig - repeat near a true sequence end, stronger recombination candidate")
            return pd.Series({"distance_to_contig_end": dist, "putative_recomb_repeat": flagged, "note": note})

        flags = df.apply(compute, axis=1)
        df = pd.concat([df, flags], axis=1)
        df.to_csv(path, sep="\t", index=False)
        n_processed += 1

    print(f"[info] flag_recomb_repeats: processed {n_processed} per-species files", file=sys.stderr)


if __name__ == "__main__":
    main()
