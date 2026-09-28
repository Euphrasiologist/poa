#!/usr/bin/env python3
"""Select which species x organelle should get the heavier gfatk-resolve
fallback tier (read-path evidence via GraphAligner), as opposed to relying
on the cheap graph-only 01_gfatk_linear.py alone.

Scope, by design, is narrow: whole-genome raw PacBio HiFi read sets for
this dataset total ~24TB compressed (verified directly, du across
meta/file.txt's ~3250 files) - running GraphAligner against the full read
set for every species would be a colossal, unjustified compute cost for a
teaching resource that's explicitly meant to stay lightweight. Read-path
evidence is only worth that cost where oatk's own Pathfinder produced
*nothing at all* to work with (no_resolved_ctg_fasta) - gfatk linear (tier
1) is a genuine independent alternative there, but it's a coverage-greedy
heuristic with no read evidence behind it either, so a species that still
has no usable linearization after tier 1 is exactly where the extra cost
of read recruitment + GraphAligner + gfatk resolve is justified.

Reads qc_summary.tsv (status_reasons contains "no_resolved_ctg_fasta") and
writes analysis/linearize/work/resolve_targets.tsv: species, organelle, gfa.
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

# Defaults to the old same-tree layout (code and data siblings under one
# repo root); set PLANT_ORGANELLE_DATA_ROOT to point at a separate data
# checkout instead (e.g. when this code is installed from its own repo).
ROOT_DIR = Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT", str(Path(__file__).resolve().parents[3])))
ANALYSIS_DIR = ROOT_DIR / "analysis"
sys.path.insert(0, str(ANALYSIS_DIR / "common"))

import pandas as pd  # noqa: E402

import species_discovery as sd  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--qc-reason", default="no_resolved_ctg_fasta",
                     help="only species x organelle whose status_reasons contains this substring (default: no_resolved_ctg_fasta)")
    ap.add_argument("--species-list", default=None)
    ap.add_argument("--out", default=str(ANALYSIS_DIR / "linearize" / "work" / "resolve_targets.tsv"))
    args = ap.parse_args()

    qc_path = ANALYSIS_DIR / "qc_basic_stats" / "results" / "qc_summary.tsv"
    qc = pd.read_csv(qc_path, sep="\t")
    qc["status_reasons"] = qc["status_reasons"].fillna("")
    hits = qc[qc["status_reasons"].str.contains(args.qc_reason, regex=False)]

    species_filter = sd.load_species_list(Path(args.species_list)) if args.species_list else None

    rows = []
    n_no_gfa = 0
    for _, row in hits.iterrows():
        species, organelle = row["species"], row["organelle"]
        if species_filter is not None and species not in species_filter:
            continue
        data_root = sd.repo_data_root(ROOT_DIR, organelle)
        sf = sd.resolve_species(data_root / species, organelle)
        if sf.gfa is None:
            n_no_gfa += 1
            print(f"[warn] {species} ({organelle}): matched '{args.qc_reason}' but has no GFA - skipped",
                  file=sys.stderr)
            continue
        rows.append({"species": species, "organelle": organelle, "gfa": sf.gfa})

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=["species", "organelle", "gfa"]).to_csv(out_path, sep="\t", index=False)
    print(f"[info] select_resolve_targets: matched={len(hits)} usable={len(rows)} no_gfa={n_no_gfa} "
          f"-> {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
