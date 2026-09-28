#!/usr/bin/env python3
"""Concatenate per-gene alignments into one supermatrix + partition file via AMAS.

Partitioning (not naive concatenation) is the scientifically correct way to
combine genes with different evolutionary rates/models - this is the second
layer of protection around the multi-exon caveat (the first being that
those genes were excluded from the marker set entirely in
01_select_core_genes.py). AMAS also transparently pads any taxon missing
from a given gene's alignment with gap characters, so genes don't all need
identical taxon sampling.

Output: analysis/phylogeny/results/concat/{mito,pltd}_concat.fasta
        analysis/phylogeny/results/concat/{mito,pltd}_partitions.txt (RAxML-style, IQ-TREE compatible)
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import os
import sys
from pathlib import Path

# Defaults to the old same-tree layout (code and data siblings under one
# repo root); set PLANT_ORGANELLE_DATA_ROOT to point at a separate data
# checkout instead (e.g. when this code is installed from its own repo).
ROOT_DIR = Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT", str(Path(__file__).resolve().parents[3])))
ANALYSIS_DIR = ROOT_DIR / "analysis"
AMAS = shutil.which("AMAS.py") or "/software/team301/AMAS/amas/AMAS.py"
PYTHON3 = shutil.which("python3") or "python3"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organelle", default="both", choices=["mito", "pltd", "both"])
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    organelles = ["mito", "pltd"] if args.organelle == "both" else [args.organelle]
    out_dir = ANALYSIS_DIR / "phylogeny" / "results" / "concat"
    out_dir.mkdir(parents=True, exist_ok=True)

    for organelle in organelles:
        aln_dir = ANALYSIS_DIR / "phylogeny" / "results" / "alignments" / organelle
        aln_files = sorted(aln_dir.glob("*.aln.fasta"))
        if not aln_files:
            print(f"[warn] {organelle}: no alignments found in {aln_dir}, skipping", file=sys.stderr)
            continue

        concat_out = out_dir / f"{organelle}_concat.fasta"
        part_out = out_dir / f"{organelle}_partitions.txt"
        if not args.force and concat_out.exists() and part_out.exists() and \
                concat_out.stat().st_mtime >= max(f.stat().st_mtime for f in aln_files):
            print(f"[info] concat_alignment: {organelle}: up to date, skipping", file=sys.stderr)
            continue

        cmd = [PYTHON3, AMAS, "concat", "-i", *[str(f) for f in aln_files], "-f", "fasta", "-d", "dna",
               "-t", str(concat_out), "-p", str(part_out), "-y", "raxml"]
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        if proc.returncode != 0:
            print(f"[err] AMAS concat failed for {organelle}: {proc.stderr.strip()}", file=sys.stderr)
            continue
        print(f"[info] concat_alignment: {organelle}: {len(aln_files)} genes -> {concat_out}", file=sys.stderr)


if __name__ == "__main__":
    main()
