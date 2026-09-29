#!/usr/bin/env python3
"""Align each marker gene's per-species FASTA with mafft --auto.

Cheap: only ~15-30 short organelle genes per organelle, no LSF needed.

Output: analysis/phylogeny/results/alignments/{mito,pltd}/<gene>.aln.fasta
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

# CODE_DIR is this checkout's analysis/ (shared modules, bundled reference
# data); data/ and every module's results/ and work/ live under the data
# root - PLANT_ORGANELLE_DATA_ROOT, or this checkout if unset.
CODE_DIR = Path(__file__).resolve().parents[2]
ROOT_DIR = Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT", str(CODE_DIR.parent)))
ANALYSIS_DIR = ROOT_DIR / "analysis"
# The top-level /software/team301/mafft wrapper is broken (hardcoded to a
# nonexistent /usr/local/libexec/mafft). The working install is
# mafft-7.525-with-extensions/core/mafft, with MAFFT_BINARIES pointed at
# that same core/ dir - verified directly (see analysis/common/tool_paths.sh).
MAFFT = os.environ.get("MAFFT") or shutil.which("mafft") or "/software/team301/mafft-7.525-with-extensions/core/mafft"
MAFFT_ENV = {**os.environ, "MAFFT_BINARIES": os.environ.get(
    "MAFFT_BINARIES", "/software/team301/mafft-7.525-with-extensions/core")}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organelle", default="both", choices=["mito", "pltd", "both"])
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    organelles = ["mito", "pltd"] if args.organelle == "both" else [args.organelle]
    n_aligned = n_skipped = n_failed = 0
    for organelle in organelles:
        marker_path = ANALYSIS_DIR / "phylogeny" / "results" / f"marker_genes_{organelle}.txt"
        if not marker_path.exists():
            print(f"[err] missing {marker_path} - run 01_select_core_genes.py first", file=sys.stderr)
            continue
        genes = [g for g in marker_path.read_text().splitlines() if g]

        genes_dir = ANALYSIS_DIR / "annotation" / "results" / "genes" / organelle
        out_dir = ANALYSIS_DIR / "phylogeny" / "results" / "alignments" / organelle
        out_dir.mkdir(parents=True, exist_ok=True)

        for gene in genes:
            in_fasta = genes_dir / f"{gene}.fasta"
            out_fasta = out_dir / f"{gene}.aln.fasta"
            if not in_fasta.exists():
                print(f"[warn] {organelle}/{gene}: missing {in_fasta}, skipping", file=sys.stderr)
                continue
            if not args.force and out_fasta.exists() and out_fasta.stat().st_mtime >= in_fasta.stat().st_mtime:
                n_skipped += 1
                continue
            proc = subprocess.run([MAFFT, "--auto", "--quiet", str(in_fasta)],
                                   capture_output=True, text=True, timeout=600, env=MAFFT_ENV)
            if proc.returncode != 0 or not proc.stdout.strip():
                print(f"[warn] mafft failed for {organelle}/{gene}: {proc.stderr.strip()[:200]}", file=sys.stderr)
                n_failed += 1
                continue
            out_fasta.write_text(proc.stdout)
            n_aligned += 1

    print(f"[info] align_genes: aligned={n_aligned} skipped={n_skipped} failed={n_failed}", file=sys.stderr)


if __name__ == "__main__":
    main()
