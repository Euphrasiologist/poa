#!/usr/bin/env python3
"""Build a partitioned ML tree with iqtree2, one run per organelle (not per
species) - genuinely cheap even at ~1200 taxa since the alignment is a
handful of short organelle marker genes, not whole genomes. No treePL
dating step (needs external fossil calibration - out of scope for a
lightweight teaching resource); ultrafast bootstrap gives support values
for free.

Output: analysis/phylogeny/results/trees/{mito,pltd}.treefile (+ iqtree2 side files)
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
IQTREE2 = shutil.which("iqtree2") or "/software/team301/iqtree-2.4.0-Linux-intel/bin/iqtree2"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organelle", default="both", choices=["mito", "pltd", "both"])
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--threads", default="2",
                     help="passed to iqtree2 -nt (default 2 - deliberately conservative since this "
                          "is often run interactively on a shared login node; pass -nt AUTO or a "
                          "higher number yourself inside an LSF job)")
    args = ap.parse_args()

    organelles = ["mito", "pltd"] if args.organelle == "both" else [args.organelle]
    concat_dir = ANALYSIS_DIR / "phylogeny" / "results" / "concat"
    trees_dir = ANALYSIS_DIR / "phylogeny" / "results" / "trees"
    trees_dir.mkdir(parents=True, exist_ok=True)

    for organelle in organelles:
        concat_fasta = concat_dir / f"{organelle}_concat.fasta"
        partitions = concat_dir / f"{organelle}_partitions.txt"
        if not concat_fasta.exists() or not partitions.exists():
            print(f"[err] missing concat alignment/partitions for {organelle} - run 03_concat_alignment.py first",
                  file=sys.stderr)
            continue

        prefix = trees_dir / organelle
        treefile = prefix.with_suffix(".treefile")
        log_path = prefix.with_suffix(".log")
        ckp_path = prefix.with_suffix(".ckp.gz")

        # A run that finished cleanly writes "Analysis results written to:" near the
        # end of its .log; a run killed externally (e.g. LSF runlimit) does not, even
        # though it may have left a .treefile behind - that file is just the last
        # snapshot of the best tree seen so far, not a finished result. Don't trust
        # treefile-mtime alone to mean "done".
        completed = (log_path.exists() and "Analysis results written to:" in log_path.read_text()
                     and treefile.exists() and treefile.stat().st_mtime >= concat_fasta.stat().st_mtime)
        if not args.force and completed:
            print(f"[info] build_tree: {organelle}: up to date, skipping", file=sys.stderr)
            continue

        # Resume from an interrupted-but-still-valid run (checkpoint newer than the
        # current alignment - i.e. it was started against this same alignment, not a
        # stale one from before the marker-gene set changed) rather than discarding
        # real progress: only pass -redo when there's nothing valid to resume from.
        # Verified necessary: a previous version always passed -redo unconditionally,
        # which would have silently thrown away a real, mostly-converged tree search
        # (killed at iteration 120 by an LSF runlimit, not a crash) had this not been fixed.
        can_resume = (not args.force and ckp_path.exists()
                      and ckp_path.stat().st_mtime >= concat_fasta.stat().st_mtime)

        cmd = [IQTREE2, "-s", str(concat_fasta), "-p", str(partitions), "-m", "MFP",
               "-bb", "1000", "-nt", args.threads, "-pre", str(prefix)]
        if not can_resume:
            cmd.append("-redo")
        else:
            print(f"[info] build_tree: {organelle}: resuming from checkpoint {ckp_path}", file=sys.stderr)

        # No subprocess-level timeout: MFP model selection + 1000 ultrafast bootstrap
        # replicates on the full dataset (hundreds of taxa, dozens of partitions) can
        # legitimately run for hours - a fixed timeout here previously killed a real,
        # still-progressing run and crashed the whole pipeline. The enclosing LSF job's
        # own runlimit is the correct backstop, not a number guessed inside this script.
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            print(f"[err] iqtree2 failed for {organelle}: {proc.stderr.strip()[-2000:]}", file=sys.stderr)
            continue
        print(f"[info] build_tree: {organelle} -> {treefile}", file=sys.stderr)


if __name__ == "__main__":
    main()
