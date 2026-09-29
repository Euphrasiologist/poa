#!/usr/bin/env python3
"""Build a summary table from completed 03_recruit_and_resolve.sh outputs.

Reads every analysis/linearize/results/resolve/<organelle>/<species>.resolve.fasta
and parses gfatk resolve's own FASTA headers. Two header shapes seen on real
output:
  - a resolved circuit:  `>Species.k1001.s31.c80.pltd_circuit1:n_segments=4:bp=153046`
  - a bubble arm (an alternative path at a branch point that read evidence
    didn't select for the main circuit, kept for reference with its own read
    support - see `gfatk resolve --help`'s --no-bubble-fasta):
    `>bubble_arm:u50:between=u59-_to_u49-:support_reads=116`
Only circuit headers count toward n_circuits/length_bp; bubble arms are
counted separately (n_bubble_arms) as an informational complexity signal,
not treated as unparsed/unexpected.

Cross-references against work/resolve_targets.tsv so species that were
submitted but produced no resolved circuit at all (empty/missing output -
GraphAligner found nothing alignable, or gfatk resolve's solver couldn't
close a single-circuit constraint) show up as a row with n_circuits=0 rather
than silently vanishing.

Output: results/resolve_summary.tsv
"""
from __future__ import annotations

import os
import argparse
import re
import sys
from pathlib import Path

# CODE_DIR is the installed poa/pipeline/ (shared modules, bundled reference
# data); data/ and every module's results/ and work/ live under the data
# root - PLANT_ORGANELLE_DATA_ROOT, or the current directory if unset.
CODE_DIR = Path(__file__).resolve().parents[2]
ROOT_DIR = Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT") or os.getcwd())
ANALYSIS_DIR = ROOT_DIR / "analysis"

CIRCUIT_RE = re.compile(r"^.+_circuit(?P<circuit>\d+):n_segments=(?P<n_segments>\d+):bp=(?P<bp>\d+)$")
BUBBLE_ARM_RE = re.compile(r"^bubble_arm:")

COLUMNS = ["species", "organelle", "n_circuits", "n_bubble_arms", "circuit_id", "n_segments", "length_bp"]


def parse_resolve_fasta(path: Path) -> tuple[list[dict], int]:
    """Return (circuit rows, n_bubble_arms)."""
    circuits, n_bubble_arms = [], 0
    if not path.exists() or path.stat().st_size == 0:
        return circuits, n_bubble_arms
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            if BUBBLE_ARM_RE.match(header):
                n_bubble_arms += 1
                continue
            m = CIRCUIT_RE.match(header)
            if not m:
                print(f"[warn] unparsed gfatk resolve header in {path}: {header!r}", file=sys.stderr)
                continue
            circuits.append({"circuit_id": m.group("circuit"), "n_segments": int(m.group("n_segments")),
                              "length_bp": int(m.group("bp"))})
    return circuits, n_bubble_arms


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--targets", default=str(ANALYSIS_DIR / "linearize" / "work" / "resolve_targets.tsv"))
    ap.add_argument("--out", default=str(ANALYSIS_DIR / "linearize" / "results" / "resolve_summary.tsv"))
    args = ap.parse_args()

    targets_path = Path(args.targets)
    if not targets_path.exists():
        print(f"[err] targets not found: {targets_path} (run 02_select_resolve_targets.py first)", file=sys.stderr)
        sys.exit(1)

    resolve_dir = ANALYSIS_DIR / "linearize" / "results" / "resolve"
    rows = []
    n_resolved = n_unresolved = n_missing = 0

    with open(targets_path) as fh:
        next(fh)  # header
        for line in fh:
            species, organelle, _gfa = line.rstrip("\n").split("\t")
            fasta_path = resolve_dir / organelle / f"{species}.resolve.fasta"
            if not fasta_path.exists():
                n_missing += 1
                continue
            circuits, n_bubble_arms = parse_resolve_fasta(fasta_path)
            if not circuits:
                rows.append({"species": species, "organelle": organelle, "n_circuits": 0,
                             "n_bubble_arms": n_bubble_arms, "circuit_id": None,
                             "n_segments": None, "length_bp": None})
                n_unresolved += 1
                continue
            for c in circuits:
                rows.append({"species": species, "organelle": organelle, "n_circuits": len(circuits),
                             "n_bubble_arms": n_bubble_arms, **c})
            n_resolved += 1

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    import pandas as pd
    df = pd.DataFrame(rows, columns=COLUMNS)
    df.to_csv(out_path, sep="\t", index=False)

    print(f"[info] summarize_resolve: resolved={n_resolved} unresolved={n_unresolved} "
          f"not_yet_run={n_missing} -> {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
