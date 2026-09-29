#!/usr/bin/env python3
"""Stage 2 of the reference-profile bootstrap (see `00_build_exon_profiles.py`
and `../README.md`): rebuild each gene's whole-gene profile from THIS
dataset's own high-confidence reconstructions instead of the 5-7 GenBank
genomes alone.

Why this matters concretely, not just in principle: `rps3`'s first-pass
whole-gene profile (5 GenBank references, one fast-evolving/size-variable
gene per the literature review) gave strongly *negative* whole-gene scores
even for reconstructions whose own per-exon evidence looked solid (both
exons independently high-scoring, correct assignment, no internal stops) -
a real quality issue in a small reference-only profile, not (necessarily)
a wrong reconstruction. Once `01_reconstruct.py` has run once, this script
uses its own successful output as a much larger, dataset-scale training set
for a second pass.

Quality gate for "successful" (deliberately conservative - this feeds back
into a profile other reconstructions will be scored against, so a wrong
inclusion here would propagate): `complete=True` (every exon slot filled)
AND `whole_gene_score > 0` (positively identifiable as this gene, not just
"the least bad option"). Gap characters ('-', i.e. profile columns the
first-pass alignment couldn't fill) are stripped before use - they're
alignment bookkeeping, not real residues.

Output: overwrites `../reference/profiles/<gene>/whole_gene.pssm` in place
(per-exon-slot profiles are untouched - the GenBank exon boundaries/order
they encode are structural ground truth, not something to re-derive from
this tool's own output). Re-run `01_reconstruct.py` afterwards to benefit.
"""
from __future__ import annotations

import os
import argparse
import math
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path

import pandas as pd

# CODE_DIR is this checkout's analysis/ (shared modules, bundled reference
# data); data/ and every module's results/ and work/ live under the data
# root - PLANT_ORGANELLE_DATA_ROOT, or this checkout if unset.
CODE_DIR = Path(__file__).resolve().parents[2]
ROOT_DIR = Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT", str(CODE_DIR.parent)))
ANALYSIS_DIR = ROOT_DIR / "analysis"
REFERENCE_DIR = CODE_DIR / "trans_splicing" / "reference"
PROFILES_DIR = REFERENCE_DIR / "profiles"
WORK_DIR = ANALYSIS_DIR / "trans_splicing" / "work"
RESULTS_DIR = ANALYSIS_DIR / "trans_splicing" / "results"

MAFFT = shutil.which("mafft") or "/software/team301/mafft-7.525-with-extensions/core/mafft"
MAFFT_BINARIES_DIR = os.environ.get("MAFFT_BINARIES", "/software/team301/mafft-7.525-with-extensions/core")

AA_ORDER = "ARNDCQEGHILKMFPSTWYV"
PSEUDOCOUNT = 0.5
GAP_COLUMN_THRESHOLD = 0.9
MIN_TRAINING_SEQS = 5


def run_mafft(protein_fasta: Path, out_fasta: Path) -> bool:
    env = ({"MAFFT_BINARIES": MAFFT_BINARIES_DIR, "PATH": "/usr/bin:/bin"}
           if Path(MAFFT_BINARIES_DIR).is_dir() else dict(os.environ))
    proc = subprocess.run(
        [MAFFT, "--auto", "--quiet", str(protein_fasta)],
        capture_output=True, text=True, timeout=600, env=env,
    )
    if proc.returncode != 0 or not proc.stdout.strip():
        print(f"[warn] mafft failed for {protein_fasta.name}: {proc.stderr.strip()[:300]}", file=sys.stderr)
        return False
    out_fasta.write_text(proc.stdout)
    return True


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    name = None
    chunks: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(chunks)
            name = line[1:].strip()
            chunks = []
        elif line.strip():
            chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def build_profile(alignment: dict[str, str]) -> list[dict[str, float]] | None:
    n = len(alignment)
    if n == 0:
        return None
    aln_len = len(next(iter(alignment.values())))
    bg_counts = Counter()
    for seq in alignment.values():
        bg_counts.update(c for c in seq if c in AA_ORDER)
    bg_total = sum(bg_counts.values()) or 1
    background = {aa: (bg_counts.get(aa, 0) + PSEUDOCOUNT) / (bg_total + PSEUDOCOUNT * 20) for aa in AA_ORDER}
    columns = []
    for i in range(aln_len):
        col_chars = [seq[i] for seq in alignment.values()]
        gap_frac = col_chars.count("-") / n
        if gap_frac > GAP_COLUMN_THRESHOLD:
            continue
        counts = Counter(c for c in col_chars if c in AA_ORDER)
        total = sum(counts.values()) or 1
        scores = {}
        for aa in AA_ORDER:
            freq = (counts.get(aa, 0) + PSEUDOCOUNT) / (total + PSEUDOCOUNT * 20)
            scores[aa] = math.log2(freq / background[aa])
        columns.append(scores)
    return columns or None


def write_profile(columns: list[dict[str, float]], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(AA_ORDER)]
    for col in columns:
        lines.append("\t".join(f"{col[aa]:.4f}" for aa in AA_ORDER))
    out_path.write_text("\n".join(lines) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--min-score", type=float, default=0.0,
                     help="whole_gene_score floor for inclusion in the refined training set")
    args = ap.parse_args()

    genes_path = RESULTS_DIR / "reconstructed_genes.tsv"
    if not genes_path.exists():
        print(f"[err] missing {genes_path} - run 01_reconstruct.py first", file=sys.stderr)
        sys.exit(1)
    df = pd.read_csv(genes_path, sep="\t")

    for gene, sub in df.groupby("gene"):
        good = sub[(sub["complete"]) & (sub["whole_gene_score"] > args.min_score)]
        if len(good) < MIN_TRAINING_SEQS:
            print(f"[warn] {gene}: only {len(good)} complete, score>{args.min_score} reconstructions "
                  f"(need >= {MIN_TRAINING_SEQS}) - keeping existing whole_gene.pssm", file=sys.stderr)
            continue

        protein_seqs = {row.species: row.protein.replace("-", "") for row in good.itertuples()}
        label = f"{gene}.whole_gene_refined"
        protein_fasta = WORK_DIR / f"{label}.protein.fasta"
        protein_fasta.write_text("".join(f">{sp}\n{seq}\n" for sp, seq in protein_seqs.items()))
        aligned_fasta = WORK_DIR / f"{label}.aligned.fasta"
        if not run_mafft(protein_fasta, aligned_fasta):
            continue
        alignment = read_fasta(aligned_fasta)
        columns = build_profile(alignment)
        if columns is None:
            print(f"[warn] {gene}: no usable alignment columns from refined training set - skipping", file=sys.stderr)
            continue
        out_path = PROFILES_DIR / gene / "whole_gene.pssm"
        write_profile(columns, out_path)
        print(f"[info] {gene}: whole_gene.pssm refined from {len(protein_seqs)} dataset-derived sequences "
              f"(was: 5-7 GenBank references) -> {out_path}", file=sys.stderr)

    print("[info] bootstrap_refine: done - re-run 01_reconstruct.py to use the refined profiles", file=sys.stderr)


if __name__ == "__main__":
    main()
