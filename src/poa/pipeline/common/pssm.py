"""Amino-acid profiles (PSSMs) for orfedit and transsplice: translate,
align with mafft, turn the alignment into per-column log-odds, write it.

One copy of this maths for every profile poa builds - editing's
dataset-built and lineage profiles, trans_splicing's exon-slot and
whole-gene templates - so all of them are made the same way.
"""

from __future__ import annotations

import logging
import math
import os
import shutil
import subprocess
from collections import Counter
from pathlib import Path

log = logging.getLogger(__name__)

MAFFT = shutil.which("mafft") or "/software/team301/mafft-7.525-with-extensions/core/mafft"
MAFFT_BINARIES_DIR = os.environ.get("MAFFT_BINARIES", "/software/team301/mafft-7.525-with-extensions/core")

AA_ORDER = "ARNDCQEGHILKMFPSTWYV"
GAP_COLUMN_THRESHOLD = 0.9  # drop alignment columns gappier than this
# 0.9, not the more obvious 0.5: verified empirically (see editing/
# README.md) that a 0.5 cutoff systematically truncated almost every profile
# by a few C-/N-terminal columns and made ~98% of zero-edit calls look
# "changed" by ~1-2 codons - not real editing, just the training alignment
# itself containing plenty of boundary-truncated raw calls (the exact
# problem this tool exists to fix), so a naive per-column gap-majority
# filter throws away genuine terminal signal, not just rare lineage-specific
# insertions. 0.9 keeps columns unless they're gap in nearly everyone.
PSEUDOCOUNT = 0.5

CODON_TABLE = {}
_bases = "TCAG"
_aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, b1 in enumerate(_bases):
    for j, b2 in enumerate(_bases):
        for k, b3 in enumerate(_bases):
            CODON_TABLE[b1 + b2 + b3] = _aas[i * 16 + j * 4 + k]


def translate(seq: str) -> str:
    seq = seq.upper()
    codons = [seq[i : i + 3] for i in range(0, len(seq) - len(seq) % 3, 3)]
    return "".join(CODON_TABLE.get(c, "X") for c in codons)


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


def run_mafft(protein_fasta: Path, out_fasta: Path) -> bool:
    env = (
        {"MAFFT_BINARIES": MAFFT_BINARIES_DIR, "PATH": "/usr/bin:/bin"}
        if Path(MAFFT_BINARIES_DIR).is_dir()
        else dict(os.environ)
    )
    if os.environ.get("TMPDIR"):  # mafft's mktemp otherwise needs a writable /tmp
        env["TMPDIR"] = os.environ["TMPDIR"]
    proc = subprocess.run(
        [MAFFT, "--auto", "--quiet", str(protein_fasta)],
        capture_output=True,
        text=True,
        timeout=600,
        env=env,
    )
    if proc.returncode != 0 or not proc.stdout.strip():
        log.warning(f"mafft failed for {protein_fasta.name}: {proc.stderr.strip()[:300]}")
        return False
    out_fasta.write_text(proc.stdout)
    return True


def composition(seqs) -> dict[str, float]:
    """Amino acid frequencies of `seqs`, for build_profile's `background`."""
    counts = Counter()
    for seq in seqs:
        counts.update(c for c in seq if c in AA_ORDER)
    total = sum(counts.values())
    return {aa: (counts.get(aa, 0) + 1) / (total + 20) for aa in AA_ORDER}


def build_profile(
    alignment: dict[str, str],
    gap_threshold: float = GAP_COLUMN_THRESHOLD,
    background: dict[str, float] | None = None,
    pseudocount_mass: float | None = None,
) -> list[dict[str, float]] | None:
    """background, pseudocount_mass: score against these amino acid
    frequencies, with `pseudocount_mass` pseudo-observations per column
    shared out in proportion to them. Without them the background is the
    alignment's own composition and every residue gets PSEUDOCOUNT, which
    only works when the alignment is long enough to hold every residue at
    something like its real frequency: in a short one (a 7-codon exon) a
    residue seen nowhere has a background far below its pseudocount and
    scores as a match in every column."""
    n = len(alignment)
    if n == 0:
        return None
    aln_len = len(next(iter(alignment.values())))

    if background is None:
        # Background = overall amino acid frequency across all non-gap residues.
        bg_counts = Counter()
        for seq in alignment.values():
            bg_counts.update(c for c in seq if c in AA_ORDER)
        bg_total = sum(bg_counts.values()) or 1
        background = {aa: (bg_counts.get(aa, 0) + PSEUDOCOUNT) / (bg_total + PSEUDOCOUNT * 20) for aa in AA_ORDER}

    columns = []
    for i in range(aln_len):
        col_chars = [seq[i] for seq in alignment.values()]
        gap_frac = col_chars.count("-") / n
        if gap_frac > gap_threshold:
            continue
        counts = Counter(c for c in col_chars if c in AA_ORDER)
        total = sum(counts.values()) or 1
        scores = {}
        for aa in AA_ORDER:
            if pseudocount_mass is None:
                freq = (counts.get(aa, 0) + PSEUDOCOUNT) / (total + PSEUDOCOUNT * 20)
            else:
                freq = (counts.get(aa, 0) + pseudocount_mass * background[aa]) / (total + pseudocount_mass)
            scores[aa] = math.log2(freq / background[aa])
        columns.append(scores)
    return columns or None


def write_profile(columns: list[dict[str, float]], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(AA_ORDER)]
    for col in columns:
        lines.append("\t".join(f"{col[aa]:.4f}" for aa in AA_ORDER))
    out_path.write_text("\n".join(lines) + "\n")


def align_and_build(
    proteins: dict[str, str],
    work_dir: Path,
    label: str,
    gap_threshold: float = GAP_COLUMN_THRESHOLD,
    background: dict[str, float] | None = None,
    pseudocount_mass: float | None = None,
) -> list[dict[str, float]] | None:
    """Write `proteins` to work_dir/<label>.protein.fasta, align them to
    <label>.aligned.fasta and return the profile columns (None if mafft
    failed or no column survived). Columns gappier than gap_threshold are
    dropped."""
    protein_fasta = work_dir / f"{label}.protein.fasta"
    protein_fasta.write_text("".join(f">{name}\n{seq}\n" for name, seq in proteins.items()))
    aligned_fasta = work_dir / f"{label}.aligned.fasta"
    if not run_mafft(protein_fasta, aligned_fasta):
        return None
    return build_profile(read_fasta(aligned_fasta), gap_threshold, background, pseudocount_mass)
