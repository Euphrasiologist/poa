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
# pseudo-observations per profile column, shared out in proportion to the
# background (build_profile)
PSEUDOCOUNT_MASS = PSEUDOCOUNT * 20
# Editing sites. A reference protein can say which of its residues RNA
# editing made, in lower case (mark_edited: an annotated edited protein
# against its gene's own translation; a protein derived from RNA-seq). A
# profile built from such proteins gets a 21st column, `edit`: the share of
# the groups (genera) of proteins present at that column that are edited
# there, 0 unless at least EDIT_SITE_MIN_GROUPS are. orfedit and transsplice
# charge --novel-edit-penalty, not --edit-penalty, for a C-to-U edit at a
# column with 0. Two genera, not one: public RNA-seq is often from another
# plant than the assembly, and where the two differ by a C/T the reads look
# edited; such sites fall anywhere, and seldom at the same column twice.
EDIT_SITE_MIN_GROUPS = 2
EDIT_COLUMN = "edit"
# build_profile's trim_ends: an end column goes if more than this fraction
# of sequences lack it
END_GAP_THRESHOLD = 0.5

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
        counts.update(c for c in seq.upper() if c in AA_ORDER)
    total = sum(counts.values())
    return {aa: (counts.get(aa, 0) + 1) / (total + 20) for aa in AA_ORDER}


def mark_edited(dna: str, protein: str) -> str:
    """`protein` with the residues RNA editing made in lower case: those
    that differ from the genome's own translation of its coding sequence
    `dna`. Unmarked when the two do not line up residue for residue."""
    dna = dna.upper()
    codons = [dna[i : i + 3] for i in range(0, len(dna) - len(dna) % 3, 3)]
    if codons and CODON_TABLE.get(codons[-1]) == "*":
        codons.pop()
    protein = protein.upper()
    if len(codons) != len(protein):
        return protein
    return "".join(p.lower() if p != CODON_TABLE.get(c, "X") and p in c_to_u_readings(c) else p for p, c in zip(protein, codons))


def c_to_u_readings(codon: str) -> set[str]:
    """What `codon` can be read as with one or more of its Cs edited to U."""
    cs = [i for i, b in enumerate(codon) if b == "C"]
    out = set()
    for mask in range(1, 1 << len(cs)):
        c = list(codon)
        for k, i in enumerate(cs):
            if mask >> k & 1:
                c[i] = "T"
        out.add(CODON_TABLE.get("".join(c), "X"))
    return out


def edited_residues(protein: str) -> set[int]:
    """Indices of the residues marked (lower case) as made by editing."""
    return {i for i, c in enumerate(protein) if c.islower()}


def build_profile(
    alignment: dict[str, str],
    gap_threshold: float = GAP_COLUMN_THRESHOLD,
    background: dict[str, float] | None = None,
    pseudocount_mass: float | None = None,
    trim_ends: tuple[bool, bool] = (False, False),
    edited: dict[str, set[int]] | None = None,
    groups: dict[str, str] | None = None,
) -> list[dict[str, float]] | None:
    """edited: per sequence name, the indices (ungapped) of the residues
    RNA editing made; each column then also gets EDIT_COLUMN, the share of
    the groups present there with such a residue (0 below
    EDIT_SITE_MIN_GROUPS). groups: sequence name -> group (genus), each
    sequence its own without it.

    trim_ends (start, end): at that end, also drop the run of columns
    most sequences lack (END_GAP_THRESHOLD). For an exon slot's intron-side
    ends: where one reference annotates the intron boundary a few codons
    out (mito ccmFc exon 2 starts 7 codons early in 2 of 10 genomes), its
    extra columns otherwise stay in the profile, and every other species'
    exon is aligned that far into its intron.

    background, pseudocount_mass: score against these amino acid
    frequencies, with `pseudocount_mass` pseudo-observations per column
    shared out in proportion to them. Without them the background is the
    alignment's own composition and the mass PSEUDOCOUNT_MASS.

    (Until 0.3.0 the default gave every residue the same PSEUDOCOUNT in
    every column. A residue the proteins hardly use - Cys or Trp in a small
    ribosomal protein - then had a background far below its pseudocount and
    scored as a match in every column (+5 to +9 bits), which invited edits
    to it: Arg to Cys in all three arginine codons of plastid rps19 in 44 of
    46 angiosperms, none of them in their RNA. Shared out in proportion, a
    residue a column does not have scores the same, below zero, whatever
    its background.)"""
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

    mass = PSEUDOCOUNT_MASS if pseudocount_mass is None else pseudocount_mass
    # per sequence, the residue index at each alignment column (None at a gap)
    residue_at: dict[str, list[int | None]] = {}
    if edited is not None:
        for name, seq in alignment.items():
            n_seen, row = 0, []
            for ch in seq:
                row.append(None if ch == "-" else n_seen)
                n_seen += ch != "-"
            residue_at[name] = row
    gap_fracs = [sum(seq[i] == "-" for seq in alignment.values()) / n for i in range(aln_len)]
    first, last = 0, aln_len
    if trim_ends[0]:
        while first < last and gap_fracs[first] > END_GAP_THRESHOLD:
            first += 1
    if trim_ends[1]:
        while last > first and gap_fracs[last - 1] > END_GAP_THRESHOLD:
            last -= 1

    columns = []
    for i in range(first, last):
        col_chars = [seq[i] for seq in alignment.values()]
        gap_frac = gap_fracs[i]
        if gap_frac > gap_threshold:
            continue
        counts = Counter(c for c in col_chars if c in AA_ORDER)
        total = sum(counts.values()) or 1
        scores = {}
        for aa in AA_ORDER:
            freq = (counts.get(aa, 0) + mass * background[aa]) / (total + mass)
            scores[aa] = math.log2(freq / background[aa])
        if edited is not None:
            present, hit = set(), set()
            for name, row in residue_at.items():
                if row[i] is None:
                    continue
                group = (groups or {}).get(name, name)
                present.add(group)
                if row[i] in edited.get(name, ()):
                    hit.add(group)
            scores[EDIT_COLUMN] = len(hit) / len(present) if present and len(hit) >= EDIT_SITE_MIN_GROUPS else 0.0
        columns.append(scores)
    return columns or None


def write_profile(columns: list[dict[str, float]], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    sites = any(EDIT_COLUMN in col for col in columns)
    lines = ["\t".join(AA_ORDER) + (f"\t{EDIT_COLUMN}" if sites else "")]
    for col in columns:
        lines.append("\t".join(f"{col[aa]:.4f}" for aa in AA_ORDER) + (f"\t{col.get(EDIT_COLUMN, 0.0):.3f}" if sites else ""))
    out_path.write_text("\n".join(lines) + "\n")


def align_and_build(
    proteins: dict[str, str],
    work_dir: Path,
    label: str,
    gap_threshold: float = GAP_COLUMN_THRESHOLD,
    background: dict[str, float] | None = None,
    pseudocount_mass: float | None = None,
    trim_ends: tuple[bool, bool] = (False, False),
    groups: dict[str, str] | None = None,
    sites: bool | None = None,
) -> list[dict[str, float]] | None:
    """Write `proteins` to work_dir/<label>.protein.fasta, align them to
    <label>.aligned.fasta and return the profile columns (None if mafft
    failed or no column survived). Columns gappier than gap_threshold are
    dropped. Proteins with residues in lower case (mark_edited) give the
    profile its editing sites, counted by `groups` (build_profile). sites:
    True for proteins that mark their edited residues even if none of
    these has any (the gene is not edited: the profile says so, with 0 in
    every column), False for none, None to go by the lower case."""
    edited = {name: edited_residues(seq) for name, seq in proteins.items()}
    if sites is False or (sites is None and not any(edited.values())):
        edited = None
    protein_fasta = work_dir / f"{label}.protein.fasta"
    protein_fasta.write_text("".join(f">{name}\n{seq.upper()}\n" for name, seq in proteins.items()))
    aligned_fasta = work_dir / f"{label}.aligned.fasta"
    if len(proteins) == 1:  # nothing to align
        aligned_fasta.write_text(protein_fasta.read_text())
    elif not run_mafft(protein_fasta, aligned_fasta):
        return None
    alignment = {name: seq.upper() for name, seq in read_fasta(aligned_fasta).items()}
    return build_profile(alignment, gap_threshold, background, pseudocount_mass, trim_ends, edited, groups)
