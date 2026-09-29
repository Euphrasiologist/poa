#!/usr/bin/env python3
"""Build one amino-acid profile (PSSM) per gene, for `orfedit` to align against.

Reuses this dataset's own already-annotated orthologs as the homology
signal - no external reference genome/database needed. Source:
`annotation/results/genes/{mito,pltd}/<gene>.fasta` (already strand-corrected,
one sequence per species - see `annotation/README.md`). Only genes that
already have such a FASTA are covered (i.e. genes single-copy in >=80% of
species where present - the phylogeny module's marker-gene set). Multi-exon/
trans-spliced genes (nad1/2/4/5/7, cox2, ccmFc/Fn, rps3) have no such FASTA
and are explicitly out of scope here - see `../README.md`.

Caveat worth remembering: if the *majority* of the dataset shares the same
boundary truncation (e.g. everyone's oatk call stops one codon short of the
true edited stop), this profile bakes that truncation in as "the consensus"
too, since it's built from those same calls. This dataset-internal approach
finds genes that are truncated in the *minority* of species relative to
the majority (a real, useful correction), not universal, lineage-wide
truncation. Worth spot-checking against literature-known cases before
trusting a gene family's corrected calls wholesale.

Output: analysis/editing/results/profiles/<gene>.pssm
        (TSV: header = 20 amino acid one-letter codes, one row per
        alignment column of log-odds scores)
"""
from __future__ import annotations

import os
import argparse
import shutil
import subprocess
import sys
from pathlib import Path

# CODE_DIR is the installed poa/pipeline/ (shared modules, bundled reference
# data); data/ and every module's results/ and work/ live under the data
# root - PLANT_ORGANELLE_DATA_ROOT, or the current directory if unset.
CODE_DIR = Path(__file__).resolve().parents[2]
ROOT_DIR = Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT") or os.getcwd())
ANALYSIS_DIR = ROOT_DIR / "analysis"
sys.path.insert(0, str(CODE_DIR / "common"))

import species_discovery as sd  # noqa: E402

MAFFT = shutil.which("mafft") or "/software/team301/mafft-7.525-with-extensions/core/mafft"
MAFFT_BINARIES_DIR = os.environ.get("MAFFT_BINARIES", "/software/team301/mafft-7.525-with-extensions/core")

AA_ORDER = "ARNDCQEGHILKMFPSTWYV"
MIN_SEQUENCES = 5  # below this, not enough signal for a meaningful profile
MIN_PROTEIN_LEN = 20  # aa
MAX_INTERNAL_STOPS = 0  # a stop anywhere but the very last codon = broken translation
GAP_COLUMN_THRESHOLD = 0.9  # drop alignment columns gappier than this
# 0.9, not the more obvious 0.5: verified empirically (see analysis/editing/
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


def best_frame_translation(nuc: str) -> str:
    """`.ctg.bed` intervals are HMM domain-hit boundaries, not guaranteed to
    start on a codon boundary (confirmed empirically: nad9 extractions came
    out 569/572/573nt across species - not a consistent multiple-of-3
    offset - and frame 0 translations were full of spurious stops that
    vanished in the correct frame). Try all 3 and keep whichever has fewest
    internal stops, exactly like `orfedit` itself does at scan time."""
    candidates = [translate(nuc[f:]) for f in range(3)]
    return min(candidates, key=lambda p: (p[:-1] if p.endswith("*") else p).count("*"))


def clean_proteins(nuc_seqs: dict[str, str]) -> dict[str, str]:
    """Translate (best of 3 frames) and drop obviously-broken entries (this
    is a filter, not a correction pass - orfedit is what fixes boundary
    truncation; this step just keeps the training set from being poisoned
    by clearly-wrong calls)."""
    out = {}
    for species, nuc in nuc_seqs.items():
        prot = best_frame_translation(nuc)
        body = prot[:-1] if prot.endswith("*") else prot
        if len(body) < MIN_PROTEIN_LEN:
            continue
        if body.count("*") > MAX_INTERNAL_STOPS:
            continue
        out[species] = body
    return out


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


def build_profile(alignment: dict[str, str]) -> list[dict[str, float]] | None:
    n = len(alignment)
    if n == 0:
        return None
    aln_len = len(next(iter(alignment.values())))

    # Background = overall amino acid frequency across all non-gap residues.
    from collections import Counter
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
            scores[aa] = __import__("math").log2(freq / background[aa])
        columns.append(scores)
    return columns or None


def write_profile(columns: list[dict[str, float]], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(AA_ORDER)]
    for col in columns:
        lines.append("\t".join(f"{col[aa]:.4f}" for aa in AA_ORDER))
    out_path.write_text("\n".join(lines) + "\n")


def genes_for_organelle(organelle: str) -> list[str]:
    """Marker genes with an extracted FASTA, restricted to protein-coding
    ones - the marker-gene set (annotation/README.md) also includes tRNA
    (`trn*`) and rRNA (`rrn*`) genes, which aren't ORFs and would produce a
    meaningless "protein" profile if translated."""
    genes_dir = ANALYSIS_DIR / "annotation" / "results" / "genes" / organelle
    if not genes_dir.exists():
        print(f"[err] missing {genes_dir} - run annotation/src/02_gene_fasta_extract.py first", file=sys.stderr)
        return []
    genes = sorted(p.stem for p in genes_dir.glob("*.fasta"))
    return [g for g in genes if not (g.startswith("trn") or g.startswith("rrn"))]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organelle", default="both", choices=["mito", "pltd", "both"])
    ap.add_argument("--force", action="store_true", help="accepted for CLI consistency; this script "
                     "always does a full recompute (cheap - a couple dozen short genes)")
    args = ap.parse_args()

    organelles = ["mito", "pltd"] if args.organelle == "both" else [args.organelle]
    work_dir = ANALYSIS_DIR / "editing" / "work"
    work_dir.mkdir(parents=True, exist_ok=True)
    profiles_dir = ANALYSIS_DIR / "editing" / "results" / "profiles"

    n_built = 0
    for organelle in organelles:
        genes = genes_for_organelle(organelle)
        for gene in genes:
            gene_fasta = ANALYSIS_DIR / "annotation" / "results" / "genes" / organelle / f"{gene}.fasta"
            nuc_seqs = read_fasta(gene_fasta)
            proteins = clean_proteins(nuc_seqs)
            if len(proteins) < MIN_SEQUENCES:
                print(f"[warn] {organelle}/{gene}: only {len(proteins)} clean sequences (need >= "
                      f"{MIN_SEQUENCES}) - skipping", file=sys.stderr)
                continue

            protein_fasta = work_dir / f"{organelle}.{gene}.protein.fasta"
            protein_fasta.write_text("".join(f">{sp}\n{seq}\n" for sp, seq in proteins.items()))
            aligned_fasta = work_dir / f"{organelle}.{gene}.aligned.fasta"
            if not run_mafft(protein_fasta, aligned_fasta):
                continue

            alignment = read_fasta(aligned_fasta)
            columns = build_profile(alignment)
            if columns is None:
                print(f"[warn] {organelle}/{gene}: no usable alignment columns - skipping", file=sys.stderr)
                continue

            write_profile(columns, profiles_dir / organelle / f"{gene}.pssm")
            n_built += 1
            print(f"[info] {organelle}/{gene}: profile from {len(proteins)} sequences, "
                  f"{len(columns)} columns -> profiles/{organelle}/{gene}.pssm", file=sys.stderr)

    print(f"[info] build_reference_profiles: {n_built} profiles -> {profiles_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
