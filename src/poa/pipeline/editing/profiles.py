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

The profile maths (alignment, log-odds columns) is common/pssm.py's,
shared with every other profile poa builds. src/00_build_reference_profiles.py
is the command-line wrapper (main() here).
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

from poa.pipeline.common.errors import StageError
from poa.pipeline.common.pssm import align_and_build, read_fasta, translate, write_profile

MIN_SEQUENCES = 5  # below this, not enough signal for a meaningful profile
MIN_PROTEIN_LEN = 20  # aa
MAX_INTERNAL_STOPS = 0  # a stop anywhere but the very last codon = broken translation


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


def genes_for_organelle(analysis_dir: Path, organelle: str) -> list[str]:
    """Marker genes with an extracted FASTA, restricted to protein-coding
    ones - the marker-gene set (annotation/README.md) also includes tRNA
    (`trn*`) and rRNA (`rrn*`) genes, which aren't ORFs and would produce a
    meaningless "protein" profile if translated."""
    genes_dir = analysis_dir / "annotation" / "results" / "genes" / organelle
    if not genes_dir.exists():
        raise StageError(f"missing {genes_dir} - run annotation/src/02_gene_fasta_extract.py first")
    genes = sorted(p.stem for p in genes_dir.glob("*.fasta"))
    return [g for g in genes if not (g.startswith("trn") or g.startswith("rrn"))]


def build_dataset_profiles(root: Path, organelles: list[str] | tuple[str, ...] = ("mito", "pltd")) -> Path:
    """One profile per marker gene, from data root `root`'s own extracted
    gene FASTAs. Returns the profiles directory."""
    analysis_dir = Path(root) / "analysis"
    work_dir = analysis_dir / "editing" / "work"
    work_dir.mkdir(parents=True, exist_ok=True)
    profiles_dir = analysis_dir / "editing" / "results" / "profiles"

    n_built = 0
    for organelle in organelles:
        genes = genes_for_organelle(analysis_dir, organelle)
        for gene in genes:
            gene_fasta = analysis_dir / "annotation" / "results" / "genes" / organelle / f"{gene}.fasta"
            nuc_seqs = read_fasta(gene_fasta)
            proteins = clean_proteins(nuc_seqs)
            if len(proteins) < MIN_SEQUENCES:
                print(
                    f"[warn] {organelle}/{gene}: only {len(proteins)} clean sequences (need >= "
                    f"{MIN_SEQUENCES}) - skipping",
                    file=sys.stderr,
                )
                continue

            columns = align_and_build(proteins, work_dir, f"{organelle}.{gene}")
            if columns is None:
                print(f"[warn] {organelle}/{gene}: no usable alignment - skipping", file=sys.stderr)
                continue

            write_profile(columns, profiles_dir / organelle / f"{gene}.pssm")
            n_built += 1
            print(
                f"[info] {organelle}/{gene}: profile from {len(proteins)} sequences, "
                f"{len(columns)} columns -> profiles/{organelle}/{gene}.pssm",
                file=sys.stderr,
            )

    print(f"[info] build_reference_profiles: {n_built} profiles -> {profiles_dir}", file=sys.stderr)
    return profiles_dir


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--organelle", default="both", choices=["mito", "pltd", "both"])
    ap.add_argument(
        "--force",
        action="store_true",
        help="accepted for CLI consistency; this script "
        "always does a full recompute (cheap - a couple dozen short genes)",
    )
    args = ap.parse_args(argv)
    root = Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT") or os.getcwd())
    try:
        build_dataset_profiles(root, ["mito", "pltd"] if args.organelle == "both" else [args.organelle])
    except StageError as e:
        print(f"[err] {e}", file=sys.stderr)
        return 1
    return 0
