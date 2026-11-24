"""Stage 2 of the reference-profile bootstrap (see exon_profiles.py
and `../README.md`): rebuild each gene's whole-gene profile from THIS
dataset's own high-confidence reconstructions instead of the 5-7 GenBank
genomes alone.

Why this matters concretely, not just in principle: `rps3`'s first-pass
whole-gene profile (5 GenBank references, one fast-evolving/size-variable
gene per the literature review) gave strongly *negative* whole-gene scores
even for reconstructions whose own per-exon evidence looked solid (both
exons independently high-scoring, correct assignment, no internal stops) -
a real quality issue in a small reference-only profile, not (necessarily)
a wrong reconstruction. Once reconstruct.py has run once, this
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
this tool's own output). Re-run reconstruct afterwards to benefit.

src/02_bootstrap_refine.py is the command-line wrapper (main() here).
"""

from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path

import pandas as pd

from poa import PIPELINE_DIR
from poa.pipeline.common import lineage as lin
from poa.pipeline.common.errors import StageError
from poa.pipeline.common.pssm import align_and_build, write_profile
from poa.pipeline.common.log import setup_logging

log = logging.getLogger(__name__)

MIN_TRAINING_SEQS = 5


def refine_whole_gene_profiles(root: Path, min_score: float = 0.0) -> list[str]:
    """Rebuild the package's angiosperm whole_gene.pssm for each gene with
    enough complete reconstructions (whole_gene_score > min_score) under
    data root `root`. Returns the genes refined."""
    analysis_dir = Path(root) / "analysis"
    work_dir = analysis_dir / "trans_splicing" / "work"
    profiles_dir = lin.trans_splicing_reference(PIPELINE_DIR, lin.DEFAULT_LINEAGE) / "profiles"
    genes_path = analysis_dir / "trans_splicing" / "results" / "reconstructed_genes.tsv"
    if not genes_path.exists():
        raise StageError(f"missing {genes_path} - run trans_splicing/src/01_reconstruct.py first")
    df = pd.read_csv(genes_path, sep="\t")
    # these are the angiosperm templates: a moss reconstruction (built on the
    # moss ones, common/lineage.py) must not leak into them
    lineage_of = lin.load(analysis_dir) or {}
    df = df[df["species"].map(lineage_of).fillna(lin.DEFAULT_LINEAGE) == lin.DEFAULT_LINEAGE]

    refined = []
    for gene, sub in df.groupby("gene"):
        # (a protein with a stop inside it, `*`, is no reference)
        good = sub[(sub["complete"]) & (sub["whole_gene_score"] > min_score) & ~sub["protein"].str.contains("*", regex=False)]
        if len(good) < MIN_TRAINING_SEQS:
            log.warning(
                f"{gene}: only {len(good)} complete, score>{min_score} reconstructions "
                f"(need >= {MIN_TRAINING_SEQS}) - keeping existing whole_gene.pssm"
            )
            continue

        protein_seqs = {row.species: row.protein.replace("-", "") for row in good.itertuples()}
        columns = align_and_build(protein_seqs, work_dir, f"{gene}.whole_gene_refined")
        if columns is None:
            log.warning(f"{gene}: no usable alignment from refined training set - skipping")
            continue
        out_path = profiles_dir / gene / "whole_gene.pssm"
        write_profile(columns, out_path)
        refined.append(gene)
        log.info(
            f"{gene}: whole_gene.pssm refined from {len(protein_seqs)} dataset-derived sequences "
            f"(was: 5-7 GenBank references) -> {out_path}"
        )

    log.info("bootstrap_refine: done - re-run 01_reconstruct.py to use the refined profiles")
    return refined


def main(argv: list[str] | None = None) -> int:
    setup_logging()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "--min-score", type=float, default=0.0, help="whole_gene_score floor for inclusion in the refined training set"
    )
    args = ap.parse_args(argv)
    try:
        refine_whole_gene_profiles(Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT") or os.getcwd()), args.min_score)
    except StageError as e:
        log.error(e)
        return 1
    return 0
