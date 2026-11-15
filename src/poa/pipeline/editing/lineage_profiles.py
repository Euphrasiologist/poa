"""Build the bundled editing profiles for a non-angiosperm lineage
(common/lineage.py), from that lineage's RefSeq mitogenomes rather than
from a dataset.

profiles.py (src/00_build_reference_profiles.py) builds profiles from a dataset's own gene
calls, so on a mostly-angiosperm dataset they are angiosperm profiles. A
moss aligned against them gets "corrected" towards the angiosperm protein:
98% of moss mito calls changed, against 37% for dicots, although mosses
edit very little (Physcomitrella has ~11 mito sites). These profiles are
built from RefSeq's curated proteins instead, for the genes the bundled
mito profiles cover, so a moss is compared against mosses.

Per gene: every single-exon CDS in the lineage's mitogenomes
(trans_splicing/reference_genes.LINEAGE_MITOGENOMES; GenBank cached in
trans_splicing/reference/genbank/) whose genomic translation verifies
against RefSeq's protein (common/refseq_mito.verify). A gene that is
multi-exon in most of the lineage (moss nad5, nad9) gets no profile: a
single-locus correction is the wrong model for it, and scan.py
leaves it uncorrected for that lineage rather than borrow the angiosperm
one.

Output: editing/reference/lineages/<lineage>/mito/<gene>.pssm (bundled)

src/00b_build_lineage_profiles.py is the command-line wrapper (main() here).
"""

from __future__ import annotations

import argparse
import logging
import os
from collections import Counter
from pathlib import Path

from Bio import SeqIO

from poa import PIPELINE_DIR
from poa.pipeline.common import lineage as lin
from poa.pipeline.common.errors import StageError
from poa.pipeline.common.refseq_mito import cds_exons, fetch, verify
from poa.pipeline.common import pssm
from poa.pipeline.editing import profiles
from poa.pipeline.trans_splicing import reference_genes as rg
from poa.pipeline.common.log import setup_logging

log = logging.getLogger(__name__)

GB_DIR = lin.trans_splicing_reference(PIPELINE_DIR, lin.DEFAULT_LINEAGE) / "genbank"
BUNDLED_MITO = PIPELINE_DIR / "editing" / "reference" / "profiles" / "mito"
NON_DEFAULT = [x for x in lin.LINEAGES if x != lin.DEFAULT_LINEAGE]


def build_lineage_profiles(lineage: str, work_root: Path) -> Path:
    """(Re)build `lineage`'s bundled mito profiles in the installed package;
    alignments go under work_root/analysis/editing/work/<lineage>. Returns
    the profiles directory."""
    if lineage not in NON_DEFAULT:
        raise ValueError(f"lineage must be one of {NON_DEFAULT}, not {lineage!r}")
    synonyms = {**rg.GENE_SYNONYMS, "tatc": "mttB", "mttb": "mttB"}
    genes = sorted(p.stem for p in BUNDLED_MITO.glob("*.pssm"))
    wanted = {g.lower(): g for g in genes}
    proteins: dict[str, dict[str, str]] = {g: {} for g in genes}
    n_exons: dict[str, Counter] = {g: Counter() for g in genes}
    for acc, organism in rg.LINEAGE_MITOGENOMES[lineage].items():
        rec = SeqIO.read(fetch(acc, GB_DIR), "genbank")
        for f in rec.features:
            if f.type != "CDS" or "gene" not in f.qualifiers or "translation" not in f.qualifiers:
                continue
            name = f.qualifiers["gene"][0].lower()
            gene = wanted.get(name) or wanted.get(synonyms.get(name, "").lower())
            if gene is None or acc in proteins[gene]:
                continue
            parts, exons = cds_exons(rec, f)
            n_exons[gene][len(parts)] += 1
            if len(parts) != 1:
                continue
            ok, note = verify(exons[0], f.qualifiers["translation"][0])
            if not ok:
                log.warning(f"{acc} {organism} {gene}: {note} - not used")
                continue
            proteins[gene][acc] = f.qualifiers["translation"][0]

    work_dir = Path(work_root) / "analysis" / "editing" / "work" / lineage
    work_dir.mkdir(parents=True, exist_ok=True)
    out_dir = lin.editing_profiles(PIPELINE_DIR, lineage) / "mito"
    for gene in genes:
        modal = n_exons[gene].most_common(1)[0][0] if n_exons[gene] else None
        stale = out_dir / f"{gene}.pssm"
        if modal != 1 or len(proteins[gene]) < profiles.MIN_SEQUENCES:
            log.info(
                f"{gene}: no {lineage} profile (exon counts {dict(n_exons[gene])}, "
                f"{len(proteins[gene])} single-exon verified loci)"
            )
            stale.unlink(missing_ok=True)
            continue
        columns = pssm.align_and_build(proteins[gene], work_dir, gene)
        if columns is None:
            raise StageError(f"{gene}: no profile from {work_dir / gene}.protein.fasta (mafft failed or no columns)")
        pssm.write_profile(columns, out_dir / f"{gene}.pssm")
        log.info(
            f"{gene}: {lineage} profile from {len(proteins[gene])} RefSeq proteins, "
            f"{len(columns)} columns -> {out_dir / gene}.pssm"
        )

    return out_dir


def main(argv: list[str] | None = None) -> int:
    setup_logging()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lineage", required=True, choices=NON_DEFAULT)
    args = ap.parse_args(argv)
    try:
        build_lineage_profiles(args.lineage, Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT") or os.getcwd()))
    except StageError as e:
        log.error(e)
        return 1
    return 0
