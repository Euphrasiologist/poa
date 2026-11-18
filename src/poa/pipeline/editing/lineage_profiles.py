"""Build bundled mito editing profiles from RefSeq mitogenomes rather than
from a dataset: per non-angiosperm lineage (common/lineage.py), and for
angiosperm the genes the dataset-built profiles don't cover.

profiles.py (src/00_build_reference_profiles.py) builds profiles from a dataset's own gene
calls, so on a mostly-angiosperm dataset they are angiosperm profiles. A
moss aligned against them gets "corrected" towards the angiosperm protein:
98% of moss mito calls changed, against 37% for dicots, although mosses
edit very little (Physcomitrella has ~11 mito sites). These profiles are
built from RefSeq's curated proteins instead, for the genes the bundled
mito profiles cover, so a moss is compared against mosses.

profiles.py also only covers genes single-copy in >=80% of species, which
in angiosperm mitochondria leaves 7 (atp4 ccmB ccmC mttB nad5 nad9
rps14): atp6, atp9, cox1, rps12 and the rest are too often duplicated or
missing. `--lineage angiosperm` builds RefSeq profiles for every
single-exon mito protein gene (REFSEQ_GENES) into
editing/reference/profiles_refseq/mito/; scan.py uses one for a gene only
when the dataset-built (or bundled dataset) profiles have none, so the 7
are unchanged. RefSeq's /translation is the edited protein, which is what
orfedit should align towards.

Per gene: every single-exon CDS in the lineage's mitogenomes
(trans_splicing/reference_genes.LINEAGE_MITOGENOMES; GenBank cached in
trans_splicing/reference/genbank/) whose genomic translation verifies
against RefSeq's protein (common/refseq_mito.verify), at least
profiles.MIN_SEQUENCES of them. A gene that is multi-exon in most of the
lineage (moss nad5, nad9) gets no profile: a single-locus correction is
the wrong model for it, and scan.py leaves it uncorrected for that
lineage rather than borrow the angiosperm one.

Split ccmFn: `--lineage angiosperm` also builds ccmFn1 and ccmFn2
profiles, one per piece of the gene that Brassicaceae and Allium carry as
two (SPLIT_CCMFN_MITOGENOMES; see denovo_annotation/gene_calls.py).

`--exclude ACC` leaves mitogenomes out and `--out-dir` writes elsewhere,
for held-out validation (the Arabidopsis smoke check compares against
NC_037304.1, which is one of the angiosperm references).

Output: editing/reference/lineages/<lineage>/mito/<gene>.pssm (moss)
        editing/reference/profiles_refseq/mito/<gene>.pssm (angiosperm)

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
# Moss profiles cover the genes the bundled dataset profiles covered when the
# lineage call was calibrated (common/lineage.py: 5 shared genes, nad5/nad9
# multi-exon in mosses). Pinned: adding genes here changes lineage calls.
MOSS_GENES = ["atp4", "ccmB", "ccmC", "mttB", "nad5", "nad9", "rps14"]
# every oatkDB mito protein family that is single-exon in angiosperms
# (trans_splicing handles the spliced ones); those with fewer than
# profiles.MIN_SEQUENCES verified RefSeq proteins (rpl10 rps2 rps11 sdh3
# sdh4, often lost to the nucleus) are skipped at build time
REFSEQ_GENES = [
    "atp1", "atp4", "atp6", "atp8", "atp9", "ccmB", "ccmC", "ccmFn", "cob", "cox1", "cox3",
    "matR", "mttB", "nad3", "nad4L", "nad6", "nad9", "rpl5", "rpl10", "rpl16", "rps1", "rps2",
    "rps4", "rps7", "rps11", "rps12", "rps13", "rps14", "rps19", "sdh3", "sdh4",
]  # fmt: skip
LINEAGE_GENES = {lin.DEFAULT_LINEAGE: REFSEQ_GENES, "moss": MOSS_GENES}
# ccmFn is two genes in Brassicaceae and Allium (denovo_annotation's
# gene_calls renames the pieces ccmFn1/ccmFn2), so each piece gets its own
# profile, from these RefSeq mitogenomes: a spread of Brassicaceae, not
# five near-identical Brassica, plus Allium cepa. Boechera and Capsella
# call the first piece plain "ccmFn" (380 aa).
SPLIT_CCMFN_GENES = ["ccmFn1", "ccmFn2"]
SPLIT_CCMFN_MITOGENOMES = {
    "NC_037304.1": "Arabidopsis thaliana",
    "NC_081483.1": "Arabidopsis lyrata",
    "NC_042143.1": "Boechera stricta",
    "NC_042883.1": "Capsella rubella",
    "NC_037070.1": "Arabis alpina",
    "NC_085821.1": "Descurainia sophia",
    "NC_088489.1": "Lepidium apetalum",
    "NC_008285.1": "Brassica napus",
    "NC_018551.1": "Raphanus sativus",
    "NC_031896.1": "Sinapis arvensis",
    "NC_030100.1": "Allium cepa",
}
SPLIT_CCMFN_SYNONYMS = {"ccmfn": "ccmFn1", "ccmfn1": "ccmFn1", "ccmfn2": "ccmFn2"}
# Drop alignment columns gappier than this. pssm's default (0.9) is for
# thousands of dataset calls, many boundary-truncated; with ~10 curated
# RefSeq proteins it keeps a column only one species fills (1/9 = 89% gap),
# so one species' N-terminal extension or read-through becomes part of the
# profile and orfedit over-extends towards it (held-out Arabidopsis: atp1
# 1845 bp vs RefSeq 1524 at 0.9). Moss keeps the default its lineage call
# was calibrated with.
LINEAGE_GAP_THRESHOLD = {lin.DEFAULT_LINEAGE: 0.5, "moss": pssm.GAP_COLUMN_THRESHOLD}
# GenBank gene-name spellings -> ours (compared lower-cased). Angiosperm
# records also use the pre-function names; moss keeps the set its bundled
# profiles were built with.
SYNONYMS = {**rg.GENE_SYNONYMS, "tatc": "mttB"}
LINEAGE_SYNONYMS = {
    lin.DEFAULT_LINEAGE: {
        **SYNONYMS,
        "mat-r": "matR",
        "orf25": "atp4",
        "orfb": "atp8",
        "orfx": "mttB",
        "ymf16": "mttB",
        "ccmfn1": "ccmFn",
        "yejr": "ccmFn",
    },
    "moss": SYNONYMS,
}


def refseq_profiles_dir(code_dir: Path = PIPELINE_DIR) -> Path:
    """Bundled angiosperm RefSeq profiles (<dir>/mito/<gene>.pssm)."""
    return Path(code_dir) / "editing" / "reference" / "profiles_refseq"


def default_out_dir(lineage: str) -> Path:
    if lineage == lin.DEFAULT_LINEAGE:
        return refseq_profiles_dir() / "mito"
    return lin.editing_profiles(PIPELINE_DIR, lineage) / "mito"


def collect_proteins(
    mitogenomes: dict[str, str], genes: list[str], synonyms: dict[str, str], exclude: tuple[str, ...]
) -> tuple[dict[str, dict[str, str]], dict[str, Counter]]:
    """Per gene, {accession: verified single-exon RefSeq protein}, and the
    gene's exon counts across the mitogenomes."""
    wanted = {g.lower(): g for g in genes}
    proteins: dict[str, dict[str, str]] = {g: {} for g in genes}
    n_exons: dict[str, Counter] = {g: Counter() for g in genes}
    for acc, organism in mitogenomes.items():
        if acc in exclude:
            continue
        rec = SeqIO.read(fetch(acc, GB_DIR), "genbank")
        for f in rec.features:
            if f.type != "CDS" or "gene" not in f.qualifiers or "translation" not in f.qualifiers:
                continue
            name = f.qualifiers["gene"][0].lower()
            gene = wanted.get(synonyms.get(name, "").lower()) or wanted.get(name)
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
    return proteins, n_exons


def build_lineage_profiles(
    lineage: str, work_root: Path, *, exclude: tuple[str, ...] = (), out_dir: Path | None = None
) -> Path:
    """(Re)build `lineage`'s bundled mito profiles in the installed package
    (or in out_dir), leaving out the mitogenomes in `exclude`; alignments go
    under work_root/analysis/editing/work/<lineage>. Returns the profiles
    directory."""
    if lineage not in LINEAGE_GENES:
        raise ValueError(f"lineage must be one of {sorted(LINEAGE_GENES)}, not {lineage!r}")
    known = set(rg.LINEAGE_MITOGENOMES[lineage])
    if lineage == lin.DEFAULT_LINEAGE:
        known |= set(SPLIT_CCMFN_MITOGENOMES)
    if unknown := set(exclude) - known:
        raise ValueError(f"--exclude: {sorted(unknown)} not among the {lineage} mitogenomes")
    genes = list(LINEAGE_GENES[lineage])
    proteins, n_exons = collect_proteins(rg.LINEAGE_MITOGENOMES[lineage], genes, LINEAGE_SYNONYMS[lineage], exclude)
    if lineage == lin.DEFAULT_LINEAGE:
        split, split_exons = collect_proteins(SPLIT_CCMFN_MITOGENOMES, SPLIT_CCMFN_GENES, SPLIT_CCMFN_SYNONYMS, exclude)
        genes += SPLIT_CCMFN_GENES
        proteins.update(split)
        n_exons.update(split_exons)

    work_dir = Path(work_root) / "analysis" / "editing" / "work" / lineage
    work_dir.mkdir(parents=True, exist_ok=True)
    out_dir = Path(out_dir) if out_dir else default_out_dir(lineage)
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
        columns = pssm.align_and_build(proteins[gene], work_dir, gene, LINEAGE_GAP_THRESHOLD[lineage])
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
    ap.add_argument("--lineage", required=True, choices=sorted(LINEAGE_GENES))
    ap.add_argument("--exclude", nargs="+", default=[], metavar="ACC", help="leave these RefSeq mitogenomes out")
    ap.add_argument("--out-dir", type=Path, help="write the profiles here instead of the bundled location")
    args = ap.parse_args(argv)
    try:
        build_lineage_profiles(
            args.lineage,
            Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT") or os.getcwd()),
            exclude=tuple(args.exclude),
            out_dir=args.out_dir,
        )
    except StageError as e:
        log.error(e)
        return 1
    return 0
