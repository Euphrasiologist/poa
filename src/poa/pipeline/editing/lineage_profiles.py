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
from poa.pipeline.trans_splicing import column_mask
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
# fern: every single-exon gene the references support. Never part of the
# automatic lineage call, so nothing is pinned.
# Not matR or nad9: on the 29 DToL fern assemblies their fern profiles did
# worse than the angiosperm ones. The references' nad9 runs about 55
# residues past the angiosperm protein's end, and 12 of 26 DToL fern
# copies aligned end to end with 3 frameshift calls, against 18 and none;
# the matR profile's first 150 columns are ones the dataset's ferns lack,
# and it drew 25 frameshift calls in 24 copies against 7.
# Nor atp1, cox1 and rps14, which have an intron in ferns and are
# reconstructed from exon templates (rg.LINEAGE_EXTRA_SPLICED).
FERN_GENES = [g for g in REFSEQ_GENES if g not in ("matR", "nad9", *rg.LINEAGE_EXTRA_SPLICED["fern"])]
# horsetail: the single-exon genes of Equisetum (cox1 has introns there).
# matR too, although the Equisetum reading frame is two thirds the length
# of the angiosperm protein and unedited: against the angiosperm profile
# it drew 24 edit calls a copy, none of them in the RNA.
HORSETAIL_GENES = [g for g in REFSEQ_GENES if g not in rg.LINEAGE_EXTRA_SPLICED["horsetail"]]
# lycophyte: every single-exon gene of clubmosses (rg.LINEAGE_EXTRA_SPLICED
# lists the ones with introns there), from RNA-edited proteins only.
LYCOPHYTE_GENES = [g for g in REFSEQ_GENES if g not in rg.LINEAGE_EXTRA_SPLICED["lycophyte"]]
LINEAGE_GENES = {
    lin.DEFAULT_LINEAGE: REFSEQ_GENES,
    "moss": MOSS_GENES,
    "fern": FERN_GENES,
    "horsetail": HORSETAIL_GENES,
    "lycophyte": LYCOPHYTE_GENES,
}
# fewest verified proteins a gene's profile is built from. Ferns have
# five species with fully edited proteins (rg.FERN_EDITED_MITOGENOMES).
# Horsetail is one genus (Equisetum), and lycophyte has one deeply
# sequenced transcriptome (Phylloglossum): see exon_profiles.LINEAGE_MIN_LOCI.
LINEAGE_MIN_SEQUENCES = {"fern": 4, "horsetail": 1, "lycophyte": 1}
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
LINEAGE_GAP_THRESHOLD = {
    lin.DEFAULT_LINEAGE: 0.5,
    "moss": pssm.GAP_COLUMN_THRESHOLD,
    "fern": 0.5,
    "horsetail": 0.5,
    "lycophyte": 0.5,
}
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
LINEAGE_SYNONYMS["fern"] = {**LINEAGE_SYNONYMS[lin.DEFAULT_LINEAGE], "mat-lepto": "matR"}
LINEAGE_SYNONYMS["horsetail"] = LINEAGE_SYNONYMS[lin.DEFAULT_LINEAGE]
LINEAGE_SYNONYMS["lycophyte"] = LINEAGE_SYNONYMS[lin.DEFAULT_LINEAGE]


def refseq_profiles_dir(code_dir: Path = PIPELINE_DIR) -> Path:
    """Bundled angiosperm RefSeq profiles (<dir>/mito/<gene>.pssm)."""
    return Path(code_dir) / "editing" / "reference" / "profiles_refseq"


def default_out_dir(lineage: str) -> Path:
    if lineage == lin.DEFAULT_LINEAGE:
        return refseq_profiles_dir() / "mito"
    return lin.editing_profiles(PIPELINE_DIR, lineage) / "mito"


def collect_proteins(
    mitogenomes: dict[str, str],
    genes: list[str],
    synonyms: dict[str, str],
    exclude: tuple[str, ...],
    reverse: bool = False,
) -> tuple[dict[str, dict[str, str]], dict[str, Counter]]:
    """Per gene, {accession: verified single-exon RefSeq protein}, and the
    gene's exon counts across the mitogenomes. reverse: the lineage edits U
    to C too, which verification then allows for."""
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
            ok, note = verify(exons[0], f.qualifiers["translation"][0], reverse)
            if not ok:
                log.warning(f"{acc} {organism} {gene}: {note} - not used")
                continue
            proteins[gene][acc] = f.qualifiers["translation"][0]
    return proteins, n_exons


# Proteins that are not GenBank records: single-exon mito genes of DToL
# assemblies with the editing their own RNA-seq shows applied
# (reference/lineages/<lineage>/README.md), headers `species|gene`.
RNA_PROTEINS_NAME = "rna_proteins.fasta"
# Pseudo-observations per profile column, shared out by the proteins'
# composition (pssm.build_profile), for lineages whose profiles are built
# that way; the others keep the flat pseudocount their calls were
# calibrated with.
# Horsetail and lycophyte: one to four proteins per gene, which 2
# pseudo-observations would outweigh (exon_profiles.LINEAGE_PSEUDOCOUNT).
LINEAGE_PSEUDOCOUNT: dict[str, float] = {"fern": 2.0, "horsetail": 0.25, "lycophyte": 0.25}


# the same for plastid genes, for a lineage with plastid profiles of its
# own (rg.LINEAGE_PLASTID)
RNA_PLASTID_PROTEINS_NAME = "rna_proteins.pltd.fasta"


def rna_proteins(lineage: str, exclude: tuple[str, ...] = (), organelle: str = "mito") -> dict[str, dict[str, str]]:
    """{gene: {"RNA_<species>": protein}} from the lineage's bundled
    RNA-edited proteins; `exclude` holds such RNA_<species> names."""
    name = RNA_PROTEINS_NAME if organelle == "mito" else RNA_PLASTID_PROTEINS_NAME
    path = lin.editing_profiles(PIPELINE_DIR, lineage) / name
    out: dict[str, dict[str, str]] = {}
    if path.exists():
        for header, seq in pssm.read_fasta(path).items():
            species, gene = header.split("|")[:2]
            if f"RNA_{species}" not in exclude:
                out.setdefault(gene, {})[f"RNA_{species}"] = seq
    return out


def one_per_species(proteins: dict[str, dict[str, str]], mitogenomes: dict[str, str]) -> dict[str, dict[str, str]]:
    """`proteins` with a species' identical copies of a gene counted once.
    A mitogenome in several records (Haplopteris is nine circles, with nad9
    on six of them) would otherwise outvote the other species: columns only
    it has would pass the gap threshold and its residues would dominate."""
    out: dict[str, dict[str, str]] = {}
    for gene, by_acc in proteins.items():
        seen: set[tuple[str, str]] = set()
        out[gene] = {}
        for acc, seq in by_acc.items():
            key = (mitogenomes.get(acc, acc), seq)
            if key not in seen:
                seen.add(key)
                out[gene][acc] = seq
    return out


def build_plastid_profiles(
    lineage: str, work_root: Path, *, exclude: tuple[str, ...] = (), out_dir: Path | None = None
) -> Path:
    """(Re)build `lineage`'s bundled plastid editing profiles
    (<lineage>/pltd/<gene>.pssm) from its RNA-edited plastid proteins, one
    per gene with enough of them, leaving out the RNA_<species> in
    `exclude`. Returns the profiles directory."""
    if lineage not in rg.LINEAGE_PLASTID:
        raise ValueError(f"lineage must be one of {sorted(rg.LINEAGE_PLASTID)}, not {lineage!r}")
    proteins = rna_proteins(lineage, tuple(exclude), "pltd")
    min_sequences = LINEAGE_MIN_SEQUENCES.get(lineage, profiles.MIN_SEQUENCES)
    work_dir = Path(work_root) / "analysis" / "editing" / "work" / f"{lineage}_pltd"
    work_dir.mkdir(parents=True, exist_ok=True)
    out_dir = Path(out_dir) if out_dir else lin.editing_profiles(PIPELINE_DIR, lineage) / "pltd"
    out_dir.mkdir(parents=True, exist_ok=True)
    for stale in out_dir.glob("*.pssm"):
        stale.unlink()
    for gene, by_name in sorted(proteins.items()):
        if len(by_name) < min_sequences:
            log.info(f"{gene}: no {lineage} plastid profile ({len(by_name)} proteins)")
            continue
        sharp = {}
        if lineage in LINEAGE_PSEUDOCOUNT:
            sharp = dict(background=pssm.composition(by_name.values()), pseudocount_mass=LINEAGE_PSEUDOCOUNT[lineage])
        columns = pssm.align_and_build(by_name, work_dir, gene, LINEAGE_GAP_THRESHOLD[lineage], **sharp)
        if columns is None:
            raise StageError(f"{gene}: no profile from {work_dir / gene}.protein.fasta (mafft failed or no columns)")
        pssm.write_profile(columns, out_dir / f"{gene}.pssm")
        log.info(f"{gene}: {lineage} plastid profile from {len(by_name)} proteins, {len(columns)} columns")
    return out_dir


def build_lineage_profiles(
    lineage: str,
    work_root: Path,
    *,
    exclude: tuple[str, ...] = (),
    out_dir: Path | None = None,
    use_mask: bool = True,
) -> Path:
    """(Re)build `lineage`'s bundled mito profiles in the installed package
    (or in out_dir), leaving out the mitogenomes in `exclude`; alignments go
    under work_root/analysis/editing/work/<lineage>. Returns the profiles
    directory."""
    if lineage not in LINEAGE_GENES:
        raise ValueError(f"lineage must be one of {sorted(LINEAGE_GENES)}, not {lineage!r}")
    mitogenomes = rg.LINEAGE_EDITED_MITOGENOMES.get(lineage, rg.LINEAGE_MITOGENOMES[lineage])
    known = set(mitogenomes)
    if lineage == lin.DEFAULT_LINEAGE:
        known |= set(SPLIT_CCMFN_MITOGENOMES)
    known |= {name for by_name in rna_proteins(lineage).values() for name in by_name}
    if unknown := set(exclude) - known:
        raise ValueError(f"--exclude: {sorted(unknown)} not among the {lineage} mitogenomes")
    genes = list(LINEAGE_GENES[lineage])
    proteins, n_exons = collect_proteins(
        mitogenomes,
        genes,
        LINEAGE_SYNONYMS[lineage],
        exclude,
        reverse=lin.stop_edit_penalty(lineage) is not None,
    )
    proteins = one_per_species(proteins, mitogenomes)
    for gene, by_name in rna_proteins(lineage, tuple(exclude)).items():
        if gene in proteins:
            proteins[gene].update(by_name)
    min_sequences = LINEAGE_MIN_SEQUENCES.get(lineage, profiles.MIN_SEQUENCES)
    # columns the lineage's own species don't have (trans_splicing/column_mask.py)
    mask = column_mask.load_mask(lineage, editing=True) if use_mask and lineage != lin.DEFAULT_LINEAGE else {}
    if lineage == lin.DEFAULT_LINEAGE:
        split, split_exons = collect_proteins(SPLIT_CCMFN_MITOGENOMES, SPLIT_CCMFN_GENES, SPLIT_CCMFN_SYNONYMS, exclude)
        genes += SPLIT_CCMFN_GENES
        proteins.update(split)
        n_exons.update(split_exons)

    work_dir = Path(work_root) / "analysis" / "editing" / "work" / lineage
    work_dir.mkdir(parents=True, exist_ok=True)
    out_dir = Path(out_dir) if out_dir else default_out_dir(lineage)
    for gene in genes:
        # a lineage with no GenBank mitogenome at all (horsetail) has only
        # its RNA-derived proteins, single-exon by construction
        modal = n_exons[gene].most_common(1)[0][0] if n_exons[gene] else (None if mitogenomes else 1)
        stale = out_dir / f"{gene}.pssm"
        if modal != 1 or len(proteins[gene]) < min_sequences:
            log.info(
                f"{gene}: no {lineage} profile (exon counts {dict(n_exons[gene])}, "
                f"{len(proteins[gene])} single-exon verified loci)"
            )
            stale.unlink(missing_ok=True)
            continue
        sharp = {}
        if lineage in LINEAGE_PSEUDOCOUNT:
            sharp = dict(
                background=pssm.composition(proteins[gene].values()), pseudocount_mass=LINEAGE_PSEUDOCOUNT[lineage]
            )
        columns = pssm.align_and_build(proteins[gene], work_dir, gene, LINEAGE_GAP_THRESHOLD[lineage], **sharp)
        if columns is None:
            raise StageError(f"{gene}: no profile from {work_dir / gene}.protein.fasta (mafft failed or no columns)")
        pssm.write_profile(columns, out_dir / f"{gene}.pssm")
        if mask.get(gene):
            n = column_mask.apply_mask(out_dir / f"{gene}.pssm", mask[gene])
            log.info(f"{gene}: dropped {n} masked column(s)")
        log.info(
            f"{gene}: {lineage} profile from {len(proteins[gene])} proteins, "
            f"{len(columns)} columns -> {out_dir / gene}.pssm"
        )

    return out_dir


def main(argv: list[str] | None = None) -> int:
    setup_logging()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lineage", required=True, choices=sorted(LINEAGE_GENES))
    ap.add_argument(
        "--organelle",
        default="mito",
        choices=["mito", "pltd"],
        help="pltd: the lineage's plastid profiles, from its RNA-edited plastid proteins (fern)",
    )
    ap.add_argument("--exclude", nargs="+", default=[], metavar="ACC", help="leave these RefSeq mitogenomes out")
    ap.add_argument("--out-dir", type=Path, help="write the profiles here instead of the bundled location")
    ap.add_argument(
        "--no-mask",
        action="store_true",
        help="keep every profile column, ignoring the lineage's profile_mask.tsv (do this before deriving a new mask)",
    )
    args = ap.parse_args(argv)
    work_root = Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT") or os.getcwd())
    if args.organelle == "pltd":
        build_plastid_profiles(args.lineage, work_root, exclude=tuple(args.exclude), out_dir=args.out_dir)
        return 0
    try:
        build_lineage_profiles(
            args.lineage,
            Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT") or os.getcwd()),
            exclude=tuple(args.exclude),
            out_dir=args.out_dir,
            use_mask=not args.no_mask,
        )
    except StageError as e:
        log.error(e)
        return 1
    return 0
