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
from poa.pipeline.common.refseq_mito import cds_exons, fetch, n_edited, unedited, verify
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
# liverwort: every single-exon gene of liverworts (rg.LINEAGE_EXTRA_SPLICED
# lists the ones with introns there; nad1, rps3, rps10 and ccmFc, single-exon
# in liverworts, are one-slot exon templates). Not matR, which liverworts
# lack. rpl6 and rps8, which they keep, have no gene model to call them.
LIVERWORT_GENES = [g for g in REFSEQ_GENES if g not in ("matR", *rg.LINEAGE_EXTRA_SPLICED["liverwort"])]
LINEAGE_GENES = {
    lin.DEFAULT_LINEAGE: REFSEQ_GENES,
    "moss": MOSS_GENES,
    "liverwort": LIVERWORT_GENES,
    "fern": FERN_GENES,
    "horsetail": HORSETAIL_GENES,
    "lycophyte": LYCOPHYTE_GENES,
}
# fewest verified proteins a gene's profile is built from. Ferns have
# five species with fully edited proteins (rg.FERN_EDITED_MITOGENOMES).
# Horsetail is one genus (Equisetum), and lycophyte has one deeply
# sequenced transcriptome (Phylloglossum): see exon_profiles.LINEAGE_MIN_LOCI.
# Angiosperm: 14 RefSeq mitogenomes are annotated with edited proteins, and
# rps19 and sdh3 are single-exon and edited in four of them.
LINEAGE_MIN_SEQUENCES = {"angiosperm": 4, "fern": 4, "horsetail": 1, "lycophyte": 1}
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
    "liverwort": 0.5,
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
LINEAGE_SYNONYMS["liverwort"] = LINEAGE_SYNONYMS[lin.DEFAULT_LINEAGE]
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
    edited_only: bool = False,
) -> tuple[dict[str, dict[str, str]], dict[str, Counter]]:
    """Per gene, {accession: verified single-exon RefSeq protein}, and the
    gene's exon counts across the mitogenomes. reverse: the lineage edits U
    to C too, which verification then allows for. edited_only: leave out a
    protein annotated unedited where the other records' are edited
    (refseq_mito.unedited)."""
    wanted = {g.lower(): g for g in genes}
    proteins: dict[str, dict[str, str]] = {g: {} for g in genes}
    n_exons: dict[str, Counter] = {g: Counter() for g in genes}
    edits: dict[str, dict[str, int]] = {g: {} for g in genes}
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
            edits[gene][acc] = n_edited(exons[0], f.qualifiers["translation"][0])
            if edited_only:  # annotated edited proteins (lin.SITE_LINEAGES): which residues editing made
                proteins[gene][acc] = pssm.mark_edited(exons[0], proteins[gene][acc])
    if edited_only:
        for gene in genes:
            for acc in sorted(unedited(edits[gene])):
                log.info(f"{acc} {gene}: annotated unedited ({edits[gene][acc]} edits; others {sorted(edits[gene].values())}) - not used")
                del proteins[gene][acc]
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


# Plastid proteins from GenBank, for a lineage whose plastid references are
# RefSeq plastomes (rg.LINEAGE_PLASTOMES: moss, liverwort): the single-exon
# genes the bundled land-plant plastid profiles cover, as annotated,
# headers `accession|gene`. Written by fetch_plastid_proteins and shipped,
# so building the profiles needs no GenBank record.
REFSEQ_PLASTID_PROTEINS_NAME = "refseq_proteins.pltd.fasta"
PLASTID_SYNONYMS = {"pafii": "ycf4", "pbf1": "psbN", "lhba": "psbZ", "ycf9": "psbZ", "ycf10": "cemA"}


# Plastid genes with a first exon of 6 to 9 bases. The gene call and the
# default profiles cover the second exon, so a RefSeq protein gives that
# part: without it moss and liverwort had no profile of their own for these
# three, and every moss took an Ala-to-Val "edit" in petD from the
# angiosperm profile.
SHORT_FIRST_EXON = ("petB", "petD", "rpl16")


def plastid_profile_genes() -> list[str]:
    """The plastid genes with a bundled land-plant editing profile."""
    return sorted(p.stem for p in (PIPELINE_DIR / "editing" / "reference" / "profiles" / "pltd").glob("*.pssm"))


def fetch_plastid_proteins(lineage: str) -> Path:
    """Write `lineage`'s REFSEQ_PLASTID_PROTEINS_NAME from its RefSeq
    plastomes: every single-exon CDS of a profiled gene whose genomic
    translation is the annotated protein, bar C-to-U edits, and the second
    exon's part of the SHORT_FIRST_EXON genes. GenBank records are cached
    beside the lineage's plastid exon sets."""
    if lineage not in rg.LINEAGE_PLASTOMES:
        raise ValueError(f"lineage must be one of {sorted(rg.LINEAGE_PLASTOMES)}, not {lineage!r}")
    gb_dir = rg.reference_dir(PIPELINE_DIR, "pltd", lineage) / "genbank"
    wanted = {g.lower(): g for g in plastid_profile_genes()}
    lines, n_bad = [], 0
    for acc, organism in rg.LINEAGE_PLASTOMES[lineage].items():
        rec = SeqIO.read(fetch(acc, gb_dir), "genbank")
        seen = set()
        for f in rec.features:
            if f.type != "CDS" or "gene" not in f.qualifiers or "translation" not in f.qualifiers:
                continue
            name = f.qualifiers["gene"][0].lower()
            gene = wanted.get(PLASTID_SYNONYMS.get(name, name).lower())
            if gene is None or gene in seen:
                continue
            parts, exons = cds_exons(rec, f)
            dna, protein = exons[0], f.qualifiers["translation"][0]
            if len(parts) == 2 and gene in SHORT_FIRST_EXON and len(exons[0]) <= 9:
                # the second exon from its first whole codon, and the residues of those codons
                dna, protein = exons[1][-len(exons[0]) % 3 :], protein[-(-len(exons[0]) // 3) :]
            elif len(parts) != 1:
                continue
            ok, note = verify(dna, protein)
            if not ok:
                n_bad += 1
                log.warning(f"{acc} {organism} {gene}: {note} - not used")
                continue
            seen.add(gene)
            lines.append(f">{acc}|{gene}\n{protein}\n")
    path = lin.editing_profiles(PIPELINE_DIR, lineage) / REFSEQ_PLASTID_PROTEINS_NAME
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(lines))
    log.info(f"{lineage}: {len(lines)} plastid proteins from {len(rg.LINEAGE_PLASTOMES[lineage])} plastomes, {n_bad} not verified -> {path}")
    return path


def refseq_plastid_proteins(lineage: str, exclude: tuple[str, ...] = ()) -> dict[str, dict[str, str]]:
    """{gene: {accession: protein}} from the lineage's bundled RefSeq
    plastid proteins, without the accessions in `exclude`."""
    path = lin.editing_profiles(PIPELINE_DIR, lineage) / REFSEQ_PLASTID_PROTEINS_NAME
    out: dict[str, dict[str, str]] = {}
    if path.exists():
        for header, seq in pssm.read_fasta(path).items():
            acc, gene = header.split("|")[:2]
            if acc not in exclude:
                out.setdefault(gene, {})[acc] = seq
    return out


def rna_proteins(lineage: str, exclude: tuple[str, ...] = (), organelle: str = "mito") -> dict[str, dict[str, str]]:
    """{gene: {"RNA_<species>": protein}} from the lineage's bundled
    RNA-edited proteins; `exclude` holds such RNA_<species> names. For a
    lineage in lin.SITE_LINEAGES the residues editing made stay in lower case."""
    name = RNA_PROTEINS_NAME if organelle == "mito" else RNA_PLASTID_PROTEINS_NAME
    path = lin.editing_profiles(PIPELINE_DIR, lineage) / name
    out: dict[str, dict[str, str]] = {}
    if path.exists():
        for header, seq in pssm.read_fasta(path).items():
            species, gene = header.split("|")[:2]
            if f"RNA_{species}" not in exclude:
                out.setdefault(gene, {})[f"RNA_{species}"] = seq if lineage in lin.SITE_LINEAGES else seq.upper()
    return out


def genus_groups(names, organisms: dict[str, str]) -> dict[str, str]:
    """Reference protein name -> genus, for counting editing sites
    (pssm.build_profile): RNA_<Genus>_<species>, or an accession of
    `organisms` ({accession: "Genus species"})."""
    return {n: n.split("_")[1] if n.startswith("RNA_") else organisms.get(n, n).split()[0] for n in names}


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
    (<lineage>/pltd/<gene>.pssm) from its RNA-edited plastid proteins and
    its RefSeq ones (moss, liverwort), one per gene with enough of them,
    leaving out the RNA_<species> and accessions in `exclude`. Returns the
    profiles directory."""
    # (angiosperm: the default plastid profiles, from RNA-edited proteins of
    # DToL angiosperms, in the same place as a lineage's)
    if lineage not in rg.LINEAGE_PLASTID and lineage != lin.DEFAULT_LINEAGE:
        raise ValueError(f"lineage must be one of {sorted((lin.DEFAULT_LINEAGE, *rg.LINEAGE_PLASTID))}, not {lineage!r}")
    proteins = rna_proteins(lineage, tuple(exclude), "pltd")
    for gene, by_acc in refseq_plastid_proteins(lineage, tuple(exclude)).items():
        proteins.setdefault(gene, {}).update(by_acc)
    # RefSeq plastomes: a column must be in half of them (as for the RefSeq
    # mito profiles, LINEAGE_GAP_THRESHOLD)
    gap_threshold = 0.5 if lineage in rg.LINEAGE_PLASTOMES else LINEAGE_GAP_THRESHOLD[lineage]
    min_sequences = LINEAGE_MIN_SEQUENCES.get(lineage, profiles.MIN_SEQUENCES)
    if lineage == lin.DEFAULT_LINEAGE:
        min_sequences = profiles.MIN_SEQUENCES
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
        columns = pssm.align_and_build(
            by_name, work_dir, gene, gap_threshold, groups=genus_groups(by_name, {}), sites=lineage in lin.SITE_LINEAGES, **sharp
        )
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
        edited_only=lineage == lin.DEFAULT_LINEAGE,
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
        columns = pssm.align_and_build(
            proteins[gene],
            work_dir,
            gene,
            LINEAGE_GAP_THRESHOLD[lineage],
            groups=genus_groups(proteins[gene], {**mitogenomes, **SPLIT_CCMFN_MITOGENOMES}),
            # (not the split ccmFn halves, whose proteins are as annotated, unmarked)
            sites=lineage in lin.SITE_LINEAGES and gene not in SPLIT_CCMFN_GENES,
            **sharp,
        )
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
        help="pltd: the lineage's plastid profiles, from its RNA-edited or RefSeq plastid proteins",
    )
    ap.add_argument(
        "--fetch",
        action="store_true",
        help="with --organelle pltd: first rewrite the lineage's RefSeq plastid proteins from GenBank (moss, liverwort)",
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
        if args.fetch:
            fetch_plastid_proteins(args.lineage)
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
