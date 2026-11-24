"""Which reference lineage a species' mito annotation uses.

poa's editing profiles and spliced-gene templates were built from
angiosperms. Mosses differ in both: their mito genes carry a different,
phylum-wide intron set (nad1 3 exons, not 5 - see
trans_splicing/src/reference_genes.py), and they barely edit, so
angiosperm-built profiles under-fill every spliced gene and over-correct
the rest. Each lineage other than angiosperm has its own templates and
editing profiles, built from that lineage's RefSeq mitogenomes.

The call is made once, in editing/src/01_scan_editing.py: every mito gene
with a profile in a lineage is aligned (orfedit) against it and against the
angiosperm profiles, and the species takes the lineage that beats the
angiosperm profiles by the most per gene, provided it reaches that
lineage's MIN_MARGIN. Anything less clear-cut, or with fewer than
MIN_SHARED_GENES genes scored both ways, stays angiosperm - the references
everything else in poa was validated on. The result is
analysis/editing/results/lineage.tsv, which trans_splicing reads.

Liverworts have a lineage of this kind too: an intron set of their own
(reference_genes.LINEAGE_EXTRA_SPLICED), none of it trans-spliced, and
little editing (none in the complex thalloids). Its templates and
profiles come from 27 RefSeq mitogenomes.

"fern" is a lineage of a different kind: it sits on top of the angiosperm
references and changes the editing model. It is looked for, with horsetail
and lycophyte, among the species the call above leaves on the angiosperm
references (DETECTED, choose), or forced with --lineage fern. Ferns edit U
to C as well as C to U, and a fern organelle gene has genomic in-frame
stop codons that editing turns into Gln or Arg - about nine per
mitochondrial gene in Polypodiales. With the lineage set, orfedit and
transsplice may read such a stop through at STOP_EDIT_PENALTY. It has
exon templates of its own (trans_splicing/reference/lineages/fern/) and
editing profiles for the single-exon mito genes its references cover
(editing/reference/lineages/fern/), on which other U-to-C edits are
called too, at REVERSE_EDIT_PENALTY; the genes it has no profile for
borrow the angiosperm references. It also has plastid exon templates and
editing profiles of its own (reference_genes.LINEAGE_PLASTID), used with
the stop read-through but without the other U-to-C edits. Forced, it also
suits quillworts (Isoetes) and hornworts, which edit the same way; they
are not detected.

"horsetail" and "lycophyte" are chosen the same way and likewise sit on
top of the angiosperm references. Horsetail
(Equisetum, the one living genus) edits C to U only, so it changes
nothing in the editing model; its exon templates and editing profiles
come from Equisetum RNA-seq, there being no annotated horsetail
mitogenome. Lycophyte is for clubmosses (Lycopodiaceae): an intron set of
their own, and a handful of in-frame stop codons edited away, so it reads
stops through like fern. The four RefSeq clubmoss mitogenomes give the
exon structures but unedited proteins, so its references are the RefSeq
annotation of Phylloglossum with the editing its total-RNA library shows
applied, plus loci of three Darwin Tree of Life clubmosses. Selaginella
and Isoetes mitogenomes are unlike either and are not covered.

Genes with introns in one of these lineages only (fern atp1, cox1, rps14;
horsetail cox1; lycophyte atp6, atp9, cob, cox1, cox3, nad3, rps14) are
reconstructed from that lineage's exon templates, and editing leaves them
alone (reference_genes.LINEAGE_EXTRA_SPLICED).
"""

from __future__ import annotations

import csv
from pathlib import Path

DEFAULT_LINEAGE = "angiosperm"
LINEAGES = (DEFAULT_LINEAGE, "moss", "liverwort")
# a lineage that borrows another's references: {name: the lineage whose
# templates and profiles it uses}
BORROWS = {"fern": DEFAULT_LINEAGE, "horsetail": DEFAULT_LINEAGE, "lycophyte": DEFAULT_LINEAGE}
FERN = "fern"
CHOICES = ("auto", *LINEAGES, *BORROWS)
# Cost of reading an in-frame TAA/TAG as Gln or TGA as Arg (a U-to-C
# edit), per lineage; the same as a C-to-U edit. Calibrated on 43 DToL
# fern assemblies against 60 angiosperm ones: at 2.5, fern mito genes
# aligned over 99% of their profile go from 152 to 292 of 702, and U-to-C
# comes to 32% of mito and 6% of plastid edits (33-35% and 1-10% measured
# in fern transcriptomes); 4 and 6 lose alignments, and allowing other
# U-to-C edits as well (at 8 or 12) changes nothing for ferns.
# Lycophyte (Lycopodiaceae): the four RefSeq mitogenomes annotate the same
# few in-frame stops as edited (cob codon 46 TAG, nad1 177 and nad5 544
# TAA to Gln, rpl2 36 TGA to Arg), so their genes need the read-through
# too, at the same cost.
STOP_EDIT_PENALTY = {"fern": 2.5, "lycophyte": 2.5}
# Cost of any other U-to-C edit (one that changes an amino acid, not a
# stop), for a lineage's mitochondrial genes aligned to its own profiles
# and exon templates only. Never for genes still on the borrowed
# (angiosperm) references or for plastid genes: against profiles that are
# not the lineage's own edited proteins a U-to-C edit is too easy to call
# (2.1 false ones per angiosperm mito gene at a cost of 2.5).
REVERSE_EDIT_PENALTY = {"fern": 5.0}
# Cost of a C-to-U edit for a lineage's genes on its own profiles and exon
# templates, where the default (2.5, editing/scan.py) is too cheap.
# Moss: mosses edit little, and at 2.5 a variable residue is "edited"
# towards the other mosses'. Against RNA-seq of 16 Darwin Tree of Life
# mosses, 579 edits were called in 626 plastid genes and 127 of them were
# in the RNA; four mosses whose RNA shows 1 to 5 plastid edits drew 15 to
# 23 calls each (held out by genus), one or none right. At 6 those are
# gone (0 to 2 a species; at 4, 4 to 6 are left) and half of what is
# called is right. The cost is the real edits of the mosses that do edit
# their plastids: Brachythecium rutabulum and Fissidens dubius, about 80
# in the RNA of each, 60 found at 2.5 and 10 to 15 at 6 - with as many
# false calls as right ones at 2.5. A false edit is a wrong residue in the
# protein, a missed one the genome's own translation. Mitochondrial genes
# likewise: 53 edits called in 78 genes of 13 mosses at 2.5, 11 of them in
# the RNA.
# Not liverworts: the leafy ones edit, and a higher cost loses their edits
# as fast as their false calls.
EDIT_PENALTY = {"moss": 6.0}
# Lineages whose reference proteins say which residues editing made (lower
# case, pssm.mark_edited), so that their profiles carry editing sites
# (pssm.EDIT_COLUMN). All of a lineage's proteins or none: a profile with
# sites from some of them only would charge the others' edits as novel.
SITE_LINEAGES = (DEFAULT_LINEAGE,)
# Costs of a C-to-U edit for a profile with editing sites: at a column
# edited in its references, and at any other. A profile without sites
# takes the ordinary cost (editing/scan.DEFAULT_EDIT_PENALTY, or
# EDIT_PENALTY) at every column.
# Held out by genus on 109 DToL angiosperms with RNA-seq, at the ordinary
# cost of 2.5 everywhere 40% of the edits called in single-exon plastid
# genes were in the RNA, and 94% of the mitochondrial ones; at 10 away
# from the known columns, 93% and 97%. A known column can then cost less
# than 2.5: at 1.5 the mitochondrial sites found go from 71% to 77%, and
# the edits annotated in 14 RefSeq mitogenomes from 84% to 89%.
KNOWN_EDIT_PENALTY = 1.5
NOVEL_EDIT_PENALTY = 10.0
# ...for a species that edits where the references do. The default
# references also serve gymnosperms, hornworts, quillworts and spikemosses,
# which edit heavily at columns no angiosperm does: at these costs Ginkgo's
# plastid calls went from 299 to 23 and a hornwort's from 616 to 22. So a
# species is first corrected at the ordinary cost everywhere, and takes the
# two costs only if enough of the C-to-U edits then called in its genes
# fall at known columns: SITE_MIN_SHARE of its mitochondrial ones, given
# SITE_MIN_CALLS of them, or else of its plastid ones (a plastid alone has
# many more calls away from the known columns). Of 89 angiosperms (70
# DToL assemblies and 19 RefSeq mitogenomes; references from 111 species)
# 85 take the two costs, with
# 67 to 99% of their mitochondrial calls at known columns (median 95%);
# four sedges, whose mitochondrial genes have drifted furthest, have 58
# to 62% and keep the one cost. Ginkgo has 42%, a yew 38%, a hornwort
# 32%, a spikemoss 25%, a quillwort 25%; ferns, mosses and liverworts 14
# to 41%. Plastid: angiosperms 14 to 59% (median 36%; one rock-rose 4%),
# the others 5% or less.
SITE_MIN_SHARE = {"mito": 0.65, "pltd": 0.10}
SITE_MIN_CALLS = {"mito": 20, "pltd": 10}
# Bases transsplice may slide an exon junction either way (its
# --junction-slide), for a lineage's mito genes on its own templates.
# Fern: nad5's 22 bp exon 3 sat a base or two off in all ten ferns with
# intron-spanning RNA reads to check it against; with the slide 135 of 144
# junctions agree with the reads, against 126. Not for the others: held
# out on the angiosperm RefSeq genomes it moves two plastid exons off.
JUNCTION_SLIDE = {"fern": 6}
# How a species' lineage is chosen (--lineage auto). Its mito genes are
# aligned (orfedit, plain C-to-U editing) to the default profiles and to
# each other lineage's, and a lineage's margin is the mean score advantage
# of its profiles over the default ones, over the genes scored both ways.
# A species takes the lineage with the largest margin among those that
# reach their MIN_MARGIN over at least MIN_SHARED_GENES genes, and stays
# angiosperm - the references everything else in poa was validated on -
# when none does.
#
# A lineage call needs MIN_SHARED_GENES genes: ferns and lycophytes have
# lost ccmB and ccmC, and on the three moss-profiled genes they have left,
# 16 of them once preferred moss.
#
# The margins, calibrated 2026-11-24 on 1,003 species: fresh poa runs of
# the dataset's 64 liverworts, 30 ferns, 28 mosses, 9 horsetails and 4
# clubmosses, and the dataset's stored gene calls for the rest
# (dev/rna_editome/detect_cal.py). Every bryophyte and lycophyte fits
# another's profiles better than the angiosperm ones, so each lineage has
# a margin of its own, set midway between its lowest member and the
# highest of the others that it has to keep out:
# - moss: 175 mosses 257-324; liverworts up to 179, horsetails 75.
# - liverwort: liverworts 335-368; mosses 203-221 (they take moss, the
#   larger margin), a hornwort 160, clubmosses 153, ferns up to 101.
# - horsetail: horsetails 475-503; the nearest other, a fern, 42.
# - lycophyte: clubmosses 380-433; liverworts up to 144, mosses 137, a
#   hornwort 105, Selaginella 74.
# - fern: 28 ferns 173-288, Cryptogramma 140, Ophioglossum 87; the
#   horsetails reach 87 too and take horsetail; of the species that can
#   be left on the default references, mosses reach 22, a hornwort 16 and
#   angiosperms 12.
# The highest angiosperm (716) on any lineage's profiles is 50. These
# moved with the angiosperm profiles: on the first ones, built from
# unedited proteins, angiosperms themselves scored lower, moss needed 50
# and fern 150, and Cryptogramma and Ophioglossum were left out.
MIN_SHARED_GENES = 4
MIN_MARGIN = {"moss": 215.0, "liverwort": 275.0, "fern": 60.0, "horsetail": 250.0, "lycophyte": 250.0}
# The lineages that sit on the default references and have mito editing
# profiles of their own: looked for among the species left on the default
# ones (editing/scan.py own_profile_margins), which are then corrected
# again as that lineage.
DETECTED = ("fern", "horsetail", "lycophyte")
MARGIN_COLUMNS = [f"{x}_margin" for x in (*LINEAGES[1:], *DETECTED)]
# edit_sites: yes for a species corrected at the editing-site costs
# (follows_sites); <organelle>_site_calls: its C-to-U calls at the ordinary
# cost in genes whose profile has sites, and how many were at known columns
SITE_COLUMNS = ["edit_sites", *[f"{o}_{x}" for o in ("mito", "pltd") for x in ("site_calls", "known_site_calls")]]
COLUMNS = ["species", "lineage", "decided_by", "n_genes", *MARGIN_COLUMNS, *SITE_COLUMNS]


def reference(lineage: str) -> str:
    """The lineage whose templates and profiles `lineage` uses."""
    return BORROWS.get(lineage, lineage)


def stop_edit_penalty(lineage: str) -> float | None:
    """Cost of a stop read-through for `lineage`, or None if it has none."""
    return STOP_EDIT_PENALTY.get(lineage)


def edit_penalty(lineage: str) -> float | None:
    """Cost of a C-to-U edit for `lineage`'s genes on its own references,
    or None for the default."""
    return EDIT_PENALTY.get(lineage)


def junction_slide(lineage: str) -> int | None:
    """How far transsplice may slide an exon junction for `lineage`'s mito
    genes on its own templates, or None for no slide."""
    return JUNCTION_SLIDE.get(lineage)


def reverse_edit_penalty(lineage: str) -> float | None:
    """Cost of any other U-to-C edit for `lineage`'s mito genes on its own
    profiles and templates, or None if it calls none."""
    return REVERSE_EDIT_PENALTY.get(lineage)


def table_path(analysis_dir: Path) -> Path:
    return Path(analysis_dir) / "editing" / "results" / "lineage.tsv"


def trans_splicing_reference(code_dir: Path, lineage: str) -> Path:
    """reference/ for angiosperm (where it always was), reference/lineages/<x>/ otherwise."""
    base = Path(code_dir) / "trans_splicing" / "reference"
    return base if lineage == DEFAULT_LINEAGE else base / "lineages" / lineage


def editing_profiles(code_dir: Path, lineage: str) -> Path:
    """Bundled editing profiles for a non-default lineage (<dir>/mito/<gene>.pssm)."""
    return Path(code_dir) / "editing" / "reference" / "lineages" / lineage


def margins(scores: dict[str, dict[str, float]]) -> dict[str, tuple[float, int]]:
    """scores: {lineage: {gene: orfedit score}} for one species, the default
    lineage among them. Returns {other lineage: (mean score advantage over
    the default profiles, genes scored both ways)}."""
    base = scores.get(DEFAULT_LINEAGE) or {}
    out = {}
    for name, genes in scores.items():
        shared = set(genes) & set(base)
        if name != DEFAULT_LINEAGE and shared:
            out[name] = (sum(genes[g] - base[g] for g in shared) / len(shared), len(shared))
    return out


def choose(margins: dict[str, tuple[float, int]], among: tuple[str, ...] | None = None) -> str:
    """The lineage a species takes, from its `margins` ({lineage: (margin,
    genes compared)}): the largest margin among the lineages (of `among`,
    if given) that reach their MIN_MARGIN over MIN_SHARED_GENES genes, else
    the default lineage."""
    ok = {
        x: m
        for x, (m, n) in margins.items()
        if (among is None or x in among) and x in MIN_MARGIN and n >= MIN_SHARED_GENES and m >= MIN_MARGIN[x]
    }
    return max(ok, key=ok.get) if ok else DEFAULT_LINEAGE


def follows_sites(calls: dict[str, tuple[int, int]]) -> bool:
    """Whether a species edits where the default references do, from its
    C-to-U calls at the ordinary cost, {organelle: (at known columns, all)}:
    judged on its mitochondrial calls if there are SITE_MIN_CALLS, else on
    its plastid ones, and not at all (False) with too few of either."""
    for organelle in ("mito", "pltd"):
        known, n = calls.get(organelle, (0, 0))
        if n >= SITE_MIN_CALLS[organelle]:
            return known / n >= SITE_MIN_SHARE[organelle]
    return False


def load(analysis_dir: Path) -> dict[str, str] | None:
    """{species: lineage} from lineage.tsv, or None if editing hasn't written one."""
    path = table_path(analysis_dir)
    if not path.exists():
        return None
    with open(path) as fh:
        return {r["species"]: r["lineage"] for r in csv.DictReader(fh, delimiter="\t")}


def load_site_species(analysis_dir: Path) -> set[str]:
    """The species lineage.tsv marks as corrected at the editing-site costs."""
    path = table_path(analysis_dir)
    if not path.exists():
        return set()
    with open(path) as fh:
        return {r["species"] for r in csv.DictReader(fh, delimiter="\t") if r.get("edit_sites") == "yes"}
