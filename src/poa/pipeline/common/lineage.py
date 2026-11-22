"""Which reference lineage a species' mito annotation uses.

poa's editing profiles and spliced-gene templates were built from
angiosperms. Mosses differ in both: their mito genes carry a different,
phylum-wide intron set (nad1 3 exons, not 5 - see
trans_splicing/src/reference_genes.py), and they barely edit, so
angiosperm-built profiles under-fill every spliced gene and over-correct
the rest. Each lineage other than angiosperm has its own templates and
editing profiles, built from that lineage's RefSeq mitogenomes.

The call is made once, in editing/src/01_scan_editing.py: every mito gene
with a profile in each lineage is aligned (orfedit) against both, and the
species takes the lineage its genes score best against in total, provided
it wins by at least MIN_MARGIN per gene. Anything less clear-cut, or with
fewer than MIN_SHARED_GENES genes scored both ways, stays angiosperm - the
references everything else in poa was validated on. The result is
analysis/editing/results/lineage.tsv, which trans_splicing reads.

"fern" is a lineage of a different kind: it is never chosen automatically,
only with --lineage fern, and it changes the editing model. Ferns edit U
to C as well as C to U, and a fern organelle gene has genomic in-frame
stop codons that editing turns into Gln or Arg - about nine per
mitochondrial gene in Polypodiales. With the lineage set, orfedit and
transsplice may read such a stop through at STOP_EDIT_PENALTY. It has
exon templates of its own (trans_splicing/reference/lineages/fern/) and
editing profiles for the single-exon mito genes its references cover
(editing/reference/lineages/fern/), on which other U-to-C edits are
called too, at REVERSE_EDIT_PENALTY; the genes it has no profile for, and
its plastid genes, borrow the angiosperm references. It is for ferns
proper, quillworts (Isoetes) and hornworts; horsetails, clubmosses and
Selaginella do not edit this way and should stay on the default.
"""

from __future__ import annotations

import csv
from pathlib import Path

DEFAULT_LINEAGE = "angiosperm"
LINEAGES = (DEFAULT_LINEAGE, "moss")
# a lineage that borrows another's references: {name: the lineage whose
# templates and profiles it uses}
BORROWS = {"fern": DEFAULT_LINEAGE}
FERN = "fern"
CHOICES = ("auto", *LINEAGES, *BORROWS)
# Cost of reading an in-frame TAA/TAG as Gln or TGA as Arg (a U-to-C
# edit), per lineage; the same as a C-to-U edit. Calibrated on 43 DToL
# fern assemblies against 60 angiosperm ones: at 2.5, fern mito genes
# aligned over 99% of their profile go from 152 to 292 of 702, and U-to-C
# comes to 32% of mito and 6% of plastid edits (33-35% and 1-10% measured
# in fern transcriptomes); 4 and 6 lose alignments, and allowing other
# U-to-C edits as well (at 8 or 12) changes nothing for ferns.
STOP_EDIT_PENALTY = {"fern": 2.5}
# Cost of any other U-to-C edit (one that changes an amino acid, not a
# stop), for a lineage's mitochondrial genes aligned to its own profiles
# and exon templates only. Never for genes still on the borrowed
# (angiosperm) references or for plastid genes: against profiles that are
# not the lineage's own edited proteins a U-to-C edit is too easy to call
# (2.1 false ones per angiosperm mito gene at a cost of 2.5).
REVERSE_EDIT_PENALTY = {"fern": 5.0}
# Bases transsplice may slide an exon junction either way (its
# --junction-slide), for a lineage's mito genes on its own templates.
# Fern: nad5's 22 bp exon 3 sat a base or two off in all ten ferns with
# intron-spanning RNA reads to check it against; with the slide 135 of 144
# junctions agree with the reads, against 126. Not for the others: held
# out on the angiosperm RefSeq genomes it moves two plastid exons off.
JUNCTION_SLIDE = {"fern": 6}
# Calibrated on the 956-species mito dataset (ToL-ID clades as truth). Every
# moss has all 5 moss-profiled genes (atp4 ccmB ccmC mttB rps14) and wins by
# 88-150 per gene; dicots and monocots lose by 57-430, liverworts by -25 to
# +6. Ferns and lycophytes are the hard case - they have lost ccmB/ccmC, and
# on the 3 genes left 16 of them preferred moss by 55-93 per gene - so a
# lineage call needs 4 genes. With that, 174/174 mosses are called moss and
# nothing else is.
MIN_SHARED_GENES = 4
# mean per-gene score advantage (orfedit score) a lineage needs over
# angiosperm before it is used
MIN_MARGIN = 50.0
COLUMNS = ["species", "lineage", "decided_by", "n_genes", *[f"{x}_score" for x in LINEAGES]]


def reference(lineage: str) -> str:
    """The lineage whose templates and profiles `lineage` uses."""
    return BORROWS.get(lineage, lineage)


def stop_edit_penalty(lineage: str) -> float | None:
    """Cost of a stop read-through for `lineage`, or None if it has none."""
    return STOP_EDIT_PENALTY.get(lineage)


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


def choose(scores: dict[str, dict[str, float]]) -> tuple[str, int, dict[str, float]]:
    """scores: {lineage: {gene: orfedit score}} for one species.
    Returns (lineage, number of genes compared, {lineage: total over those genes})."""
    present = [x for x in LINEAGES if scores.get(x)]
    shared = set.intersection(*(set(scores[x]) for x in present)) if present else set()
    totals = {x: sum(scores[x][g] for g in shared) for x in present}
    if DEFAULT_LINEAGE not in totals or len(shared) < MIN_SHARED_GENES:
        return DEFAULT_LINEAGE, len(shared), totals
    best = max(totals, key=totals.get)
    if best != DEFAULT_LINEAGE and (totals[best] - totals[DEFAULT_LINEAGE]) / len(shared) < MIN_MARGIN:
        best = DEFAULT_LINEAGE
    return best, len(shared), totals


def load(analysis_dir: Path) -> dict[str, str] | None:
    """{species: lineage} from lineage.tsv, or None if editing hasn't written one."""
    path = table_path(analysis_dir)
    if not path.exists():
        return None
    with open(path) as fh:
        return {r["species"]: r["lineage"] for r in csv.DictReader(fh, delimiter="\t")}
