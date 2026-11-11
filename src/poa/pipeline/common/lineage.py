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
"""
from __future__ import annotations

import csv
from pathlib import Path

DEFAULT_LINEAGE = "angiosperm"
LINEAGES = (DEFAULT_LINEAGE, "moss")
CHOICES = ("auto", *LINEAGES)
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
