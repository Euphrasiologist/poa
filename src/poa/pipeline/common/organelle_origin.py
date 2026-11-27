"""Which organelle a gene call belongs to, for the genes both organelles
have (ribosomal proteins: rpl2, rpl16, rps3, rps4, rps12, rps14, rps19).

A call's protein is scored against the same gene's references in both
organelles - the profile poa corrects that gene with in each (editing's
single-exon profiles, or trans_splicing's whole-gene template for a spliced
gene) - by a local alignment to the profile, as orfedit scores a gene. A
call scoring at least FOREIGN_MIN_RATIO times as well, and at least
FOREIGN_MIN_SCORE, against the other organelle's profile is that
organelle's DNA: plastid DNA in a mitogenome (a "MTPT"), or mitochondrial
DNA in a plastome. A gene's native copy scores far better against its own
organelle's profile (the other's scores 0.02 to 0.37 as well, at the median
per gene, over the DToL dataset).

Over the dataset, 148 mito calls in 131 species score more than twice as
well against the plastid profile, and none less than 81 against it: rpl2
84 (mostly grasses and Apiaceae, reported as an incomplete mito rpl2), rps14
25, rpl16 18, rps12 12, rps19 5, rps4 4. Between 1 and 2 times as well,
the calls are alignments that failed against both profiles (scores of 9 to
16), but one. No plastid call scores better against the mito profile with
a score of any size. The references are the default (angiosperm) ones in
both organelles, for every lineage: what decides is which organelle a gene
looks like, which they tell as well for a moss or a fern."""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

from poa import PIPELINE_DIR

ORIGIN_GENES = ("rpl2", "rpl16", "rps3", "rps4", "rps12", "rps14", "rps19")
FOREIGN_MIN_RATIO = 2.0
FOREIGN_MIN_SCORE = 50.0
# cost of a gap position in the local alignment (one residue or column left out)
GAP = -4.0


def profile_paths(organelle: str, gene: str) -> list[Path]:
    """Where `gene`'s whole-protein profile for `organelle` ("mito" or
    "pltd") may be, first choice first."""
    ed, ts = PIPELINE_DIR / "editing" / "reference", PIPELINE_DIR / "trans_splicing" / "reference"
    if organelle == "mito":
        return [ed / "profiles_refseq" / "mito" / f"{gene}.pssm", ts / "profiles" / gene / "whole_gene.pssm", ed / "profiles" / "mito" / f"{gene}.pssm"]
    return [
        ed / "lineages" / "angiosperm" / "pltd" / f"{gene}.pssm",
        ts / "plastid" / "profiles" / gene / "whole_gene.pssm",
        ed / "profiles" / "pltd" / f"{gene}.pssm",
    ]


@lru_cache(maxsize=None)
def load_profile(organelle: str, gene: str) -> tuple[dict[str, int], np.ndarray] | None:
    """(residue -> column of the score matrix, columns x 20 scores) of the
    first profile there is, or None."""
    for path in profile_paths(organelle, gene):
        if path.exists():
            rows = [line.rstrip("\n").split("\t") for line in open(path) if line.strip()]
            residues = rows[0][:20]
            return {a: i for i, a in enumerate(residues)}, np.array([[float(x) for x in r[:20]] for r in rows[1:]])
    return None


def local_score(protein: str, profile: tuple[dict[str, int], np.ndarray]) -> float:
    """Best local alignment score of `protein` (deleted columns `-` and
    stops left out, case ignored) to `profile`, linear gaps of GAP."""
    index, scores = profile
    residues = [index[c] for c in str(protein).upper() if c in index]
    n = scores.shape[0]
    j = np.arange(n + 1)
    row = np.zeros(n + 1)
    best = 0.0
    for a in residues:
        diag = np.concatenate(([0.0], row[:-1] + scores[:, a]))
        here = np.maximum(np.maximum(diag, row + GAP), 0.0)
        # a run of profile columns left out: row[j] = max over k <= j of here[k] + GAP * (j - k)
        row = np.maximum(np.maximum.accumulate(here - GAP * j) + GAP * j, 0.0)
        best = max(best, float(row.max()))
    return best


def call_key(gene, contig=None, start=None, end=None) -> str:
    """What a call is known by: an edited copy (several can share a gene's
    name) by its location too, `gene@contig:start-end`; a reconstructed
    gene by its name."""
    if contig is None or start is None or pd.isna(start):
        return str(gene)
    return f"{gene}@{contig}:{int(start)}-{int(end)}"


def foreign_calls(calls: pd.DataFrame) -> dict[tuple[str, str, str], str]:
    """(species, organelle, call_key) -> the other organelle ("pltd" or
    "mito") for each call in `calls` (columns species, organelle, gene,
    protein, and for edited copies contig, corrected_start, corrected_end)
    whose protein is that organelle's."""
    out: dict[tuple[str, str, str], str] = {}
    if calls.empty or "protein" not in calls.columns:
        return out
    located = {"contig", "corrected_start", "corrected_end"} <= set(calls.columns)
    for r in calls.dropna(subset=["protein"]).itertuples():
        name = str(r.gene).split("|")[0]
        if name not in ORIGIN_GENES or r.organelle not in ("mito", "pltd"):
            continue
        other = "pltd" if r.organelle == "mito" else "mito"
        own_profile, other_profile = load_profile(r.organelle, name), load_profile(other, name)
        if own_profile is None or other_profile is None:
            continue
        own, there = local_score(r.protein, own_profile), local_score(r.protein, other_profile)
        if there >= FOREIGN_MIN_SCORE and there >= FOREIGN_MIN_RATIO * own:
            key = call_key(r.gene, r.contig, r.corrected_start, r.corrected_end) if located else call_key(r.gene)
            out[(r.species, r.organelle, key)] = other
    return out
