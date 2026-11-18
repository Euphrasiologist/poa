"""Reshape `03_build_combined_gff.py`'s one-GFF-per-species output into a
table with the SAME columns as `annotation/results/gene_calls.tsv` - a
deliberate schema match so this is a drop-in alternative anchor source for
`editing`/`trans_splicing` (`--gene-calls` flag on both), not a new format
they need to learn. The combined GFF is the canonical, inspectable output
of this module; this table is a derived convenience for those two
consumers, not a second source of truth.

Caveat worth remembering: `score` units differ by source column (`oatkDB`:
HMMER bitscore, higher better; `tRNAscan-SE`: its own confidence score,
higher better; `barrnap:*`: an E-value, *lower* better) - never compare
`score` across sources, only within one. The combined GFF is not
deduplicated across sources (see `03_build_combined_gff.py`), so a species
can show e.g. both an `oatkDB` and a `tRNAscan-SE` call for the same tRNA -
both are kept here too, as independent rows.

Beyond the shared 8 columns, four more are appended (readers select by
name, so this stays a drop-in): `source` (oatkDB|tRNAscan-SE|barrnap),
and for oatkDB hits `hmm_from`/`hmm_to`/`model_len` - where along the
gene's HMM the hit lies, joined back from nhmmscan's tblout (the GFF
carries only genomic coordinates). That's what lets `gff_export` tell a
gene's successive exons (consecutive model ranges) from a second copy
(the same model range again).

Split ccmFn: in Brassicaceae (and Allium, Galeopsis) mitochondrial ccmFn
is two genes, ccmFN1 and ccmFN2 in RefSeq, kilobases apart. nhmmscan
finds them as two hits that meet in the HMM (around position 1360 of
2137) with no full-length hit. relabel_split_ccmfn() renames those hits
ccmFn1/ccmFn2, so gff_export doesn't chain them into one cis-spliced
gene and editing doesn't correct either piece against a full-length
profile (on Arabidopsis that over-extended ccmFN1: 1380 bp vs RefSeq's
1149). A species is split when, among its mito ccmFn oatkDB hits, none
covers SPLIT_MAX_SPAN of the HMM, the best N-terminal piece and the best
C-terminal piece meet within SPLIT_MAX_HMM_GAP HMM positions, and they
lie on different contigs or at least SPLIT_MIN_DISTANCE apart. Two
hits with a gap of 160-800 positions in HMM coverage between them
(Butomus, Papaver, Plantago, Sagina, Cerastium, Eriocaulon) are one gene
the HMM matches poorly in the middle, and stay ccmFn. On the DToL
dataset 37 species are split: 31 Brassicaceae, 3 Allium, Galeopsis
tetrahit, Lunaria annua and Matthiola sinuata.

Output: analysis/denovo_annotation/results/gene_calls.tsv

src/04_build_gene_calls.py is the command-line wrapper (main() here).
"""

from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path
from typing import Iterable

from poa.pipeline.common import io_utils
from poa.pipeline.common import species_discovery as sd
from poa.pipeline.common.errors import StageError
from poa.pipeline.denovo_annotation import targets
from poa.pipeline.common.log import setup_logging

log = logging.getLogger(__name__)

COLUMNS = [
    "species",
    "organelle",
    "contig_id",
    "start",
    "end",
    "gene",
    "score",
    "strand",
    "source",
    "hmm_from",
    "hmm_to",
    "model_len",
]


def read_tblout_models(path: Path) -> dict[tuple[str, str, int, int], tuple[int, int, int]]:
    """(contig, target, env_lo, env_hi) -> (hmm_from, hmm_to, model_len), 1-based
    env coords as nhmmscan reports them (the GFF's start/end are the envelope)."""
    out = {}
    if not path.exists():
        return out
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        f = line.split()
        target, contig = f[0], f[2]
        hmm_from, hmm_to, env_a, env_b, model_len = int(f[4]), int(f[5]), int(f[8]), int(f[9]), int(f[10])
        out[(contig, target, min(env_a, env_b), max(env_a, env_b))] = (hmm_from, hmm_to, model_len)
    return out


AA_3TO1 = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Gln": "Q",
    "Glu": "E",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
    "SeC": "U",
    "Sup": "X",
    "Pseudo": "X",
}

RRNA_NAME = {
    "mito": {"16S_rRNA": "rrn18", "23S_rRNA": "rrn26", "5S_rRNA": "rrn5", "4.5S_rRNA": "rrn4"},
    "pltd": {"16S_rRNA": "rrn16", "23S_rRNA": "rrn23", "5S_rRNA": "rrn5", "4.5S_rRNA": "rrn4.5"},
}


def parse_attrs(attr_str: str) -> dict[str, str]:
    out = {}
    for kv in attr_str.strip().rstrip(";").split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
    return out


def gff_to_rows(path: Path, species: str, organelle: str, nhmmscan_dir: Path) -> list[dict]:
    if not path.exists():
        return []
    models = read_tblout_models(nhmmscan_dir / f"{species}.{organelle}.filtered.tblout")
    rows = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 9:
            continue
        contig, source, feature, start, end, score_s, strand = f[0], f[1], f[2], f[3], f[4], f[5], f[6]
        attrs = parse_attrs(f[8])

        if source == "oatkDB":
            gene = attrs.get("target_name")
            if gene is None:
                continue
        elif source == "tRNAscan-SE":
            if feature != "tRNA":
                continue  # skip the redundant "exon" sub-feature rows
            isotype = attrs.get("isotype", "Undet")
            anticodon = attrs.get("anticodon", "NNN")
            if isotype in ("Undet", "Sup", "Pseudo") or "N" in anticodon:
                continue
            aa1 = AA_3TO1.get(isotype)
            if aa1 is None:
                continue
            gene = f"trn{aa1}-{anticodon.replace('T', 'U')}"
        elif source.startswith("barrnap"):
            if feature != "rRNA":
                continue
            gene = RRNA_NAME[organelle].get(attrs.get("Name", ""))
            if gene is None:
                continue
        else:
            continue  # unrecognised source - skip rather than guess

        try:
            score = float(score_s)
        except ValueError:
            score = float("nan")
        # GFF3 is 1-based inclusive; every other gene_calls.tsv-shaped table
        # in this repo (annotation/'s own, and what the Rust tools assume)
        # is 0-based half-open (BED-style) - convert here, once, at the
        # only place raw GFF coordinates enter this table, rather than
        # leaving a silent off-by-one between "compatible" tables.
        hmm_from = hmm_to = model_len = ""
        if source == "oatkDB":
            m = models.get((contig, gene, int(start), int(end)))
            if m:
                hmm_from, hmm_to, model_len = m
        rows.append(
            {
                "species": species,
                "organelle": organelle,
                "contig_id": contig,
                "start": int(start) - 1,
                "end": int(end),
                "gene": gene,
                "score": score,
                "strand": strand,
                "source": source.split(":")[0],
                "hmm_from": hmm_from,
                "hmm_to": hmm_to,
                "model_len": model_len,
            }
        )
    return rows


SPLIT_GENE = "ccmFn"
SPLIT_NAMES = ("ccmFn1", "ccmFn2")
SPLIT_MAX_SPAN = 0.8  # a hit covering this much of the HMM is a whole gene
SPLIT_MAX_HMM_GAP = 100  # nt of HMM between where piece 1 ends and piece 2 starts
# bp between the pieces on one contig: Brassicaceae pieces are 0.9 kb
# (Cardamine amara/pratensis/raphanifolia) to 240 kb apart; this only
# rules out one gene nhmmscan reported as two abutting hits
SPLIT_MIN_DISTANCE = 500
SPLIT_MIN_SCORE = 200.0


def relabel_split_ccmfn(rows: list[dict]) -> int:
    """Rename one species' mito ccmFn hits ccmFn1/ccmFn2 in place when the
    gene is split (module docstring). Returns the number of hits renamed."""
    hits = [
        r
        for r in rows
        if r["organelle"] == "mito" and r["gene"] == SPLIT_GENE and r["source"] == "oatkDB" and r["model_len"] != ""
    ]
    if not hits:
        return 0
    if any((r["hmm_to"] - r["hmm_from"] + 1) >= SPLIT_MAX_SPAN * r["model_len"] for r in hits):
        return 0
    strong = [r for r in hits if r["score"] >= SPLIT_MIN_SCORE]
    n_side = [r for r in strong if r["hmm_from"] < 0.2 * r["model_len"] and r["hmm_to"] < 0.75 * r["model_len"]]
    c_side = [r for r in strong if r["hmm_from"] > 0.5 * r["model_len"]]
    if not n_side or not c_side:
        return 0
    n, c = max(n_side, key=lambda r: r["score"]), max(c_side, key=lambda r: r["score"])
    if abs(c["hmm_from"] - n["hmm_to"]) > SPLIT_MAX_HMM_GAP:
        return 0
    apart = (
        n["contig_id"] != c["contig_id"] or max(n["start"], c["start"]) - min(n["end"], c["end"]) >= SPLIT_MIN_DISTANCE
    )
    if not apart:
        return 0
    breakpoint = (n["hmm_to"] + c["hmm_from"]) / 2
    for r in hits:
        r["gene"] = SPLIT_NAMES[0] if (r["hmm_from"] + r["hmm_to"]) / 2 < breakpoint else SPLIT_NAMES[1]
    return len(hits)


def build_gene_calls(
    root: Path, organelles: Iterable[str] = targets.ORGANELLES, *, species: set[str] | None = None
) -> Path:
    """gene_calls.tsv from the combined GFFs under data root `root`. With
    `species`, those species' rows are merged into the existing table
    rather than replacing it. Returns the table's path."""
    root = Path(root)
    organelles = targets.check_organelles(organelles)
    analysis_dir = root / "analysis"
    work_dir = analysis_dir / "denovo_annotation" / "work"
    gff_dir = analysis_dir / "denovo_annotation" / "results" / "gff"

    all_rows: list[dict] = []
    for organelle in organelles:
        n_species = n_split = 0
        for r in targets.assemblies(root, organelle, species):
            rows = gff_to_rows(gff_dir / f"{r.species}.{organelle}.gff", r.species, organelle, work_dir / "nhmmscan")
            if relabel_split_ccmfn(rows):
                n_split += 1
            if rows:
                all_rows.extend(rows)
                n_species += 1
        log.info(f"build_gene_calls: {organelle}: {n_species} species with >=1 call")
        if n_split:
            log.info(f"build_gene_calls: {organelle}: {n_split} species with ccmFn split into ccmFn1/ccmFn2")

    out_path = analysis_dir / "denovo_annotation" / "results" / "gene_calls.tsv"
    # Always write what this run actually computed to its own path first;
    # merge_species_scoped_tsv then decides how that relates to what's
    # already in out_path (replace entirely if unscoped, species-scoped
    # merge otherwise - see its docstring for why that distinction matters).
    new_path = work_dir / "gene_calls.new.tsv"
    new_path.parent.mkdir(parents=True, exist_ok=True)
    with open(new_path, "w") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        for row in all_rows:
            fh.write("\t".join(str(row[c]) for c in COLUMNS) + "\n")
    io_utils.merge_species_scoped_tsv(out_path, new_path, species)
    log.info(
        f"build_gene_calls: {len(all_rows)} rows "
        f"({'species-scoped merge' if species else 'full rewrite'}) -> {out_path}"
    )
    return out_path


def main(argv: list[str] | None = None) -> int:
    setup_logging()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--organelle", default="both", choices=["mito", "pltd", "both"])
    ap.add_argument("--species-list", default=None)
    args = ap.parse_args(argv)
    root = Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT") or os.getcwd())
    try:
        build_gene_calls(
            root,
            list(targets.ORGANELLES) if args.organelle == "both" else [args.organelle],
            species=sd.load_species_list(Path(args.species_list)) if args.species_list else None,
        )
    except StageError as e:
        log.error(e)
        return 1
    return 0
