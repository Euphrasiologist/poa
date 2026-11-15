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

Output: analysis/denovo_annotation/results/gene_calls.tsv

src/04_build_gene_calls.py is the command-line wrapper (main() here).
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from typing import Iterable

from poa.pipeline.common import io_utils
from poa.pipeline.common import species_discovery as sd
from poa.pipeline.common.errors import StageError
from poa.pipeline.denovo_annotation import targets

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
        n_species = 0
        for r in targets.assemblies(root, organelle, species):
            rows = gff_to_rows(gff_dir / f"{r.species}.{organelle}.gff", r.species, organelle, work_dir / "nhmmscan")
            if rows:
                all_rows.extend(rows)
                n_species += 1
        print(f"[info] build_gene_calls: {organelle}: {n_species} species with >=1 call", file=sys.stderr)

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
    print(
        f"[info] build_gene_calls: {len(all_rows)} rows "
        f"({'species-scoped merge' if species else 'full rewrite'}) -> {out_path}",
        file=sys.stderr,
    )
    return out_path


def main(argv: list[str] | None = None) -> int:
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
        print(f"[err] {e}", file=sys.stderr)
        return 1
    return 0
