"""Combine the three independent annotation sources - `nhmmscan` (protein-
coding, converted to GFF by `nhmmscan_gff.py`), `tRNAscan-SE`, and `barrnap` -
into ONE per-species GFF3, sorted by contig then start. Same append/strip-
header/sort pattern as `mito_structural_variation/annotation/src/5_get_rna_genes.bash`
and `6_get_rrna_genes.bash` use to build up one combined GFF per species.

Deliberately not deduplicated: oatkDB's own `.fam` database includes tRNA/
rRNA gene-family models too (built with `oatkdb`'s default `--no-trna`/
`--no-rrna` flags *not* set), so `nhmmscan`'s own hits can already include
some `trn*`/`rrn*` calls alongside `tRNAscan-SE`/`barrnap`'s independent
ones for the same genes. Keeping all three, tagged by source (column 2),
is a feature - agreement across independent methods is itself useful
signal, and dropping rows would hide disagreement rather than surface it.
`04_build_gene_calls.py` reshapes this combined GFF into the
`gene_calls.tsv`-shaped table `editing`/`trans_splicing` consume.

Output: analysis/denovo_annotation/results/gff/<species>.<organelle>.gff

src/03_build_combined_gff.py is the command-line wrapper (main() here).
"""

from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path
from typing import Iterable

from poa.pipeline.common import species_discovery as sd
from poa.pipeline.common.errors import StageError
from poa.pipeline.denovo_annotation import targets
from poa.pipeline.common.log import setup_logging

log = logging.getLogger(__name__)

SOURCES = ("nhmmscan", "trnascan", "barrnap")  # work/ subdirectories


def gff_data_lines(path: Path) -> list[str]:
    if not path.exists():
        return []
    return [line for line in path.read_text().splitlines() if line.strip() and not line.startswith("#")]


def sort_key(line: str) -> tuple:
    f = line.split("\t")
    return (f[0], int(f[3])) if len(f) > 3 and f[3].isdigit() else (f[0], 0)


def build_combined_gff(
    root: Path, organelles: Iterable[str] = targets.ORGANELLES, *, species: set[str] | None = None, force: bool = False
) -> Path:
    """One combined GFF per assembly under data root `root`, from the three
    scans' GFFs. Returns the output directory."""
    root = Path(root)
    organelles = targets.check_organelles(organelles)
    work_dir = root / "analysis" / "denovo_annotation" / "work"
    out_dir = root / "analysis" / "denovo_annotation" / "results" / "gff"
    out_dir.mkdir(parents=True, exist_ok=True)

    for organelle in organelles:
        n_written = 0
        for r in targets.assemblies(root, organelle, species):
            out_path = out_dir / f"{r.species}.{organelle}.gff"
            if not force and out_path.exists() and out_path.stat().st_size > 0:
                continue

            lines = []
            for source in SOURCES:
                lines += gff_data_lines(work_dir / source / f"{r.species}.{organelle}.gff")
            if not lines:
                continue
            lines.sort(key=sort_key)
            out_path.write_text("##gff-version 3\n" + "\n".join(lines) + "\n")
            n_written += 1
        log.info(f"build_combined_gff: {organelle}: {n_written} species -> {out_dir}")
    return out_dir


def main(argv: list[str] | None = None) -> int:
    setup_logging()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--organelle", default="both", choices=["mito", "pltd", "both"])
    ap.add_argument("--species-list", default=None)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args(argv)
    root = Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT") or os.getcwd())
    try:
        build_combined_gff(
            root,
            list(targets.ORGANELLES) if args.organelle == "both" else [args.organelle],
            species=sd.load_species_list(Path(args.species_list)) if args.species_list else None,
            force=args.force,
        )
    except StageError as e:
        log.error(e)
        return 1
    return 0
