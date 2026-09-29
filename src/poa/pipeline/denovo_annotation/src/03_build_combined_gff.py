#!/usr/bin/env python3
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
"""
from __future__ import annotations

import os
import argparse
import sys
from pathlib import Path

# CODE_DIR is the installed poa/pipeline/ (shared modules, bundled reference
# data); data/ and every module's results/ and work/ live under the data
# root - PLANT_ORGANELLE_DATA_ROOT, or the current directory if unset.
CODE_DIR = Path(__file__).resolve().parents[2]
ROOT_DIR = Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT") or os.getcwd())
ANALYSIS_DIR = ROOT_DIR / "analysis"
sys.path.insert(0, str(CODE_DIR / "common"))

import species_discovery as sd  # noqa: E402


def gff_data_lines(path: Path) -> list[str]:
    if not path.exists():
        return []
    return [line for line in path.read_text().splitlines() if line.strip() and not line.startswith("#")]


def sort_key(line: str) -> tuple:
    f = line.split("\t")
    return (f[0], int(f[3])) if len(f) > 3 and f[3].isdigit() else (f[0], 0)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organelle", default="both", choices=["mito", "pltd", "both"])
    ap.add_argument("--species-list", default=None)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    organelles = ["mito", "pltd"] if args.organelle == "both" else [args.organelle]
    species_filter = sd.load_species_list(Path(args.species_list)) if args.species_list else None
    work_dir = ANALYSIS_DIR / "denovo_annotation" / "work"
    out_dir = ANALYSIS_DIR / "denovo_annotation" / "results" / "gff"
    out_dir.mkdir(parents=True, exist_ok=True)

    for organelle in organelles:
        data_root = sd.repo_data_root(ROOT_DIR, organelle)
        resolved = sd.discover_all(data_root, organelle, species_filter=species_filter)
        n_written = 0
        for r in resolved:
            if r.status != "ok":
                continue
            out_path = out_dir / f"{r.species}.{organelle}.gff"
            if not args.force and out_path.exists() and out_path.stat().st_size > 0:
                continue

            lines = []
            lines += gff_data_lines(work_dir / "nhmmscan" / f"{r.species}.{organelle}.gff")
            lines += gff_data_lines(work_dir / "trnascan" / f"{r.species}.{organelle}.gff")
            lines += gff_data_lines(work_dir / "barrnap" / f"{r.species}.{organelle}.gff")
            if not lines:
                continue
            lines.sort(key=sort_key)
            out_path.write_text("##gff-version 3\n" + "\n".join(lines) + "\n")
            n_written += 1
        print(f"[info] build_combined_gff: {organelle}: {n_written} species -> {out_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
