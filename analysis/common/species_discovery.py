#!/usr/bin/env python3
"""Resolve per-species organelle assembly files robustly.

Handles two confirmed real-world quirks of data/{mito,plastid}/<Species>/:

1. Filenames may or may not carry a "<Species>.k1001.s31.c<N>." parameter
   infix before the ".mito."/".pltd." suffix (982 species with the infix,
   263 without, in the mito set). Discovery must never assume the infix is
   present or absent.

2. A small number of species have MORE THAN ONE candidate run (e.g. a bare
   run plus a later revision run with a parameter infix). Naively preferring
   "the file with an infix" or "the newest file" is WRONG:
     - Lathyrus_aphaca: the infixed run (newer mtime) is a completely EMPTY
       (0-byte) failed rerun; the bare run is the valid one.
     - Solanum_nigrum: the bare run is an orphaned .mito.gfa with no sibling
       .ctg.fasta/.ctg.bed/.annot_mito.txt; the infixed run is the complete,
       usable one.
   The correct rule is: prefer whichever candidate run has the most complete
   set of non-empty sibling files (gfa, ctg.fasta, ctg.bed, annot), breaking
   remaining ties by larger ctg.fasta size, then newest mtime. Every
   multi-candidate species is logged so the choice is auditable.

This module is intentionally only about *file resolution* - it does not
judge assembly quality (empty graphs, dead-ends, etc.); that's
qc_basic_stats' job, using the paths this module resolves.
"""
from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass, fields
from pathlib import Path
from typing import Optional

ORGANELLE_SUFFIX = {
    "mito": {"gfa": ".mito.gfa", "ctg_bed": ".mito.ctg.bed", "ctg_fasta": ".mito.ctg.fasta",
             "bed": ".mito.bed", "annot": ".annot_mito.txt"},
    "pltd": {"gfa": ".pltd.gfa", "ctg_bed": ".pltd.ctg.bed", "ctg_fasta": ".pltd.ctg.fasta",
             "bed": ".pltd.bed", "annot": ".annot_pltd.txt"},
}

REQUIRED_FOR_COMPLETENESS = ("gfa", "ctg_fasta", "ctg_bed", "annot")


@dataclass
class SpeciesFiles:
    species: str
    organelle: str
    gfa: Optional[str] = None
    ctg_bed: Optional[str] = None
    ctg_fasta: Optional[str] = None
    bed: Optional[str] = None
    annot: Optional[str] = None
    run_prefix: Optional[str] = None
    n_candidates: int = 0
    status: str = "missing_gfa"  # missing_gfa | no_ctg_fasta | ok

    def as_row(self) -> dict:
        return {f.name: getattr(self, f.name) for f in fields(self)}


def _nonempty(path: Optional[Path]) -> bool:
    return path is not None and path.exists() and path.stat().st_size > 0


def _candidate_prefixes(species_dir: Path, gfa_suffix: str) -> list[str]:
    """All run prefixes present, i.e. every '<prefix>' for '<prefix><gfa_suffix>'."""
    prefixes = []
    for p in species_dir.glob(f"*{gfa_suffix}"):
        prefix = p.name[: -len(gfa_suffix)]
        prefixes.append(prefix)
    return sorted(prefixes)


def _completeness_score(species_dir: Path, prefix: str, suffixes: dict) -> tuple:
    """Return a sortable tuple: (n_complete_nonempty_files, ctg_fasta_size, mtime)."""
    paths = {key: species_dir / f"{prefix}{suf}" for key, suf in suffixes.items()}
    n_complete = sum(1 for k in REQUIRED_FOR_COMPLETENESS if _nonempty(paths.get(k)))
    fasta_size = paths["ctg_fasta"].stat().st_size if _nonempty(paths.get("ctg_fasta")) else 0
    gfa_path = paths.get("gfa")
    mtime = gfa_path.stat().st_mtime if gfa_path and gfa_path.exists() else 0
    return (n_complete, fasta_size, mtime)


def resolve_species(species_dir: Path, organelle: str, logger=None) -> SpeciesFiles:
    species = species_dir.name
    suffixes = ORGANELLE_SUFFIX[organelle]
    prefixes = _candidate_prefixes(species_dir, suffixes["gfa"])

    if not prefixes:
        return SpeciesFiles(species=species, organelle=organelle, status="missing_gfa", n_candidates=0)

    if len(prefixes) > 1 and logger is not None:
        logger(f"{species} ({organelle}): {len(prefixes)} candidate runs found: {prefixes}")

    best_prefix = max(prefixes, key=lambda pfx: _completeness_score(species_dir, pfx, suffixes))

    if len(prefixes) > 1 and logger is not None:
        logger(f"{species} ({organelle}): selected run prefix '{best_prefix}'")

    paths = {key: species_dir / f"{best_prefix}{suf}" for key, suf in suffixes.items()}
    resolved = {k: (str(v) if _nonempty(v) else None) for k, v in paths.items()}

    if resolved["gfa"] is None:
        status = "missing_gfa"
    elif resolved["ctg_fasta"] is None:
        status = "no_ctg_fasta"
    else:
        status = "ok"

    return SpeciesFiles(
        species=species,
        organelle=organelle,
        gfa=resolved["gfa"],
        ctg_bed=resolved["ctg_bed"],
        ctg_fasta=resolved["ctg_fasta"],
        bed=resolved["bed"],
        annot=resolved["annot"],
        run_prefix=best_prefix,
        n_candidates=len(prefixes),
        status=status,
    )


def discover_all(data_root: Path, organelle: str, species_filter: Optional[set] = None,
                  logger=None) -> list[SpeciesFiles]:
    results = []
    for species_dir in sorted(data_root.iterdir()):
        if not species_dir.is_dir():
            continue
        if species_filter is not None and species_dir.name not in species_filter:
            continue
        results.append(resolve_species(species_dir, organelle, logger=logger))
    return results


def load_species_list(path: Path) -> set:
    with open(path) as fh:
        return {line.strip() for line in fh if line.strip() and not line.startswith("#")}


def repo_data_root(root_dir: Path, organelle: str) -> Path:
    return root_dir / "data" / ("mito" if organelle == "mito" else "plastid")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organelle", required=True, choices=["mito", "pltd"])
    ap.add_argument("--root-dir", default=None, help="repo root (default: 3 levels up from this file)")
    ap.add_argument("--species-list", default=None)
    ap.add_argument("--print-table", action="store_true", help="emit TSV to stdout")
    args = ap.parse_args()

    root_dir = Path(args.root_dir) if args.root_dir else Path(__file__).resolve().parents[2]
    data_root = repo_data_root(root_dir, args.organelle)
    species_filter = load_species_list(Path(args.species_list)) if args.species_list else None

    def logger(msg):
        print(f"[warn] {msg}", file=sys.stderr)

    results = discover_all(data_root, args.organelle, species_filter=species_filter, logger=logger)

    if args.print_table:
        writer = csv.DictWriter(sys.stdout, fieldnames=[f.name for f in fields(SpeciesFiles)], delimiter="\t")
        writer.writeheader()
        for r in results:
            writer.writerow(r.as_row())

    n_ok = sum(1 for r in results if r.status == "ok")
    n_no_fasta = sum(1 for r in results if r.status == "no_ctg_fasta")
    n_missing = sum(1 for r in results if r.status == "missing_gfa")
    print(f"[info] {args.organelle}: {len(results)} species, ok={n_ok} no_ctg_fasta={n_no_fasta} missing_gfa={n_missing}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
