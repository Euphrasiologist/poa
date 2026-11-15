"""Which assemblies the de novo scans run on (the QC gate itself is
common/qc_gate.py), and running a scan over them.

Shared by nhmmscan, trnascan and rrna: several scanner processes
concurrently (`jobs`), each still multi-threaded through its own option -
these scan steps had no parallelism at all until full-dataset-scale runs
made that a real problem (orf_scan's README documents the identical issue
for its own heavy step).
"""

from __future__ import annotations

import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Callable, Iterable

from poa.pipeline.common import species_discovery as sd

log = logging.getLogger(__name__)

ORGANELLES = ("mito", "pltd")


def check_organelles(organelles: Iterable[str]) -> list[str]:
    organelles = list(organelles)
    if unknown := set(organelles) - set(ORGANELLES):
        raise ValueError(f"unknown organelle(s) {sorted(unknown)}; expected {ORGANELLES}")
    return organelles


def assemblies(root: Path, organelle: str, species: set[str] | None = None, qc_eligible=None) -> list:
    """The resolved assemblies (species_discovery.SpeciesFiles) to scan."""
    resolved = sd.discover_all(sd.repo_data_root(Path(root), organelle), organelle, species_filter=species)
    return [r for r in resolved if r.status == "ok" and (qc_eligible is None or (r.species, organelle) in qc_eligible)]


def scan_all(tasks: list, scan_one: Callable[[object], str], jobs: int) -> dict[str, int]:
    """scan_one(task) -> "scanned" | "skipped" | "failed", over `jobs` threads."""
    counts = {"scanned": 0, "skipped": 0, "failed": 0}
    with ThreadPoolExecutor(max_workers=max(1, jobs)) as pool:
        for fut in as_completed([pool.submit(scan_one, t) for t in tasks]):
            counts[fut.result()] += 1
    return counts


def report(label: str, organelle: str, counts: dict[str, int], work_dir: Path) -> None:
    log.info(
        f"{label}: {organelle}: scanned={counts['scanned']} skipped={counts['skipped']} "
        f"failed={counts['failed']} -> {work_dir}"
    )
