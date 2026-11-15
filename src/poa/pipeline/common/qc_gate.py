"""The QC gate: which assemblies a dataset-scale stage runs on.

qc_basic_stats 01-05 judge each assembly against the rest of the dataset
(qc_summary.tsv); the per-assembly stages skip the ones whose status isn't
in the caller's list (default: pass only). With no qc_summary.tsv - `poa
run`, or a dataset before QC - nothing is filtered.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path
from typing import Iterable


def qc_eligible_species(root: Path, qc_status: Iterable[str]) -> set[tuple[str, str]] | None:
    """(species, organelle) pairs whose qc_summary.tsv status is in
    `qc_status` - the same gate `annotation/src/01_gene_table.py` applies
    by default. Confirmed necessary here too: a real "fail"-status assembly
    (`Stellaria_graminea`, fragmented/non-circular/<50% core genes)
    otherwise still gets the full ~4min tRNAscan-SE treatment for nothing,
    which adds up at full-dataset scale. None (meaning "don't filter") if
    the QC table doesn't exist yet, rather than silently filtering
    everything out."""
    qc_path = Path(root) / "analysis" / "qc_basic_stats" / "results" / "qc_summary.tsv"
    if not qc_path.exists():
        print(f"[warn] missing {qc_path} - not filtering by QC status", file=sys.stderr)
        return None
    qc_status = set(qc_status)
    with open(qc_path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return {(r["species"], r["organelle"]) for r in rows if r["status"] in qc_status}
