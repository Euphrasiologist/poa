"""rRNA calling on the assembled contigs - the rRNA half of the de novo
annotation pass. Originally `barrnap --kingdom plant` (a database added to
one cluster's barrnap); now `rrna_nhmmer.py`, a byte-identical Python port
that needs only nhmmer and poa's bundled plant.hmm.

Several nhmmer processes run concurrently (`jobs`), each with `threads`
threads - see targets.py (this step alone is fast even single-threaded,
but kept consistent with the other two scans at full-dataset scale).

Output: analysis/denovo_annotation/work/barrnap/<species>.<organelle>.gff

src/02_run_barrnap.py is the command-line wrapper (main() here).
"""

from __future__ import annotations

import argparse
import logging
import os
import subprocess
from pathlib import Path
from typing import Iterable

from poa.pipeline.common import species_discovery as sd
from poa.pipeline.common.errors import StageError
from poa.pipeline.common.qc_gate import qc_eligible_species
from poa.pipeline.denovo_annotation import targets
from poa.pipeline.denovo_annotation.rrna_nhmmer import call_rrna
from poa.pipeline.common.log import setup_logging

log = logging.getLogger(__name__)

DEFAULT_EVALUE = "1e-5"


def scan_one(
    species: str, organelle: str, ctg_fasta: str, evalue: str, threads: int | str, work_dir: Path, force: bool
) -> str:
    out_gff = work_dir / f"{species}.{organelle}.gff"
    if not force and out_gff.exists() and out_gff.stat().st_size > 0:
        return "skipped"
    try:
        gff = call_rrna(ctg_fasta, evalue, threads)
    except (OSError, ValueError, subprocess.SubprocessError) as e:
        err = getattr(e, "stderr", None) or str(e)
        log.warning(f"rRNA calling failed for {species} ({organelle}): {err.strip()[:300]}")
        return "failed"
    out_gff.write_text(gff)
    return "scanned"


def run_rrna(
    root: Path,
    organelles: Iterable[str] = targets.ORGANELLES,
    *,
    species: set[str] | None = None,
    evalue: str = DEFAULT_EVALUE,
    threads: int = 1,
    jobs: int = 1,
    qc_status: Iterable[str] = ("pass",),
    force: bool = False,
) -> dict[str, dict[str, int]]:
    """Call rRNAs in every assembly under data root `root`. Returns
    {organelle: {"scanned"|"skipped"|"failed": count}}."""
    root = Path(root)
    organelles = targets.check_organelles(organelles)
    qc_eligible = qc_eligible_species(root, qc_status)
    work_dir = root / "analysis" / "denovo_annotation" / "work" / "barrnap"
    work_dir.mkdir(parents=True, exist_ok=True)

    out = {}
    for organelle in organelles:
        tasks = targets.assemblies(root, organelle, species, qc_eligible)
        out[organelle] = counts = targets.scan_all(
            tasks, lambda r: scan_one(r.species, organelle, r.ctg_fasta, evalue, threads, work_dir, force), jobs
        )
        targets.report("run_barrnap", organelle, counts, work_dir)
    return out


def main(argv: list[str] | None = None) -> int:
    setup_logging()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--organelle", default="both", choices=["mito", "pltd", "both"])
    ap.add_argument("--species-list", default=None)
    ap.add_argument("--evalue", default=DEFAULT_EVALUE)
    ap.add_argument("--threads", type=int, default=1, help="threads per nhmmer process (default 1)")
    ap.add_argument("--jobs", type=int, default=1, help="number of nhmmer processes to run concurrently (default 1)")
    ap.add_argument("--qc-status", default="pass", help="comma-separated subset of pass,flag,fail (default: pass)")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args(argv)
    root = Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT") or os.getcwd())
    try:
        run_rrna(
            root,
            list(targets.ORGANELLES) if args.organelle == "both" else [args.organelle],
            species=sd.load_species_list(Path(args.species_list)) if args.species_list else None,
            evalue=args.evalue,
            threads=args.threads,
            jobs=args.jobs,
            qc_status=args.qc_status.split(","),
            force=args.force,
        )
    except StageError as e:
        log.error(e)
        return 1
    return 0
