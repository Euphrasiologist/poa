"""Run `tRNAscan-SE` (organelle mode) on oatk's assembled contigs - the
tRNA half of a de novo, oatk-independent annotation pass. Same tool/flags
already proven in `mito_structural_variation/annotation/src/5_get_rna_genes.bash`.

Several tRNAscan-SE processes run concurrently (`jobs`), each with
`threads` threads - see targets.py.

Output: analysis/denovo_annotation/work/trnascan/<species>.<organelle>.gff

src/01_run_trnascan.py is the command-line wrapper (main() here).
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Iterable

from poa.pipeline.common import species_discovery as sd
from poa.pipeline.common.errors import StageError
from poa.pipeline.common.qc_gate import qc_eligible_species
from poa.pipeline.denovo_annotation import targets


def trnascan_exe() -> str:
    return shutil.which("tRNAscan-SE") or "/software/team301/tRNAscan-SE/tRNAscan-SE"


def trnascan_conf() -> str:
    return os.environ.get("TRNASCAN_CONF", "/software/team301/tRNAscan-SE/tRNAscan-SE.conf")


def scan_one(species: str, organelle: str, ctg_fasta: str, threads: int | str, work_dir: Path, force: bool) -> str:
    out_gff = work_dir / f"{species}.{organelle}.gff"
    if not force and out_gff.exists() and out_gff.stat().st_size > 0:
        return "skipped"
    conf = trnascan_conf()
    proc = subprocess.run(
        # -Q/--forceow: tRNAscan-SE itself prompts interactively (and hangs
        # forever under subprocess.run, which never answers) if its own
        # output file already exists - confirmed directly (a `--force`
        # re-run on already-scanned species hung two workers indefinitely
        # and drove the job over its memory limit). Always pass -Q; this
        # module's own `force`/skip-if-exists check above is what actually
        # controls re-scanning, not tRNAscan-SE's own overwrite prompt.
        [
            trnascan_exe(),
            *(["-c", conf] if Path(conf).exists() else []),
            "-O",
            "-I",
            "-Q",
            "--thread",
            str(threads),
            ctg_fasta,
            "-j",
            str(out_gff),
        ],
        capture_output=True,
        text=True,
        timeout=1800,
    )
    if proc.returncode != 0 or not out_gff.exists():
        print(f"[warn] tRNAscan-SE failed for {species} ({organelle}): {proc.stderr.strip()[:300]}", file=sys.stderr)
        return "failed"
    return "scanned"


def run_trnascan(
    root: Path,
    organelles: Iterable[str] = targets.ORGANELLES,
    *,
    species: set[str] | None = None,
    threads: int = 2,
    jobs: int = 1,
    qc_status: Iterable[str] = ("pass",),
    force: bool = False,
) -> dict[str, dict[str, int]]:
    """tRNAscan-SE every assembly under data root `root`. Returns
    {organelle: {"scanned"|"skipped"|"failed": count}}."""
    root = Path(root)
    organelles = targets.check_organelles(organelles)
    qc_eligible = qc_eligible_species(root, qc_status)
    work_dir = root / "analysis" / "denovo_annotation" / "work" / "trnascan"
    work_dir.mkdir(parents=True, exist_ok=True)

    out = {}
    for organelle in organelles:
        tasks = targets.assemblies(root, organelle, species, qc_eligible)
        out[organelle] = counts = targets.scan_all(
            tasks, lambda r: scan_one(r.species, organelle, r.ctg_fasta, threads, work_dir, force), jobs
        )
        targets.report("run_trnascan", organelle, counts, work_dir)
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--organelle", default="both", choices=["mito", "pltd", "both"])
    ap.add_argument("--species-list", default=None)
    ap.add_argument("--thread", type=int, default=2, help="threads per tRNAscan-SE process (default 2)")
    ap.add_argument(
        "--jobs", type=int, default=1, help="number of tRNAscan-SE processes to run concurrently (default 1)"
    )
    ap.add_argument("--qc-status", default="pass", help="comma-separated subset of pass,flag,fail (default: pass)")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args(argv)
    root = Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT") or os.getcwd())
    try:
        run_trnascan(
            root,
            list(targets.ORGANELLES) if args.organelle == "both" else [args.organelle],
            species=sd.load_species_list(Path(args.species_list)) if args.species_list else None,
            threads=args.thread,
            jobs=args.jobs,
            qc_status=args.qc_status.split(","),
            force=args.force,
        )
    except StageError as e:
        print(f"[err] {e}", file=sys.stderr)
        return 1
    return 0
