"""Run `nhmmscan` directly against oatkDB's own gene-family databases, on
oatk's assembled contigs - decoupled from oatk's own bundled annotation
step. See `README.md` for why (short version: this is oatk's own
targeted gene search, just run independently, and it surfaces materially
more of a trans-spliced gene's real exon structure than oatk's bundled
per-species calling does - see the `nad5` comparison there).

Databases: the oatkDB builds bundled with poa (viridiplantae_mito_v20250217,
viridiplantae_pltd_v20260928), prepared on first use, or your own via
$OATKDB_MITO_FAM / $OATKDB_PLTD_FAM - see common/oatkdb.py.

Several `nhmmscan` processes run concurrently (`jobs`), each still
multi-threaded via `cpu` - see targets.py.

Output: analysis/denovo_annotation/work/nhmmscan/<species>.<organelle>.tblout
        analysis/denovo_annotation/work/nhmmscan/<species>.<organelle>.filtered.tblout
        analysis/denovo_annotation/work/nhmmscan/<species>.<organelle>.gff

src/00_run_nhmmscan.py is the command-line wrapper (main() here).
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Iterable

from poa.pipeline.common import oatkdb
from poa.pipeline.common import species_discovery as sd
from poa.pipeline.common.errors import StageError
from poa.pipeline.common.qc_gate import qc_eligible_species
from poa.pipeline.denovo_annotation import targets
from poa.pipeline.denovo_annotation.nhmmscan_gff import filter_and_convert

DEFAULT_EVALUE = "1e-5"


def nhmmscan_exe() -> str:
    return shutil.which("nhmmscan") or "/software/team301/hmmer-3.4/src/nhmmscan"


def scan_one(
    species: str, organelle: str, ctg_fasta: str, fam: str, evalue: str, cpu: int | str, work_dir: Path, force: bool
) -> str:
    tblout = work_dir / f"{species}.{organelle}.tblout"
    filtered = work_dir / f"{species}.{organelle}.filtered.tblout"
    gff = work_dir / f"{species}.{organelle}.gff"
    if not force and gff.exists() and gff.stat().st_size > 0:
        return "skipped"

    proc = subprocess.run(
        [nhmmscan_exe(), "--cpu", str(cpu), "--tblout", str(tblout), fam, ctg_fasta],
        capture_output=True,
        text=True,
        timeout=1800,
    )
    if proc.returncode != 0:
        print(f"[warn] nhmmscan failed for {species} ({organelle}): {proc.stderr.strip()[:300]}", file=sys.stderr)
        return "failed"

    try:
        filtered_text, gff_text = filter_and_convert(tblout, Path(fam), evalue)
    except (ValueError, KeyError, IndexError) as e:
        print(f"[warn] tblout -> GFF failed for {species} ({organelle}): {e!r}", file=sys.stderr)
        return "failed"
    filtered.write_text(filtered_text)
    gff.write_text(gff_text)
    return "scanned"


def run_nhmmscan(
    root: Path,
    organelles: Iterable[str] = targets.ORGANELLES,
    *,
    species: set[str] | None = None,
    evalue: str = DEFAULT_EVALUE,
    cpu: int = 2,
    jobs: int = 1,
    qc_status: Iterable[str] = ("pass",),
    force: bool = False,
) -> dict[str, dict[str, int]]:
    """nhmmscan every assembly under data root `root`. Returns
    {organelle: {"scanned"|"skipped"|"failed": count}}; a species that
    fails is warned about and counted, not raised."""
    root = Path(root)
    organelles = targets.check_organelles(organelles)
    qc_eligible = qc_eligible_species(root, qc_status)
    work_dir = root / "analysis" / "denovo_annotation" / "work" / "nhmmscan"
    work_dir.mkdir(parents=True, exist_ok=True)

    out = {}
    for organelle in organelles:
        try:
            fam = oatkdb.fam_path(organelle)
        except (OSError, subprocess.CalledProcessError) as e:
            raise StageError(f"{organelle}: can't prepare the oatkDB database: {e}") from e
        if not fam.exists():
            print(f"[warn] {organelle}: no oatkDB database at {fam} - skipping", file=sys.stderr)
            continue
        tasks = targets.assemblies(root, organelle, species, qc_eligible)
        out[organelle] = counts = targets.scan_all(
            tasks,
            lambda r: scan_one(r.species, organelle, r.ctg_fasta, str(fam), evalue, cpu, work_dir, force),
            jobs,
        )
        targets.report("run_nhmmscan", organelle, counts, work_dir)
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--organelle", default="both", choices=["mito", "pltd", "both"])
    ap.add_argument("--species-list", default=None)
    ap.add_argument("--evalue", default=DEFAULT_EVALUE)
    ap.add_argument("--cpu", type=int, default=2, help="threads per nhmmscan process (default 2)")
    ap.add_argument("--jobs", type=int, default=1, help="number of nhmmscan processes to run concurrently (default 1)")
    ap.add_argument(
        "--qc-status",
        default="pass",
        help="comma-separated subset of pass,flag,fail (default: pass) - "
        "same gate annotation/src/01_gene_table.py applies, skipped only if qc_summary.tsv is missing",
    )
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args(argv)
    root = Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT") or os.getcwd())
    try:
        run_nhmmscan(
            root,
            list(targets.ORGANELLES) if args.organelle == "both" else [args.organelle],
            species=sd.load_species_list(Path(args.species_list)) if args.species_list else None,
            evalue=args.evalue,
            cpu=args.cpu,
            jobs=args.jobs,
            qc_status=args.qc_status.split(","),
            force=args.force,
        )
    except StageError as e:
        print(f"[err] {e}", file=sys.stderr)
        return 1
    return 0
