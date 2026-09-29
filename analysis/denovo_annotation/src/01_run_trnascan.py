#!/usr/bin/env python3
"""Run `tRNAscan-SE` (organelle mode) on oatk's assembled contigs - the
tRNA half of a de novo, oatk-independent annotation pass. Same tool/flags
already proven in `mito_structural_variation/annotation/src/5_get_rna_genes.bash`.

Same `--jobs`/`--thread` pattern as `00_run_nhmmscan.py`/
`orf_scan/src/02_scan_pfam.py` - see that script's docstring for why.

Output: analysis/denovo_annotation/work/trnascan/<species>.<organelle>.gff
"""
from __future__ import annotations

import os
import argparse
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

# CODE_DIR is this checkout's analysis/ (shared modules, bundled reference
# data); data/ and every module's results/ and work/ live under the data
# root - PLANT_ORGANELLE_DATA_ROOT, or this checkout if unset.
CODE_DIR = Path(__file__).resolve().parents[2]
ROOT_DIR = Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT", str(CODE_DIR.parent)))
ANALYSIS_DIR = ROOT_DIR / "analysis"
sys.path.insert(0, str(CODE_DIR / "common"))

import species_discovery as sd  # noqa: E402

TRNASCAN = shutil.which("tRNAscan-SE") or "/software/team301/tRNAscan-SE/tRNAscan-SE"
TRNASCAN_CONF = os.environ.get("TRNASCAN_CONF", "/software/team301/tRNAscan-SE/tRNAscan-SE.conf")


def qc_eligible_species(qc_status: list[str]) -> set[tuple[str, str]] | None:
    """Same gate `annotation/src/01_gene_table.py` applies by default -
    confirmed necessary here too: a real "fail"-status assembly
    (`Stellaria_graminea`, fragmented/non-circular/<50% core genes) still
    got the full ~4min tRNAscan-SE treatment for nothing before this,
    which adds up at full-dataset scale. Returns None (meaning "don't
    filter") if the QC table doesn't exist yet."""
    qc_path = ANALYSIS_DIR / "qc_basic_stats" / "results" / "qc_summary.tsv"
    if not qc_path.exists():
        print(f"[warn] missing {qc_path} - not filtering by QC status", file=sys.stderr)
        return None
    import csv
    with open(qc_path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return {(r["species"], r["organelle"]) for r in rows if r["status"] in qc_status}


def scan_one(species: str, organelle: str, ctg_fasta: str, threads: str, work_dir: Path, force: bool) -> str:
    out_gff = work_dir / f"{species}.{organelle}.gff"
    if not force and out_gff.exists() and out_gff.stat().st_size > 0:
        return "skipped"
    proc = subprocess.run(
        # -Q/--forceow: tRNAscan-SE itself prompts interactively (and hangs
        # forever under subprocess.run, which never answers) if its own
        # output file already exists - confirmed directly (a `--force`
        # re-run on already-scanned species hung two workers indefinitely
        # and drove the job over its memory limit). Always pass -Q; this
        # script's own `force`/skip-if-exists check above is what actually
        # controls re-scanning, not tRNAscan-SE's own overwrite prompt.
        [TRNASCAN, *(["-c", TRNASCAN_CONF] if Path(TRNASCAN_CONF).exists() else []), "-O", "-I", "-Q", "--thread", threads, ctg_fasta, "-j", str(out_gff)],
        capture_output=True, text=True, timeout=1800,
    )
    if proc.returncode != 0 or not out_gff.exists():
        print(f"[warn] tRNAscan-SE failed for {species} ({organelle}): {proc.stderr.strip()[:300]}", file=sys.stderr)
        return "failed"
    return "scanned"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organelle", default="both", choices=["mito", "pltd", "both"])
    ap.add_argument("--species-list", default=None)
    ap.add_argument("--thread", default="2", help="threads per tRNAscan-SE process (default 2)")
    ap.add_argument("--jobs", type=int, default=1, help="number of tRNAscan-SE processes to run concurrently (default 1)")
    ap.add_argument("--qc-status", default="pass", help="comma-separated subset of pass,flag,fail (default: pass)")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    organelles = ["mito", "pltd"] if args.organelle == "both" else [args.organelle]
    species_filter = sd.load_species_list(Path(args.species_list)) if args.species_list else None
    qc_eligible = qc_eligible_species(args.qc_status.split(","))
    work_dir = ANALYSIS_DIR / "denovo_annotation" / "work" / "trnascan"
    work_dir.mkdir(parents=True, exist_ok=True)

    for organelle in organelles:
        data_root = sd.repo_data_root(ROOT_DIR, organelle)
        resolved = sd.discover_all(data_root, organelle, species_filter=species_filter)
        tasks = [r for r in resolved if r.status == "ok"
                 and (qc_eligible is None or (r.species, organelle) in qc_eligible)]

        counts = {"scanned": 0, "skipped": 0, "failed": 0}
        with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
            futures = {
                pool.submit(scan_one, r.species, organelle, r.ctg_fasta, args.thread, work_dir, args.force): r.species
                for r in tasks
            }
            for fut in as_completed(futures):
                counts[fut.result()] += 1
        print(f"[info] run_trnascan: {organelle}: scanned={counts['scanned']} skipped={counts['skipped']} "
              f"failed={counts['failed']} -> {work_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
