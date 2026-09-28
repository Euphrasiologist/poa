#!/usr/bin/env python3
"""Run `barrnap` (plant kingdom model) on oatk's assembled contigs - the
rRNA half of a de novo, oatk-independent annotation pass. Same tool/flags
already proven in `mito_structural_variation/annotation/src/6_get_rrna_genes.bash`.

Same `--jobs`/`--threads` pattern as `00_run_nhmmscan.py`/
`orf_scan/src/02_scan_pfam.py` - see that script's docstring for why (this
step alone is fast even single-threaded, but kept consistent with the
other two `denovo_annotation` scan scripts rather than being the odd one
out at full-dataset scale).

Output: analysis/denovo_annotation/work/barrnap/<species>.<organelle>.gff
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import os
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

# Defaults to the old same-tree layout (code and data siblings under one
# repo root); set PLANT_ORGANELLE_DATA_ROOT to point at a separate data
# checkout instead (e.g. when this code is installed from its own repo).
ROOT_DIR = Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT", str(Path(__file__).resolve().parents[3])))
ANALYSIS_DIR = ROOT_DIR / "analysis"
sys.path.insert(0, str(ANALYSIS_DIR / "common"))

import species_discovery as sd  # noqa: E402

BARRNAP = shutil.which("barrnap") or "/software/team301/barrnap/bin/barrnap"


def qc_eligible_species(qc_status: list[str]) -> set[tuple[str, str]] | None:
    """Same gate `annotation/src/01_gene_table.py` applies by default - see
    `00_run_nhmmscan.py`'s identical helper for why (a real "fail"-status
    assembly otherwise still gets scanned for nothing). Returns None
    (meaning "don't filter") if the QC table doesn't exist yet."""
    qc_path = ANALYSIS_DIR / "qc_basic_stats" / "results" / "qc_summary.tsv"
    if not qc_path.exists():
        print(f"[warn] missing {qc_path} - not filtering by QC status", file=sys.stderr)
        return None
    import csv
    with open(qc_path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return {(r["species"], r["organelle"]) for r in rows if r["status"] in qc_status}


def scan_one(species: str, organelle: str, ctg_fasta: str, evalue: str, threads: str,
             work_dir: Path, force: bool) -> str:
    out_gff = work_dir / f"{species}.{organelle}.gff"
    if not force and out_gff.exists() and out_gff.stat().st_size > 0:
        return "skipped"
    proc = subprocess.run(
        [BARRNAP, "--evalue", evalue, "--kingdom", "plant", "--threads", threads, ctg_fasta],
        capture_output=True, text=True, timeout=600,
    )
    if proc.returncode != 0:
        print(f"[warn] barrnap failed for {species} ({organelle}): {proc.stderr.strip()[:300]}", file=sys.stderr)
        return "failed"
    out_gff.write_text(proc.stdout)
    return "scanned"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organelle", default="both", choices=["mito", "pltd", "both"])
    ap.add_argument("--species-list", default=None)
    ap.add_argument("--evalue", default="1e-5")
    ap.add_argument("--threads", default="1", help="threads per barrnap process (default 1)")
    ap.add_argument("--jobs", type=int, default=1, help="number of barrnap processes to run concurrently (default 1)")
    ap.add_argument("--qc-status", default="pass", help="comma-separated subset of pass,flag,fail (default: pass)")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    organelles = ["mito", "pltd"] if args.organelle == "both" else [args.organelle]
    species_filter = sd.load_species_list(Path(args.species_list)) if args.species_list else None
    qc_eligible = qc_eligible_species(args.qc_status.split(","))
    work_dir = ANALYSIS_DIR / "denovo_annotation" / "work" / "barrnap"
    work_dir.mkdir(parents=True, exist_ok=True)

    for organelle in organelles:
        data_root = sd.repo_data_root(ROOT_DIR, organelle)
        resolved = sd.discover_all(data_root, organelle, species_filter=species_filter)
        tasks = [r for r in resolved if r.status == "ok"
                 and (qc_eligible is None or (r.species, organelle) in qc_eligible)]

        counts = {"scanned": 0, "skipped": 0, "failed": 0}
        with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
            futures = {
                pool.submit(scan_one, r.species, organelle, r.ctg_fasta, args.evalue, args.threads,
                            work_dir, args.force): r.species
                for r in tasks
            }
            for fut in as_completed(futures):
                counts[fut.result()] += 1
        print(f"[info] run_barrnap: {organelle}: scanned={counts['scanned']} skipped={counts['skipped']} "
              f"failed={counts['failed']} -> {work_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
