#!/usr/bin/env python3
"""Run `nhmmscan` directly against oatkDB's own gene-family databases, on
oatk's assembled contigs - decoupled from oatk's own bundled annotation
step. See `../README.md` for why (short version: this is oatk's own
targeted gene search, just run independently, and it surfaces materially
more of a trans-spliced gene's real exon structure than oatk's bundled
per-species calling does - see the `nad5` comparison there).

Mito: `/software/team301/OatkDB/viridiplantae_mito_v20250217.fam` (already
built, hmmpress'd). Plastid: `viridiplantae_pltd_v20260928.fam`, built via
`oatkdb` on the `long` queue - see `../README.md`'s plastid-database
section for the `nquire --http1.0` bug that blocked every earlier attempt
and how it was fixed.

Same `--jobs`/`--cpu` pattern as `orf_scan/src/02_scan_pfam.py` (several
`nhmmscan` processes concurrently, each still multi-threaded via `--cpu`) -
this and the other two `denovo_annotation` scan scripts had no
parallelism at all until full-dataset-scale runs made that a real
problem (orf_scan's README already documents the identical issue for its
own heavy step).

Output: analysis/denovo_annotation/work/nhmmscan/<species>.<organelle>.tblout
        analysis/denovo_annotation/work/nhmmscan/<species>.<organelle>.filtered.tblout
        analysis/denovo_annotation/work/nhmmscan/<species>.<organelle>.gff
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

NHMMSCAN = shutil.which("nhmmscan") or "/software/team301/hmmer-3.4/src/nhmmscan"
FILTER_TBLOUT = shutil.which("filter_tblout") or str(Path.home() / ".cargo" / "bin" / "filter_tblout")
HMM_TO_GFF = shutil.which("hmm_to_gff") or str(Path.home() / ".cargo" / "bin" / "hmm_to_gff")

FAM_PATHS = {
    "mito": os.environ.get("OATKDB_MITO_FAM", "/software/team301/OatkDB/viridiplantae_mito_v20250217.fam"),
    "pltd": os.environ.get("OATKDB_PLTD_FAM", "/software/team301/OatkDB/viridiplantae_pltd_v20260928.fam"),
}

DEFAULT_EVALUE = "1e-5"


def qc_eligible_species(qc_status: list[str]) -> set[tuple[str, str]] | None:
    """Same gate `annotation/src/01_gene_table.py` applies by default -
    confirmed necessary here too: a real "fail"-status assembly
    (`Stellaria_graminea`, fragmented/non-circular/<50% core genes)
    otherwise still gets the full ~4min tRNAscan-SE treatment for nothing,
    which adds up at full-dataset scale across however many species fail
    QC. Returns None (meaning "don't filter") if the QC table doesn't
    exist yet, rather than silently filtering everything out."""
    qc_path = ANALYSIS_DIR / "qc_basic_stats" / "results" / "qc_summary.tsv"
    if not qc_path.exists():
        print(f"[warn] missing {qc_path} - not filtering by QC status", file=sys.stderr)
        return None
    import csv
    with open(qc_path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return {(r["species"], r["organelle"]) for r in rows if r["status"] in qc_status}


def scan_one(species: str, organelle: str, ctg_fasta: str, fam: str, evalue: str, cpu: str,
             work_dir: Path, force: bool) -> str:
    tblout = work_dir / f"{species}.{organelle}.tblout"
    filtered = work_dir / f"{species}.{organelle}.filtered.tblout"
    gff = work_dir / f"{species}.{organelle}.gff"
    if not force and gff.exists() and gff.stat().st_size > 0:
        return "skipped"

    proc = subprocess.run(
        [NHMMSCAN, "--cpu", cpu, "--tblout", str(tblout), fam, ctg_fasta],
        capture_output=True, text=True, timeout=1800,
    )
    if proc.returncode != 0:
        print(f"[warn] nhmmscan failed for {species} ({organelle}): {proc.stderr.strip()[:300]}", file=sys.stderr)
        return "failed"

    filt = subprocess.run([FILTER_TBLOUT, str(tblout), evalue], capture_output=True, text=True, timeout=120)
    if filt.returncode != 0:
        print(f"[warn] filter_tblout failed for {species} ({organelle}): {filt.stderr.strip()[:300]}", file=sys.stderr)
        return "failed"
    filtered.write_text(filt.stdout)

    togff = subprocess.run([HMM_TO_GFF, str(filtered), "nhmmscan", organelle], capture_output=True, text=True, timeout=120)
    if togff.returncode != 0:
        print(f"[warn] hmm_to_gff failed for {species} ({organelle}): {togff.stderr.strip()[:300]}", file=sys.stderr)
        return "failed"
    gff.write_text(togff.stdout)
    return "scanned"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organelle", default="both", choices=["mito", "pltd", "both"])
    ap.add_argument("--species-list", default=None)
    ap.add_argument("--evalue", default=DEFAULT_EVALUE)
    ap.add_argument("--cpu", default="2", help="threads per nhmmscan process (default 2)")
    ap.add_argument("--jobs", type=int, default=1, help="number of nhmmscan processes to run concurrently (default 1)")
    ap.add_argument("--qc-status", default="pass", help="comma-separated subset of pass,flag,fail (default: pass) - "
                     "same gate annotation/src/01_gene_table.py applies, skipped only if qc_summary.tsv is missing")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    organelles = ["mito", "pltd"] if args.organelle == "both" else [args.organelle]
    species_filter = sd.load_species_list(Path(args.species_list)) if args.species_list else None
    qc_eligible = qc_eligible_species(args.qc_status.split(","))
    work_dir = ANALYSIS_DIR / "denovo_annotation" / "work" / "nhmmscan"
    work_dir.mkdir(parents=True, exist_ok=True)

    for organelle in organelles:
        fam = Path(FAM_PATHS[organelle])
        if not fam.exists():
            print(f"[warn] {organelle}: no oatkDB database at {fam} - skipping "
                  f"(see ../README.md's plastid-database caveat)", file=sys.stderr)
            continue

        data_root = sd.repo_data_root(ROOT_DIR, organelle)
        resolved = sd.discover_all(data_root, organelle, species_filter=species_filter)
        tasks = [r for r in resolved if r.status == "ok"
                 and (qc_eligible is None or (r.species, organelle) in qc_eligible)]

        counts = {"scanned": 0, "skipped": 0, "failed": 0}
        with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
            futures = {
                pool.submit(scan_one, r.species, organelle, r.ctg_fasta, str(fam), args.evalue, args.cpu,
                            work_dir, args.force): r.species
                for r in tasks
            }
            for fut in as_completed(futures):
                counts[fut.result()] += 1
        print(f"[info] run_nhmmscan: {organelle}: scanned={counts['scanned']} skipped={counts['skipped']} "
              f"failed={counts['failed']} -> {work_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
