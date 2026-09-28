#!/usr/bin/env python3
"""Find ORFs on each species' resolved contigs and keep only the ones NOT
already explained by oatk's own core-gene annotation ("non-core" ORFs) -
candidates for transposable elements, mitovirus-derived sequences, or other
orphan content that a targeted core-gene HMM search wouldn't catch.

Uses NCBI's ORFfinder directly (protein output, outfmt 0), genetic code 1
(standard - correct for plant mito/plastid, unlike animal mitochondria
which use code 2), `-ml 300` (100 aa minimum - verified empirically that
the default 75nt cutoff finds ~4000 ORFs on a single ~770kb mitogenome,
almost all spurious background noise; 300nt/100aa cuts that to ~200 while
still comfortably catching TE/viral domains, which are typically hundreds
of aa), and `-n true` (ignore ORFs nested entirely within another).

An ORF is dropped as "core" if it overlaps any .ctg.bed gene interval by
more than 50% of the ORF's own length - verified empirically that
ORFfinder independently re-finds real annotated genes (NADH dehydrogenase,
COX3, ribosomal proteins, ...) as ORFs, so without this filter the scan
would mostly just rediscover known genes rather than surface novel content.

Output:
  analysis/orf_scan/work/<Species>.<organelle>.noncore_orfs.faa (scratch - protein FASTA, input to 02_scan_pfam.py)
  analysis/orf_scan/results/per_species/<Species>.<organelle>.orfs.tsv
    Columns: orf_id contig_id start end strand length_aa
"""
from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import os
import sys
from pathlib import Path

import pandas as pd

# Defaults to the old same-tree layout (code and data siblings under one
# repo root); set PLANT_ORGANELLE_DATA_ROOT to point at a separate data
# checkout instead (e.g. when this code is installed from its own repo).
ROOT_DIR = Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT", str(Path(__file__).resolve().parents[3])))
ANALYSIS_DIR = ROOT_DIR / "analysis"
sys.path.insert(0, str(ANALYSIS_DIR / "common"))

import io_utils  # noqa: E402
import species_discovery as sd  # noqa: E402

ORFFINDER = shutil.which("ORFfinder") or "/software/team301/ORFfinder"
MIN_LENGTH_NT = 300  # -ml: 100 aa minimum, see module docstring
OVERLAP_FRACTION_THRESHOLD = 0.5  # drop ORF if >=50% of its length overlaps a core gene call

HEADER_RE = re.compile(r"^>lcl\|(ORF\d+)_(.+):(\d+):(\d+)")
COLUMNS = ["orf_id", "contig_id", "start", "end", "strand", "length_aa"]


def run_orffinder(fasta_path: str) -> list[tuple[str, str, int, int, str, str]]:
    """Returns list of (orf_id, contig_id, start, end, strand, protein_seq), 1-based inclusive coords."""
    cmd = [ORFFINDER, "-in", fasta_path, "-g", "1", "-s", "0", "-n", "true",
           "-strand", "both", "-ml", str(MIN_LENGTH_NT), "-outfmt", "0"]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    if proc.returncode != 0:
        return []

    orfs = []
    header = None
    seq_lines: list[str] = []
    for line in proc.stdout.splitlines():
        if line.startswith(">"):
            if header is not None:
                orfs.append(header + ("".join(seq_lines),))
            m = HEADER_RE.match(line)
            if m:
                orf_id, contig_id, a, b = m.groups()
                a, b = int(a), int(b)
                if a <= b:
                    header = (orf_id, contig_id, a, b, "+")
                else:
                    header = (orf_id, contig_id, b, a, "-")
            else:
                header = None
            seq_lines = []
        else:
            seq_lines.append(line.strip())
    if header is not None:
        orfs.append(header + ("".join(seq_lines),))
    return orfs


def load_core_intervals(ctg_bed: str) -> dict:
    """contig_id -> list of (start, end), 0-based half-open as in .ctg.bed."""
    df = pd.read_csv(ctg_bed, sep="\t", comment="#", names=["contig_id", "start", "end", "gene", "score", "strand"])
    out: dict = {}
    for contig_id, group in df.groupby("contig_id"):
        out[contig_id] = list(zip(group["start"], group["end"]))
    return out


def overlap_fraction(orf_start: int, orf_end: int, intervals: list[tuple[int, int]]) -> float:
    """Max fraction of the (1-based inclusive) ORF interval covered by any single core-gene interval."""
    orf_len = orf_end - orf_start + 1
    if orf_len <= 0 or not intervals:
        return 0.0
    best = 0.0
    for cs, ce in intervals:  # 0-based half-open
        cs1, ce1 = cs + 1, ce  # convert to 1-based inclusive
        overlap = max(0, min(orf_end, ce1) - max(orf_start, cs1) + 1)
        best = max(best, overlap / orf_len)
    return best


def process_species(species: str, organelle: str, fasta_path: str, ctg_bed: str) -> tuple[list[dict], list[str]]:
    orfs = run_orffinder(fasta_path)
    core_intervals = load_core_intervals(ctg_bed)

    rows, fasta_records = [], []
    for orf_id, contig_id, start, end, strand, protein in orfs:
        frac = overlap_fraction(start, end, core_intervals.get(contig_id, []))
        if frac >= OVERLAP_FRACTION_THRESHOLD:
            continue
        rows.append({"orf_id": orf_id, "contig_id": contig_id, "start": start, "end": end,
                      "strand": strand, "length_aa": len(protein)})
        fasta_records.append(f">{species}|{organelle}|{orf_id}\n{protein}\n")
    return rows, fasta_records


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organelle", default="both", choices=["mito", "pltd", "both"])
    ap.add_argument("--species-list", default=None)
    ap.add_argument("--qc-status", default="pass")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    qc_path = ANALYSIS_DIR / "qc_basic_stats" / "results" / "qc_summary.tsv"
    qc = pd.read_csv(qc_path, sep="\t")
    qc_status = args.qc_status.split(",")
    eligible = set(zip(qc.loc[qc["status"].isin(qc_status), "species"], qc.loc[qc["status"].isin(qc_status), "organelle"]))

    organelles = ["mito", "pltd"] if args.organelle == "both" else [args.organelle]
    species_filter = sd.load_species_list(Path(args.species_list)) if args.species_list else None

    out_dir = ANALYSIS_DIR / "orf_scan" / "results" / "per_species"
    work_dir = ANALYSIS_DIR / "orf_scan" / "work"
    out_dir.mkdir(parents=True, exist_ok=True)
    work_dir.mkdir(parents=True, exist_ok=True)

    n_computed = n_skipped = n_no_orfs = 0
    for organelle in organelles:
        data_root = sd.repo_data_root(ROOT_DIR, organelle)
        resolved = sd.discover_all(data_root, organelle, species_filter=species_filter)
        for r in resolved:
            if (r.species, organelle) not in eligible or r.status != "ok" or not r.ctg_bed:
                continue
            out_tsv = out_dir / f"{r.species}.{organelle}.orfs.tsv"
            faa_path = work_dir / f"{r.species}.{organelle}.noncore_orfs.faa"
            if not args.force and out_tsv.exists() and Path(r.ctg_fasta).stat().st_mtime <= out_tsv.stat().st_mtime:
                n_skipped += 1
                continue

            rows, fasta_records = process_species(r.species, organelle, r.ctg_fasta, r.ctg_bed)
            io_utils.atomic_write_tsv(out_tsv, rows, COLUMNS)
            faa_path.write_text("".join(fasta_records))
            if rows:
                n_computed += 1
            else:
                n_no_orfs += 1

    print(f"[info] find_orfs: computed={n_computed} skipped={n_skipped} no_noncore_orfs={n_no_orfs} -> {out_dir}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
