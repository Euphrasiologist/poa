#!/usr/bin/env python3
"""Per-contig length/circularity/GC%, parsed from .ctg.fasta headers + seqkit.

oatk .ctg.fasta headers look like:
  >ctg000001c length=769913 wlength=... nv=12 circular=true path=u75+,...
`length=`/`nv=`/`circular=` are read straight from the header (no sequence
scan needed for those); GC% comes from `seqkit fx2tab -g` (real sequence
composition, not the graph-level value gfatk reports, so the two can be
cross-checked).

A handful of species (59, all plastid, confirmed by direct inspection) have
a .pltd.gfa + .pltd.ctg.bed but no .pltd.ctg.fasta. For those, fall back to
`gfatools gfa2fa` on the GFA (unitig sequences, NOT the resolved circularised
contig path) so basic stats/QC still get *something* to measure size/GC
from, clearly labelled contig_source=gfa_fallback rather than silently
treated as equivalent to a real resolved contig.

Output: analysis/qc_basic_stats/results/contig_stats.tsv
Columns: species organelle contig_id length nv circular gc_pct contig_source
"""
from __future__ import annotations

import os
import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

# CODE_DIR is this checkout's analysis/ (shared modules, bundled reference
# data); data/ and every module's results/ and work/ live under the data
# root - PLANT_ORGANELLE_DATA_ROOT, or this checkout if unset.
CODE_DIR = Path(__file__).resolve().parents[2]
ROOT_DIR = Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT", str(CODE_DIR.parent)))
ANALYSIS_DIR = ROOT_DIR / "analysis"
sys.path.insert(0, str(CODE_DIR / "common"))

import io_utils  # noqa: E402
import species_discovery as sd  # noqa: E402

COLUMNS = ["species", "organelle", "contig_id", "length", "nv", "circular", "gc_pct", "contig_source"]
KEY_COLS = ["species", "organelle"]

SEQKIT = shutil.which("seqkit") or "/software/team301/seqkit"
GFATOOLS = shutil.which("gfatools") or "/software/team301/gfatools/gfatools"

HEADER_RE = re.compile(r"^>(\S+).*?length=(\d+).*?nv=(\d+).*?circular=(true|false)")

WORK_DIR = ANALYSIS_DIR / "qc_basic_stats" / "work"


def parse_fasta_headers(fasta_path: Path) -> dict:
    """contig_id -> {length, nv, circular} straight from header attributes, if present."""
    out = {}
    with open(fasta_path) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            m = HEADER_RE.match(line)
            if m:
                cid, length, nv, circular = m.groups()
                out[cid] = {"length": int(length), "nv": int(nv), "circular": circular}
            else:
                cid = line[1:].split()[0]
                out[cid] = {"length": None, "nv": None, "circular": None}
    return out


def seqkit_gc(fasta_path: Path) -> dict:
    """contig_id -> gc_pct via seqkit fx2tab."""
    proc = subprocess.run([SEQKIT, "fx2tab", "-n", "-i", "-l", "-g", str(fasta_path)],
                           capture_output=True, text=True, timeout=120)
    out = {}
    if proc.returncode != 0:
        return out
    for line in proc.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) >= 3:
            out[parts[0]] = float(parts[2])
    return out


def gfa_fallback_fasta(gfa_path: Path, species: str, organelle: str) -> Path | None:
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    out_path = WORK_DIR / f"{species}.{organelle}.gfa_fallback.fasta"
    try:
        with open(out_path, "w") as fh:
            proc = subprocess.run([GFATOOLS, "gfa2fa", str(gfa_path)], stdout=fh, stderr=subprocess.PIPE,
                                   text=True, timeout=120)
        if proc.returncode != 0 or out_path.stat().st_size == 0:
            return None
        return out_path
    except Exception:  # noqa: BLE001
        return None


def rows_for_species(species: str, organelle: str, ctg_fasta: str | None, gfa: str | None) -> list[dict]:
    base = {"species": species, "organelle": organelle}
    fasta_path = None
    contig_source = "no_assembly"

    if ctg_fasta:
        fasta_path = Path(ctg_fasta)
        contig_source = "ctg_fasta"
    elif gfa:
        fasta_path = gfa_fallback_fasta(Path(gfa), species, organelle)
        contig_source = "gfa_fallback" if fasta_path else "no_assembly"

    if fasta_path is None:
        return [{**base, "contig_id": None, "length": None, "nv": None, "circular": None,
                 "gc_pct": None, "contig_source": "no_assembly"}]

    headers = parse_fasta_headers(fasta_path)
    gc = seqkit_gc(fasta_path)
    if not headers:
        return [{**base, "contig_id": None, "length": None, "nv": None, "circular": None,
                 "gc_pct": None, "contig_source": "no_assembly"}]

    rows = []
    for cid, attrs in headers.items():
        rows.append({**base, "contig_id": cid, "length": attrs["length"], "nv": attrs["nv"],
                     "circular": attrs["circular"], "gc_pct": gc.get(cid), "contig_source": contig_source})
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organelle", default="both", choices=["mito", "pltd", "both"])
    ap.add_argument("--species-list", default=None)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    out_path = ANALYSIS_DIR / "qc_basic_stats" / "results" / "contig_stats.tsv"
    existing_df, existing_keys = io_utils.load_existing(out_path, KEY_COLS)

    organelles = ["mito", "pltd"] if args.organelle == "both" else [args.organelle]
    species_filter = sd.load_species_list(Path(args.species_list)) if args.species_list else None

    new_rows = []
    n_computed = n_skipped = 0
    for organelle in organelles:
        data_root = sd.repo_data_root(ROOT_DIR, organelle)
        resolved = sd.discover_all(data_root, organelle, species_filter=species_filter)
        for r in resolved:
            key = (r.species, organelle)
            input_path = Path(r.ctg_fasta or r.gfa) if (r.ctg_fasta or r.gfa) else data_root / r.species
            if io_utils.needs_recompute(key, existing_keys, out_path, input_path, args.force):
                new_rows.extend(rows_for_species(r.species, organelle, r.ctg_fasta, r.gfa))
                n_computed += 1
            else:
                n_skipped += 1

    merged = io_utils.merge_rows(existing_df, new_rows, KEY_COLS, COLUMNS)
    io_utils.atomic_write_tsv(out_path, merged.to_dict("records"), COLUMNS)
    print(f"[info] contig_stats: computed={n_computed} skipped={n_skipped} total_rows={len(merged)} -> {out_path}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
