#!/usr/bin/env python3
"""Run `gfatk stats` per species x organelle and record per-subgraph metrics.

A genome can legitimately have >1 subgraph (multipartite plant mitochondria
are real biology, not an assembly error), so this writes one row per
(species, organelle, subgraph_id), plus a single summary row (subgraph_id
NA) for species with a missing/empty/error graph. Aggregation to one row
per species (total length, worst dead-end count, any-circular) is done
downstream by 05_qc_summary.py, not duplicated here.

Output: analysis/qc_basic_stats/results/gfa_stats.tsv
Columns: species organelle subgraph_id n_subgraphs n_nodes n_edges circular
         dead_end_nodes total_length gc_content avg_coverage parse_status
"""
from __future__ import annotations

import argparse
import re
import subprocess
import os
import sys
from pathlib import Path

# Defaults to the old same-tree layout (code and data siblings under one
# repo root); set PLANT_ORGANELLE_DATA_ROOT to point at a separate data
# checkout instead (e.g. when this code is installed from its own repo).
ROOT_DIR = Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT", str(Path(__file__).resolve().parents[3])))
ANALYSIS_DIR = ROOT_DIR / "analysis"
sys.path.insert(0, str(ANALYSIS_DIR / "common"))

import io_utils  # noqa: E402
import species_discovery as sd  # noqa: E402

COLUMNS = ["species", "organelle", "subgraph_id", "n_subgraphs", "n_nodes", "n_edges",
           "circular", "dead_end_nodes", "total_length", "gc_content", "avg_coverage", "parse_status"]
KEY_COLS = ["species", "organelle"]

import shutil  # noqa: E402
GFATK = shutil.which("gfatk") or str(Path.home() / ".cargo" / "bin" / "gfatk")


def parse_gfatk_stats(text: str):
    total_match = re.search(r"Total number of subgraphs:\s*(\d+)", text)
    n_subgraphs = int(total_match.group(1)) if total_match else None
    blocks = re.split(r"Subgraph \d+:", text)[1:]
    subgraphs = []
    for block in blocks:
        def grab(pattern):
            m = re.search(pattern, block)
            return m.group(1) if m else None
        subgraphs.append({
            "n_nodes": grab(r"Number of nodes/segments:\s*(\d+)"),
            "n_edges": grab(r"Number of edges/links:\s*(\d+)"),
            "circular": grab(r"Circular:\s*(\w+)"),
            "dead_end_nodes": grab(r"Dead-end nodes:\s*(\d+)"),
            "total_length": grab(r"Total sequence length:\s*(\d+)"),
            "gc_content": grab(r"GC content of total sequence:\s*([\d.eE+-]+)"),
            "avg_coverage": grab(r"Average coverage of total segments:\s*([\d.eE+-]+)"),
        })
    return n_subgraphs, subgraphs


def rows_for_species(species: str, organelle: str, gfa_path: str | None) -> list[dict]:
    base = {"species": species, "organelle": organelle}
    if not gfa_path:
        return [{**base, "subgraph_id": None, "n_subgraphs": None, "n_nodes": None, "n_edges": None,
                 "circular": None, "dead_end_nodes": None, "total_length": None, "gc_content": None,
                 "avg_coverage": None, "parse_status": "missing_gfa"}]
    try:
        proc = subprocess.run([GFATK, "stats", gfa_path], capture_output=True, text=True, timeout=120)
    except Exception as exc:  # noqa: BLE001
        return [{**base, "subgraph_id": None, "n_subgraphs": None, "n_nodes": None, "n_edges": None,
                 "circular": None, "dead_end_nodes": None, "total_length": None, "gc_content": None,
                 "avg_coverage": None, "parse_status": f"gfatk_error:{exc}"}]
    if proc.returncode != 0:
        return [{**base, "subgraph_id": None, "n_subgraphs": None, "n_nodes": None, "n_edges": None,
                 "circular": None, "dead_end_nodes": None, "total_length": None, "gc_content": None,
                 "avg_coverage": None, "parse_status": "gfatk_error"}]

    n_subgraphs, subgraphs = parse_gfatk_stats(proc.stdout)
    if n_subgraphs is None:
        return [{**base, "subgraph_id": None, "n_subgraphs": None, "n_nodes": None, "n_edges": None,
                 "circular": None, "dead_end_nodes": None, "total_length": None, "gc_content": None,
                 "avg_coverage": None, "parse_status": "gfatk_error"}]
    if n_subgraphs == 0:
        return [{**base, "subgraph_id": None, "n_subgraphs": 0, "n_nodes": None, "n_edges": None,
                 "circular": None, "dead_end_nodes": None, "total_length": None, "gc_content": None,
                 "avg_coverage": None, "parse_status": "empty"}]

    rows = []
    for i, sg in enumerate(subgraphs, start=1):
        rows.append({**base, "subgraph_id": i, "n_subgraphs": n_subgraphs, "parse_status": "ok", **sg})
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organelle", default="both", choices=["mito", "pltd", "both"])
    ap.add_argument("--species-list", default=None)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    out_path = ANALYSIS_DIR / "qc_basic_stats" / "results" / "gfa_stats.tsv"
    existing_df, existing_keys = io_utils.load_existing(out_path, KEY_COLS)

    organelles = ["mito", "pltd"] if args.organelle == "both" else [args.organelle]
    species_filter = sd.load_species_list(Path(args.species_list)) if args.species_list else None

    new_rows = []
    n_computed = n_skipped = 0
    for organelle in organelles:
        data_root = sd.repo_data_root(ROOT_DIR, organelle)
        resolved = sd.discover_all(data_root, organelle, species_filter=species_filter,
                                    logger=lambda m: print(f"[warn] {m}", file=sys.stderr))
        for r in resolved:
            key = (r.species, organelle)
            input_path = Path(r.gfa) if r.gfa else data_root / r.species
            if io_utils.needs_recompute(key, existing_keys, out_path, input_path, args.force):
                new_rows.extend(rows_for_species(r.species, organelle, r.gfa))
                n_computed += 1
            else:
                n_skipped += 1

    merged = io_utils.merge_rows(existing_df, new_rows, KEY_COLS, COLUMNS)
    io_utils.atomic_write_tsv(out_path, merged.to_dict("records"), COLUMNS)
    print(f"[info] gfa_stats: computed={n_computed} skipped={n_skipped} total_rows={len(merged)} -> {out_path}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
