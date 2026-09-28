#!/usr/bin/env python3
"""Render a matplotlib dotplot for every PAF under synteny/results/{mtpt,paf}/.

Each alignment block is drawn as a line segment from (target_start,
query_start) to (target_end, query_end) - or mirrored for a '-' strand
block - colored by strand, so inversions are visually obvious (a classic,
simple way to read structural rearrangements from a PAF without extra
tooling).

Output: analysis/synteny/results/dotplots/<paf_stem>.png
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# Defaults to the old same-tree layout (code and data siblings under one
# repo root); set PLANT_ORGANELLE_DATA_ROOT to point at a separate data
# checkout instead (e.g. when this code is installed from its own repo).
ROOT_DIR = Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT", str(Path(__file__).resolve().parents[3])))
ANALYSIS_DIR = ROOT_DIR / "analysis"
sys.path.insert(0, str(ANALYSIS_DIR / "common"))
import minimap_utils  # noqa: E402


def load_paf(path: Path) -> list[dict]:
    records = []
    for line in path.read_text().splitlines():
        fields = line.split("\t")
        if len(fields) < 12:
            continue
        rec = dict(zip(minimap_utils.PAF_COLS, fields[:12]))
        for k in ("query_length", "query_start", "query_end", "target_length", "target_start", "target_end"):
            rec[k] = int(rec[k])
        records.append(rec)
    return records


def plot_paf(records: list[dict], title: str, out_path: Path):
    fig, ax = plt.subplots(figsize=(6, 6))
    for r in records:
        color = "#2b6cb0" if r["strand"] == "+" else "#c53030"
        if r["strand"] == "+":
            xs, ys = [r["target_start"], r["target_end"]], [r["query_start"], r["query_end"]]
        else:
            xs, ys = [r["target_start"], r["target_end"]], [r["query_end"], r["query_start"]]
        ax.plot(xs, ys, color=color, linewidth=1.5, solid_capstyle="round")

    target_name = records[0]["target_name"] if records else "target"
    query_name = records[0]["query_name"] if records else "query"
    ax.set_xlabel(f"target: {target_name}")
    ax.set_ylabel(f"query: {query_name}")
    ax.set_title(title, fontsize=10)
    ax.plot([], [], color="#2b6cb0", label="+ strand")
    ax.plot([], [], color="#c53030", label="- strand")
    ax.legend(loc="upper left", fontsize=8)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=120)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    out_dir = ANALYSIS_DIR / "synteny" / "results" / "dotplots"
    paf_dirs = [ANALYSIS_DIR / "synteny" / "results" / "mtpt", ANALYSIS_DIR / "synteny" / "results" / "paf"]

    n_plotted = n_skipped = n_empty = 0
    for paf_dir in paf_dirs:
        if not paf_dir.exists():
            continue
        for paf_path in sorted(paf_dir.glob("*.paf")):
            out_path = out_dir / f"{paf_path.stem}.png"
            if not args.force and out_path.exists() and out_path.stat().st_mtime >= paf_path.stat().st_mtime:
                n_skipped += 1
                continue
            records = load_paf(paf_path)
            if not records:
                n_empty += 1
                continue
            plot_paf(records, paf_path.stem, out_path)
            n_plotted += 1

    print(f"[info] dotplot: plotted={n_plotted} skipped={n_skipped} empty={n_empty} -> {out_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
