#!/usr/bin/env python3
"""Render a gene-labeled Bandage assembly-graph topology plot (genes
colored/labeled via Bandage's own internal BLAST search against the gene
query FASTA from 01_make_gene_queries.py). Only the gene-labeled version is
rendered - an earlier version also rendered a plain (no gene labels)
variant, but that's a strict subset of what the gene-labeled one shows
(same topology, plus gene positions), so keeping both was pure duplication.

Reproduces the approach in the sibling mito_structural_variation repo's
figures/fig1b_genome_topology_examples/ (see that dir's README) against
this repo's data layout - no manual editing, so see the caveat below.

Output format is PNG by default (smaller files, viewable in any image
viewer with no SVG renderer needed) - Bandage's own `image` command
supports .png/.jpg/.svg output natively, so no external SVG->PNG
conversion step (e.g. resvg) is needed either way. Pass --format svg if
you want vector output instead (matches the reference figures' own
_genes.svg naming).

Known limitation (not fixed here - inherent to Bandage): its force-directed
graph layout is unseeded, with no way to pass a fixed layout back in. For
a genuinely tangled multi-node single-component graph, two renders of the
same GFA can look meaningfully different. Single-node and clearly-separate
multipartite topologies aren't affected (nothing to tangle). If a specific
plot looks bad, just re-render that one species with --force - you may get
a cleaner layout on the next attempt.

Output: analysis/topology_plots/results/plots/<organelle>/<species>.genes.<ext>
"""
from __future__ import annotations

import argparse
import os
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

import pandas as pd  # noqa: E402
import species_discovery as sd  # noqa: E402

BANDAGE_ENV = os.environ.get("BANDAGE_ENV", "/nfs/users/nfs_m/mb39/miniconda3/envs/bandage")
BANDAGE = os.environ.get("BANDAGE", str(Path(BANDAGE_ENV) / "bin" / "Bandage")
                          )
BANDAGE_RUN_ENV = {**os.environ, "PATH": f"{BANDAGE_ENV}/bin:{os.environ.get('PATH', '')}",
                    "QT_QPA_PLATFORM": os.environ.get("QT_QPA_PLATFORM", "offscreen")}


def render(gfa: str, out_path: Path, query_fasta: str | None, height: int, timeout: int) -> bool:
    cmd = [BANDAGE, "image", gfa, str(out_path), "--height", str(height), "--names", "--lengths", "--fontsize", "12"]
    if query_fasta:
        cmd += ["--query", query_fasta, "--colour", "blastsolid", "--blasthits"]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, env=BANDAGE_RUN_ENV)
    if proc.returncode != 0 or not out_path.exists():
        print(f"[warn] Bandage failed for {out_path.name}: {proc.stderr.strip()[:300]}", file=sys.stderr)
        return False
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organelle", default="both", choices=["mito", "pltd", "both"])
    ap.add_argument("--species-list", default=None)
    ap.add_argument("--qc-status", default="fail,flag",
                     help="comma-separated subset of pass,flag,fail (default fail,flag)")
    ap.add_argument("--format", default="png", choices=["png", "svg", "jpg"],
                     help="output image format, passed straight to Bandage (default png)")
    ap.add_argument("--height", type=int, default=900)
    ap.add_argument("--timeout", type=int, default=120, help="per-render timeout in seconds (default 120)")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    if not Path(BANDAGE).exists():
        print(f"[err] Bandage not found at {BANDAGE} - see analysis/common/tool_paths.sh", file=sys.stderr)
        sys.exit(1)

    qc_path = ANALYSIS_DIR / "qc_basic_stats" / "results" / "qc_summary.tsv"
    qc = pd.read_csv(qc_path, sep="\t")
    qc_status = args.qc_status.split(",")
    eligible = set(zip(qc.loc[qc["status"].isin(qc_status), "species"],
                        qc.loc[qc["status"].isin(qc_status), "organelle"]))

    organelles = ["mito", "pltd"] if args.organelle == "both" else [args.organelle]
    species_filter = sd.load_species_list(Path(args.species_list)) if args.species_list else None
    work_dir = ANALYSIS_DIR / "topology_plots" / "work"
    plots_dir = ANALYSIS_DIR / "topology_plots" / "results" / "plots"
    ext = args.format

    n_rendered = n_skipped = n_no_query = 0
    for organelle in organelles:
        data_root = sd.repo_data_root(ROOT_DIR, organelle)
        resolved = sd.discover_all(data_root, organelle, species_filter=species_filter)
        out_dir = plots_dir / organelle
        for r in resolved:
            if (r.species, organelle) not in eligible or r.status != "ok":
                continue

            genes_out = out_dir / f"{r.species}.genes.{ext}"
            query_fasta = work_dir / f"{r.species}.{organelle}.genes.fasta"

            if not args.force and genes_out.exists() and \
                    genes_out.stat().st_mtime >= Path(r.gfa).stat().st_mtime:
                n_skipped += 1
                continue

            if query_fasta.exists() and query_fasta.stat().st_size > 0:
                render(r.gfa, genes_out, str(query_fasta), args.height, args.timeout)
                n_rendered += 1
            else:
                n_no_query += 1
                print(f"[warn] {r.species}.{organelle}: no gene query fasta "
                      f"(run 01_make_gene_queries.py first) - skipped", file=sys.stderr)

    print(f"[info] render_topology: rendered={n_rendered} skipped={n_skipped} "
          f"no_query={n_no_query} -> {plots_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
