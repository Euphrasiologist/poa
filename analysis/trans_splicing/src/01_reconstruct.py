#!/usr/bin/env python3
"""Run `transsplice` over every species' calls for the genes in
reference_genes.GENES - trans-spliced (nad1/nad2/nad5/rps3) and, since
2026-09-29, cis-spliced ones too (transsplice classifies each junction
cis/trans per species) - reconstructing full-length genes from their per-exon
`.ctg.bed` hits.

Unlike `editing/01_scan_editing.py` (best-scoring hit only per gene), every
raw `.ctg.bed` row for these genes is a candidate exon fragment here - all
of them go into the manifest, and `transsplice` itself does the
merging/assignment (see its README).

Output: analysis/trans_splicing/results/reconstructed_genes.tsv
        analysis/trans_splicing/results/reconstructed_exons.tsv
        analysis/trans_splicing/results/reconstructed_junctions.tsv
"""
from __future__ import annotations

import os
import argparse
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

# CODE_DIR is this checkout's analysis/ (shared modules, bundled reference
# data); data/ and every module's results/ and work/ live under the data
# root - PLANT_ORGANELLE_DATA_ROOT, or this checkout if unset.
CODE_DIR = Path(__file__).resolve().parents[2]
ROOT_DIR = Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT", str(CODE_DIR.parent)))
ANALYSIS_DIR = ROOT_DIR / "analysis"
sys.path.insert(0, str(CODE_DIR / "common"))

import species_discovery as sd  # noqa: E402

TRANSSPLICE = shutil.which("transsplice") or str(Path.home() / ".cargo" / "bin" / "transsplice")
sys.path.insert(0, str(Path(__file__).resolve().parent))
from reference_genes import GENES  # noqa: E402  (shared with 00a/01/gff_export)


def build_manifest(gene_calls: pd.DataFrame, species_filter: set | None, profiles_dir: Path) -> list[dict]:
    rows = []
    available_genes = [g for g in GENES if (profiles_dir / g / "whole_gene.pssm").exists()]
    missing = set(GENES) - set(available_genes)
    if missing:
        print(f"[warn] no template for {sorted(missing)} - run 00_build_exon_profiles.py first", file=sys.stderr)

    sub = gene_calls[(gene_calls["organelle"] == "mito") & (gene_calls["gene"].isin(available_genes))]
    data_root = sd.repo_data_root(ROOT_DIR, "mito")
    fasta_cache: dict[str, str | None] = {}
    for row in sub.itertuples():
        if species_filter is not None and row.species not in species_filter:
            continue
        if row.species not in fasta_cache:
            r = sd.resolve_species(data_root / row.species, "mito")
            fasta_cache[row.species] = r.ctg_fasta if r.status == "ok" else None
        fasta_path = fasta_cache[row.species]
        if fasta_path is None:
            continue
        rows.append({
            "species": row.species, "organelle": "mito", "gene": row.gene,
            "fasta": fasta_path, "template_dir": str(profiles_dir / row.gene),
            "contig": row.contig_id, "start": int(row.start), "end": int(row.end),
            "strand": row.strand,
        })
    return rows


def write_manifest(rows: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    cols = ["species", "organelle", "gene", "fasta", "template_dir", "contig", "start", "end", "strand"]
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--species-list", default=None)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--merge-distance", type=int, default=5000)
    ap.add_argument("--junction-distance-threshold", type=int, default=10000)
    ap.add_argument("--exon-flank", type=int, default=60)
    ap.add_argument("--gene-calls", default=None, help="path to a gene_calls.tsv-shaped table to anchor "
                     "on (default: annotation/results/gene_calls.tsv, oatk's own calls) - e.g. point this "
                     "at denovo_annotation/results/gene_calls.tsv to compare/use the oatk-independent "
                     "nhmmscan-based calls instead")
    ap.add_argument("--force", action="store_true", help="accepted for CLI consistency; transsplice "
                     "always does a full recompute over the manifest (no incremental skip yet)")
    args = ap.parse_args()

    if not Path(TRANSSPLICE).exists():
        print(f"[err] transsplice binary not found at {TRANSSPLICE} (`cargo install --path .` in its repo)",
              file=sys.stderr)
        sys.exit(1)

    gene_calls_path = Path(args.gene_calls) if args.gene_calls else ANALYSIS_DIR / "annotation" / "results" / "gene_calls.tsv"
    if not gene_calls_path.exists():
        print(f"[err] missing {gene_calls_path} - run annotation/src/01_gene_table.py first "
              f"(or pass --gene-calls)", file=sys.stderr)
        sys.exit(1)
    gene_calls = pd.read_csv(gene_calls_path, sep="\t")

    species_filter = sd.load_species_list(Path(args.species_list)) if args.species_list else None
    profiles_dir = CODE_DIR / "trans_splicing" / "reference" / "profiles"

    rows = build_manifest(gene_calls, species_filter, profiles_dir)
    if not rows:
        print("[err] empty manifest - nothing to reconstruct", file=sys.stderr)
        sys.exit(1)

    work_dir = ANALYSIS_DIR / "trans_splicing" / "work"
    manifest_path = work_dir / "manifest.tsv"
    write_manifest(rows, manifest_path)
    print(f"[info] reconstruct: manifest has {len(rows)} fragment rows -> {manifest_path}", file=sys.stderr)

    results_dir = ANALYSIS_DIR / "trans_splicing" / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        TRANSSPLICE, "scan-batch",
        "--manifest", str(manifest_path),
        "--genes-out", str(results_dir / "reconstructed_genes.tsv"),
        "--exons-out", str(results_dir / "reconstructed_exons.tsv"),
        "--junctions-out", str(results_dir / "reconstructed_junctions.tsv"),
        "--edits-out", str(results_dir / "reconstructed_edits.tsv"),
        "--threads", str(args.threads),
        "--merge-distance", str(args.merge_distance),
        "--junction-distance-threshold", str(args.junction_distance_threshold),
        "--exon-flank", str(args.exon_flank),
    ]
    proc = subprocess.run(cmd)
    if proc.returncode != 0:
        print(f"[err] transsplice scan-batch failed (exit {proc.returncode})", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
