#!/usr/bin/env python3
"""Build a per-species BLAST query FASTA of gene sequences, for Bandage's
--query/--blasthits gene-labeling of assembly graph plots.

One record per gene name (the longest interval across the whole assembly
if a gene has multiple hits/copies - Bandage's own internal BLAST search
against the graph will still find every location a gene's sequence occurs,
regardless of which single copy was used as the query), reverse-complemented
for '-' strand genes so all queries are in a consistent orientation. tRNAs
are excluded (too small and numerous to label legibly on a graph plot).

This mirrors the approach in the sibling mito_structural_variation repo's
figures/fig1b_genome_topology_examples/extract_gene_queries.py, reimplemented
against this repo's file layout (oatk's own .ctg.bed/.ctg.fasta - the same
file shape, just a different root path).

Output: analysis/topology_plots/work/<Species>.<organelle>.genes.fasta
"""
from __future__ import annotations

import os
import argparse
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

BED_COLS = ["contig_id", "start", "end", "gene", "score", "strand"]
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def load_fasta(path: str) -> dict:
    seqs, cur_id, cur_lines = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur_id is not None:
                    seqs[cur_id] = "".join(cur_lines)
                cur_id = line[1:].split()[0]
                cur_lines = []
            else:
                cur_lines.append(line.strip())
    if cur_id is not None:
        seqs[cur_id] = "".join(cur_lines)
    return seqs


def build_query_fasta(ctg_bed: str, ctg_fasta: str, out_path: Path) -> int:
    df = pd.read_csv(ctg_bed, sep="\t", comment="#", names=BED_COLS)
    df = df[~df["gene"].str.startswith("trn")]
    if df.empty:
        return 0
    df = df.copy()
    df["length"] = df["end"] - df["start"]
    best = df.loc[df.groupby("gene")["length"].idxmax()]

    seqs = load_fasta(ctg_fasta)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with open(out_path, "w") as fh:
        for row in best.itertuples():
            if row.contig_id not in seqs:
                continue
            seq = seqs[row.contig_id][row.start:row.end]
            if row.strand == "-":
                seq = revcomp(seq)
            fh.write(f">{row.gene}\n{seq}\n")
            n += 1
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organelle", default="both", choices=["mito", "pltd", "both"])
    ap.add_argument("--species-list", default=None)
    ap.add_argument("--qc-status", default="fail,flag",
                     help="comma-separated subset of pass,flag,fail (default fail,flag - "
                          "this tool exists to visually triage the species QC already flagged)")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    qc_path = ANALYSIS_DIR / "qc_basic_stats" / "results" / "qc_summary.tsv"
    qc = pd.read_csv(qc_path, sep="\t")
    qc_status = args.qc_status.split(",")
    eligible = set(zip(qc.loc[qc["status"].isin(qc_status), "species"],
                        qc.loc[qc["status"].isin(qc_status), "organelle"]))

    organelles = ["mito", "pltd"] if args.organelle == "both" else [args.organelle]
    species_filter = sd.load_species_list(Path(args.species_list)) if args.species_list else None

    work_dir = ANALYSIS_DIR / "topology_plots" / "work"
    n_built = n_skipped = n_no_genes = 0
    for organelle in organelles:
        data_root = sd.repo_data_root(ROOT_DIR, organelle)
        resolved = sd.discover_all(data_root, organelle, species_filter=species_filter)
        for r in resolved:
            if (r.species, organelle) not in eligible or r.status != "ok":
                continue
            out_path = work_dir / f"{r.species}.{organelle}.genes.fasta"
            if not args.force and out_path.exists() and Path(r.ctg_bed).stat().st_mtime <= out_path.stat().st_mtime:
                n_skipped += 1
                continue
            n_genes = build_query_fasta(r.ctg_bed, r.ctg_fasta, out_path)
            if n_genes == 0:
                n_no_genes += 1
            else:
                n_built += 1

    print(f"[info] make_gene_queries: built={n_built} skipped={n_skipped} no_genes={n_no_genes} -> {work_dir}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
