#!/usr/bin/env python3
"""Extract one FASTA per phylogeny marker gene, one sequence per species.

Marker genes = qc_basic_stats/results/core_genes_{organelle}.tsv entries
with pct_present>=0.90 AND pct_single_copy>=0.90 (computed there, not
hardcoded - this mechanically excludes multi-exon/trans-spliced genes).
For each gene, the single best-scoring .ctg.bed hit per species is
extracted from that species' .ctg.fasta via `samtools faidx`, one
subprocess call per species covering all needed genes at once (not one
call per gene) to keep this cheap at ~1200 species.

Output: analysis/annotation/results/genes/{mito,pltd}/<gene>.fasta
        (headers are just the species name)
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

SAMTOOLS = shutil.which("samtools") or "/software/team301/samtools/samtools"
MARKER_PRESENT_THRESHOLD = 0.90
# 0.80, not 0.90: verified at full-dataset scale (~1200 species) that mito genes
# top out around 89% single-copy (ccmB, the best, is 89.2%) - real biology (plant
# mitochondrial genomes duplicate genes across multipartite structure far more than
# plastids, which comfortably clear 93%+), not a bug. A 90% bar excluded every mito
# gene; 80% captures the genuinely-mostly-single-copy ones (ccmB, nad5, trnI-CAU,
# mttB, trnE-UUC, ccmC, trnY-GUA, trnQ-UUG, atp4, ...) without being arbitrary.
MARKER_SINGLE_COPY_THRESHOLD = 0.80
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def marker_genes(organelle: str) -> list[str]:
    path = ANALYSIS_DIR / "qc_basic_stats" / "results" / f"core_genes_{organelle}.tsv"
    if not path.exists():
        print(f"[err] missing {path} - run qc_basic_stats/src/04_core_gene_list.py first", file=sys.stderr)
        return []
    df = pd.read_csv(path, sep="\t")
    sub = df[(df["pct_present"] >= MARKER_PRESENT_THRESHOLD) & (df["pct_single_copy"] >= MARKER_SINGLE_COPY_THRESHOLD)]
    return sorted(sub["gene"])


def best_hits_per_species(gene_calls: pd.DataFrame, organelle: str, genes: list[str]) -> pd.DataFrame:
    sub = gene_calls[(gene_calls["organelle"] == organelle) & (gene_calls["gene"].isin(genes))]
    if sub.empty:
        return sub
    idx = sub.groupby(["species", "gene"])["score"].idxmax()
    return sub.loc[idx]


def extract_for_species(species: str, fasta_path: str, hits: pd.DataFrame) -> dict[str, str]:
    """Return {gene: sequence} for this species' marker-gene hits.

    `samtools faidx fasta region1 region2 ...` guarantees output records
    are in the same order as the requested regions, so records are matched
    to genes by position rather than by re-parsing headers (simpler and
    avoids off-by-one bugs when hand-rolling a streaming FASTA parser).
    """
    rows = list(hits.itertuples())
    regions = [f"{row.contig_id}:{row.start + 1}-{row.end}" for row in rows]
    if not regions:
        return {}
    subprocess.run([SAMTOOLS, "faidx", fasta_path], capture_output=True, timeout=60)  # ensure .fai exists
    proc = subprocess.run([SAMTOOLS, "faidx", fasta_path, *regions], capture_output=True, text=True, timeout=60)
    if proc.returncode != 0:
        print(f"[warn] samtools faidx failed for {species}: {proc.stderr.strip()}", file=sys.stderr)
        return {}

    records = [chunk.split("\n", 1)[1].replace("\n", "") for chunk in proc.stdout.split(">")[1:]]
    if len(records) != len(rows):
        print(f"[warn] {species}: expected {len(rows)} faidx records, got {len(records)} - skipping", file=sys.stderr)
        return {}

    seqs = {}
    for row, seq in zip(rows, records):
        seqs[row.gene] = revcomp(seq) if row.strand == "-" else seq
    return seqs


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organelle", default="both", choices=["mito", "pltd", "both"])
    ap.add_argument("--species-list", default=None)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    gene_calls_path = ANALYSIS_DIR / "annotation" / "results" / "gene_calls.tsv"
    if not gene_calls_path.exists():
        print(f"[err] missing {gene_calls_path} - run 01_gene_table.py first", file=sys.stderr)
        sys.exit(1)
    gene_calls = pd.read_csv(gene_calls_path, sep="\t")

    species_filter = sd.load_species_list(Path(args.species_list)) if args.species_list else None
    organelles = ["mito", "pltd"] if args.organelle == "both" else [args.organelle]

    for organelle in organelles:
        genes = marker_genes(organelle)
        if not genes:
            print(f"[warn] {organelle}: no marker genes found, skipping", file=sys.stderr)
            continue
        out_dir = ANALYSIS_DIR / "annotation" / "results" / "genes" / organelle
        out_dir.mkdir(parents=True, exist_ok=True)

        # skip whole-organelle rebuild if not forced and every gene file already has every eligible species
        eligible_species = sorted(gene_calls.loc[gene_calls["organelle"] == organelle, "species"].unique())
        if species_filter is not None:
            eligible_species = [s for s in eligible_species if s in species_filter]
        if not args.force:
            def has_all_species(gene_file: Path) -> bool:
                if not gene_file.exists():
                    return False
                headers = {l[1:].strip() for l in gene_file.read_text().splitlines() if l.startswith(">")}
                return set(eligible_species).issubset(headers)
            if all(has_all_species(out_dir / f"{g}.fasta") for g in genes):
                print(f"[info] gene_fasta_extract: {organelle}: up to date, skipping ({len(genes)} genes)",
                      file=sys.stderr)
                continue

        hits = best_hits_per_species(gene_calls, organelle, genes)
        data_root = sd.repo_data_root(ROOT_DIR, organelle)
        gene_seqs: dict[str, dict[str, str]] = {g: {} for g in genes}

        for species in eligible_species:
            sp_hits = hits[hits["species"] == species]
            if sp_hits.empty:
                continue
            r = sd.resolve_species(data_root / species, organelle)
            if r.status != "ok":
                continue
            seqs = extract_for_species(species, r.ctg_fasta, sp_hits)
            for gene, seq in seqs.items():
                gene_seqs[gene][species] = seq

        for gene in genes:
            out_path = out_dir / f"{gene}.fasta"
            with open(out_path, "w") as fh:
                for species, seq in sorted(gene_seqs[gene].items()):
                    fh.write(f">{species}\n{seq}\n")
        n_seqs = sum(len(v) for v in gene_seqs.values())
        print(f"[info] gene_fasta_extract: {organelle}: {len(genes)} genes, {n_seqs} sequences -> {out_dir}",
              file=sys.stderr)


if __name__ == "__main__":
    main()
