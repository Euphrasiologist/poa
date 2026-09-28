#!/usr/bin/env python3
"""Build one standardized, hierarchical GFF3 per species, merging the
three annotation modules by precedence:

1. `trans_splicing` reconstruction (`nad1`/`nad2`/`nad5`/`rps3` only, if a
   reconstruction with >=1 filled exon exists - partial is still more
   informative than a single raw fragment).
2. `editing`'s RNA-editing-corrected single-exon call, if this gene has one.
3. `denovo_annotation`'s raw best-scoring hit, otherwise (also the only
   tier for tRNA/rRNA genes - neither of the other two modules touches
   those gene classes).

Which tier won is recorded per gene as `annotation_tier=reconstructed|
edited|raw` - an explicit confidence signal, not just implicit in which
script you happened to run.

GFF3 mechanics used here (see `../README.md` for the full explanation):
- Protein-coding genes are `gene` -> `mRNA` -> one or more `CDS` rows
  sharing one `ID` (that's the spec's actual mechanism for a discontinuous/
  multi-exon feature - it does not require those rows to share a strand or
  seqid, which is exactly what a trans-spliced gene needs).
- tRNA/rRNA are `gene` -> `tRNA`/`rRNA` directly - no `CDS`/`mRNA` layer,
  they're not protein-coding.
- Predicted RNA-editing sites are child `sequence_alteration` features
  under the relevant `CDS` `ID` (a real, general Sequence Ontology term
  for "this position differs from a reference" - a pragmatic choice, not
  a claim that this is THE canonical way to encode RNA editing in GFF3).
- `phase` (GFF3 column 8) is written as `0` for every `CDS` - a verified
  guarantee for `editing`/`trans_splicing`-derived exons (their DP only
  ever emits codon-aligned boundaries by construction) but an unverified
  *assumption* for tier-3 raw `denovo_annotation`-only calls (HMMER
  envelope coordinates aren't codon-boundary-enforced) - see README.

Output: analysis/gff_export/results/<species>.<organelle>.gff
"""
from __future__ import annotations

import argparse
import os
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

# Defaults to the old same-tree layout (code and data siblings under one
# repo root); set PLANT_ORGANELLE_DATA_ROOT to point at a separate data
# checkout instead (e.g. when this code is installed from its own repo).
ROOT_DIR = Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT", str(Path(__file__).resolve().parents[3])))
ANALYSIS_DIR = ROOT_DIR / "analysis"
sys.path.insert(0, str(ANALYSIS_DIR / "common"))

import species_discovery as sd  # noqa: E402

TRANS_SPLICING_GENES = {"nad1", "nad2", "nad5", "rps3"}


def best_hit_per_gene(gene_calls: pd.DataFrame) -> pd.DataFrame:
    if gene_calls.empty:
        return gene_calls
    idx = gene_calls.groupby(["species", "organelle", "gene"])["score"].idxmax()
    return gene_calls.loc[idx]


def gff_attrs(d: dict) -> str:
    return ";".join(f"{k}={v}" for k, v in d.items() if v is not None and v != "")


def emit_rna_gene(lines: list[str], contig: str, start: int, end: int, strand: str,
                   score: float, gene: str, feature: str) -> None:
    gid = f"gene-{gene}"
    lines.append("\t".join([contig, "gff_export", "gene", str(start + 1), str(end), f"{score:.3f}", strand, ".",
                             gff_attrs({"ID": gid, "Name": gene, "annotation_tier": "raw"})]))
    lines.append("\t".join([contig, "gff_export", feature, str(start + 1), str(end), f"{score:.3f}", strand, ".",
                             gff_attrs({"ID": f"{feature.lower()}-{gene}", "Parent": gid, "Name": gene})]))


def emit_protein_coding_gene(lines: list[str], gene: str, tier: str, parts: list[dict],
                              edits: list[dict], gene_score: float, n_edits_summary: int) -> None:
    """parts: list of {contig, start, end, strand, score, exon_number} in exon order.
    edits: list of {genomic_pos, resulting_aa} to attach as sequence_alteration children."""
    contig = parts[0]["contig"]
    gene_start = min(p["start"] for p in parts)
    gene_end = max(p["end"] for p in parts)
    gene_strand = parts[0]["strand"] if len(parts) == 1 else "."  # mixed-strand for real trans-spliced genes
    gid, mid, cid = f"gene-{gene}", f"mRNA-{gene}", f"cds-{gene}"

    lines.append("\t".join([contig, "gff_export", "gene", str(gene_start + 1), str(gene_end), f"{gene_score:.3f}",
                             gene_strand, ".",
                             gff_attrs({"ID": gid, "Name": gene, "annotation_tier": tier})]))
    lines.append("\t".join([contig, "gff_export", "mRNA", str(gene_start + 1), str(gene_end), f"{gene_score:.3f}",
                             gene_strand, ".",
                             gff_attrs({"ID": mid, "Parent": gid, "Name": gene,
                                        "n_exons": len(parts), "n_edits": n_edits_summary})]))
    for p in parts:
        lines.append("\t".join([p["contig"], "gff_export", "CDS", str(p["start"] + 1), str(p["end"]),
                                 f"{p['score']:.3f}", p["strand"], "0",
                                 gff_attrs({"ID": cid, "Parent": mid, "exon_number": p["exon_number"]})]))
    for i, e in enumerate(edits, start=1):
        pos = int(e["genomic_pos"])
        lines.append("\t".join([contig, "gff_export", "sequence_alteration", str(pos + 1), str(pos + 1), ".",
                                 gene_strand if gene_strand != "." else ".", ".",
                                 gff_attrs({"ID": f"sa-{gene}-{i}", "Parent": cid,
                                            "edited_from": "C", "edited_to": "T",
                                            "resulting_aa": e.get("resulting_aa", "")})]))


def build_species_gff(species: str, organelle: str, raw_best: pd.DataFrame,
                       editing_calls: pd.DataFrame, editing_edits: pd.DataFrame,
                       ts_genes: pd.DataFrame, ts_exons: pd.DataFrame, ts_edits: pd.DataFrame) -> list[str]:
    lines: list[str] = []
    sub = raw_best[(raw_best.species == species) & (raw_best.organelle == organelle)]

    for row in sub.itertuples():
        gene = row.gene
        if gene.startswith("trn"):
            emit_rna_gene(lines, row.contig_id, row.start, row.end, row.strand, row.score, gene, "tRNA")
            continue
        if gene.startswith("rrn"):
            emit_rna_gene(lines, row.contig_id, row.start, row.end, row.strand, row.score, gene, "rRNA")
            continue

        if gene in TRANS_SPLICING_GENES:
            g = ts_genes[(ts_genes.species == species) & (ts_genes.gene == gene)]
            if not g.empty and int(g.iloc[0].n_filled) >= 1:
                grow = g.iloc[0]
                exons = ts_exons[(ts_exons.species == species) & (ts_exons.gene == gene)].sort_values("slot")
                parts = [{"contig": e.contig, "start": int(e.start), "end": int(e.end), "strand": e.strand,
                          "score": e.score, "exon_number": int(e.slot)} for e in exons.itertuples()]
                edits = ts_edits[(ts_edits.species == species) & (ts_edits.gene == gene)].to_dict("records")
                emit_protein_coding_gene(lines, gene, "reconstructed", parts, edits,
                                          grow.whole_gene_score, int(grow.whole_gene_n_edits))
                continue

        e = editing_calls[(editing_calls.species == species) & (editing_calls.gene == gene)]
        if not e.empty:
            erow = e.iloc[0]
            part = {"contig": row.contig_id, "start": int(erow.corrected_start), "end": int(erow.corrected_end),
                    "strand": row.strand, "score": erow.score, "exon_number": 1}
            edits = editing_edits[(editing_edits.species == species) & (editing_edits.gene == gene)].to_dict("records")
            emit_protein_coding_gene(lines, gene, "edited", [part], edits, erow.score, int(erow.n_edits))
            continue

        part = {"contig": row.contig_id, "start": int(row.start), "end": int(row.end),
                "strand": row.strand, "score": row.score, "exon_number": 1}
        emit_protein_coding_gene(lines, gene, "raw", [part], [], row.score, 0)

    def sort_key(line: str):
        f = line.split("\t")
        return (f[0], int(f[3]))
    lines.sort(key=sort_key)
    return lines


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organelle", default="both", choices=["mito", "pltd", "both"])
    ap.add_argument("--species-list", default=None)
    ap.add_argument("--gene-calls", default=None, help="raw gene_calls.tsv to use for tier-3 fallback + "
                     "tRNA/rRNA + contig/strand lookups (default: denovo_annotation/results/gene_calls.tsv)")
    args = ap.parse_args()

    species_filter = sd.load_species_list(Path(args.species_list)) if args.species_list else None

    raw_path = Path(args.gene_calls) if args.gene_calls else ANALYSIS_DIR / "denovo_annotation" / "results" / "gene_calls.tsv"
    if not raw_path.exists():
        print(f"[err] missing {raw_path} - run denovo_annotation/src/04_build_gene_calls.py first "
              f"(or pass --gene-calls)", file=sys.stderr)
        sys.exit(1)
    raw = pd.read_csv(raw_path, sep="\t")
    if species_filter is not None:
        raw = raw[raw.species.isin(species_filter)]
    raw_best = best_hit_per_gene(raw)

    def load(path: Path, columns: list[str]) -> pd.DataFrame:
        # Always returns a DataFrame with `columns` present (empty if the
        # upstream file doesn't exist yet, e.g. plastid before editing/
        # trans_splicing ever cover it) so every `.species`/`.gene` filter
        # downstream is safe regardless of which modules have been run.
        if not path.exists():
            return pd.DataFrame(columns=columns)
        df = pd.read_csv(path, sep="\t")
        return df[df.species.isin(species_filter)] if species_filter is not None and not df.empty else df

    editing_calls = load(ANALYSIS_DIR / "editing" / "results" / "edited_gene_calls.tsv",
                          ["species", "organelle", "gene", "corrected_start", "corrected_end", "score", "n_edits"])
    editing_edits = load(ANALYSIS_DIR / "editing" / "results" / "predicted_edits.tsv",
                          ["species", "organelle", "gene", "genomic_pos", "resulting_aa"])
    ts_genes = load(ANALYSIS_DIR / "trans_splicing" / "results" / "reconstructed_genes.tsv",
                     ["species", "organelle", "gene", "n_filled", "whole_gene_score", "whole_gene_n_edits"])
    ts_exons = load(ANALYSIS_DIR / "trans_splicing" / "results" / "reconstructed_exons.tsv",
                     ["species", "organelle", "gene", "slot", "contig", "start", "end", "strand", "score"])
    ts_edits = load(ANALYSIS_DIR / "trans_splicing" / "results" / "reconstructed_edits.tsv",
                     ["species", "organelle", "gene", "slot", "genomic_pos", "resulting_aa"])

    organelles = ["mito", "pltd"] if args.organelle == "both" else [args.organelle]
    out_dir = ANALYSIS_DIR / "gff_export" / "results"
    out_dir.mkdir(parents=True, exist_ok=True)

    n_written = 0
    for organelle in organelles:
        species_here = sorted(raw_best[raw_best.organelle == organelle].species.unique())
        for species in species_here:
            lines = build_species_gff(species, organelle, raw_best, editing_calls, editing_edits,
                                       ts_genes, ts_exons, ts_edits)
            if not lines:
                continue
            out_path = out_dir / f"{species}.{organelle}.gff"
            out_path.write_text("##gff-version 3\n" + "\n".join(lines) + "\n")
            n_written += 1
    print(f"[info] build_gff: {n_written} species -> {out_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
