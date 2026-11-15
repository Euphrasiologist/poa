"""Fetch RefSeq mitogenomes and extract verified exon sets, in the format
exon_profiles.py reads (<raw>/<acc>_<gene>[_<locus>].fasta +
junctions.tsv rows).

--lineage angiosperm (the default) writes reference/raw/, and only the
genes in reference_genes.CIS_SPLICED_MITO: the nad1/nad2/nad5/rps3
references were fetched ad hoc before this script existed and are left
exactly as they are. Any other lineage writes every gene in
reference_genes.GENES to reference/lineages/<lineage>/raw/, from
reference_genes.LINEAGE_MITOGENOMES[lineage]. Either way an accession+gene
that already has a file is skipped unless --force.

Verification, per locus: translate the joined exons from the genome and
compare with RefSeq's own curated /translation (common/refseq_mito.verify).
A locus that fails is still written but marked UNRESOLVED in
junctions.tsv, which exon_profiles.py excludes from training.

GenBank records are cached in reference/genbank/.

src/00a_fetch_references.py is the command-line wrapper (main() here).
"""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path

from Bio import SeqIO

from poa import PIPELINE_DIR
from poa.pipeline.common import lineage as lin
from poa.pipeline.common.refseq_mito import cds_exons, fetch, verify
from poa.pipeline.trans_splicing.reference_genes import CIS_SPLICED_MITO, GENE_SYNONYMS, GENES, LINEAGE_MITOGENOMES
from poa.pipeline.common.log import setup_logging

log = logging.getLogger(__name__)

GB_DIR = lin.trans_splicing_reference(PIPELINE_DIR, lin.DEFAULT_LINEAGE) / "genbank"
CIS_MAX_GAP = 10_000
JUNCTION_COLS = [
    "accession",
    "gene",
    "locus_tag",
    "exon_num",
    "n_exons",
    "start",
    "end",
    "strand",
    "junction_to_next",
    "gap_bp",
    "notes",
]


def fetch_references(lineage: str = lin.DEFAULT_LINEAGE, force: bool = False) -> Path:
    """Write `lineage`'s verified exon sets into the package's reference
    raw/ directory (returned). force: rewrite accession+gene files that
    already exist."""
    if lineage not in LINEAGE_MITOGENOMES:
        raise ValueError(f"lineage must be one of {sorted(LINEAGE_MITOGENOMES)}, not {lineage!r}")
    raw_dir = lin.trans_splicing_reference(PIPELINE_DIR, lineage) / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    genes = CIS_SPLICED_MITO if lineage == lin.DEFAULT_LINEAGE else GENES
    junctions_path = raw_dir / "junctions.tsv"
    existing_rows = []
    if junctions_path.exists():
        with open(junctions_path) as fh:
            existing_rows = list(csv.DictReader(fh, delimiter="\t"))
    done = {(r["accession"], r["gene"]) for r in existing_rows}
    wanted = {g.lower(): g for g in genes}
    new_rows, n_ok, n_bad = [], 0, 0

    for acc, organism in LINEAGE_MITOGENOMES[lineage].items():
        rec = SeqIO.read(fetch(acc, GB_DIR), "genbank")
        seen: dict[str, int] = {}
        for f in rec.features:
            if f.type != "CDS" or "gene" not in f.qualifiers or "translation" not in f.qualifiers:
                continue
            name = f.qualifiers["gene"][0].lower()
            gene = wanted.get(name) or (GENE_SYNONYMS.get(name) if GENE_SYNONYMS.get(name) in genes else None)
            if gene is None or ((acc, gene) in done and not force):
                continue
            parts, exons = cds_exons(rec, f)
            ok, note = verify("".join(exons), f.qualifiers["translation"][0])
            locus = f.qualifiers.get("locus_tag", [f"{gene}_{seen.get(gene, 0) + 1}"])[0]
            seen[gene] = seen.get(gene, 0) + 1
            n_ok += ok
            n_bad += not ok

            k = len(parts)
            suffix = (
                f"_{locus}"
                if seen[gene] > 1
                or sum(1 for g in rec.features if g.type == "CDS" and g.qualifiers.get("gene", [""])[0].lower() == name)
                > 1
                else ""
            )
            lines = []
            for n, (p, seq) in enumerate(zip(parts, exons), start=1):
                strand = "-" if p.strand == -1 else "+"
                lines.append(f">{acc}|{gene}|{locus}|exon{n}of{k}|{int(p.start) + 1}-{int(p.end)}|{strand}\n{seq}\n")
                if n < k:
                    q = parts[n]
                    same = q.strand == p.strand
                    gap = (int(q.start) - int(p.end)) if p.strand != -1 else (int(p.start) - int(q.end))
                    kind = "cis" if same and 0 < gap <= CIS_MAX_GAP else "trans"
                else:
                    gap, kind = "", "none"
                new_rows.append(
                    {
                        "accession": acc,
                        "gene": gene,
                        "locus_tag": locus,
                        "exon_num": n,
                        "n_exons": k,
                        "start": int(p.start) + 1,
                        "end": int(p.end),
                        "strand": strand,
                        "junction_to_next": kind,
                        "gap_bp": gap,
                        "notes": f"{organism}: {note}",
                    }
                )
            lines.append(
                f">{acc}|{gene}|{locus}|translation_REFSEQ_ANNOTATED|len={len(f.qualifiers['translation'][0])}\n"
                f"{f.qualifiers['translation'][0]}\n"
            )
            (raw_dir / f"{acc}_{gene}{suffix}.fasta").write_text("".join(lines))
            (log.info if ok else log.warning)(
                f"{'verified' if ok else 'UNRESOLVED'}: {acc} {organism}: {gene} ({locus}) {k} exon(s) - {note}"
            )

    if force:
        refreshed = {(r["accession"], r["gene"]) for r in new_rows}
        existing_rows = [r for r in existing_rows if (r["accession"], r["gene"]) not in refreshed]
    with open(junctions_path, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=JUNCTION_COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(existing_rows + new_rows)
    log.info(f"fetch_references: {n_ok} verified, {n_bad} UNRESOLVED loci written -> {raw_dir}")
    return raw_dir


def main(argv: list[str] | None = None) -> int:
    setup_logging()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lineage", default=lin.DEFAULT_LINEAGE, choices=sorted(LINEAGE_MITOGENOMES))
    ap.add_argument("--force", action="store_true", help="rewrite accession+gene files that already exist")
    args = ap.parse_args(argv)
    fetch_references(args.lineage, args.force)
    return 0
