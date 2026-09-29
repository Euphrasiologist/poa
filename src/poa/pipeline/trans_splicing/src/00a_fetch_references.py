#!/usr/bin/env python3
"""Fetch RefSeq mitogenomes and extract verified exon sets for the genes in
reference_genes.CIS_SPLICED_MITO, in the format 00_build_exon_profiles.py
reads (reference/raw/<acc>_<gene>[_<locus>].fasta + junctions.tsv rows).

The nad1/nad2/nad5/rps3 references were fetched ad hoc before this script
existed; they're left exactly as they are (this only writes the cis genes,
and skips any accession+gene that already has a file unless --force).

Verification, per locus: translate the joined exons from the genome and
compare with RefSeq's own curated /translation. The two must be the same
length, and every differing residue must be explainable by C->U RNA
editing of the genomic codon (the RefSeq protein is the edited one). A
locus that fails is still written but marked UNRESOLVED in junctions.tsv,
which 00_build_exon_profiles.py excludes from training.

GenBank records are cached in reference/genbank/.
"""
from __future__ import annotations

import argparse
import csv
import itertools
import sys
import time
import urllib.request
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

sys.path.insert(0, str(Path(__file__).resolve().parent))
from reference_genes import CIS_SPLICED_MITO, GENE_SYNONYMS, REFERENCE_MITOGENOMES  # noqa: E402

REFERENCE_DIR = Path(__file__).resolve().parents[1] / "reference"
RAW_DIR = REFERENCE_DIR / "raw"
GB_DIR = REFERENCE_DIR / "genbank"
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gbwithparts&retmode=text&id={}"
CIS_MAX_GAP = 10_000
JUNCTION_COLS = ["accession", "gene", "locus_tag", "exon_num", "n_exons", "start", "end", "strand",
                 "junction_to_next", "gap_bp", "notes"]


def fetch(acc: str) -> Path:
    GB_DIR.mkdir(parents=True, exist_ok=True)
    path = GB_DIR / f"{acc}.gb"
    if not path.exists() or path.stat().st_size == 0:
        with urllib.request.urlopen(EFETCH.format(acc), timeout=120) as r:
            path.write_bytes(r.read())
        time.sleep(0.4)  # stay under NCBI's 3 requests/s without an API key
    return path


def editing_explains(genomic_codon: str, edited_aa: str) -> bool:
    """Can C->T changes at some subset of the codon's C positions give edited_aa?"""
    cs = [i for i, b in enumerate(genomic_codon) if b == "C"]
    for k in range(1, len(cs) + 1):
        for subset in itertools.combinations(cs, k):
            c = list(genomic_codon)
            for i in subset:
                c[i] = "T"
            if str(Seq("".join(c)).translate()) == edited_aa:
                return True
    return False


def verify(dna: str, refseq_protein: str) -> tuple[bool, str]:
    dna = dna.upper()
    codons = [dna[i:i + 3] for i in range(0, len(dna) - len(dna) % 3, 3)]
    naive = str(Seq("".join(codons)).translate())
    naive_core = naive[:-1] if naive.endswith("*") else naive
    if len(naive_core) == len(refseq_protein) + 1 and editing_explains(codons[len(refseq_protein)], "*"):
        # the stop codon is only created by editing (genomic CGA/CAA/CAG ->
        # UGA/UAA/UAG), so the unedited genome reads one residue further
        naive_core = naive_core[:-1]
    if len(naive_core) != len(refseq_protein):
        return False, f"UNRESOLVED length naive {len(naive_core)} vs RefSeq {len(refseq_protein)}"
    unexplained, edited = 0, 0
    for i, (a, b) in enumerate(zip(naive_core, refseq_protein)):
        if a == b:
            continue
        if editing_explains(codons[i], b) or (i == 0 and b == "M"):  # ACG->AUG starts included
            edited += 1
        else:
            unexplained += 1
    if unexplained:
        return False, f"UNRESOLVED {unexplained} residue(s) not explained by C-to-U editing"
    return True, f"naive vs RefSeq differ only at {edited} C-to-U editing-consistent site(s)"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--force", action="store_true", help="rewrite accession+gene files that already exist")
    args = ap.parse_args()

    junctions_path = RAW_DIR / "junctions.tsv"
    with open(junctions_path) as fh:
        existing_rows = list(csv.DictReader(fh, delimiter="\t"))
    done = {(r["accession"], r["gene"]) for r in existing_rows}
    wanted = {g.lower(): g for g in CIS_SPLICED_MITO}
    new_rows, n_ok, n_bad = [], 0, 0

    for acc, organism in REFERENCE_MITOGENOMES.items():
        rec = SeqIO.read(fetch(acc), "genbank")
        seen: dict[str, int] = {}
        for f in rec.features:
            if f.type != "CDS" or "gene" not in f.qualifiers or "translation" not in f.qualifiers:
                continue
            name = f.qualifiers["gene"][0].lower()
            gene = wanted.get(name) or (GENE_SYNONYMS.get(name) if GENE_SYNONYMS.get(name) in CIS_SPLICED_MITO else None)
            if gene is None or ((acc, gene) in done and not args.force):
                continue
            parts = list(f.location.parts)
            if "".join(str(p.extract(rec.seq)) for p in parts) != str(f.extract(rec.seq)):
                parts = parts[::-1]  # put exons in transcript order
            exons = [str(p.extract(rec.seq)).upper() for p in parts]
            ok, note = verify("".join(exons), f.qualifiers["translation"][0])
            locus = f.qualifiers.get("locus_tag", [f"{gene}_{seen.get(gene, 0) + 1}"])[0]
            seen[gene] = seen.get(gene, 0) + 1
            n_ok += ok
            n_bad += not ok

            k = len(parts)
            suffix = f"_{locus}" if seen[gene] > 1 or sum(
                1 for g in rec.features if g.type == "CDS" and g.qualifiers.get("gene", [""])[0].lower() == name) > 1 else ""
            lines = []
            for n, (p, seq) in enumerate(zip(parts, exons), start=1):
                strand = "-" if p.strand == -1 else "+"
                lines.append(f">{acc}|{gene}|{locus}|exon{n}of{k}|{int(p.start) + 1}-{int(p.end)}|{strand}\n{seq}\n")
                if n < k:
                    q = parts[n]
                    same = (q.strand == p.strand)
                    gap = (int(q.start) - int(p.end)) if p.strand != -1 else (int(p.start) - int(q.end))
                    kind = "cis" if same and 0 < gap <= CIS_MAX_GAP else "trans"
                else:
                    gap, kind = "", "none"
                new_rows.append({"accession": acc, "gene": gene, "locus_tag": locus, "exon_num": n, "n_exons": k,
                                 "start": int(p.start) + 1, "end": int(p.end), "strand": strand,
                                 "junction_to_next": kind, "gap_bp": gap, "notes": f"{organism}: {note}"})
            lines.append(f">{acc}|{gene}|{locus}|translation_REFSEQ_ANNOTATED|len={len(f.qualifiers['translation'][0])}\n"
                         f"{f.qualifiers['translation'][0]}\n")
            (RAW_DIR / f"{acc}_{gene}{suffix}.fasta").write_text("".join(lines))
            print(f"[{'ok' if ok else 'UNRESOLVED'}] {acc} {organism}: {gene} ({locus}) {k} exon(s) - {note}",
                  file=sys.stderr)

    if args.force:
        refreshed = {(r["accession"], r["gene"]) for r in new_rows}
        existing_rows = [r for r in existing_rows if (r["accession"], r["gene"]) not in refreshed]
    with open(junctions_path, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=JUNCTION_COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(existing_rows + new_rows)
    print(f"[info] fetch_references: {n_ok} verified, {n_bad} UNRESOLVED loci written -> {RAW_DIR}", file=sys.stderr)


if __name__ == "__main__":
    main()
