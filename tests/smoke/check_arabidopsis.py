#!/usr/bin/env python3
"""Check an Arabidopsis thaliana mito GFF from gff_export against RefSeq.

    check_arabidopsis.py --gff <data_root>/analysis/gff_export/results/Arabidopsis_thaliana.mito.gff \
                         --fasta <data_root>/data/mito/Arabidopsis_thaliana/<run>.mito.ctg.fasta

The assembly is a different Arabidopsis accession from RefSeq NC_037304.1
(bundled at trans_splicing/reference/genbank/), so coordinates are not
comparable; gene content, exon structure and protein sequence are.

  - gene content: every RefSeq CDS/tRNA/rRNA gene name (normalised:
    tRNA anticodons in RNA alphabet, ccmFN1/ccmFN2 -> ccmFn) present as a
    GFF gene.
  - per protein gene: the GFF model's CDS (exons in exon_number order,
    each on its own strand) is extracted from the assembly and translated;
    exon count and per-exon length are compared with RefSeq's, and the
    protein is globally aligned (BLOSUM62) to RefSeq's curated /translation.
    That protein is the edited one, so a mismatch counts as explained when
    C-to-U editing of the genomic codon gives RefSeq's residue (the rule
    00a_fetch_references.py verifies references with); `identity` counts
    those as matches, `raw_identity` doesn't.
    `annotation_tier=raw` models are nhmmscan spans with approximate
    boundaries, not reading frames, so for those only the best of three
    frames is scored (right locus, not right model) and no exon/identity
    threshold applies. Where the GFF has several copies of a gene, the
    best-matching one counts.

Writes a per-gene TSV (--out) and exits 1 if any EXPECT threshold fails.
"""
from __future__ import annotations

import argparse
import csv
import itertools
import re
import sys
from collections import defaultdict
from pathlib import Path

from Bio import Align, SeqIO
from Bio.Align import substitution_matrices
from Bio.Seq import Seq

CODE_DIR = Path(__file__).resolve().parents[2] / "analysis"
REFSEQ_GB = CODE_DIR / "trans_splicing" / "reference" / "genbank" / "NC_037304.1.gb"

# Thresholds sit just under the dataset repo's results at 9d627214, so
# regressions fail and improvements pass. min_identity is editing-aware
# protein identity vs RefSeq; max_offset bounds the largest per-exon length
# difference from RefSeq (1 = a split-codon boundary placed one base off).
MIN_GENES_FOUND = 51  # of 51 normalised RefSeq names
MIN_EDITED_IDENTITY = 0.93  # median over edited-tier (orfedit-corrected) single-exon genes
EXPECT = {
    # trans-spliced (reconstructed tier since stage 1)
    "nad1": {"n_exons": 5, "min_identity": 0.90},
    "nad2": {"n_exons": 5, "max_offset": 1, "min_identity": 0.97},
    "nad5": {"n_exons": 5, "max_offset": 1, "min_identity": 0.97},
    "rps3": {"n_exons": 2, "max_offset": 0, "min_identity": 0.98},
    # cis-spliced (stage 2a, transsplice 0.2.1): base-exact per trans_splicing/README.md
    "nad4": {"n_exons": 4, "max_offset": 0, "min_identity": 0.97},
    "nad7": {"n_exons": 5, "max_offset": 0, "min_identity": 0.97},
}


def norm(name: str) -> str:
    m = re.match(r"^(?:trn|tRNA)([A-Za-z]+)\(([ACGTU]{3})\)$", name)
    if m:
        aa = "fM" if m.group(1) == "fM" else m.group(1)
        return f"trn{aa}-{m.group(2).replace('T', 'U')}".lower()
    if re.match(r"^ccmFN\d$", name, re.I):
        return "ccmfn"
    return name.lower()


def editing_explains(codon: str, aa: str) -> bool:
    """Can C->T changes at some of the codon's C positions give aa?"""
    cs = [i for i, b in enumerate(codon) if b == "C"]
    for k in range(1, len(cs) + 1):
        for subset in itertools.combinations(cs, k):
            c = list(codon)
            for i in subset:
                c[i] = "T"
            if str(Seq("".join(c)).translate()) == aa:
                return True
    return False


def refseq_genes(gb: Path) -> tuple[set[str], dict[str, dict]]:
    rec = SeqIO.read(gb, "genbank")
    names, proteins = set(), {}
    for f in rec.features:
        if f.type not in ("CDS", "tRNA", "rRNA") or "gene" not in f.qualifiers:
            continue
        name = norm(f.qualifiers["gene"][0])
        names.add(name)
        if f.type == "CDS" and "translation" in f.qualifiers:
            parts = list(f.location.parts)
            if "".join(str(p.extract(rec.seq)) for p in parts) != str(f.extract(rec.seq)):
                parts = parts[::-1]  # transcript order
            proteins[name] = {"exon_lengths": [len(p) for p in parts],
                              "protein": f.qualifiers["translation"][0]}
    return names, proteins


def gff_models(gff: Path) -> tuple[set[str], dict[str, list[tuple[str, list[tuple]]]]]:
    """Gene names, and per gene name a list of (tier, CDS rows) models."""
    names, gene_tier, mrna_gene, cds = set(), {}, {}, defaultdict(list)
    for line in open(gff):
        if line.startswith("#") or not line.strip():
            continue
        seqid, _, ftype, start, end, _, strand, _, attrs = line.rstrip("\n").split("\t")
        a = dict(kv.split("=", 1) for kv in attrs.split(";") if "=" in kv)
        if ftype == "gene":
            names.add(norm(a["Name"]))
            gene_tier[a["ID"]] = a.get("annotation_tier", "raw")
        elif ftype == "mRNA":
            mrna_gene[a["ID"]] = (norm(a["Name"]), gene_tier.get(a["Parent"], "raw"))
        elif ftype == "CDS":
            cds[a["Parent"]].append((int(a.get("exon_number", 1)), seqid, int(start), int(end), strand))
    models = defaultdict(list)
    for mrna, rows in cds.items():
        if mrna in mrna_gene:
            name, tier = mrna_gene[mrna]
            models[name].append((tier, sorted(rows)))
    return names, models


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gff", required=True, type=Path)
    ap.add_argument("--fasta", required=True, type=Path, help="the assembly the GFF annotates (.ctg.fasta)")
    ap.add_argument("--genbank", type=Path, default=REFSEQ_GB)
    ap.add_argument("--out", type=Path, help="per-gene TSV")
    args = ap.parse_args()

    genome = {r.id: r.seq for r in SeqIO.parse(args.fasta, "fasta")}
    ref_names, ref_prot = refseq_genes(args.genbank)
    gff_names, models = gff_models(args.gff)

    aligner = Align.PairwiseAligner(mode="global", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                                    open_gap_score=-10, extend_gap_score=-0.5)

    def identity(cds: str, ref: str) -> tuple[float, float]:
        """(editing-aware, raw) identity of cds's translation vs ref."""
        cds = cds[: len(cds) - len(cds) % 3].upper()
        protein = str(Seq(cds).translate()).rstrip("*")
        if not protein:
            return 0.0, 0.0
        aln = aligner.align(protein, ref)[0]
        raw = edited = 0
        qi = 0
        for x, y in zip(aln[0], aln[1]):
            if x != "-":
                if x == y:
                    raw += 1
                    edited += 1
                elif y != "-" and editing_explains(cds[3 * qi:3 * qi + 3], y):
                    edited += 1
                qi += 1
        n = max(len(protein), len(ref))
        return edited / n, raw / n

    rows, failures = [], []
    for gene, ref in sorted(ref_prot.items()):
        best = None
        for tier, model in models.get(gene, []):
            seq = Seq("")
            for _, seqid, start, end, strand in model:
                piece = genome[seqid][start - 1:end]
                seq += piece.reverse_complement() if strand == "-" else piece
            frames = [0] if tier != "raw" else [0, 1, 2]
            ident, raw_ident = max(identity(str(seq[f:]), ref["protein"]) for f in frames)
            lengths = [end - start + 1 for _, _, start, end, _ in model]
            ref_lengths = list(ref["exon_lengths"])
            ref_lengths[-1] -= 3  # RefSeq CDS includes the stop codon; compare without it
            # the last exon may or may not carry the stop codon
            offsets = [abs(x - y) if i < len(lengths) - 1 else min(abs(x - y), abs(x - y - 3))
                       for i, (x, y) in enumerate(zip(lengths, ref_lengths))]
            same_count = len(lengths) == len(ref_lengths)
            exact = same_count and max(offsets) == 0
            cand = {"gene": gene, "tier": tier, "n_models": len(models[gene]), "gff_n_exons": len(lengths),
                    "refseq_n_exons": len(ref_lengths), "gff_exon_lengths": ",".join(map(str, lengths)),
                    "refseq_exon_lengths": ",".join(map(str, ref["exon_lengths"])),
                    "exact_exons": exact, "max_exon_offset": max(offsets) if same_count else "",
                    "protein_identity": round(ident, 4), "raw_identity": round(raw_ident, 4)}
            if best is None or (cand["exact_exons"], ident) > (best["exact_exons"], best["protein_identity"]):
                best = cand
        rows.append(best or {"gene": gene, "tier": "missing", "n_models": 0, "gff_n_exons": 0,
                             "refseq_n_exons": len(ref["exon_lengths"]), "gff_exon_lengths": "",
                             "refseq_exon_lengths": ",".join(map(str, ref["exon_lengths"])),
                             "exact_exons": False, "max_exon_offset": "", "protein_identity": 0.0,
                             "raw_identity": 0.0})

    missing = sorted(ref_names - gff_names)
    found = len(ref_names) - len(missing)
    edited = sorted(r["protein_identity"] for r in rows if r["refseq_n_exons"] == 1 and r["tier"] == "edited")
    median_edited = edited[len(edited) // 2] if edited else 0.0

    by_gene = {r["gene"]: r for r in rows}
    if found < MIN_GENES_FOUND:
        failures.append(f"gene content: {found}/{len(ref_names)} found, expected >= {MIN_GENES_FOUND}")
    if median_edited < MIN_EDITED_IDENTITY:
        failures.append(f"edited-tier protein identity median {median_edited:.3f} < {MIN_EDITED_IDENTITY}")
    for gene, want in EXPECT.items():
        r = by_gene.get(gene)
        if r is None:
            failures.append(f"{gene}: not in RefSeq record")
            continue
        if "max_offset" in want and (r["max_exon_offset"] == "" or r["max_exon_offset"] > want["max_offset"]):
            failures.append(f"{gene}: exons {r['gff_exon_lengths'] or '-'} vs RefSeq {r['refseq_exon_lengths']}")
        if "n_exons" in want and r["gff_n_exons"] != want["n_exons"]:
            failures.append(f"{gene}: {r['gff_n_exons']} exons, expected {want['n_exons']}")
        if "min_identity" in want and r["protein_identity"] < want["min_identity"]:
            failures.append(f"{gene}: protein identity {r['protein_identity']:.3f} < {want['min_identity']}")

    if args.out:
        with open(args.out, "w") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)

    print(f"gene content: {found}/{len(ref_names)} RefSeq genes found" + (f"; missing: {', '.join(missing)}" if missing else ""))
    print(f"edited-tier single-exon genes: median identity {median_edited:.3f} over {len(edited)}")
    print(f"{'gene':8} {'tier':>13} {'exons gff/ref':>13} {'offset':>6} {'identity':>9} {'raw':>6}")
    for r in rows:
        print(f"{r['gene']:8} {r['tier']:>13} {r['gff_n_exons']:>6}/{r['refseq_n_exons']:<6} "
              f"{str(r['max_exon_offset']):>6} {r['protein_identity']:>9.3f} {r['raw_identity']:>6.3f}")
    if failures:
        print("FAIL\n  " + "\n  ".join(failures))
        sys.exit(1)
    print("PASS")


if __name__ == "__main__":
    main()
