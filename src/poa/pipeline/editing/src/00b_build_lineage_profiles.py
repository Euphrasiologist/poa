#!/usr/bin/env python3
"""Build the bundled editing profiles for a non-angiosperm lineage
(common/lineage.py), from that lineage's RefSeq mitogenomes rather than
from a dataset.

00_build_reference_profiles.py builds profiles from a dataset's own gene
calls, so on a mostly-angiosperm dataset they are angiosperm profiles. A
moss aligned against them gets "corrected" towards the angiosperm protein:
98% of moss mito calls changed, against 37% for dicots, although mosses
edit very little (Physcomitrella has ~11 mito sites). These profiles are
built from RefSeq's curated proteins instead, for the genes the bundled
mito profiles cover, so a moss is compared against mosses.

Per gene: every single-exon CDS in the lineage's mitogenomes
(trans_splicing/src/reference_genes.LINEAGE_MITOGENOMES; GenBank cached in
trans_splicing/reference/genbank/) whose genomic translation verifies
against RefSeq's protein (common/refseq_mito.verify). A gene that is
multi-exon in most of the lineage (moss nad5, nad9) gets no profile: a
single-locus correction is the wrong model for it, and 01_scan_editing.py
leaves it uncorrected for that lineage rather than borrow the angiosperm
one.

Output: editing/reference/lineages/<lineage>/mito/<gene>.pssm (bundled)
"""
from __future__ import annotations

import argparse
import importlib.util
import os
import sys
from collections import Counter
from pathlib import Path

from Bio import SeqIO

CODE_DIR = Path(__file__).resolve().parents[2]
ROOT_DIR = Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT") or os.getcwd())
sys.path.insert(0, str(CODE_DIR / "common"))
sys.path.insert(0, str(CODE_DIR / "trans_splicing" / "src"))
import lineage as lin  # noqa: E402
from refseq_mito import cds_exons, fetch, verify  # noqa: E402
from reference_genes import GENE_SYNONYMS, LINEAGE_MITOGENOMES  # noqa: E402

# the profile maths is 00_build_reference_profiles.py's, so the two kinds of
# profile are built identically
_spec = importlib.util.spec_from_file_location("build_reference_profiles",
                                               Path(__file__).with_name("00_build_reference_profiles.py"))
brp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(brp)

GB_DIR = lin.trans_splicing_reference(CODE_DIR, lin.DEFAULT_LINEAGE) / "genbank"
BUNDLED_MITO = CODE_DIR / "editing" / "reference" / "profiles" / "mito"
SYNONYMS = {**GENE_SYNONYMS, "tatc": "mttB", "mttb": "mttB"}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lineage", required=True, choices=[x for x in lin.LINEAGES if x != lin.DEFAULT_LINEAGE])
    args = ap.parse_args()

    genes = sorted(p.stem for p in BUNDLED_MITO.glob("*.pssm"))
    wanted = {g.lower(): g for g in genes}
    proteins: dict[str, dict[str, str]] = {g: {} for g in genes}
    n_exons: dict[str, Counter] = {g: Counter() for g in genes}
    for acc, organism in LINEAGE_MITOGENOMES[args.lineage].items():
        rec = SeqIO.read(fetch(acc, GB_DIR), "genbank")
        for f in rec.features:
            if f.type != "CDS" or "gene" not in f.qualifiers or "translation" not in f.qualifiers:
                continue
            name = f.qualifiers["gene"][0].lower()
            gene = wanted.get(name) or wanted.get(SYNONYMS.get(name, "").lower())
            if gene is None or acc in proteins[gene]:
                continue
            parts, exons = cds_exons(rec, f)
            n_exons[gene][len(parts)] += 1
            if len(parts) != 1:
                continue
            ok, note = verify(exons[0], f.qualifiers["translation"][0])
            if not ok:
                print(f"[warn] {acc} {organism} {gene}: {note} - not used", file=sys.stderr)
                continue
            proteins[gene][acc] = f.qualifiers["translation"][0]

    work_dir = ROOT_DIR / "analysis" / "editing" / "work" / args.lineage
    work_dir.mkdir(parents=True, exist_ok=True)
    out_dir = lin.editing_profiles(CODE_DIR, args.lineage) / "mito"
    for gene in genes:
        modal = n_exons[gene].most_common(1)[0][0] if n_exons[gene] else None
        stale = out_dir / f"{gene}.pssm"
        if modal != 1 or len(proteins[gene]) < brp.MIN_SEQUENCES:
            print(f"[info] {gene}: no {args.lineage} profile (exon counts {dict(n_exons[gene])}, "
                  f"{len(proteins[gene])} single-exon verified loci)", file=sys.stderr)
            stale.unlink(missing_ok=True)
            continue
        protein_fasta = work_dir / f"{gene}.protein.fasta"
        protein_fasta.write_text("".join(f">{acc}\n{seq}\n" for acc, seq in proteins[gene].items()))
        aligned = work_dir / f"{gene}.aligned.fasta"
        if not brp.run_mafft(protein_fasta, aligned):
            sys.exit(f"[err] mafft failed on {protein_fasta}")
        columns = brp.build_profile(brp.read_fasta(aligned))
        brp.write_profile(columns, out_dir / f"{gene}.pssm")
        print(f"[info] {gene}: {args.lineage} profile from {len(proteins[gene])} RefSeq proteins, "
              f"{len(columns)} columns -> {out_dir / gene}.pssm", file=sys.stderr)


if __name__ == "__main__":
    main()
