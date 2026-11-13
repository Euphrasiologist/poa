"""RefSeq mitogenome helpers shared by the reference builders
(trans_splicing/src/00a_fetch_references.py,
editing/src/00b_build_lineage_profiles.py): fetching a GenBank record into
a cache, pulling out a gene's CDS exons in transcript order, and checking
the genomic translation against RefSeq's curated protein.

Verification: the joined exons' translation must be the same length as
RefSeq's /translation, and every differing residue must be explainable by
C->U RNA editing of the genomic codon (the RefSeq protein is the edited
one).
"""

from __future__ import annotations

import itertools
import time
import urllib.request
from pathlib import Path

from Bio.Seq import Seq

EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gbwithparts&retmode=text&id={}"


def fetch(acc: str, cache_dir: Path) -> Path:
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / f"{acc}.gb"
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
    codons = [dna[i : i + 3] for i in range(0, len(dna) - len(dna) % 3, 3)]
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


def cds_exons(rec, feature) -> tuple[list, list[str]]:
    """A CDS feature's location parts and exon sequences, in transcript order."""
    parts = list(feature.location.parts)
    if "".join(str(p.extract(rec.seq)) for p in parts) != str(feature.extract(rec.seq)):
        parts = parts[::-1]
    return parts, [str(p.extract(rec.seq)).upper() for p in parts]
