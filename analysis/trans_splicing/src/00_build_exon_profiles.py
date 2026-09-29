#!/usr/bin/env python3
"""Build per-exon-slot and whole-gene reference profiles for `transsplice`.

Source: `../reference/raw/*.fasta` + `junctions.tsv` - real GenBank RefSeq
mitochondrial genomes (Arabidopsis, Beta vulgaris, maize, rice, Vitis),
fetched and independently re-verified (every exon re-translated from its
own raw coordinates and checked against the curated RefSeq protein - see
`../README.md`), plus, for the cis-spliced genes (reference_genes.py),
six more angiosperm RefSeq mitogenomes fetched and verified the same way
by 00a_fetch_references.py. Genes: reference_genes.GENES.

Two kinds of profile per gene, both written under `../reference/profiles/
<gene>/`:
- `exon{N}.pssm` + `exon{N}.threshold`: one per exon slot (N = 1..K),
  built from the loci that have a canonical (modal) exon count - a locus
  missing an intron entirely (e.g. Beta vulgaris `rps3`) or flagged
  UNRESOLVED (rice `rps3` - a likely assembly artifact, confirmed not a
  transcription error on our end) is simply absent from a slot's training
  set rather than distorting it. `.threshold` is FLOOR_FRACTION of the
  minimum score of each slot's own training sequences run back through
  `orfedit scan-batch` against the profile just built from them - the
  same scoring code path `transsplice` uses at runtime.
- `exon{N}.frame`: `lead<TAB>trail` - how many bases at the exon's 5'
  end finish a codon begun in the previous exon, and how many at its 3'
  end start one the next exon finishes (the modal value across training
  loci). Exon profiles are built from the exon's *full codons only*,
  translated in their true frame: translating every exon from its own
  first base put every exon after a phase-1/2 junction out of frame
  (nad1 exons 2/5, nad2 exons 3-5, nad5 exons 2/4, rps3 exon 2), and
  joining exons on codon boundaries then dropped the split codon's bases,
  frame-shifting the rest of the gene - transsplice adds lead/trail back.
- `whole_gene.pssm`: one full-length profile per gene, from every
  included locus's curated RefSeq protein (regardless of its own exon
  count - Beta vulgaris' intron-less `rps3` protein is still a real,
  correct full-length `rps3`, just not useful for the per-exon slots).
"""
from __future__ import annotations

import csv
import math
import os
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

# CODE_DIR is this checkout's analysis/ (shared modules, bundled reference
# data); data/ and every module's results/ and work/ live under the data
# root - PLANT_ORGANELLE_DATA_ROOT, or this checkout if unset.
CODE_DIR = Path(__file__).resolve().parents[2]
ROOT_DIR = Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT", str(CODE_DIR.parent)))
ANALYSIS_DIR = ROOT_DIR / "analysis"
REFERENCE_DIR = CODE_DIR / "trans_splicing" / "reference"
RAW_DIR = REFERENCE_DIR / "raw"
PROFILES_DIR = REFERENCE_DIR / "profiles"
WORK_DIR = ANALYSIS_DIR / "trans_splicing" / "work"

MAFFT = shutil.which("mafft") or "/software/team301/mafft-7.525-with-extensions/core/mafft"
MAFFT_BINARIES_DIR = os.environ.get("MAFFT_BINARIES", "/software/team301/mafft-7.525-with-extensions/core")
ORFEDIT = shutil.which("orfedit") or str(Path.home() / ".cargo" / "bin" / "orfedit")

AA_ORDER = "ARNDCQEGHILKMFPSTWYV"
PSEUDOCOUNT = 0.5
GAP_COLUMN_THRESHOLD = 0.9  # same reasoning/value as editing/00_build_reference_profiles.py
sys.path.insert(0, str(Path(__file__).resolve().parent))
from reference_genes import GENES  # noqa: E402  (shared with 00a/01/gff_export)
# Slot floor = this fraction of the lowest training self-score. With exons
# translated in frame, self-scores are high (nad5 exon 2: 1013), and a
# floor at the bare minimum rejected even Arabidopsis's own nad5 exons 1-2
# (a training genome) whenever the candidate window didn't span the whole
# exon. 0.5 recovered them with no wrong-locus assignment it didn't
# already have. Negative self-scores (tiny exons) are kept as-is.
FLOOR_FRACTION = 0.5

CODON_TABLE = {}
_bases = "TCAG"
_aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, b1 in enumerate(_bases):
    for j, b2 in enumerate(_bases):
        for k, b3 in enumerate(_bases):
            CODON_TABLE[b1 + b2 + b3] = _aas[i * 16 + j * 4 + k]


def translate(seq: str) -> str:
    seq = seq.upper()
    codons = [seq[i : i + 3] for i in range(0, len(seq) - len(seq) % 3, 3)]
    return "".join(CODON_TABLE.get(c, "X") for c in codons)


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    name = None
    chunks: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(chunks)
            name = line[1:].strip()
            chunks = []
        elif line.strip():
            chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def load_all_raw_data() -> dict[tuple[str, str, str], dict]:
    """{(accession, gene, locus_tag): {exon_num_int: dna_seq, "refseq": protein_seq}}"""
    data: dict[tuple[str, str, str], dict] = defaultdict(dict)
    for fasta_path in RAW_DIR.glob("*.fasta"):
        for header, seq in read_fasta(fasta_path).items():
            parts = header.split("|")
            accession, gene, locus_tag, record_type = parts[0], parts[1], parts[2], parts[3]
            key = (accession, gene, locus_tag)
            if record_type.startswith("exon"):
                n = int(record_type[len("exon") :].split("of")[0])
                data[key][n] = seq
            elif record_type == "translation_REFSEQ_ANNOTATED":
                data[key]["refseq"] = seq
    return data


def load_junctions() -> dict[tuple[str, str, str], dict]:
    """{(accession, gene, locus_tag): {"n_exons": int, "excluded": bool, "reason": str}}"""
    out: dict[tuple[str, str, str], dict] = {}
    with open(REFERENCE_DIR / "raw" / "junctions.tsv") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            key = (row["accession"], row["gene"], row["locus_tag"])
            if key in out:
                continue
            excluded = "UNRESOLVED" in row["notes"]
            out[key] = {"n_exons": int(row["n_exons"]), "excluded": excluded}
    return out


def run_mafft(protein_fasta: Path, out_fasta: Path) -> bool:
    # a minimal env, but keep TMPDIR: mafft's mktemp otherwise needs a writable /tmp
    env = ({"MAFFT_BINARIES": MAFFT_BINARIES_DIR, "PATH": "/usr/bin:/bin"}
           if Path(MAFFT_BINARIES_DIR).is_dir() else dict(os.environ))
    if os.environ.get("TMPDIR"):
        env["TMPDIR"] = os.environ["TMPDIR"]
    proc = subprocess.run(
        [MAFFT, "--auto", "--quiet", str(protein_fasta)],
        capture_output=True, text=True, timeout=120, env=env,
    )
    if proc.returncode != 0 or not proc.stdout.strip():
        print(f"[warn] mafft failed for {protein_fasta.name}: {proc.stderr.strip()[:300]}", file=sys.stderr)
        return False
    out_fasta.write_text(proc.stdout)
    return True


def build_profile(alignment: dict[str, str]) -> list[dict[str, float]] | None:
    n = len(alignment)
    if n == 0:
        return None
    aln_len = len(next(iter(alignment.values())))
    bg_counts = Counter()
    for seq in alignment.values():
        bg_counts.update(c for c in seq if c in AA_ORDER)
    bg_total = sum(bg_counts.values()) or 1
    background = {aa: (bg_counts.get(aa, 0) + PSEUDOCOUNT) / (bg_total + PSEUDOCOUNT * 20) for aa in AA_ORDER}

    columns = []
    for i in range(aln_len):
        col_chars = [seq[i] for seq in alignment.values()]
        gap_frac = col_chars.count("-") / n
        if gap_frac > GAP_COLUMN_THRESHOLD:
            continue
        counts = Counter(c for c in col_chars if c in AA_ORDER)
        total = sum(counts.values()) or 1
        scores = {}
        for aa in AA_ORDER:
            freq = (counts.get(aa, 0) + PSEUDOCOUNT) / (total + PSEUDOCOUNT * 20)
            scores[aa] = math.log2(freq / background[aa])
        columns.append(scores)
    return columns or None


def write_profile(columns: list[dict[str, float]], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(AA_ORDER)]
    for col in columns:
        lines.append("\t".join(f"{col[aa]:.4f}" for aa in AA_ORDER))
    out_path.write_text("\n".join(lines) + "\n")


def build_and_write_profile(protein_seqs: dict[str, str], out_path: Path, label: str) -> bool:
    if len(protein_seqs) < 2:
        print(f"[warn] {label}: only {len(protein_seqs)} sequences, need >= 2 - skipping", file=sys.stderr)
        return False
    protein_fasta = WORK_DIR / f"{label}.protein.fasta"
    protein_fasta.write_text("".join(f">{name}\n{seq}\n" for name, seq in protein_seqs.items()))
    aligned_fasta = WORK_DIR / f"{label}.aligned.fasta"
    if not run_mafft(protein_fasta, aligned_fasta):
        return False
    alignment = read_fasta(aligned_fasta)
    columns = build_profile(alignment)
    if columns is None:
        print(f"[warn] {label}: no usable alignment columns - skipping", file=sys.stderr)
        return False
    write_profile(columns, out_path)
    return True


def compute_min_self_score(gene: str, slot: int, dna_seqs: dict[str, str], profile_path: Path) -> float | None:
    """Score every training sequence for this slot back against the profile
    just built from it, via a real `orfedit scan-batch` call (same code
    path `transsplice` uses at runtime) - the floor is the minimum of these,
    not an invented cutoff."""
    manifest_rows = []
    for name, seq in dna_seqs.items():
        fasta_path = WORK_DIR / f"selfscore.{gene}.exon{slot}.{name}.fasta"
        fasta_path.write_text(f">self\n{seq}\n")
        manifest_rows.append({
            "species": name, "organelle": "mito", "gene": f"{gene}_exon{slot}",
            "fasta": str(fasta_path), "profile": str(profile_path),
            "contig": "self", "start": 0, "end": len(seq), "strand": "+",
        })
    manifest_path = WORK_DIR / f"selfscore.{gene}.exon{slot}.manifest.tsv"
    cols = ["species", "organelle", "gene", "fasta", "profile", "contig", "start", "end", "strand"]
    with open(manifest_path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in manifest_rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    out_path = WORK_DIR / f"selfscore.{gene}.exon{slot}.out.tsv"
    edits_path = WORK_DIR / f"selfscore.{gene}.exon{slot}.edits.tsv"
    proc = subprocess.run([
        ORFEDIT, "scan-batch", "--manifest", str(manifest_path),
        "--gene-calls-out", str(out_path), "--edits-out", str(edits_path),
        "--threads", "1", "--flank", "0",
    ], capture_output=True, text=True, timeout=120)
    if proc.returncode != 0:
        print(f"[warn] {gene} exon{slot}: self-scoring orfedit call failed: {proc.stderr.strip()[:300]}", file=sys.stderr)
        return None
    scores = []
    with open(out_path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            scores.append(float(row["score"]))
    return min(scores) if scores else None


def main() -> None:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes", nargs="+", default=GENES,
                    help="build only these (default: all). Rebuilding a gene overwrites its whole_gene.pssm, "
                         "including one 02_bootstrap_refine.py refined from this dataset")
    args = ap.parse_args()
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    PROFILES_DIR.mkdir(parents=True, exist_ok=True)

    raw = load_all_raw_data()
    junctions = load_junctions()

    for gene in args.genes:
        loci = [key for key in raw if key[1] == gene]
        included = [key for key in loci if not junctions.get(key, {}).get("excluded", False)]
        excluded = [key for key in loci if key not in included]
        if excluded:
            print(f"[info] {gene}: excluding {[k[2] for k in excluded]} (flagged in junctions.tsv)", file=sys.stderr)

        gene_dir = PROFILES_DIR / gene

        # Whole-gene profile: every included locus's curated RefSeq protein,
        # regardless of its own exon count.
        whole_seqs = {}
        for key in included:
            acc, _, locus = key
            if "refseq" in raw[key]:
                whole_seqs[f"{acc}_{locus}"] = raw[key]["refseq"]
        build_and_write_profile(whole_seqs, gene_dir / "whole_gene.pssm", f"{gene}.whole_gene")

        # Per-exon-slot profiles: only loci with the canonical (modal) exon
        # count contribute a given slot - naturally excludes Beta vulgaris'
        # intron-less rps3 from exon2 (it simply has no exon2 record).
        n_exons_seen = Counter(junctions[key]["n_exons"] for key in included if key in junctions)
        if not n_exons_seen:
            print(f"[warn] {gene}: no junction data for included loci - skipping exon slots", file=sys.stderr)
            continue
        k = n_exons_seen.most_common(1)[0][0]
        print(f"[info] {gene}: {len(included)} loci, canonical exon count K={k} "
              f"(seen: {dict(n_exons_seen)})", file=sys.stderr)

        # lead/trail per (locus, slot) from the locus's own cumulative exon
        # lengths - the junction phase, which the translation must respect
        frames: dict[tuple, dict[int, tuple[int, int]]] = {}
        for key in included:
            if junctions.get(key, {}).get("n_exons") != k or not all(n in raw[key] for n in range(1, k + 1)):
                continue
            cum, per = 0, {}
            for n in range(1, k + 1):
                length = len(raw[key][n])
                lead = (3 - cum % 3) % 3
                trail = (cum + length) % 3
                per[n] = (lead, trail)
                cum += length
            frames[key] = per

        for slot in range(1, k + 1):
            dna_seqs = {}
            for key in included:
                acc, _, locus = key
                # Require this locus's OWN exon count to match the
                # canonical K, not just "has a record at this slot index" -
                # an intron-less locus (e.g. Beta vulgaris rps3, n_exons=1)
                # stores its single whole-CDS record as "exon1of1", which
                # would otherwise corrupt exon1's profile with a ~551aa
                # full-length sequence instead of the real ~24aa leader
                # exon (caught empirically: min_self_score was -483 before
                # this check).
                if junctions.get(key, {}).get("n_exons") != k:
                    continue
                if slot in raw[key] and key in frames:
                    lead, trail = frames[key][slot]
                    seq = raw[key][slot]
                    dna_seqs[f"{acc}_{locus}"] = seq[lead:len(seq) - trail]
            if len(dna_seqs) < 2:
                print(f"[warn] {gene} exon{slot}: only {len(dna_seqs)} sequences - skipping", file=sys.stderr)
                continue
            protein_seqs = {name: translate(seq).rstrip("*") for name, seq in dna_seqs.items()}
            profile_path = gene_dir / f"exon{slot}.pssm"
            if not build_and_write_profile(protein_seqs, profile_path, f"{gene}.exon{slot}"):
                continue
            min_score = compute_min_self_score(gene, slot, dna_seqs, profile_path)
            if min_score is None:
                print(f"[warn] {gene} exon{slot}: could not compute min_self_score - skipping", file=sys.stderr)
                continue
            floor = min_score * FLOOR_FRACTION if min_score > 0 else min_score
            (gene_dir / f"exon{slot}.threshold").write_text(f"{floor:.4f}\n")
            lead, trail = Counter(frames[key][slot] for key in frames).most_common(1)[0][0]
            (gene_dir / f"exon{slot}.frame").write_text(f"{lead}\t{trail}\n")
            print(f"[info] {gene} exon{slot}: {len(dna_seqs)} training seqs, "
                  f"min_self_score={min_score:.3f} -> {profile_path.name}", file=sys.stderr)


if __name__ == "__main__":
    main()
