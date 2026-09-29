#!/usr/bin/env python3
"""Self-vs-self blastn on each species' resolved contigs to find repeat copies,
clustered into repeat families with cd-hit-est.

Reimplemented fresh (not ported) from the approach used in the sibling
mito_structural_variation repo's repeats/src/old/recomb_repeats.bash: a
clean, per-species, organelle-only template (self-blastn -> filter trivial/
mirrored hits -> cd-hit-est cluster). Organelle contigs (tens of kb to
~1.5 Mb) make full self-blastn cheap - no genome-scale masking/indexing
needed.

Output: analysis/repeats/results/per_species/<Species>.<organelle>.repeats.tsv
Columns: contig_id start end repeat_family copy_length pct_identity strand
(strand is always '+' here - it describes the interval on the contig's own
coordinate system, not orientation relative to a partner copy; a segment's
orientation relative to any one partner can be direct or inverted, which is
exactly what makes flip-flop/recombinogenic repeats interesting, but isn't
a property of the segment in isolation.)

02_flag_recomb_repeats.py adds contig-end-distance/putative_recomb_repeat
columns on top of this.
"""
from __future__ import annotations

import os
import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

# CODE_DIR is this checkout's analysis/ (shared modules, bundled reference
# data); data/ and every module's results/ and work/ live under the data
# root - PLANT_ORGANELLE_DATA_ROOT, or this checkout if unset.
CODE_DIR = Path(__file__).resolve().parents[2]
ROOT_DIR = Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT", str(CODE_DIR.parent)))
ANALYSIS_DIR = ROOT_DIR / "analysis"
sys.path.insert(0, str(CODE_DIR / "common"))

import io_utils  # noqa: E402
import species_discovery as sd  # noqa: E402

BLASTN = shutil.which("blastn") or "/software/team301/ncbi-blast-2.16.0+/bin/blastn"
SAMTOOLS = shutil.which("samtools") or "/software/team301/samtools/samtools"
CDHIT_EST = shutil.which("cd-hit-est") or "/software/team301/cdhit/cd-hit-est"

COLUMNS = ["contig_id", "start", "end", "repeat_family", "copy_length", "pct_identity", "strand"]
BLAST_FIELDS = ["qseqid", "sseqid", "pident", "length", "qstart", "qend", "sstart", "send"]


def self_blastn(fasta_path: str) -> pd.DataFrame:
    cmd = [BLASTN, "-query", fasta_path, "-subject", fasta_path,
           "-word_size", "50", "-reward", "1", "-penalty", "-4",
           "-evalue", "1e-3", "-dust", "no", "-soft_masking", "false",
           "-outfmt", "6 " + " ".join(BLAST_FIELDS)]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    if proc.returncode != 0 or not proc.stdout.strip():
        return pd.DataFrame(columns=BLAST_FIELDS)
    df = pd.read_csv(pd.io.common.StringIO(proc.stdout), sep="\t", names=BLAST_FIELDS)
    for c in ("pident", "length", "qstart", "qend", "sstart", "send"):
        df[c] = pd.to_numeric(df[c])
    return df


def normalize_segment(seqid: str, a: int, b: int) -> tuple[str, int, int]:
    return (seqid, min(a, b), max(a, b))


def filter_hits(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df
    q = df.apply(lambda r: normalize_segment(r.qseqid, r.qstart, r.qend), axis=1)
    s = df.apply(lambda r: normalize_segment(r.sseqid, r.sstart, r.send), axis=1)
    df = df.assign(q_seg=q, s_seg=s)
    # drop the trivial whole-sequence-vs-itself hit
    df = df[df["q_seg"] != df["s_seg"]]
    if df.empty:
        return df
    # dedupe mirrored pairs (A vs B and B vs A): canonical order by sorting the two segment tuples
    df["pair_key"] = df.apply(lambda r: tuple(sorted([r.q_seg, r.s_seg])), axis=1)
    df = df.drop_duplicates(subset="pair_key")
    return df


def cluster_segments(fasta_path: str, segments: list[tuple[str, int, int]], work_dir: Path) -> dict:
    """Extract each segment, cluster with cd-hit-est, return {segment: family_id}."""
    if not segments:
        return {}
    ids = [f"{c}:{s}-{e}" for c, s, e in segments]
    regions = [f"{c}:{s + 1}-{e}" for c, s, e in segments]  # 1-based inclusive for samtools
    subprocess.run([SAMTOOLS, "faidx", fasta_path], capture_output=True, timeout=60)
    proc = subprocess.run([SAMTOOLS, "faidx", fasta_path, *regions], capture_output=True, text=True, timeout=120)
    if proc.returncode != 0:
        return {}
    records = [chunk.split("\n", 1)[1].replace("\n", "") for chunk in proc.stdout.split(">")[1:]]
    if len(records) != len(segments):
        return {}

    with tempfile.TemporaryDirectory(dir=work_dir) as tmpdir:
        in_fasta = Path(tmpdir) / "segments.fasta"
        out_prefix = Path(tmpdir) / "clustered"
        with open(in_fasta, "w") as fh:
            for seg_id, seq in zip(ids, records):
                fh.write(f">{seg_id}\n{seq}\n")
        cmd = [CDHIT_EST, "-i", str(in_fasta), "-o", str(out_prefix),
               "-c", "0.90", "-n", "8", "-d", "0", "-M", "0", "-T", "1", "-l", "20"]
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        if proc.returncode != 0:
            return {}
        clstr_path = Path(f"{out_prefix}.clstr")
        if not clstr_path.exists():
            return {}
        seg_to_family = {}
        cluster_id = -1
        for line in clstr_path.read_text().splitlines():
            if line.startswith(">Cluster"):
                cluster_id += 1
            else:
                seg_id = line.split(">")[1].split("...")[0]
                seg_to_family[seg_id] = cluster_id
    return seg_to_family


def merge_overlapping_segments(segments: list[tuple[str, int, int]]) -> list[tuple[str, int, int]]:
    """Collapse overlapping/nested repeat-copy intervals per contig into their
    maximal union block. blastn reports many overlapping/re-extended HSPs
    for what is really one repeat locus; without this a single repeat
    region shows up as dozens of near-duplicate rows differing by a few bp.
    """
    from collections import defaultdict
    by_contig: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for c, s, e in segments:
        by_contig[c].append((s, e))
    merged = []
    for c, intervals in by_contig.items():
        intervals.sort()
        cur_s, cur_e = intervals[0]
        for s, e in intervals[1:]:
            if s <= cur_e:
                cur_e = max(cur_e, e)
            else:
                merged.append((c, cur_s, cur_e))
                cur_s, cur_e = s, e
        merged.append((c, cur_s, cur_e))
    return merged


def rows_for_species(species: str, organelle: str, fasta_path: str, work_dir: Path) -> list[dict]:
    hits = filter_hits(self_blastn(fasta_path))
    if hits.empty:
        return []

    raw_segments = sorted(set(hits["q_seg"]) | set(hits["s_seg"]))
    segments = merge_overlapping_segments(raw_segments)
    seg_to_family = cluster_segments(fasta_path, segments, work_dir)
    if not seg_to_family:
        return []

    # map each merged block -> max pident among the raw hit segments it now contains
    best_pident: dict[tuple, float] = {seg: 0.0 for seg in segments}
    for _, r in hits.iterrows():
        for raw_seg in (r.q_seg, r.s_seg):
            for merged_seg in segments:
                if raw_seg[0] == merged_seg[0] and raw_seg[1] >= merged_seg[1] and raw_seg[2] <= merged_seg[2]:
                    best_pident[merged_seg] = max(best_pident[merged_seg], r.pident)
                    break

    rows = []
    for seg in segments:
        seg_id = f"{seg[0]}:{seg[1]}-{seg[2]}"
        family = seg_to_family.get(seg_id)
        if family is None:
            continue
        rows.append({
            "contig_id": seg[0], "start": seg[1], "end": seg[2],
            "repeat_family": f"fam{family}", "copy_length": seg[2] - seg[1],
            "pct_identity": round(best_pident.get(seg, 0.0), 2), "strand": "+",
        })
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organelle", default="both", choices=["mito", "pltd", "both"])
    ap.add_argument("--species-list", default=None)
    ap.add_argument("--qc-status", default="pass")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    qc_path = ANALYSIS_DIR / "qc_basic_stats" / "results" / "qc_summary.tsv"
    qc = pd.read_csv(qc_path, sep="\t")
    qc_status = args.qc_status.split(",")
    eligible = set(zip(qc.loc[qc["status"].isin(qc_status), "species"], qc.loc[qc["status"].isin(qc_status), "organelle"]))

    organelles = ["mito", "pltd"] if args.organelle == "both" else [args.organelle]
    species_filter = sd.load_species_list(Path(args.species_list)) if args.species_list else None

    out_dir = ANALYSIS_DIR / "repeats" / "results" / "per_species"
    out_dir.mkdir(parents=True, exist_ok=True)
    work_dir = ANALYSIS_DIR / "repeats" / "work"
    work_dir.mkdir(parents=True, exist_ok=True)

    n_computed = n_skipped = 0
    for organelle in organelles:
        data_root = sd.repo_data_root(ROOT_DIR, organelle)
        resolved = sd.discover_all(data_root, organelle, species_filter=species_filter)
        for r in resolved:
            if (r.species, organelle) not in eligible or r.status != "ok":
                continue
            out_path = out_dir / f"{r.species}.{organelle}.repeats.tsv"
            if not args.force and out_path.exists() and Path(r.ctg_fasta).stat().st_mtime <= out_path.stat().st_mtime:
                n_skipped += 1
                continue
            rows = rows_for_species(r.species, organelle, r.ctg_fasta, work_dir)
            io_utils.atomic_write_tsv(out_path, rows, COLUMNS)
            n_computed += 1

    print(f"[info] find_repeats: computed={n_computed} skipped={n_skipped} -> {out_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
