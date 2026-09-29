"""Small shared minimap2/PAF helpers, used by qc_basic_stats (a coarse
mito-vs-plastid contamination pre-screen) and synteny (the full MTPT/dotplot
module). Keeping this in one place means both modules agree on what
"total aligned length" means.
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

MINIMAP2 = shutil.which("minimap2") or "/software/team301/minimap2/minimap2"

PAF_COLS = ["query_name", "query_length", "query_start", "query_end", "strand",
            "target_name", "target_length", "target_start", "target_end",
            "n_matches", "aln_block_length", "mapq"]


def run_asm5(ref_fasta: str | Path, query_fasta: str | Path, timeout: int = 120) -> list[dict]:
    """Run `minimap2 -x asm5 ref query` and return parsed PAF records.

    Matches minimap2's own positional convention: first fasta = reference/
    target (indexed), second = query. PAF "query_*" columns describe
    query_fasta, "target_*" columns describe ref_fasta.
    """
    proc = subprocess.run([MINIMAP2, "-x", "asm5", str(ref_fasta), str(query_fasta)],
                           capture_output=True, text=True, timeout=timeout)
    if proc.returncode != 0:
        return []
    records = []
    for line in proc.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 12:
            continue
        rec = dict(zip(PAF_COLS, fields[:12]))
        for k in ("query_length", "query_start", "query_end", "target_length",
                  "target_start", "target_end", "n_matches", "aln_block_length", "mapq"):
            rec[k] = int(rec[k])
        records.append(rec)
    return records


def total_aligned_length(records: list[dict], side: str = "query") -> int:
    """Length of the query or target sequence actually covered by at least
    one alignment block, merging overlapping/duplicate blocks per sequence
    name first.

    Naively summing block lengths without merging was verified to produce
    nonsensical "fraction aligned" values above 1.0 for real species in
    this dataset (e.g. a repeated element in one organelle matching several
    loci in the other, each reported as a separate PAF block covering the
    same query region) - i.e. this isn't a rare edge case at real dataset
    scale, it was inflating exactly the high-alignment cases that matter
    most for the contamination check that uses this function.
    """
    name_key = "query_name" if side == "query" else "target_name"
    end_key = "query_end" if side == "query" else "target_end"
    start_key = "query_start" if side == "query" else "target_start"

    by_name: dict[str, list[tuple[int, int]]] = {}
    for r in records:
        by_name.setdefault(r[name_key], []).append((r[start_key], r[end_key]))

    total = 0
    for intervals in by_name.values():
        intervals.sort()
        cur_start, cur_end = intervals[0]
        for start, end in intervals[1:]:
            if start <= cur_end:
                cur_end = max(cur_end, end)
            else:
                total += cur_end - cur_start
                cur_start, cur_end = start, end
        total += cur_end - cur_start
    return total
