#!/usr/bin/env python3
"""Force a linear representation of every assembly with `gfatk linear`.

This is a *universal, cheap, graph-only* linearization pass - runs on every
species x organelle that has a GFA at all (not just QC-eligible ones),
takes seconds each, needs no raw reads. It's a genuinely independent second
opinion on path resolution from oatk's own Pathfinder: `gfatk linear`
greedily picks the highest-cumulative-coverage path through each subgraph,
with no gene-annotation weighting, so it can and does disagree with oatk
(verified on Acaena_ovalifolia: oatk calls it circular via its own
annotation-aware Pathfinder, gfatk linear picks a *linear* 2-node path and
leaves a 112kb second segment as an orphan - not a bug, just a cruder
heuristic). Treat this as a second data point to compare against
data/*/*.ctg.fasta, not a replacement for it.

For `no_resolved_ctg_fasta` species specifically - where oatk's Pathfinder
produced nothing at all - this is often the *only* linearization available
short of the much heavier read-evidence-based 02/03 (gfatk resolve) tier.

Output parsing note (verified against real multi-subgraph output): within
each subgraph, `gfatk linear -e` emits exactly one `gfatk_linear:path=...`
record (its chosen best path) plus zero or more bare `>u<N>` records - every
other node in that subgraph that didn't make it into the chosen path,
dumped individually. A subgraph with only one node ever (no internal choice
to make) has no `path=` record at all; its lone node record *is* the
linearization for that subgraph. This script keeps exactly one output
sequence per subgraph, applying that rule, and reports how much of each
subgraph's sequence got left out (n_orphan_nodes/orphan_bp) as a QC signal
on how much this heuristic actually explains.

Output:
  results/linear/<organelle>/<species>.linear.fasta   - one record per subgraph
  results/linear_summary.tsv                          - one row per species x organelle x subgraph
"""
from __future__ import annotations

import os
import argparse
import re
import subprocess
import sys
from pathlib import Path

# CODE_DIR is this checkout's analysis/ (shared modules, bundled reference
# data); data/ and every module's results/ and work/ live under the data
# root - PLANT_ORGANELLE_DATA_ROOT, or this checkout if unset.
CODE_DIR = Path(__file__).resolve().parents[2]
ROOT_DIR = Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT", str(CODE_DIR.parent)))
ANALYSIS_DIR = ROOT_DIR / "analysis"
sys.path.insert(0, str(CODE_DIR / "common"))

import pandas as pd  # noqa: E402

import io_utils as iou  # noqa: E402
import species_discovery as sd  # noqa: E402

GFATK = __import__("shutil").which("gfatk") or str(Path.home() / ".cargo" / "bin" / "gfatk")

HEADER_RE = re.compile(
    r"^(?:gfatk_linear:path=(?P<path>[^:]+):coverage=(?P<coverage>\d+)\s+)?"
    r"(?P<node>\S+)?"
    r"(?:\s+subgraph-(?P<subgraph>\d+):is_circular-(?P<circular>true|false))?$"
)

COLUMNS = ["species", "organelle", "subgraph_id", "path", "n_segments", "coverage",
           "is_circular", "length_bp", "n_orphan_nodes", "orphan_bp", "pct_subgraph_explained"]


def run_gfatk_linear(gfa: str, node_threshold: int, timeout: int) -> tuple[list[dict], str]:
    """Return (records, stderr). Each record: {header, seq}."""
    cmd = [GFATK, "linear", "-e", "-n", str(node_threshold), gfa]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    if proc.returncode != 0:
        return [], proc.stderr.strip()

    records = []
    header, seq_parts = None, []
    for line in proc.stdout.splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append({"header": header, "seq": "".join(seq_parts)})
            header, seq_parts = line[1:].strip(), []
        else:
            seq_parts.append(line.strip())
    if header is not None:
        records.append({"header": header, "seq": "".join(seq_parts)})
    return records, proc.stderr.strip()


def summarize_species(species: str, organelle: str, records: list[dict]) -> tuple[list[dict], list[dict]]:
    """Group records by subgraph, apply the one-resolved-record-per-subgraph
    rule, return (fasta_records_to_write, summary_rows)."""
    by_subgraph: dict[str, list[dict]] = {}
    for rec in records:
        m = HEADER_RE.match(rec["header"])
        if not m:
            print(f"[warn] {species} ({organelle}): unparsed gfatk linear header: {rec['header']!r}",
                  file=sys.stderr)
            continue
        # gfatk linear short-circuits for a GFA with only a single segment overall
        # ("Only a single segment detected. Printing sequence and exiting.") and
        # omits the "subgraph-N:is_circular-..." suffix entirely in that case -
        # trivially subgraph 1, circularity simply not reported.
        sub_id = m.group("subgraph") or "1"
        by_subgraph.setdefault(sub_id, []).append((m, rec))

    fasta_out, rows = [], []
    for sub_id, entries in sorted(by_subgraph.items(), key=lambda kv: int(kv[0])):
        resolved = [(m, rec) for m, rec in entries if m.group("path")]
        orphans = [(m, rec) for m, rec in entries if not m.group("path")]

        if resolved:
            m, rec = resolved[0]
            path = m.group("path")
            n_segments = len(path.split(","))
            coverage = int(m.group("coverage"))
        elif len(entries) == 1:
            # trivial single-node subgraph: its own node name stands in for a path.
            m, rec = entries[0]
            path = f"{m.group('node')}+"
            n_segments = 1
            coverage = None
            orphans = []
        else:
            print(f"[warn] {species} ({organelle}) subgraph {sub_id}: no resolved path among "
                  f"{len(entries)} records, skipping", file=sys.stderr)
            continue

        length_bp = len(rec["seq"])
        orphan_bp = sum(len(o_rec["seq"]) for _, o_rec in orphans)
        pct_explained = length_bp / (length_bp + orphan_bp) if (length_bp + orphan_bp) else None
        circular = m.group("circular")  # None only for the single-segment-overall special case

        fasta_header = (f"{species}.{organelle}.subgraph{sub_id} path={path} "
                         f"coverage={coverage if coverage is not None else 'NA'} "
                         f"circular={circular if circular is not None else 'unknown'} length={length_bp}")
        fasta_out.append({"header": fasta_header, "seq": rec["seq"]})
        rows.append({
            "species": species, "organelle": organelle, "subgraph_id": sub_id, "path": path,
            "n_segments": n_segments, "coverage": coverage, "is_circular": circular,
            "length_bp": length_bp, "n_orphan_nodes": len(orphans), "orphan_bp": orphan_bp,
            "pct_subgraph_explained": round(pct_explained, 4) if pct_explained is not None else None,
        })
    return fasta_out, rows


def write_fasta(path: Path, records: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w") as fh:
        for rec in records:
            fh.write(f">{rec['header']}\n")
            seq = rec["seq"]
            for i in range(0, len(seq), 70):
                fh.write(seq[i:i + 70] + "\n")
    tmp.replace(path)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--organelle", default="both", choices=["mito", "pltd", "both"])
    ap.add_argument("--species-list", default=None)
    ap.add_argument("--node-threshold", type=int, default=20000,
                     help="skip (sub)graphs with more nodes than this (default 20000)")
    ap.add_argument("--timeout", type=int, default=120)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    if not Path(GFATK).exists():
        print(f"[err] gfatk not found at {GFATK}", file=sys.stderr)
        sys.exit(1)

    organelles = ["mito", "pltd"] if args.organelle == "both" else [args.organelle]
    species_filter = sd.load_species_list(Path(args.species_list)) if args.species_list else None
    linear_dir = ANALYSIS_DIR / "linearize" / "results" / "linear"
    summary_path = ANALYSIS_DIR / "linearize" / "results" / "linear_summary.tsv"

    existing_df, existing_keys = iou.load_existing(summary_path, ["species", "organelle"])
    new_rows: list[dict] = []
    recomputed_keys: set[tuple] = set()
    n_ok = n_skipped = n_no_gfa = n_error = 0

    for organelle in organelles:
        data_root = sd.repo_data_root(ROOT_DIR, organelle)
        resolved = sd.discover_all(data_root, organelle, species_filter=species_filter)
        out_dir = linear_dir / organelle

        for r in resolved:
            if r.gfa is None:
                n_no_gfa += 1
                continue

            fasta_path = out_dir / f"{r.species}.linear.fasta"
            key = (r.species, organelle)
            if not iou.needs_recompute(key, existing_keys, fasta_path, Path(r.gfa), args.force):
                n_skipped += 1
                continue

            records, err = run_gfatk_linear(r.gfa, args.node_threshold, args.timeout)
            recomputed_keys.add(key)
            if not records:
                n_error += 1
                print(f"[warn] {r.species} ({organelle}): gfatk linear produced nothing"
                      + (f" ({err[:200]})" if err else ""), file=sys.stderr)
                continue

            fasta_out, rows = summarize_species(r.species, organelle, records)
            if not fasta_out:
                n_error += 1
                continue
            write_fasta(fasta_path, fasta_out)
            new_rows.extend(rows)
            n_ok += 1

    # A species' subgraph count can change between runs (that's the point of
    # rerunning), so we can't key-match individual rows for merge - instead,
    # drop ALL existing rows for every (species, organelle) recomputed this
    # run, then append the freshly computed rows for those. Untouched species
    # (e.g. outside --species-list scope) are left exactly as they were.
    if existing_df.empty:
        final_rows = new_rows
    else:
        recompute_idx = pd.MultiIndex.from_tuples(recomputed_keys, names=["species", "organelle"])
        keep = existing_df[~existing_df.set_index(["species", "organelle"]).index.isin(recompute_idx)]
        final_rows = keep.to_dict("records") + new_rows
    iou.atomic_write_tsv(summary_path, final_rows, COLUMNS)

    print(f"[info] gfatk_linear: ok={n_ok} skipped={n_skipped} no_gfa={n_no_gfa} error={n_error} "
          f"-> {linear_dir}, {summary_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
