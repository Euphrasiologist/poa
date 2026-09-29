#!/usr/bin/env python3
"""Quantify how completely `<species>.<organelle>.ctg.fasta` (oatk's
linearised assembly) accounts for the raw assembly graph it was built
from, using `00_build_unitig_map.py`'s per-unitig placement map.

Two independent completeness signals, both from the unitig map alone (no
gene annotation needed):
  - graph-side: what fraction of total *unique* unitig content (each
    unitig counted once, regardless of copy number) is placed anywhere in
    the ctg at all? Unitigs with zero placements are graph content the
    linearisation dropped entirely.
  - ctg-side: what fraction of the ctg's own length is actually accounted
    for by a placed unitig? (Should be close to 100% when unitig
    placements don't overlap and fully tile the ctg; a low value flags
    either overlap-trimming this exact-match approach doesn't model, or
    real sequence in the ctg that didn't come from any single whole
    unitig.)

`n_multicopy_unitigs` is reported separately, not folded into either
completeness score - a unitig placed more than once is real recombination
/dosage signal (see 00_build_unitig_map.py's docstring), not a defect.

Cross-references `qc_basic_stats/results/qc_summary.tsv`'s own
`any_circular`/`core_gene_pct`/`n_contigs` columns as corroborating
context (not blended into one score - each measures something different,
and conflating them would hide which one actually explains a low value).

Output: results/linearization_qc.tsv - one row per species x organelle.
"""
from __future__ import annotations

import os
import argparse
import csv
import sys
from pathlib import Path

# CODE_DIR is this checkout's analysis/ (shared modules, bundled reference
# data); data/ and every module's results/ and work/ live under the data
# root - PLANT_ORGANELLE_DATA_ROOT, or this checkout if unset.
CODE_DIR = Path(__file__).resolve().parents[2]
ROOT_DIR = Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT", str(CODE_DIR.parent)))
ANALYSIS_DIR = ROOT_DIR / "analysis"
sys.path.insert(0, str(CODE_DIR / "common"))

import species_discovery as sd  # noqa: E402


def fasta_lengths(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    name = None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                lengths[name] = 0
            elif name is not None:
                lengths[name] += len(line.strip())
    return lengths


def load_qc_summary() -> dict[tuple[str, str], dict]:
    qc_path = ANALYSIS_DIR / "qc_basic_stats" / "results" / "qc_summary.tsv"
    if not qc_path.exists():
        return {}
    with open(qc_path) as fh:
        return {(r["species"], r["organelle"]): r for r in csv.DictReader(fh, delimiter="\t")}


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [intervals[0]]
    for start, end in intervals[1:]:
        if start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def evaluate_one(species: str, organelle: str, ctg_fasta: str, map_path: Path,
                 unitig_fasta: Path) -> dict | None:
    if not map_path.exists() or not unitig_fasta.exists():
        return None
    ctg_lengths = fasta_lengths(Path(ctg_fasta))
    ctg_len = sum(ctg_lengths.values())

    # The unitig FASTA is the complete unitig list - the map only has rows
    # for placed unitigs, so orphans are only visible from here.
    unitig_lengths = fasta_lengths(unitig_fasta)
    placements: dict[str, list[tuple[str, int, int]]] = {}
    intervals_by_ctg: dict[str, list[tuple[int, int]]] = {}
    sources = set()
    with open(map_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            u, seqid = row["unitig"], row["ctg_seqid"]
            start, end = int(row["ctg_start"]), int(row["ctg_end"])
            placements.setdefault(u, []).append((seqid, start, end))
            # a trimmed seqmatch hit only verifies its middle (see 00_build_unitig_map.py)
            trim = int(row.get("match_trim") or 0)
            start, end = start + trim, end - trim
            sources.add(row.get("source", ""))
            # A seqmatch hit across a circular join has ctg_end > ctg length
            # - split it at the origin so coverage can't exceed the ctg.
            clen = ctg_lengths.get(seqid, 0)
            ivs = intervals_by_ctg.setdefault(seqid, [])
            if clen and end > clen:
                ivs.append((start, clen))
                ivs.append((0, min(end - clen, clen)))
            else:
                ivs.append((start, end))

    n_unitigs = len(unitig_lengths)
    total_unitig_bp = sum(unitig_lengths.values())
    placed_unitig_bp = sum(unitig_lengths.get(u, 0) for u in placements)
    n_orphan_unitigs = sum(1 for u in unitig_lengths if u not in placements)
    n_multicopy_unitigs = sum(1 for p in placements.values() if len(p) > 1)

    covered_bp = sum(e - s for ivs in intervals_by_ctg.values() for s, e in merge_intervals(ivs))

    return {
        "species": species,
        "organelle": organelle,
        "map_source": ",".join(sorted(s for s in sources if s)) or "none",
        "ctg_length_bp": ctg_len,
        "n_unitigs": n_unitigs,
        "total_unitig_bp": total_unitig_bp,
        "pct_unitig_content_placed": round(100 * placed_unitig_bp / total_unitig_bp, 2) if total_unitig_bp else 0.0,
        "n_orphan_unitigs": n_orphan_unitigs,
        "n_multicopy_unitigs": n_multicopy_unitigs,
        "pct_ctg_explained_by_unitigs": round(100 * covered_bp / ctg_len, 2) if ctg_len else 0.0,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organelle", default="both", choices=["mito", "pltd", "both"])
    ap.add_argument("--species-list", default=None)
    args = ap.parse_args()

    organelles = ["mito", "pltd"] if args.organelle == "both" else [args.organelle]
    species_filter = sd.load_species_list(Path(args.species_list)) if args.species_list else None
    qc_summary = load_qc_summary()
    map_dir = ANALYSIS_DIR / "unitig_coords" / "results" / "unitig_map"
    fasta_dir = ANALYSIS_DIR / "unitig_coords" / "results" / "unitig_fasta"

    rows = []
    for organelle in organelles:
        data_root = sd.repo_data_root(ROOT_DIR, organelle)
        resolved = sd.discover_all(data_root, organelle, species_filter=species_filter)
        for r in resolved:
            if r.status != "ok":
                continue
            map_path = map_dir / f"{r.species}.{organelle}.tsv"
            result = evaluate_one(r.species, organelle, r.ctg_fasta, map_path,
                                  fasta_dir / f"{r.species}.{organelle}.unitig.fasta")
            if result is None:
                continue
            qc = qc_summary.get((r.species, organelle), {})
            result["qc_status"] = qc.get("status", "")
            result["any_circular"] = qc.get("any_circular", "")
            result["core_gene_pct"] = qc.get("core_gene_pct", "")
            result["n_contigs"] = qc.get("n_contigs", "")
            rows.append(result)

    out_path = ANALYSIS_DIR / "unitig_coords" / "results" / "linearization_qc.tsv"
    cols = ["species", "organelle", "map_source", "ctg_length_bp", "n_unitigs", "total_unitig_bp",
            "pct_unitig_content_placed", "n_orphan_unitigs", "n_multicopy_unitigs",
            "pct_ctg_explained_by_unitigs", "qc_status", "any_circular", "core_gene_pct", "n_contigs"]
    with open(out_path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for row in rows:
            fh.write("\t".join(str(row[c]) for c in cols) + "\n")
    print(f"[info] linearization_qc: {len(rows)} rows -> {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
