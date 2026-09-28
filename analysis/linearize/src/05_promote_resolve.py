#!/usr/bin/env python3
"""Promote gfatk-resolve results (analysis/linearize/results/resolve/) into
data/{mito,plastid}/<species>/, so qc_basic_stats' has_ctg_fasta/status
actually reflects that these species now have a real resolved genome.

Every target here started as no_resolved_ctg_fasta (either a genuinely
empty ctg.fasta, or - the unjoined subtype - one contig per raw graph
segment, i.e. Pathfinder joined nothing) - see qc_basic_stats/README.md.
Any gfatk-resolve circuit is therefore an unambiguous improvement over
"nothing usable"; there's no dead-ends-vs-subgraphs trade-off to weigh
like src/09_promote_revision.py, so promotion here is simply: resolved a
circuit -> promote.

What this does NOT do, deliberately: write a new .ctg.bed/.annot_*.txt.
gfatk resolve's output is a bare sequence with no gene calls - oatk's own
.annot_*.txt (raw nhmmscan tblout) -> .bed/.ctg.bed conversion is internal
to its Rust binary, and reimplementing that translation from outside
oatk risks silently getting gene coordinates/strand/names wrong. Rather
than fabricate annotation, promoted species get a real .ctg.fasta and
correctly-updated has_ctg_fasta/status, but core_gene_pct stays unknown
(NaN, not penalized, not guessed) until a real re-annotation step exists -
every promoted row is logged with annotation_status=not_yet_annotated so
this gap stays visible, not silently papered over.

The existing .gfa is left untouched (it's the same graph gfatk resolve
itself read - still accurate) - only .ctg.fasta/.ctg.bed/.bed/.annot_*.txt
at the SAME existing prefix are touched (old ones archived to superseded/
if present, never deleted).

Usage:
  analysis/linearize/src/05_promote_resolve.py [--dry-run]
Output: appends to meta/promotion_log_<timestamp>.tsv (same convention as
src/09_promote_revision.py).
"""
from __future__ import annotations

import argparse
import re
import shutil
import os
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd

# Defaults to the old same-tree layout (code and data siblings under one
# repo root); set PLANT_ORGANELLE_DATA_ROOT to point at a separate data
# checkout instead (e.g. when this code is installed from its own repo).
ROOT_DIR = Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT", str(Path(__file__).resolve().parents[3])))
ANALYSIS_DIR = ROOT_DIR / "analysis"
sys.path.insert(0, str(ANALYSIS_DIR / "common"))

import species_discovery as sd  # noqa: E402

SUFFIXES_TO_ARCHIVE = {
    "mito": {"ctg_bed": ".mito.ctg.bed", "ctg_fasta": ".mito.ctg.fasta", "bed": ".mito.bed", "annot": ".annot_mito.txt"},
    "pltd": {"ctg_bed": ".pltd.ctg.bed", "ctg_fasta": ".pltd.ctg.fasta", "bed": ".pltd.bed", "annot": ".annot_pltd.txt"},
}
GFA_SUFFIX = {"mito": ".mito.gfa", "pltd": ".pltd.gfa"}
DATA_DIRNAME = {"mito": "mito", "pltd": "plastid"}

CIRCUIT_RE = re.compile(r"^.+_circuit(?P<circuit>\d+):n_segments=(?P<n_segments>\d+):bp=(?P<bp>\d+)$")
BUBBLE_ARM_RE = re.compile(r"^bubble_arm:")


def parse_circuits(fasta_path: Path) -> list[dict]:
    """Return [{header, seq, n_segments, bp}], circuits only (bubble arms excluded)."""
    if not fasta_path.exists() or fasta_path.stat().st_size == 0:
        return []
    circuits = []
    header, seq_parts, meta = None, [], None
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None and meta is not None:
                    circuits.append({"header": header, "seq": "".join(seq_parts), **meta})
                header, seq_parts = line[1:].strip(), []
                if BUBBLE_ARM_RE.match(header):
                    meta = None
                    continue
                m = CIRCUIT_RE.match(header)
                meta = {"n_segments": int(m.group("n_segments")), "bp": int(m.group("bp"))} if m else None
            else:
                seq_parts.append(line.strip())
        if header is not None and meta is not None:
            circuits.append({"header": header, "seq": "".join(seq_parts), **meta})
    return circuits


def write_promoted_fasta(path: Path, circuits: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w") as fh:
        for i, c in enumerate(circuits, start=1):
            header = (f"ctg{i:06d}c length={c['bp']} nv={c['n_segments']} circular=true "
                      f"resolver=gfatk_resolve source_header={c['header']}")
            fh.write(f">{header}\n")
            seq = c["seq"]
            for j in range(0, len(seq), 70):
                fh.write(seq[j:j + 70] + "\n")
    tmp.replace(path)


def archive_old_files(species_dir: Path, prefix: str, organelle: str) -> list[str]:
    archived = []
    archive_dir = species_dir / "superseded"
    for key, suf in SUFFIXES_TO_ARCHIVE[organelle].items():
        f = species_dir / f"{prefix}{suf}"
        if f.exists():
            archive_dir.mkdir(exist_ok=True)
            shutil.move(str(f), str(archive_dir / f.name))
            archived.append(key)
    return archived


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--resolve-dir", default=str(ANALYSIS_DIR / "linearize" / "results" / "resolve"))
    ap.add_argument("--organelle", default="both", choices=["mito", "pltd", "both"])
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    resolve_dir = Path(args.resolve_dir)
    organelles = ["mito", "pltd"] if args.organelle == "both" else [args.organelle]

    rows = []
    for organelle in organelles:
        data_root = ROOT_DIR / "data" / DATA_DIRNAME[organelle]
        org_resolve_dir = resolve_dir / organelle
        if not org_resolve_dir.is_dir():
            continue
        for fasta_path in sorted(org_resolve_dir.glob("*.resolve.fasta")):
            species = fasta_path.name[: -len(".resolve.fasta")]
            circuits = parse_circuits(fasta_path)
            if not circuits:
                rows.append({"species": species, "organelle": organelle, "promoted": False,
                             "reason": "no_circuit_resolved", "n_circuits": 0, "total_bp": None,
                             "annotation_status": None})
                continue

            species_dir = data_root / species
            gfa_suffix = GFA_SUFFIX[organelle]
            existing_gfas = list(species_dir.glob(f"*{gfa_suffix}")) if species_dir.is_dir() else []
            if not existing_gfas:
                rows.append({"species": species, "organelle": organelle, "promoted": False,
                             "reason": "no_existing_gfa_to_pair_with", "n_circuits": len(circuits),
                             "total_bp": sum(c["bp"] for c in circuits), "annotation_status": None})
                continue
            prefix = existing_gfas[0].name[: -len(gfa_suffix)]

            total_bp = sum(c["bp"] for c in circuits)
            if not args.dry_run:
                archived = archive_old_files(species_dir, prefix, organelle)
                new_ctg_fasta = species_dir / f"{prefix}{SUFFIXES_TO_ARCHIVE[organelle]['ctg_fasta']}"
                write_promoted_fasta(new_ctg_fasta, circuits)
            else:
                archived = None

            rows.append({"species": species, "organelle": organelle, "promoted": True, "reason": "circuit_resolved",
                         "n_circuits": len(circuits), "total_bp": total_bp,
                         "annotation_status": "not_yet_annotated", "archived_files": archived})

    df = pd.DataFrame(rows)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_path = ROOT_DIR / "meta" / f"promotion_log_{ts}.tsv"
    df.to_csv(log_path, sep="\t", index=False)

    n_promoted = int(df["promoted"].sum()) if not df.empty else 0
    print(df.groupby(["organelle", "promoted", "reason"]).size() if not df.empty else "no resolve outputs found",
          file=sys.stderr)
    print(f"[info] {'(dry run) ' if args.dry_run else ''}promoted={n_promoted}/{len(df)} -> {log_path}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
