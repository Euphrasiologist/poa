#!/usr/bin/env python3
"""Promote gfatk-resolve results (analysis/linearize/results/resolve/) into
data/{mito,plastid}/<species>/, so qc_basic_stats' has_ctg_fasta/status
actually reflects that these species now have a real resolved genome.

Every target here started as no_resolved_ctg_fasta (either a missing/
empty ctg.fasta, or - the unjoined subtype - one contig per raw graph
segment, i.e. Pathfinder joined nothing) - see qc_basic_stats/README.md.

Two gates before a circuit is promoted (added 2026-09-28, after "resolved
a circuit -> promote" replaced 44 plastids whose Pathfinder genome had
merely gone missing - 7 of them with a 34-55 kb piece of a 136-230 kb
genome; see 07_restore_pathfinder.py):
  1. oatk's own result first: if re-running pathfinder on the GFA + oatk's
     annotation reproduces oatk's original .ctg.bed with a resolved genome,
     don't promote - run 07_restore_pathfinder.py to put that back.
  2. coverage: the circuit must place >= GATE (plastid 90%, mito 70%) of
     the graph's unique unitig content (common/graph_coverage.py).

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

import os
import argparse
import re
import shutil
import sys
from datetime import datetime
from pathlib import Path

import importlib.util
import tempfile

import pandas as pd

# CODE_DIR is the installed poa/pipeline/ (shared modules, bundled reference
# data); data/ and every module's results/ and work/ live under the data
# root - PLANT_ORGANELLE_DATA_ROOT, or the current directory if unset.
CODE_DIR = Path(__file__).resolve().parents[2]
ROOT_DIR = Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT") or os.getcwd())
ANALYSIS_DIR = ROOT_DIR / "analysis"
sys.path.insert(0, str(CODE_DIR / "common"))

import graph_coverage  # noqa: E402
import species_discovery as sd  # noqa: E402

_spec = importlib.util.spec_from_file_location("restore_pathfinder", Path(__file__).with_name("07_restore_pathfinder.py"))
restore_pathfinder = importlib.util.module_from_spec(_spec)
sys.modules["restore_pathfinder"] = restore_pathfinder
_spec.loader.exec_module(restore_pathfinder)

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
            current = species_dir / f"{prefix}{SUFFIXES_TO_ARCHIVE[organelle]['ctg_fasta']}"
            if current.exists() and current.stat().st_size > 0:
                headers = [l for l in open(current) if l.startswith(">")]
                skip = None
                if "resolver=gfatk_resolve" in headers[0]:
                    # re-archiving would overwrite oatk's original in superseded/
                    skip = "already_promoted"
                elif any("circular=true" in h or not re.search(r"\bnv=1\b", h) for h in headers):
                    skip = "resolved_genome_present"  # e.g. restored by 07_restore_pathfinder.py
                if skip:
                    rows.append({"species": species, "organelle": organelle, "promoted": False, "reason": skip,
                                 "n_circuits": len(circuits), "total_bp": total_bp, "annotation_status": None})
                    continue
            with tempfile.TemporaryDirectory() as tmp:
                regen = restore_pathfinder.plan_one(species, organelle, existing_gfas[0], Path(tmp))
                circuit_fa = Path(tmp) / "circuits.fasta"
                circuit_fa.write_text("".join(f">c{i}\n{c['seq']}\n" for i, c in enumerate(circuits)))
                placed = graph_coverage.unitig_content_placed(existing_gfas[0], circuit_fa)
            if regen["action"] == "restore":
                rows.append({"species": species, "organelle": organelle, "promoted": False,
                             "reason": "pathfinder_genome_recoverable_run_07", "n_circuits": len(circuits),
                             "total_bp": total_bp, "annotation_status": None})
                continue
            if placed < restore_pathfinder.GATE[organelle]:
                rows.append({"species": species, "organelle": organelle, "promoted": False,
                             "reason": f"circuit_places_{placed:.0%}_of_graph", "n_circuits": len(circuits),
                             "total_bp": total_bp, "annotation_status": None})
                continue
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
    if not args.dry_run:  # a dry run changed nothing, so it isn't a promotion record
        df.to_csv(log_path, sep="\t", index=False)

    n_promoted = int(df["promoted"].sum()) if not df.empty else 0
    print(df.groupby(["organelle", "promoted", "reason"]).size() if not df.empty else "no resolve outputs found",
          file=sys.stderr)
    print(f"[info] {'(dry run) ' if args.dry_run else ''}promoted={n_promoted}/{len(df)} -> {log_path}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
