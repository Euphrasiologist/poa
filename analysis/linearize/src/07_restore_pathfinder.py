#!/usr/bin/env python3
"""Restore oatk Pathfinder's own contigs where they went missing, and undo
gfatk-resolve promotions that replaced a genome with part of one.

Background (2026-09-28): for 31 plastids promoted by `05_promote_resolve.py`
- and 2 never promoted (Salix_viminalis, Trifolium_fragiferum) - oatk had
resolved a complete circular genome, but its `.ctg.fasta` never reached
`data/` (the `.ctg.bed`/`.bed`/`annot_*.txt` did; the cause is not in this
repo's history). QC read that as "no resolved genome", and the gfatk
fallback promoted a circuit in its place - for 7 Cyperaceae/Juncaceae/
grass/Apiaceae plastids a 34-55 kb piece of a 136-230 kb genome.

Re-running `pathfinder` on the untouched GFA + oatk's own annotation
reproduces the original run exactly: its `.ctg.bed` came out byte-identical
to the archived one for every promoted plastid. That identity check is the
condition for restoring - this only ever puts back what oatk produced.

Per candidate (current contig is a gfatk_resolve circuit, or the species
has no/empty `.ctg.fasta`):
  restore  regenerated .ctg.bed identical to oatk's original, and
           resolved: plastid - >=1 circle, keep the highest-coverage one
           only (a second, disjoint circle at 8-40x lower coverage is
           another organism's plastid - logged, not kept); mito - every
           contig circular, keep them all (multichromosomal genomes are
           real).
  revert   current is a gfatk circuit placing < GATE of the graph's unitig
           content, and oatk's original (unjoined) contigs are archived:
           put those back - fragmented but complete beats complete-looking
           but partial.
  keep     otherwise (e.g. mito circuits: Pathfinder only ever gave
           unjoined pieces there, so the circuit is a real improvement).

Replaced files go to `superseded/` (a gfatk circuit as
`<prefix>.<org>.gfatk_resolve.ctg.fasta`); nothing is deleted. Dry run by
default; `--apply` to act. Log: `meta/promotion_log_<timestamp>.tsv`.
"""
from __future__ import annotations

import os
import argparse
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime
from pathlib import Path

# CODE_DIR is this checkout's analysis/ (shared modules, bundled reference
# data); data/ and every module's results/ and work/ live under the data
# root - PLANT_ORGANELLE_DATA_ROOT, or this checkout if unset.
CODE_DIR = Path(__file__).resolve().parents[2]
ROOT_DIR = Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT", str(CODE_DIR.parent)))
ANALYSIS_DIR = ROOT_DIR / "analysis"
sys.path.insert(0, str(CODE_DIR / "common"))

import graph_coverage  # noqa: E402
import species_discovery as sd  # noqa: E402

PATHFINDER = shutil.which("pathfinder") or "/software/team301/oatk/pathfinder"
# Share of the graph's unique unitig bp a gfatk circuit must place to stand.
# Mito graphs carry plastid-derived/foreign low-coverage unitigs, so less is
# expected there (current mito promotions: 75-100%; plastid: 95-100%, bar
# the three partials at 14-58%).
GATE = {"pltd": 0.90, "mito": 0.70}
ORIGINAL_SUFFIXES = {org: {"ctg_fasta": f".{org}.ctg.fasta", "ctg_bed": f".{org}.ctg.bed",
                           "bed": f".{org}.bed", "annot": f".annot_{org}.txt"} for org in ("mito", "pltd")}


def read_records(path: Path) -> list[tuple[str, list[str]]]:
    recs: list[tuple[str, list[str]]] = []
    for line in open(path):
        if line.startswith(">"):
            recs.append((line, []))
        elif recs:
            recs[-1][1].append(line)
    return recs


def header_field(header: str, key: str) -> str | None:
    for tok in header.split():
        if tok.startswith(f"{key}="):
            return tok.split("=", 1)[1]
    return None


def data_lines(path: Path | None) -> list[str] | None:
    return [l for l in open(path) if not l.startswith("#")] if path and path.exists() else None


def find(dirs: list[Path], prefix: str, suffix: str) -> Path | None:
    for d in dirs:
        p = d / f"{prefix}{suffix}"
        if p.exists() and p.stat().st_size > 0:
            return p
    return None


def archive(path: Path, dest_dir: Path, name: str) -> Path:
    dest_dir.mkdir(exist_ok=True)
    dest = dest_dir / name
    if dest.exists():
        dest = dest_dir / f"{name}.{datetime.now():%Y%m%d_%H%M%S}"
    shutil.move(str(path), str(dest))
    return dest


def plan_one(species: str, organelle: str, gfa: Path, tmp: Path) -> dict:
    d = gfa.parent
    sup = d / "superseded"
    prefix = gfa.name[: -len(f".{organelle}.gfa")]
    suf = ORIGINAL_SUFFIXES[organelle]
    cur = d / f"{prefix}{suf['ctg_fasta']}"
    cur_is_gfatk = cur.exists() and cur.stat().st_size > 0 and "resolver=gfatk_resolve" in open(cur).readline()
    cur_missing = not cur.exists() or cur.stat().st_size == 0
    row = {"species": species, "organelle": organelle, "prefix": prefix, "dir": d,
           "current": "gfatk_resolve" if cur_is_gfatk else ("missing" if cur_missing else "oatk")}
    if not (cur_is_gfatk or cur_missing):
        return {**row, "action": "not_a_candidate"}

    # oatk's originals: superseded/ first - 06_reannotate_promoted.py writes
    # the promoted circuit's own annotation into the species dir
    annot = find([sup, d], prefix, suf["annot"])
    orig_bed = find([sup, d], prefix, suf["ctg_bed"])
    if annot is None:
        return {**row, "action": "keep", "reason": "no_oatk_annotation_to_regenerate_from"}

    out = tmp / f"{species}.{organelle}"
    flag = "-p" if organelle == "pltd" else "-m"
    try:
        subprocess.run([PATHFINDER, flag, str(annot), "-o", str(out), str(gfa)],
                       capture_output=True, timeout=300, check=False)
    except subprocess.TimeoutExpired:
        return {**row, "action": "keep", "reason": "pathfinder_timeout"}
    regen_fa, regen_bed = Path(f"{out}.{organelle}.ctg.fasta"), Path(f"{out}.{organelle}.ctg.bed")
    recs = read_records(regen_fa) if regen_fa.exists() and regen_fa.stat().st_size else []
    circles = [r for r in recs if header_field(r[0], "circular") == "true"]
    reproduces = orig_bed is not None and data_lines(orig_bed) == data_lines(regen_bed)
    row.update(regen_contigs=len(recs), regen_circular=len(circles), reproduces_original=reproduces)

    # mito: only a fully circular set (a resolved multichromosomal genome) -
    # a circle among unjoined pieces is not a resolved genome (Silene_vulgaris:
    # 1 circle + 172 linear pieces)
    resolved = circles and (organelle == "pltd" or len(circles) == len(recs))
    if recs and resolved and reproduces:
        if organelle == "pltd":
            keep = [max(circles, key=lambda r: float(header_field(r[0], "wlength") or 0))]
        else:
            keep = recs
        kept_ids = {r[0][1:].split()[0] for r in keep}
        dropped = [r[0][1:].split()[0] + ":" + (header_field(r[0], "length") or "?") + "bp"
                   for r in recs if r[0][1:].split()[0] not in kept_ids]
        return {**row, "action": "restore", "keep": keep, "kept_ids": kept_ids, "regen_bed": regen_bed,
                "reason": "pathfinder_reproduces_original",
                "kept_bp": sum(int(header_field(r[0], "length") or 0) for r in keep),
                "dropped_secondary": ",".join(dropped)}

    if cur_is_gfatk:
        placed = graph_coverage.unitig_content_placed(gfa, cur)
        row["circuit_graph_placed"] = round(placed, 3)
        orig_fa = find([sup], prefix, suf["ctg_fasta"])
        if placed < GATE[organelle] and orig_fa is not None:
            return {**row, "action": "revert", "reason": f"circuit_places_<{GATE[organelle]:.0%}_of_graph"}
        return {**row, "action": "keep", "reason": "circuit_passes_gate" if placed >= GATE[organelle]
                else "circuit_below_gate_but_no_original_to_revert_to"}
    return {**row, "action": "keep", "reason": "pathfinder_gave_no_reproducible_circle"}


def apply_one(row: dict) -> None:
    d, prefix, org = row["dir"], row["prefix"], row["organelle"]
    sup = d / "superseded"
    suf = ORIGINAL_SUFFIXES[org]
    cur = d / f"{prefix}{suf['ctg_fasta']}"
    if row["current"] == "gfatk_resolve":
        archive(cur, sup, f"{prefix}.{org}.gfatk_resolve.ctg.fasta")
    elif cur.exists():  # an empty placeholder
        archive(cur, sup, f"{prefix}{suf['ctg_fasta']}.empty")
    # oatk's own files come back out of superseded/. 06_reannotate_promoted.py
    # left the circuit's gene calls as the live .ctg.bed (archive them with
    # the circuit) and a copy of oatk's annotation (identical - leave it).
    for key in ("ctg_bed", "bed", "annot") + (("ctg_fasta",) if row["action"] == "revert" else ()):
        src, dst = sup / f"{prefix}{suf[key]}", d / f"{prefix}{suf[key]}"
        if not src.exists():
            continue
        if dst.exists():
            if key != "ctg_bed" or row["current"] != "gfatk_resolve":
                continue
            archive(dst, sup, f"{prefix}.{org}.gfatk_resolve.ctg.bed")
        shutil.move(str(src), str(dst))
    if row["action"] == "restore":
        with open(cur, "w") as fh:
            for header, seq in row["keep"]:
                fh.write(header)
                fh.writelines(seq)
        # ctg.bed restricted to the kept contigs (identical to oatk's otherwise)
        bed_lines = [l for l in open(row["regen_bed"])
                     if l.startswith("#") or l.split("\t")[0] in row["kept_ids"]]
        (d / f"{prefix}{suf['ctg_bed']}").write_text("".join(bed_lines))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--organelle", default="both", choices=["mito", "pltd", "both"])
    ap.add_argument("--species-list", default=None)
    ap.add_argument("--apply", action="store_true", help="act (default: dry run - report the plan only)")
    args = ap.parse_args()

    organelles = ["mito", "pltd"] if args.organelle == "both" else [args.organelle]
    species_filter = sd.load_species_list(Path(args.species_list)) if args.species_list else None
    rows = []
    with tempfile.TemporaryDirectory() as tmp:
        for organelle in organelles:
            for r in sd.discover_all(sd.repo_data_root(ROOT_DIR, organelle), organelle, species_filter=species_filter):
                if not r.gfa:
                    continue
                row = plan_one(r.species, organelle, ROOT_DIR / r.gfa if not Path(r.gfa).is_absolute() else Path(r.gfa), Path(tmp))
                if row["action"] == "not_a_candidate":
                    continue
                if args.apply and row["action"] in ("restore", "revert"):
                    apply_one(row)
                rows.append(row)
                print(f"[{'apply' if args.apply else 'plan'}] {r.species}.{organelle}: {row['action']} "
                      f"({row.get('reason', '')})", file=sys.stderr)

    cols = ["species", "organelle", "current", "action", "reason", "regen_contigs", "regen_circular",
            "reproduces_original", "kept_bp", "dropped_secondary", "circuit_graph_placed"]
    table = "\t".join(cols) + "\n" + "".join("\t".join(str(r.get(c, "")) for c in cols) + "\n" for r in rows)
    if args.apply:
        log = ROOT_DIR / "meta" / f"promotion_log_{datetime.now():%Y%m%d_%H%M%S}.tsv"
        log.write_text(table)
        print(f"[info] restore_pathfinder: log -> {log}", file=sys.stderr)
    sys.stdout.write(table)
    counts = {}
    for r in rows:
        counts[(r["organelle"], r["action"])] = counts.get((r["organelle"], r["action"]), 0) + 1
    print(f"[info] restore_pathfinder: {dict(sorted(counts.items()))}"
          f"{'' if args.apply else ' (dry run - nothing changed)'}", file=sys.stderr)


if __name__ == "__main__":
    main()
