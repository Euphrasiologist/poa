#!/usr/bin/env python3
"""Generate .ctg.bed gene calls for the species promoted by 05_promote_resolve.py.

05_promote_resolve.py deliberately left these species without gene
annotation, flagged annotation_status=not_yet_annotated: gfatk resolve's
output is a bare sequence, and oatk's own nhmmscan-tblout -> .ctg.bed
conversion is internal to its Rust binary. This script closes that gap by
reverse-engineering and validating that conversion from real data, rather
than guessing at it.

The conversion has two parts:

1. Coordinate transform (unitig-local hit position -> final contig
   position). This is a path-walk: offset[i] = offset[i-1] + len(seg[i-1])
   - overlap(seg[i-1] -> seg[i]), reading each junction's overlap straight
   from the GFA's L-line CIGAR (M-only overlaps only - the only kind this
   dataset's GFAs use, verified). Circular genomes get one further global
   rotation: subtract the closing-junction's overlap from every offset
   (oatk's linearization starts partway into the first segment, at the
   point where the redundant closing overlap would otherwise duplicate
   sequence already at the end). A "-"-oriented segment's local
   coordinates are reflected (L - hi, L - lo) and its hit strand flips.
   VALIDATED against the whole existing (non-promoted) dataset: reproduces
   289,291/292,231 (99.0%) of real .ctg.bed rows exactly (position AND
   strand), with 701/1901 tested species matching with zero differences -
   see /dev/null (ephemeral validation script, not committed; rerun by
   comparing this script's own output against real ctg.bed files if you
   want to re-verify).

2. Inclusion filtering (which of the many raw nhmmscan hits per unitig
   oatk keeps in the final .ctg.bed). NOT successfully reverse-engineered:
   tried a global score cutoff, overlap-based non-max suppression
   (cross-gene and same-gene-only), and a per-gene empirical score floor
   learned from the real dataset - every filtering attempt REDUCED
   fidelity relative to no filtering at all (94.4%/75.6% recall vs 99.0%
   unfiltered), because real data has legitimate overlapping/tandem hits
   for the same gene that a naive suppression rule wrongly collapses.
   Per an explicit user decision: this script applies NO filtering -
   every raw hit is transformed and kept. This recovers 99% of what oatk's
   own filtering would show, at the cost of also including additional
   low-confidence hits oatk's undocumented internal filter would have
   dropped. Every output row's own score is included precisely so a
   downstream consumer can apply their own threshold if they want one.

Output: data/{mito,plastid}/<species>/<prefix>.{mito,pltd}.ctg.bed for
every promoted species, marked with a leading comment noting it was
generated this way (not oatk's own output). Also restores the archived
.annot_*.txt to the main species directory (copy, not move - the
superseded/ copy stays too) since nothing about the raw unitig-level
hits changed and species_discovery's completeness scoring expects it
there.
"""
from __future__ import annotations

import argparse
import re
import shutil
import os
import sys
from pathlib import Path

import pandas as pd

# Defaults to the old same-tree layout (code and data siblings under one
# repo root); set PLANT_ORGANELLE_DATA_ROOT to point at a separate data
# checkout instead (e.g. when this code is installed from its own repo).
ROOT_DIR = Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT", str(Path(__file__).resolve().parents[3])))
ANALYSIS_DIR = ROOT_DIR / "analysis"
sys.path.insert(0, str(ANALYSIS_DIR / "common"))

DATA_DIRNAME = {"mito": "mito", "pltd": "plastid"}
GFA_SUFFIX = {"mito": ".mito.gfa", "pltd": ".pltd.gfa"}
CTG_BED_SUFFIX = {"mito": ".mito.ctg.bed", "pltd": ".pltd.ctg.bed"}
CTG_FASTA_SUFFIX = {"mito": ".mito.ctg.fasta", "pltd": ".pltd.ctg.fasta"}

CIRCUIT_LOG_RE = re.compile(r"Circuit (\d+): (\d+) segments, (\d+) bp, path=(\S+)")
FLIP = {"+": "-", "-": "+"}


def parse_gfa(gfa_path: Path) -> tuple[dict, dict]:
    seg_len, overlaps = {}, {}
    with open(gfa_path) as fh:
        for line in fh:
            if line.startswith("S\t"):
                f = line.rstrip("\n").split("\t")
                seg_len[f[1]] = len(f[2])
            elif line.startswith("L\t"):
                f = line.rstrip("\n").split("\t")
                from_seg, from_o, to_seg, to_o, cigar = f[1], f[2], f[3], f[4], f[5]
                m = re.fullmatch(r"(\d+)M", cigar)
                if not m:
                    continue
                ov = int(m.group(1))
                overlaps[(from_seg, from_o, to_seg, to_o)] = ov
                overlaps[(to_seg, FLIP[to_o], from_seg, FLIP[from_o])] = ov
    return seg_len, overlaps


def parse_path_string(path_str: str) -> list[tuple[str, str]]:
    return [(tok[:-1], tok[-1]) for tok in path_str.split(",")]


def path_offsets(path_segs, seg_len, overlaps):
    """Circular only (every gfatk-resolve circuit is, by construction).
    Returns (list of (seg, orient, offset), rotation) or (None, None) if a
    junction's overlap isn't found (non-M CIGAR or missing link - can't
    model, caller should skip and log)."""
    offsets = []
    cur = 0
    for i, (seg, orient) in enumerate(path_segs):
        offsets.append((seg, orient, cur))
        cur += seg_len[seg]
        nxt_seg, nxt_orient = path_segs[(i + 1) % len(path_segs)]
        ov = overlaps.get((seg, orient, nxt_seg, nxt_orient))
        if ov is None:
            return None, None
        if i < len(path_segs) - 1:
            cur -= ov
    first_seg, first_orient = path_segs[0]
    last_seg, last_orient = path_segs[-1]
    rotation = overlaps.get((last_seg, last_orient, first_seg, first_orient))
    if rotation is None:
        return None, None
    return offsets, rotation


def parse_raw_hits(annot_path: Path) -> list[tuple]:
    hits = []
    with open(annot_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 14:
                continue
            gene, qseg = f[0], f[2]
            alifrom, alito = int(f[6]), int(f[7])
            strand, score = f[11], float(f[13])
            hits.append((gene, qseg, alifrom, alito, strand, score))
    return hits


def transform_circuit(path_segs, seg_len, overlaps, raw_hits) -> list[dict] | None:
    offsets, rotation = path_offsets(path_segs, seg_len, overlaps)
    if offsets is None:
        return None
    # naive (pre-rotation) total length = last offset + its segment's length;
    # the true circular genome length subtracts the closing-junction overlap once.
    last_seg, _, last_off = offsets[-1]
    total_length = last_off + seg_len[last_seg] - rotation

    rows = []
    skipped_straddling = 0
    for seg, orient, off in offsets:
        L = seg_len[seg]
        for gene, qseg, alifrom, alito, strand, score in raw_hits:
            if qseg != seg:
                continue
            lo, hi = min(alifrom, alito), max(alifrom, alito)
            if orient == "+":
                start, end, eff_strand = lo, hi, strand
            else:
                start, end, eff_strand = L - hi, L - lo, FLIP[strand]
            align_from, align_to = off + start - rotation, off + end - rotation
            # wrap into [0, total_length): only ever needed at the low end (offsets
            # are constructed to never reach total_length at the high end) - a hit
            # that straddles the origin itself (from_ <0<= to, or similar) can't be
            # expressed as one bed interval on a circular genome; skip it rather
            # than emit something wrong (rare - not seen in validation, but cheap
            # to guard against rather than assume away).
            if align_from < 0 and align_to < 0:
                align_from += total_length
                align_to += total_length
            elif align_from < 0 or align_to < 0:
                skipped_straddling += 1
                continue
            rows.append({"align_from": align_from, "align_to": align_to,
                         "gene": gene, "score": round(min(score, 1000)), "strand": eff_strand})
    if skipped_straddling:
        print(f"[warn] {skipped_straddling} hit(s) straddling the circular origin skipped "
              f"(can't be expressed as one bed interval)", file=sys.stderr)
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--promotion-log", required=True, help="the meta/promotion_log_<ts>.tsv from 05_promote_resolve.py")
    ap.add_argument("--work-dir", default=str(ANALYSIS_DIR / "linearize" / "work"))
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    promoted = pd.read_csv(args.promotion_log, sep="\t")
    promoted = promoted[promoted["promoted"]]
    work_dir = Path(args.work_dir)

    n_ok = n_no_annot = n_no_log = n_unmodeled = 0
    for _, row in promoted.iterrows():
        species, organelle = row["species"], row["organelle"]
        data_dir = ROOT_DIR / "data" / DATA_DIRNAME[organelle] / species
        superseded_dir = data_dir / "superseded"

        gfa = next(data_dir.glob(f"*{GFA_SUFFIX[organelle]}"), None)
        ctg_fasta = data_dir / f"{gfa.name[:-len(GFA_SUFFIX[organelle])]}{CTG_FASTA_SUFFIX[organelle]}" if gfa else None
        if not gfa or not ctg_fasta or not ctg_fasta.exists():
            print(f"[warn] {species} ({organelle}): missing gfa/ctg.fasta, skipped", file=sys.stderr)
            continue
        prefix = gfa.name[: -len(GFA_SUFFIX[organelle])]

        annot_candidates = list(superseded_dir.glob(f"*.annot_{organelle}.txt")) if superseded_dir.is_dir() else []
        annot = annot_candidates[0] if annot_candidates else None
        if annot is None:
            n_no_annot += 1
            print(f"[warn] {species} ({organelle}): no archived annot_*.txt (genuinely never annotated), skipped",
                  file=sys.stderr)
            continue

        log_path = work_dir / f"{species}.{organelle}.resolve.log"
        if not log_path.exists():
            n_no_log += 1
            print(f"[warn] {species} ({organelle}): no resolve.log (path info unrecoverable), skipped",
                  file=sys.stderr)
            continue
        circuit_paths = {}
        for line in open(log_path):
            m = CIRCUIT_LOG_RE.search(line)
            if m:
                circuit_paths[int(m.group(1))] = parse_path_string(m.group(4))

        seg_len, overlaps = parse_gfa(gfa)
        raw_hits = parse_raw_hits(annot)

        all_rows = []
        unmodeled = False
        for ctg_num, circuit_id in enumerate(sorted(circuit_paths), start=1):
            rows = transform_circuit(circuit_paths[circuit_id], seg_len, overlaps, raw_hits)
            if rows is None:
                unmodeled = True
                break
            for r in rows:
                r["seq_name"] = f"ctg{ctg_num:06d}c"
            all_rows.extend(rows)

        if unmodeled:
            n_unmodeled += 1
            print(f"[warn] {species} ({organelle}): unmodeled junction (non-M CIGAR or missing link), skipped",
                  file=sys.stderr)
            continue

        out_path = data_dir / f"{prefix}{CTG_BED_SUFFIX[organelle]}"
        if not args.dry_run:
            with open(out_path, "w") as fh:
                fh.write("# reannotated by analysis/linearize/src/06_reannotate_promoted.py, not oatk itself.\n"
                         "# Coordinate transform validated at 99.0% row-fidelity against the existing dataset,\n"
                         "# but NO inclusion filter applied (unlike oatk's own ctg.bed) - see this script's\n"
                         "# docstring. Every row is a real nhmmscan hit at a correctly-transformed position;\n"
                         "# some would not have made oatk's own (undocumented) cut.\n"
                         "#seq_name\talign_from\talign_to\tgene_name\tscore_capped_at_1000\tstrand\n")
                for r in sorted(all_rows, key=lambda r: (r["seq_name"], r["align_from"])):
                    fh.write(f"{r['seq_name']}\t{r['align_from']}\t{r['align_to']}\t{r['gene']}\t"
                             f"{r['score']}\t{r['strand']}\n")
            restored_annot = data_dir / annot.name
            if not restored_annot.exists():
                shutil.copy2(annot, restored_annot)
        n_ok += 1

    print(f"[info] reannotate_promoted: ok={n_ok} no_annot={n_no_annot} no_log={n_no_log} "
          f"unmodeled={n_unmodeled} / {len(promoted)} promoted rows", file=sys.stderr)


if __name__ == "__main__":
    main()
