#!/usr/bin/env python3
"""Map every GFA unitig onto the resolved `.ctg.fasta` it contributed to,
and extract the raw (non-linearised) unitig sequences as their own FASTA.

Why: a linearised contig doesn't, by itself, say which unitigs (in what
order/orientation) it was built from. Two different sources of that
sequence exist in this dataset, handled by two different methods here:

TIER 1 - oatk's own Pathfinder (the majority of species). Its ctg.fasta
header already carries the exact answer: `>ctg000001c length=368828
wlength=... nv=8 circular=true path=u66+,u68+,u69+,u64-,u66+,u67+,u69+,u65-`.
Walking that path and subtracting each junction's overlap (read from the
GFA's `L` lines, which give an exact bp overlap per link - e.g. `L u66 + u68
+ 1304M` = 1304bp) reconstructs ctg coordinates for every unitig placement
exactly, no guessing. Verified by hand on Arabidopsis_thaliana mito: the 8
segment lengths (371177bp total) minus the 7 internal + 1 circular-closing
overlap (2349bp total) equals exactly 368828 - the header's own length=.
This also naturally and correctly handles reverse-complement placement and
genuinely multi-copy unitigs (a unitig appearing twice in the path, e.g.
u66/u69 above - real recombination-mediated repeat structure, not noise).
The earlier substring-search-only version of this script missed one of
u66's two real placements; the path field catches both.

TIER 2 - `linearize`'s `gfatk resolve` fallback, promoted into data/ for
species oatk's own Pathfinder couldn't resolve at all (see
`linearize/README.md`, `05_promote_resolve.py`). These contigs' headers
use a completely different, path-free format (`..._circuit<N>:
n_segments=<N>:bp=<N>` - see 05_promote_resolve.py's CIRCUIT_RE) - gfatk's
resolve_summary.tsv records n_segments as a count only, not an
ordered/oriented walk. There is no authoritative path to parse for these,
so this falls back to exact sequence matching: each unitig's raw GFA
sequence (or its reverse complement) either appears in the ctg
byte-exact, or it doesn't. Circularity is handled by searching the ctg
doubled on itself (standard trick for matching across a circular join),
with a small symmetric trim tolerance at each end of the unitig for the
overlap a circular closure trims off (confirmed on Achillea_maritima mito:
u6's raw 202131bp is 36bp longer than the trimmed circular ctg's
202095bp).

Every row places the *whole* unitig. On a circular ctg a placement can
wrap across the origin, in which case ctg_end > ctg length (Tier 1: the
first unitig of every circular path, since oatk trims the closing overlap
off its start; Tier 2: any hit across the join). Consumers should take
positions modulo ctg length: ctg position x inside a placement is unitig
offset (x - ctg_start) mod ctg_len on `+`, counted from the unitig's end
on `-`.

Adjacent Tier 1 placements overlap by their link's bp, and within that
overlap the ctg carries only one neighbour's copy of the sequence - which
can differ from the other copy by a SNP or 1bp indel (seen in 17/1714
path-mapped species at full scale; every differing base was inside an
overlap, every placement's non-overlap core was byte-exact). So a
position inside an overlap belongs to two placements and is only
approximately placed in one of them; prefer the placement where it is
core.

Both tiers write the same map format - `source` records which produced a
given row (path|seqmatch). A Tier 1 species whose path can't be walked
(missing segment or link) drops to seqmatch rather than emit a partial
map; the run summary counts these separately so it isn't silent.

This deliberately does NOT try to resolve genes independently at unitig
level (that would fragment/lose anything spanning a join between two
placed unitigs - a real, if rare, case, and multi-copy unitigs make
"which copy" genuinely ambiguous in that direction). Translation happens
the other way: existing ctg-level results get relabelled into unitig
coordinates downstream using this map - ctg -> unitig is the unambiguous
direction.

Multi-copy unitigs are real signal, not noise - see
`01_linearization_qc.py`, and note the ctg-level GFF stays the primary
annotation output specifically because it preserves this dosage
information (a unitig-level GFF would collapse both copies onto one
sequence).

Output:
  results/unitig_map/<species>.<organelle>.tsv - one row per placement
    (a multi-copy unitig gets >1 row; an unplaced unitig gets none - the
    unitig FASTA below is the complete list): unitig, unitig_length,
    ctg_seqid, ctg_start, ctg_end, strand, source (path|seqmatch),
    match_trim - ctg_start/end are 0-based half-open, matching *.ctg.bed's
    convention elsewhere in this pipeline. match_trim (seqmatch only, 0
    for path) is how many bp at each end of the unitig were NOT verified:
    a trimmed hit is typically one of several near-identical repeat
    variants, so only [ctg_start + trim, ctg_end - trim) is trustworthy.
  results/unitig_fasta/<species>.<organelle>.unitig.fasta - the raw
    (non-linearised) unitig sequences, one record per GFA segment, named
    by their GFA segment ID (u5, u14, ...) - "the unitig reference"
"""
from __future__ import annotations

import os
import argparse
import re
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

# CODE_DIR is this checkout's analysis/ (shared modules, bundled reference
# data); data/ and every module's results/ and work/ live under the data
# root - PLANT_ORGANELLE_DATA_ROOT, or this checkout if unset.
CODE_DIR = Path(__file__).resolve().parents[2]
ROOT_DIR = Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT", str(CODE_DIR.parent)))
ANALYSIS_DIR = ROOT_DIR / "analysis"
sys.path.insert(0, str(CODE_DIR / "common"))

import species_discovery as sd  # noqa: E402

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
FLIP = {"+": "-", "-": "+"}
PATH_RE = re.compile(r"\bpath=(?P<path>\S+)")
CIRCULAR_RE = re.compile(r"\bcircular=(?P<circular>true|false)")
CIGAR_OP_RE = re.compile(r"(\d+)M")
MIN_SEQMATCH_BP = 50
MAX_TRIM_BP = 2000


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def read_gfa(gfa_path: Path) -> tuple[dict[str, str], dict[tuple[str, str, str, str], int]]:
    """Returns (segment_id -> sequence, (fromSeg,fromOrient,toSeg,toOrient) -> overlap_bp).
    The overlap lookup is populated in both directions (a GFA link is
    usable however the walk actually traverses it)."""
    segs: dict[str, str] = {}
    overlaps: dict[tuple[str, str, str, str], int] = {}
    with open(gfa_path) as fh:
        for line in fh:
            if line.startswith("S\t"):
                parts = line.rstrip("\n").split("\t")
                segs[parts[1]] = parts[2]
            elif line.startswith("L\t"):
                parts = line.rstrip("\n").split("\t")
                a, ao, b, bo, cigar = parts[1], parts[2], parts[3], parts[4], parts[5]
                bp = sum(int(n) for n in CIGAR_OP_RE.findall(cigar))
                overlaps[(a, ao, b, bo)] = bp
                overlaps[(b, FLIP[bo], a, FLIP[ao])] = bp
    return segs, overlaps


def read_ctg_fasta(fasta_path: Path) -> tuple[list[tuple[str, str | None, bool]], dict[str, str]]:
    """Returns ([(seqid, path_str_or_None, circular), ...], seqid -> sequence)."""
    headers = []
    seqs: dict[str, list[str]] = {}
    name = None
    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                header = line[1:]
                name = header.split()[0]
                m = PATH_RE.search(header)
                cm = CIRCULAR_RE.search(header)
                headers.append((name, m.group("path") if m else None,
                                cm is not None and cm.group("circular") == "true"))
                seqs[name] = []
            elif name is not None:
                seqs[name].append(line)
    return headers, {k: "".join(v) for k, v in seqs.items()}


def walk_path(path_str: str, circular: bool, segs: dict[str, str],
              overlaps: dict[tuple[str, str, str, str], int]) -> list[tuple[str, int, int, str]] | None:
    """Returns [(unitig, ctg_start, ctg_end, strand), ...] or None if any
    segment/overlap along the walk can't be resolved (caller falls back to
    seqmatch rather than emit a partial/wrong map)."""
    steps = [(tok[:-1], tok[-1]) for tok in path_str.split(",")]
    placements = []
    pos = 0
    for i, (seg, orient) in enumerate(steps):
        if seg not in segs:
            return None
        if i > 0:
            key = (*steps[i - 1], seg, orient)
            if key not in overlaps:
                return None
            pos -= overlaps[key]
        length = len(segs[seg])
        placements.append((seg, pos, pos + length, orient))
        pos += length
    if circular:
        # oatk trims the last->first closing overlap (for a one-step path,
        # the unitig's self-link) off the *start* of the first unitig, so
        # the ctg origin sits `closing` bp into it: shift everything back by
        # that, and the first placement wraps across the origin (verified
        # byte-exact on Arabidopsis_thaliana mito/pltd and
        # Achillea_maritima mito).
        key = (*steps[-1], *steps[0])
        if key not in overlaps:
            return None
        closing = overlaps[key]
        ctg_len = pos - closing
        placements = [(seg, (a - closing) % ctg_len, (a - closing) % ctg_len + (b - a), strand)
                      for seg, a, b, strand in placements]
    return placements


def _hits(haystack_doubled: str, hay_len: int, sub: str) -> list[int]:
    positions = []
    start = 0
    while True:
        p = haystack_doubled.find(sub, start)
        if p == -1 or p >= hay_len:
            return positions
        positions.append(p)
        start = p + 1


def find_all_with_tolerance(haystack_doubled: str, hay_len: int, needle: str,
                            max_trim: int) -> tuple[list[int], int]:
    """Search a self-doubled haystack for `needle`, allowing the smallest
    symmetric trim (0..max_trim) off both ends that yields a hit (handles a
    circular join's overlap trim). Hits are monotone in trim - a trimmed
    needle is a substring of every less-trimmed one - so binary search
    finds the smallest trim in ~log2(max_trim) scans. Returns (start
    positions of the *untrimmed* needle modulo hay_len, trim used) - only
    the needle's middle [trim, len - trim) is verified at those positions."""
    max_trim = min(max_trim, (len(needle) - MIN_SEQMATCH_BP) // 2)
    if max_trim < 0:
        return [], 0
    exact = _hits(haystack_doubled, hay_len, needle)
    if exact:
        return exact, 0
    if max_trim == 0 or not _hits(haystack_doubled, hay_len, needle[max_trim:len(needle) - max_trim]):
        return [], 0
    lo, hi = 0, max_trim  # lo misses, hi hits
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if _hits(haystack_doubled, hay_len, needle[mid:len(needle) - mid]):
            hi = mid
        else:
            lo = mid
    return [(p - hi) % hay_len for p in _hits(haystack_doubled, hay_len, needle[hi:len(needle) - hi])], hi


def build_one(species: str, organelle: str, gfa_path: str, ctg_fasta_path: str,
              map_dir: Path, fasta_dir: Path, force: bool) -> str:
    map_out = map_dir / f"{species}.{organelle}.tsv"
    fasta_out = fasta_dir / f"{species}.{organelle}.unitig.fasta"
    # a map older than its inputs is stale (e.g. linearize/07_restore_pathfinder
    # swapped the genome in place) and is rebuilt even without --force
    if (not force and map_out.exists() and fasta_out.exists() and map_out.stat().st_size > 0
            and map_out.stat().st_mtime >= max(Path(p).stat().st_mtime for p in (gfa_path, ctg_fasta_path))):
        return "skipped"

    segs, overlaps = read_gfa(Path(gfa_path))
    if not segs:
        return "no_segments"
    fasta_out.write_text("".join(f">{name}\n{seq}\n" for name, seq in segs.items()))

    ctg_headers, ctgs = read_ctg_fasta(Path(ctg_fasta_path))

    rows = ["unitig\tunitig_length\tctg_seqid\tctg_start\tctg_end\tstrand\tsource\tmatch_trim"]
    any_path_header = any(p is not None for _, p, _ in ctg_headers)
    any_walked = False
    for seqid, path_str, circular in ctg_headers:
        if path_str is None:
            continue
        placements = walk_path(path_str, circular, segs, overlaps)
        if placements is None:
            continue
        any_walked = True
        for uname, a, b, strand in placements:
            rows.append(f"{uname}\t{len(segs[uname])}\t{seqid}\t{a}\t{b}\t{strand}\tpath\t0")

    if not any_walked:
        doubled = {cid: cseq + cseq for cid, cseq in ctgs.items()}
        for uname, useq in segs.items():
            max_trim = min(MAX_TRIM_BP, len(useq) // 4)
            u_rc = revcomp(useq)
            for cid, cseq in ctgs.items():
                for strand, needle in (("+", useq), ("-", u_rc)):
                    positions, trim = find_all_with_tolerance(doubled[cid], len(cseq), needle, max_trim)
                    for pos in positions:
                        rows.append(f"{uname}\t{len(useq)}\t{cid}\t{pos}\t{pos + len(useq)}\t{strand}"
                                    f"\tseqmatch\t{trim}")

    map_out.write_text("\n".join(rows) + "\n")
    if any_walked:
        return "built_path"
    return "built_seqmatch_path_failed" if any_path_header else "built_seqmatch"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organelle", default="both", choices=["mito", "pltd", "both"])
    ap.add_argument("--species-list", default=None)
    ap.add_argument("--qc-status", default="pass", help="comma-separated subset of pass,flag,fail (default: pass)")
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    organelles = ["mito", "pltd"] if args.organelle == "both" else [args.organelle]
    species_filter = sd.load_species_list(Path(args.species_list)) if args.species_list else None

    qc_status = set(args.qc_status.split(","))
    qc_path = ANALYSIS_DIR / "qc_basic_stats" / "results" / "qc_summary.tsv"
    qc_eligible = None
    if qc_path.exists():
        import csv
        with open(qc_path) as fh:
            qc_eligible = {(r["species"], r["organelle"]) for r in csv.DictReader(fh, delimiter="\t")
                           if r["status"] in qc_status}
    else:
        print(f"[warn] missing {qc_path} - not filtering by QC status", file=sys.stderr)

    map_dir = ANALYSIS_DIR / "unitig_coords" / "results" / "unitig_map"
    fasta_dir = ANALYSIS_DIR / "unitig_coords" / "results" / "unitig_fasta"
    map_dir.mkdir(parents=True, exist_ok=True)
    fasta_dir.mkdir(parents=True, exist_ok=True)

    for organelle in organelles:
        data_root = sd.repo_data_root(ROOT_DIR, organelle)
        resolved = sd.discover_all(data_root, organelle, species_filter=species_filter)
        tasks = [r for r in resolved if r.status == "ok"
                 and (qc_eligible is None or (r.species, organelle) in qc_eligible)]

        counts: dict[str, int] = {}
        path_failed = []
        with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
            futures = {
                pool.submit(build_one, r.species, organelle, r.gfa, r.ctg_fasta,
                            map_dir, fasta_dir, args.force): r.species
                for r in tasks
            }
            for fut in as_completed(futures):
                outcome = fut.result()
                counts[outcome] = counts.get(outcome, 0) + 1
                if outcome == "built_seqmatch_path_failed":
                    path_failed.append(futures[fut])
        summary = " ".join(f"{k}={v}" for k, v in sorted(counts.items()))
        print(f"[info] build_unitig_map: {organelle}: {summary} -> {map_dir}", file=sys.stderr)
        for species in sorted(path_failed):
            print(f"[warn] {species}.{organelle}: path= header present but unwalkable, fell back to seqmatch",
                  file=sys.stderr)


if __name__ == "__main__":
    main()
