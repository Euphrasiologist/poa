#!/usr/bin/env python3
"""Find low-depth side paths embedded in mitochondrial assembly graphs.

Found on Silene_vulgaris (2026-09-28): a coverage-40 oatk rerun admitted
~206 kb of unitigs at ~6% of the genome's depth - 117 small pieces, each
linked into the mitochondrial graph at both ends (substoichiometric
recombination products or NUMT variants) - and a gfatk circuit then
walked the real 351 kb genome 3-8 times through them (1.25 Mb total).

What does NOT count as a problem here: gene-free sequence as such.
Multichromosomal Silene mitogenomes carry many gene-free chromosomes,
and their abundance tracks gene-bearing ones within ~2-fold (Wu et al.
2015, PNAS 112:10185 - S. noctiflora chromosomes at 12-25x; no coverage
difference between chromosomes with and without genes). So the signal
is depth and topology, never gene content:
  genome depth  = length-weighted median depth of unitigs carrying
                  mitochondrial protein genes (oatk's own nhmmscan hits)
  low depth     = below LOW_DEPTH_RATIO of genome depth (well outside the
                  ~2-fold chromosome abundance range)
  embedded      = a connected piece of low-depth unitigs (at most
                  EMBEDDED_MAX_BP) linked to the rest of the graph - a
                  side path, not a chromosome
  self-contained= a low-depth component with no link to the rest of the
                  graph, or larger than EMBEDDED_MAX_BP - could be a real
                  low-abundance chromosome; reported, never flagged
Flagged: embedded low-depth bp >= FLAG_MIN_BP or >= FLAG_MIN_FRAC of the
graph. Where unitig_coords has a map for the current assembly, also
reports how much of the assembly runs through low-depth unitigs.

Reports only - changes nothing. Output: results/low_depth_paths.tsv
"""
from __future__ import annotations

import argparse
import csv
import os
import statistics as st
import sys
from collections import defaultdict
from pathlib import Path

# CODE_DIR is this checkout's analysis/ (shared modules, bundled reference
# data); data/ and every module's results/ and work/ live under the data
# root - PLANT_ORGANELLE_DATA_ROOT, or this checkout if unset.
CODE_DIR = Path(__file__).resolve().parents[2]
ROOT_DIR = Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT", str(CODE_DIR.parent)))
ANALYSIS_DIR = ROOT_DIR / "analysis"
sys.path.insert(0, str(CODE_DIR / "common"))

import species_discovery as sd  # noqa: E402

LOW_DEPTH_RATIO = 0.2
EMBEDDED_MAX_BP = 30_000
FLAG_MIN_BP = 20_000
FLAG_MIN_FRAC = 0.10
MIN_GENE_SCORE = 300.0
MAP_DIR = ANALYSIS_DIR / "unitig_coords" / "results" / "unitig_map"


def read_graph(gfa: Path):
    depth, length, adj = {}, {}, defaultdict(set)
    for line in open(gfa):
        f = line.rstrip("\n").split("\t")
        if f[0] == "S":
            n = len(f[2])
            tags = {t[:2]: t.split(":", 2)[2] for t in f[3:] if t.count(":") >= 2}
            if "KC" in tags:
                d = float(tags["KC"]) / n
            elif "SC" in tags:
                d = float(tags["SC"])
            else:
                d = None
            depth[f[1]], length[f[1]] = d, n
        elif f[0] == "L":
            adj[f[1]].add(f[3])
            adj[f[3]].add(f[1])
    return depth, length, adj


def gene_unitigs(annot: Path) -> set[str]:
    out = set()
    for line in open(annot):
        if line.startswith("#"):
            continue
        f = line.split()
        if len(f) > 13 and not f[0].startswith(("trn", "rrn")):
            try:
                if float(f[13]) >= MIN_GENE_SCORE:
                    out.add(f[2])
            except ValueError:
                pass
    return out


def weighted_median(pairs):
    pairs = sorted(pairs)
    half, acc = sum(w for _, w in pairs) / 2, 0
    for v, w in pairs:
        acc += w
        if acc >= half:
            return v
    return None


def scan_one(species: str, gfa: Path) -> dict | None:
    d = gfa.parent
    prefix = gfa.name[: -len(".mito.gfa")]
    annot = next((p for p in (d / f"{prefix}.annot_mito.txt", d / "superseded" / f"{prefix}.annot_mito.txt")
                  if p.exists()), None)
    depth, length, adj = read_graph(gfa)
    if not depth or any(v is None for v in depth.values()):
        return {"species": species, "status": "no_depth_tags"}
    genes = gene_unitigs(annot) & depth.keys() if annot else set()
    if not genes:
        return {"species": species, "status": "no_gene_bearing_unitigs"}
    genome_depth = weighted_median([(depth[u], length[u]) for u in genes])
    low = {u for u in depth if depth[u] < LOW_DEPTH_RATIO * genome_depth}

    embedded_bp = self_bp = n_embedded = 0
    seen = set()
    for u in low:
        if u in seen:
            continue
        comp, stack = set(), [u]
        while stack:
            x = stack.pop()
            if x not in comp:
                comp.add(x)
                stack.extend(y for y in adj[x] if y in low)
        seen |= comp
        bp = sum(length[x] for x in comp)
        linked = any(y not in low for x in comp for y in adj[x])
        if linked and bp <= EMBEDDED_MAX_BP:
            embedded_bp += bp
            n_embedded += 1
        else:
            self_bp += bp

    total_bp = sum(length.values())
    row = {"species": species, "status": "ok", "run_prefix": prefix,
           "graph_bp": total_bp, "n_unitigs": len(depth), "genome_depth": round(genome_depth, 1),
           "low_depth_bp": sum(length[u] for u in low), "embedded_low_bp": embedded_bp,
           "n_embedded_pieces": n_embedded, "self_contained_low_bp": self_bp,
           "embedded_low_frac": round(embedded_bp / total_bp, 3) if total_bp else 0.0}
    row["flagged"] = embedded_bp >= FLAG_MIN_BP or row["embedded_low_frac"] >= FLAG_MIN_FRAC

    map_path = MAP_DIR / f"{species}.mito.tsv"
    if map_path.exists():
        rows = list(csv.DictReader(open(map_path), delimiter="\t"))
        placed = sum(int(r["ctg_end"]) - int(r["ctg_start"]) for r in rows)
        through_low = sum(int(r["ctg_end"]) - int(r["ctg_start"]) for r in rows if r["unitig"] in low)
        row["assembly_bp_via_low"] = through_low
        row["assembly_frac_via_low"] = round(through_low / placed, 3) if placed else ""
    return row


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species-list", default=None)
    args = ap.parse_args()
    species_filter = sd.load_species_list(Path(args.species_list)) if args.species_list else None

    rows = []
    for r in sd.discover_all(sd.repo_data_root(ROOT_DIR, "mito"), "mito", species_filter=species_filter):
        if not r.gfa:
            continue
        gfa = Path(r.gfa) if os.path.isabs(r.gfa) else ROOT_DIR / r.gfa
        row = scan_one(r.species, gfa)
        if row:
            ctg = Path(r.ctg_fasta) if r.ctg_fasta else None
            if ctg and not ctg.is_absolute():
                ctg = ROOT_DIR / ctg
            row["resolver"] = ("gfatk_resolve" if ctg and ctg.exists() and ctg.stat().st_size
                               and "resolver=gfatk_resolve" in open(ctg).readline() else "oatk_pathfinder")
            rows.append(row)

    cols = ["species", "status", "resolver", "run_prefix", "graph_bp", "n_unitigs", "genome_depth", "low_depth_bp",
            "embedded_low_bp", "n_embedded_pieces", "self_contained_low_bp", "embedded_low_frac",
            "assembly_bp_via_low", "assembly_frac_via_low", "flagged"]
    out = ANALYSIS_DIR / "qc_basic_stats" / "results" / "low_depth_paths.tsv"
    with open(out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    n_flag = sum(1 for r in rows if r.get("flagged"))
    print(f"[info] low_depth_paths: {len(rows)} mito graphs, {n_flag} flagged -> {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
