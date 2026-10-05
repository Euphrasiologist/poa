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
                  mitochondrial protein genes (oatk's own nhmmscan hits on
                  the unitigs, .annot_mito.txt; without it - poa run, or
                  any non-oatk assembly - poa's own gene calls carried onto
                  unitigs through unitig_coords' map: gene_source column)
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

# CODE_DIR is the installed poa/pipeline/ (shared modules, bundled reference
# data); data/ and every module's results/ and work/ live under the data
# root - PLANT_ORGANELLE_DATA_ROOT, or the current directory if unset.
CODE_DIR = Path(__file__).resolve().parents[2]
ROOT_DIR = Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT") or os.getcwd())
ANALYSIS_DIR = ROOT_DIR / "analysis"
sys.path.insert(0, str(CODE_DIR / "common"))

import species_discovery as sd  # noqa: E402

LOW_DEPTH_RATIO = 0.2
EMBEDDED_MAX_BP = 30_000
FLAG_MIN_BP = 20_000
FLAG_MIN_FRAC = 0.10
MIN_GENE_SCORE = 300.0
MAP_DIR = ANALYSIS_DIR / "unitig_coords" / "results" / "unitig_map"
GENE_CALLS = ANALYSIS_DIR / "denovo_annotation" / "results" / "gene_calls.tsv"
# a unitig carries a gene call if its placement covers at least this much of it
MIN_CALL_OVERLAP = 0.5


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


def load_gene_calls(path: Path) -> dict[str, list[tuple[str, int, int]]]:
    """species -> [(contig, start, end)] of confident mito protein-gene calls,
    the same filter gene_unitigs() applies to oatk's hits."""
    calls = defaultdict(list)
    if not path.exists():
        return calls
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["organelle"] != "mito" or r["gene"].startswith(("trn", "rrn")):
                continue
            try:
                if float(r["score"]) < MIN_GENE_SCORE:
                    continue
            except ValueError:
                continue
            calls[r["species"]].append((r["contig_id"], int(r["start"]), int(r["end"])))
    return calls


def fasta_lengths(path: Path) -> dict[str, int]:
    lengths, name = {}, None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                lengths[name] = 0
            elif name is not None:
                lengths[name] += len(line.strip())
    return lengths


def gene_unitigs_from_calls(calls: list[tuple[str, int, int]], map_path: Path,
                            ctg_lengths: dict[str, int]) -> set[str]:
    """Unitigs whose placement on the contigs covers >= MIN_CALL_OVERLAP of a
    gene call. A placement across a circular contig's origin runs past the
    contig's end (unitig_coords' convention) - the usual shape of a
    single-unitig circular mitogenome - so calls are also tested shifted by
    the contig length."""
    if not calls or not map_path.exists():
        return set()
    by_ctg = defaultdict(list)
    for ctg, start, end in calls:
        by_ctg[ctg].append((min(start, end), max(start, end)))
    out = set()
    with open(map_path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            trim = int(r.get("match_trim") or 0)
            u_start, u_end = int(r["ctg_start"]) + trim, int(r["ctg_end"]) - trim
            clen = ctg_lengths.get(r["ctg_seqid"], 0)
            shifts = (0, clen) if clen and u_end > clen else (0,)
            if any(min(end + k, u_end) - max(start + k, u_start) >= MIN_CALL_OVERLAP * (end - start)
                   for start, end in by_ctg.get(r["ctg_seqid"], ()) for k in shifts):
                out.add(r["unitig"])
    return out


def weighted_median(pairs):
    pairs = sorted(pairs)
    half, acc = sum(w for _, w in pairs) / 2, 0
    for v, w in pairs:
        acc += w
        if acc >= half:
            return v
    return None


def scan_one(species: str, gfa: Path, calls: list[tuple[str, int, int]], ctg: Path | None) -> dict | None:
    d = gfa.parent
    prefix = gfa.name[: -len(".mito.gfa")]
    annot = next((p for p in (d / f"{prefix}.annot_mito.txt", d / "superseded" / f"{prefix}.annot_mito.txt")
                  if p.exists()), None)
    depth, length, adj = read_graph(gfa)
    if not depth or any(v is None for v in depth.values()):
        return {"species": species, "status": "no_depth_tags"}
    map_path = MAP_DIR / f"{species}.mito.tsv"
    if annot:
        gene_source, genes = "oatk_annot", gene_unitigs(annot) & depth.keys()
    else:
        ctg_lengths = fasta_lengths(ctg) if ctg and ctg.exists() else {}
        gene_source, genes = "gene_calls", gene_unitigs_from_calls(calls, map_path, ctg_lengths) & depth.keys()
    if not genes:
        return {"species": species, "status": "no_gene_bearing_unitigs", "gene_source": gene_source}
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
    row = {"species": species, "status": "ok", "gene_source": gene_source, "run_prefix": prefix,
           "graph_bp": total_bp, "n_unitigs": len(depth), "genome_depth": round(genome_depth, 1),
           "low_depth_bp": sum(length[u] for u in low), "embedded_low_bp": embedded_bp,
           "n_embedded_pieces": n_embedded, "self_contained_low_bp": self_bp,
           "embedded_low_frac": round(embedded_bp / total_bp, 3) if total_bp else 0.0}
    row["flagged"] = embedded_bp >= FLAG_MIN_BP or row["embedded_low_frac"] >= FLAG_MIN_FRAC

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

    calls = load_gene_calls(GENE_CALLS)
    rows = []
    for r in sd.discover_all(sd.repo_data_root(ROOT_DIR, "mito"), "mito", species_filter=species_filter):
        if not r.gfa:
            continue
        gfa = Path(r.gfa) if os.path.isabs(r.gfa) else ROOT_DIR / r.gfa
        ctg = Path(r.ctg_fasta) if r.ctg_fasta else None
        if ctg and not ctg.is_absolute():
            ctg = ROOT_DIR / ctg
        row = scan_one(r.species, gfa, calls.get(r.species, []), ctg)
        if row:
            row["resolver"] = ("gfatk_resolve" if ctg and ctg.exists() and ctg.stat().st_size
                               and "resolver=gfatk_resolve" in open(ctg).readline() else "oatk_pathfinder")
            rows.append(row)

    cols = ["species", "status", "resolver", "run_prefix", "graph_bp", "n_unitigs", "genome_depth", "low_depth_bp",
            "embedded_low_bp", "n_embedded_pieces", "self_contained_low_bp", "embedded_low_frac",
            "assembly_bp_via_low", "assembly_frac_via_low", "flagged", "gene_source"]
    out = ANALYSIS_DIR / "qc_basic_stats" / "results" / "low_depth_paths.tsv"
    with open(out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    n_flag = sum(1 for r in rows if r.get("flagged"))
    print(f"[info] low_depth_paths: {len(rows)} mito graphs, {n_flag} flagged -> {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
