#!/usr/bin/env python3
"""Build one standardized, hierarchical GFF3 per species, merging the
three annotation modules.

Gene models are built from EVERY denovo_annotation hit, not just the best
per gene:
- overlapping hits of a gene collapse to one locus (antisense members -
  weak reverse-strand model hits - dropped for the strongest hit's strand);
  tRNA loci are first resolved across names, since one real tRNA is hit
  by many anticodon models (tRNAscan-SE's identity wins, oatkDB's for CAU);
- protein-coding and tRNA loci chain into cis-spliced genes (same contig
  and strand, MIN_INTRON <= intron <= MAX_INTRON, each exon continuing the
  gene model - hmm_from/hmm_to - where the last left off; a second copy
  restarts it; a shorter gap merges the pieces into one exon);
- every remaining locus/chain is a copy: `gene-X`, `gene-X-2`, ...
  Copies under PARTIAL_COVERAGE of the gene's usual model span are tagged
  `partial=true`; `copy_number` counts the rest. Extra tRNA copies need
  tRNAscan-SE and oatkDB to agree (or to be an intron-split chain).

Corrections then claim the loci they were computed on:
1. `trans_splicing` reconstruction (`nad1`/`nad2`/`nad5`/`rps3`, >=1 filled
   exon) - unless a cis-spliced model of the same loci has more exons
   (rps3 is cis-spliced in most plants).
2. `editing`'s RNA-editing-corrected call, on the single-exon copy holding
   the hit it corrected (an ORF search over one exon of a cis-spliced gene
   can't fix that gene's boundaries, so multi-exon copies stay raw).
3. Everything else is `raw` (also the only tier for tRNA/rRNA).

Which tier won is recorded per gene as `annotation_tier=reconstructed|
edited|raw` - an explicit confidence signal, not just implicit in which
script you happened to run.

GFF3 mechanics used here (see `../README.md` for the full explanation):
- Protein-coding genes are `gene` -> `mRNA` -> one or more `CDS` rows
  sharing one `ID` (that's the spec's actual mechanism for a discontinuous/
  multi-exon feature - it does not require those rows to share a strand or
  seqid, which is exactly what a trans-spliced gene needs).
- tRNA/rRNA are `gene` -> `tRNA`/`rRNA` directly - no `CDS`/`mRNA` layer,
  they're not protein-coding - with `exon` children when intron-split.
- Predicted RNA-editing sites are child `sequence_alteration` features
  under the relevant `CDS` `ID` (a real, general Sequence Ontology term
  for "this position differs from a reference" - a pragmatic choice, not
  a claim that this is THE canonical way to encode RNA editing in GFF3).
- `phase` (GFF3 column 8) is `0` for `editing`/`trans_splicing`-derived
  exons. For editing's single-exon calls that's exact; for trans_splicing
  it is NOT - see README ("Known issue"). Raw cis-spliced models take each CDS's
  phase from the exon lengths before it - consistent, but only as exact
  as HMMER envelope bounds (their mRNA says exon_boundaries=approximate).

Assembly provenance and unitig coordinates (see README):
- Each ctg gets a `##sequence-region` pragma and a `region` feature
  carrying `resolver=oatk_pathfinder|gfatk_resolve` (which linearisation
  produced it - read from the ctg.fasta header) and GFF3's `Is_circular`.
- Every feature gets `unitig_loc=<unitig>:<start>-<end>:<strand>` (1-based
  inclusive, on the raw unitig from `unitig_coords`) when it lies inside
  one unitig placement, or `unitig_span=<u>,<u>,...` when it crosses a
  junction between placements - flagged explicitly, not guessed. Needs
  `unitig_coords/results/unitig_map/`; skipped per species if absent.
- A second, unitig-level GFF places the same features on the unitig
  sequences themselves. A junction-crossing feature is split into pieces
  sharing one `ID` (GFF3's discontinuous-feature mechanism), CDS phase
  recomputed per piece. This view collapses repeat copies onto one
  sequence - the ctg-level GFF stays primary; each unitig's `region` row
  records `n_ctg_copies` so that dosage isn't silently lost.

Output:
  analysis/gff_export/results/<species>.<organelle>.gff (ctg-level, primary)
  analysis/gff_export/results/unitig/<species>.<organelle>.unitig.gff
"""
from __future__ import annotations

import os
import argparse
import csv
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

# CODE_DIR is this checkout's analysis/ (shared modules, bundled reference
# data); data/ and every module's results/ and work/ live under the data
# root - PLANT_ORGANELLE_DATA_ROOT, or this checkout if unset.
CODE_DIR = Path(__file__).resolve().parents[2]
ROOT_DIR = Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT", str(CODE_DIR.parent)))
ANALYSIS_DIR = ROOT_DIR / "analysis"
sys.path.insert(0, str(CODE_DIR / "common"))

import species_discovery as sd  # noqa: E402

sys.path.insert(0, str(CODE_DIR / "trans_splicing" / "src"))
from reference_genes import CIS_SPLICED_MITO, GENES as RECONSTRUCTED_GENES  # noqa: E402
# a cis-spliced gene's reconstruction must agree with the HMM evidence: every
# exon overlaps a raw hit for the gene and extends no further than this past it
# (a full-length template on a lineage-truncated gene - Brassicaceae rpl2 -
# otherwise aligns on into downstream sequence: 573 bp in Arabidopsis)
CIS_MAX_EXTENSION = 100
UNITIG_MAP_DIR = ANALYSIS_DIR / "unitig_coords" / "results" / "unitig_map"
RESOLVER_RE = re.compile(r"\bresolver=(\S+)")
CIRCULAR_RE = re.compile(r"\bcircular=(true|false)")
FLIP_STRAND = {"+": "-", "-": "+", ".": "."}


@dataclass
class Segment:
    """One linear stretch of a unitig placement on a ctg (a placement that
    wraps a circular ctg's origin is two segments)."""
    unitig: str
    unitig_length: int
    strand: str
    ctg_lo: int
    ctg_hi: int
    unitig_off: int  # forward-walk offset into the (oriented) unitig at ctg_lo


@dataclass
class SeqContext:
    lengths: dict[str, int] = field(default_factory=dict)
    resolver: dict[str, str] = field(default_factory=dict)
    circular: dict[str, bool] = field(default_factory=dict)
    segments: dict[str, list[Segment]] = field(default_factory=dict)  # by ctg, sorted
    tiles: dict[str, list[Segment]] = field(default_factory=dict)     # disjoint partition
    unitig_lengths: dict[str, int] = field(default_factory=dict)
    unitig_copies: dict[str, int] = field(default_factory=dict)


def load_seq_context(ctg_fasta: Path, map_path: Path) -> SeqContext:
    ctx = SeqContext()
    name = None
    with open(ctg_fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                header = line[1:].rstrip("\n")
                name = header.split()[0]
                ctx.lengths[name] = 0
                m = RESOLVER_RE.search(header)
                ctx.resolver[name] = m.group(1) if m else ("oatk_pathfinder" if "path=" in header else "unknown")
                c = CIRCULAR_RE.search(header)
                ctx.circular[name] = c is not None and c.group(1) == "true"
            elif name is not None:
                ctx.lengths[name] += len(line.strip())
    if not map_path.exists():
        return ctx

    for row in pd.read_csv(map_path, sep="\t").itertuples():
        clen = ctx.lengths.get(row.ctg_seqid)
        if not clen:
            continue
        # A trimmed seqmatch hit (near-identical repeat variant) only
        # verifies its middle, so only that is claimed - an unverified end
        # would place features on sequence the unitig doesn't carry.
        trim = int(getattr(row, "match_trim", 0) or 0)
        a, b = int(row.ctg_start) + trim, int(row.ctg_end) - trim
        segs = ctx.segments.setdefault(row.ctg_seqid, [])
        if a >= clen:  # the verified core starts past the origin
            a, b, off = a - clen, b - clen, trim
        else:
            off = trim
        segs.append(Segment(row.unitig, int(row.unitig_length), row.strand, a, min(b, clen), off))
        if b > clen:  # wraps the origin (see unitig_coords/src/00_build_unitig_map.py)
            segs.append(Segment(row.unitig, int(row.unitig_length), row.strand, 0, b - clen, off + clen - a))
        ctx.unitig_lengths[row.unitig] = int(row.unitig_length)
        ctx.unitig_copies[row.unitig] = ctx.unitig_copies.get(row.unitig, 0) + 1

    # Adjacent placements overlap by their GFA link; tiles hand each
    # overlap to the earlier placement so a split feature's pieces never
    # double-count bases.
    for ctg, segs in ctx.segments.items():
        segs.sort(key=lambda s: (s.ctg_lo, s.ctg_hi))
        tiles, covered_to = [], 0
        for s in segs:
            lo = max(s.ctg_lo, covered_to)
            if lo < s.ctg_hi:
                tiles.append(Segment(s.unitig, s.unitig_length, s.strand, lo, s.ctg_hi,
                                     s.unitig_off + (lo - s.ctg_lo)))
                covered_to = s.ctg_hi
        ctx.tiles[ctg] = tiles
    return ctx


def to_unitig(seg: Segment, s: int, e: int, strand: str) -> tuple[int, int, str]:
    """ctg [s,e) inside `seg` -> (unitig start, unitig end, strand), 0-based
    half-open on the unitig's own forward sequence."""
    o_s = seg.unitig_off + (s - seg.ctg_lo)
    o_e = o_s + (e - s)
    if seg.strand == "+":
        return o_s, o_e, strand
    return seg.unitig_length - o_e, seg.unitig_length - o_s, FLIP_STRAND[strand]


def locate(ctx: SeqContext, ctg: str, s: int, e: int) -> tuple[Segment | None, list[Segment]]:
    """(containing segment or None, tiles the feature crosses in ctg order).
    A containing segment clear of every other placement (its core) is
    preferred over one where the feature sits inside a link overlap."""
    segs = ctx.segments.get(ctg, [])
    containing = [g for g in segs if g.ctg_lo <= s and e <= g.ctg_hi]
    if containing:
        def in_core(g):
            return not any(o is not g and o.ctg_lo < e and s < o.ctg_hi for o in segs)
        return next((g for g in containing if in_core(g)), containing[0]), []
    return None, [t for t in ctx.tiles.get(ctg, []) if t.ctg_lo < e and s < t.ctg_hi]


def add_unitig_layer(lines: list[str], ctx: SeqContext) -> tuple[list[str], list[str]]:
    """Returns (ctg-level lines with unitig attributes added, unitig-level lines)."""
    ctg_out, unitig_out = [], []
    for line in lines:
        f = line.split("\t")
        ctg, s, e, strand = f[0], int(f[3]) - 1, int(f[4]), f[6]
        if ctg not in ctx.segments:
            ctg_out.append(line)
            continue
        seg, crossed = locate(ctx, ctg, s, e)
        if seg is not None:
            us, ue, ustrand = to_unitig(seg, s, e, strand)
            f[8] += f";unitig_loc={seg.unitig}:{us + 1}-{ue}:{ustrand}"
            unitig_out.append("\t".join([seg.unitig, f[1], f[2], str(us + 1), str(ue), f[5], ustrand, f[7], f[8]]))
        else:
            f[8] += f";unitig_span={','.join(t.unitig for t in crossed) or 'unplaced'}"
            for t in crossed:
                ps, pe = max(s, t.ctg_lo), min(e, t.ctg_hi)
                us, ue, ustrand = to_unitig(t, ps, pe, strand)
                phase = f[7]
                if f[2] == "CDS" and phase != ".":
                    # bases of this CDS upstream of the piece, in transcript direction
                    k = (ps - s if strand != "-" else e - pe) - int(phase)
                    phase = str((3 - k % 3) % 3)
                unitig_out.append("\t".join([t.unitig, f[1], f[2], str(us + 1), str(ue), f[5], ustrand, phase,
                                             f[8] + f";unitig_split={ctg}:{ps + 1}-{pe}"]))
        ctg_out.append("\t".join(f))
    return ctg_out, unitig_out


LOW_DEPTH_PATHS = ANALYSIS_DIR / "qc_basic_stats" / "results" / "low_depth_paths.tsv"


def load_linearisation_notes() -> dict[str, float]:
    """species -> embedded low-depth fraction, for mito graphs that
    qc_basic_stats/src/06_low_depth_paths.py flags: minor configurations
    (~5-10% depth side paths - recombination products, nuclear copies or
    artefacts; not yet told apart) mean any single linear/circular reading
    of the graph is a modelling choice. The unitig-level GFF is then the
    safer view; the contig-level one is annotated but marked uncertain."""
    if not LOW_DEPTH_PATHS.exists():
        return {}
    with open(LOW_DEPTH_PATHS) as fh:
        return {r["species"]: float(r["embedded_low_frac"]) for r in csv.DictReader(fh, delimiter="\t")
                if r.get("flagged") == "True"}


def linearisation_comment(frac: float) -> str:
    return (f"# linearisation uncertain: {frac:.0%} of the mito assembly graph is low-depth side paths "
            f"(analysis/qc_basic_stats/results/low_depth_paths.tsv) - the contig sequence is one possible "
            f"reading of the graph; prefer the unitig-level annotation (results/unitig/)")


def ctg_header(ctx: SeqContext, seqids: list[str], low_depth_frac: float | None = None) -> tuple[list[str], list[str]]:
    pragmas = [f"##sequence-region {c} 1 {ctx.lengths[c]}" for c in seqids]
    if low_depth_frac is not None:
        pragmas.append(linearisation_comment(low_depth_frac))
    note = {"linearisation": "uncertain", "linearisation_reason": "low_depth_side_paths",
            "graph_low_depth_frac": f"{low_depth_frac:.3f}"} if low_depth_frac is not None else {}
    regions = ["\t".join([c, "gff_export", "region", "1", str(ctx.lengths[c]), ".", "+", ".",
                          gff_attrs({"ID": f"region-{c}", "resolver": ctx.resolver[c],
                                     "Is_circular": "true" if ctx.circular[c] else None, **note})])
               for c in seqids]
    return pragmas, regions


def unitig_header(ctx: SeqContext, unitigs: list[str], low_depth_frac: float | None = None) -> tuple[list[str], list[str]]:
    pragmas = [f"##sequence-region {u} 1 {ctx.unitig_lengths[u]}" for u in unitigs]
    if low_depth_frac is not None:
        pragmas.append(linearisation_comment(low_depth_frac))
    regions = ["\t".join([u, "gff_export", "region", "1", str(ctx.unitig_lengths[u]), ".", "+", ".",
                          gff_attrs({"ID": f"region-{u}", "n_ctg_copies": ctx.unitig_copies[u]})])
               for u in unitigs]
    return pragmas, regions


def assemble(pragmas: list[str], regions: list[str], features: list[str]) -> str:
    """GFF3 text with each seqid's region row leading its features."""
    by_seqid: dict[str, list[str]] = defaultdict(list)
    for line in features:
        by_seqid[line.split("\t", 1)[0]].append(line)
    body = []
    for region in regions:
        seqid = region.split("\t", 1)[0]
        body.append(region)
        body.extend(sorted(by_seqid.pop(seqid, []), key=lambda l: int(l.split("\t")[3])))
    for seqid in sorted(by_seqid):  # features on a seqid with no region row (shouldn't happen)
        body.extend(by_seqid[seqid])
    return "\n".join(["##gff-version 3", *pragmas, *body]) + "\n"


MAX_INTRON = {"mito": 8000, "pltd": 4000}  # longest cis intron in the Arabidopsis mito reference: 3511bp
MIN_INTRON = 250  # a shorter gap between consecutive model pieces is one exon whose HMM hit broke
                  # (ycf1/ycf2/rpoC2 gaps of 25-200 bp) - organellar group I/II introns are longer
MODEL_OVERLAP_TOL = 60  # successive exons' HMM ranges may overlap a little at envelope edges
PARTIAL_COVERAGE = 0.8  # below this share of its gene's usual model span a copy is tagged partial -
                        # Arabidopsis mito: real copies >= 0.89, fragments/pseudogene pieces <= 0.68
TRNASCAN_SUPPORT = 35.0  # extra tRNA copies need tRNAscan-SE >= this AND an oatkDB hit; real
                         # Arabidopsis mito tRNAs score >= 41.5, plastid-derived fragments ~27


@dataclass
class Locus:
    """One genomic locus of a gene: overlapping hits collapsed (antisense
    members of the overlap - weak reverse-strand rRNA/tRNA model hits -
    dropped in favour of the strongest hit's strand)."""
    contig: str
    start: int
    end: int
    strand: str
    score: float
    hmm_from: int | None
    hmm_to: int | None
    model_len: int | None
    hits: list = field(default_factory=list)


def hit_strength(h) -> float:
    """Comparable-enough ranking across sources, used only to pick a strand
    and a representative within one overlapping cluster: bitscores as-is,
    barrnap E-values as -log10 (0.0 means below double precision)."""
    if h.source == "barrnap":
        return 300.0 if h.score <= 0 else -np.log10(h.score)
    return float(h.score)


def collapse_hits(hits: pd.DataFrame) -> list[Locus]:
    loci = []
    for contig, ch in hits.sort_values("start").groupby("contig_id", sort=False):
        cluster, cluster_end = [], -1
        for h in list(ch.itertuples()) + [None]:
            if h is not None and (not cluster or h.start < cluster_end):
                cluster.append(h)
                cluster_end = max(cluster_end, h.end)
                continue
            best = max(cluster, key=hit_strength)
            keep = [c for c in cluster if c.strand == best.strand]
            model = [c for c in keep if c.source == "oatkDB" and pd.notna(c.hmm_from)]
            loci.append(Locus(contig, min(c.start for c in keep), max(c.end for c in keep), best.strand,
                              max(float(c.score) for c in model) if model else hit_strength(best),
                              int(min(c.hmm_from for c in model)) if model else None,
                              int(max(c.hmm_to for c in model)) if model else None,
                              int(model[0].model_len) if model else None, keep))
            if h is not None:
                cluster, cluster_end = [h], h.end
    return loci


def chain_exons(loci: list[Locus], organelle: str) -> list[list[Locus]]:
    """Group a protein-coding gene's loci into cis-spliced genes: same
    contig and strand, intron <= MAX_INTRON, and each exon picking up the
    gene model where the previous one left off. A second copy restarts
    the model, so it never extends an existing chain. A gap under
    MIN_INTRON isn't an intron: the two pieces merge into one exon."""
    def tx_order(l):
        return (l.contig, l.strand, l.start if l.strand != "-" else -l.end)
    chains: list[list[Locus]] = []
    for loc in sorted(loci, key=tx_order):
        best, best_gap = None, None
        for chain in chains:
            last = chain[-1]
            if last.contig != loc.contig or last.strand != loc.strand or loc.hmm_from is None or last.hmm_to is None:
                continue
            gap = loc.start - last.end if loc.strand != "-" else last.start - loc.end
            if not (0 <= gap <= MAX_INTRON[organelle]):
                continue
            if loc.hmm_from < last.hmm_to - MODEL_OVERLAP_TOL or loc.hmm_to <= last.hmm_to:
                continue
            if best is None or gap < best_gap:
                best, best_gap = chain, gap
        if best is not None and best_gap < MIN_INTRON:
            last = best[-1]
            best[-1] = Locus(last.contig, min(last.start, loc.start), max(last.end, loc.end), last.strand,
                             max(last.score, loc.score), last.hmm_from, loc.hmm_to, last.model_len,
                             last.hits + loc.hits)
        elif best is not None:
            best.append(loc)
        else:
            chains.append([loc])
    return chains


def covered_model_bp(spans) -> int:
    return sum(e - s + 1 for s, e in merge_intervals(sorted(spans)))


def expected_model_spans(calls: pd.DataFrame) -> dict[tuple[str, str], float]:
    """(organelle, gene) -> median across species of the model span that
    species' oatkDB hits cover. oatkDB models often extend well beyond the
    CDS (rps12's is 1119 nt for a ~380 bp gene), so model_len understates
    how complete a hit is; what full-length genes actually cover does not."""
    o = calls[(calls.source == "oatkDB") & calls.hmm_from.notna()]
    per_species = o.groupby(["organelle", "gene", "species"])[["hmm_from", "hmm_to"]].apply(
        lambda d: covered_model_bp(zip(d.hmm_from.astype(int), d.hmm_to.astype(int))))
    return per_species.groupby(level=[0, 1]).median().to_dict()


def model_coverage(chain: list[Locus], expected: float | None) -> float | None:
    spans = [(l.hmm_from, l.hmm_to) for l in chain if l.hmm_from is not None]
    if not spans or not expected:
        return None
    return min(1.0, covered_model_bp(spans) / expected)


def resolve_trna_identity(trna_hits: pd.DataFrame) -> pd.DataFrame:
    """The same tRNA locus is hit by many oatkDB tRNA models (a real trnK
    also scores as trnI/trnM/trnN/trnQ/trnR/trnT/trnV/trnW...). Cluster
    overlapping tRNA hits across names and keep one identity per locus:
    tRNAscan-SE's (anticodon-based) where it called one, else the
    strongest oatkDB model's. tRNAscan-SE hits stay, relabelled, as
    support evidence; other oatkDB models' hits go - their model
    coordinates belong to a different HMM."""
    keep: dict = {}  # hit index -> resolved name
    for _, ch in trna_hits.sort_values("start").groupby("contig_id", sort=False):
        cluster, cluster_end = [], -1
        for h in list(ch.itertuples()) + [None]:
            if h is not None and (not cluster or h.start < cluster_end):
                cluster.append(h)
                cluster_end = max(cluster_end, h.end)
                continue
            ts = [c for c in cluster if c.source == "tRNAscan-SE"]
            name = max(ts or cluster, key=hit_strength).gene
            if name == "trnM-CAU":
                # tRNAscan-SE calls every CAU anticodon Met; oatkDB has
                # distinct initiator (trnfM) and lysidine (trnI-CAU) models
                cau = [c for c in cluster if c.source == "oatkDB" and c.gene.endswith("-CAU")]
                if cau:
                    name = max(cau, key=hit_strength).gene
            keep.update({c.Index: name for c in cluster if c.gene == name or c.source == "tRNAscan-SE"})
            if h is not None:
                cluster, cluster_end = [h], h.end
    out = trna_hits.loc[list(keep)]
    return out.assign(gene=[keep[i] for i in out.index])


def trna_supported(copy: dict) -> bool:
    """Both predictors agree - or it's a near-complete intron-split tRNA
    chain, which tRNAscan-SE misses (plastid trnA-UGC/trnI-GAU sit in the
    inverted repeat, so their real second copy looks exactly like this)."""
    hits = [h for l in copy["loci"] for h in l.hits]
    if len(copy["loci"]) > 1 and not copy["partial"]:
        return True
    return (any(h.source == "oatkDB" for h in hits)
            and any(h.source == "tRNAscan-SE" and h.score >= TRNASCAN_SUPPORT for h in hits))


def merge_intervals(ivs):
    out = []
    for s, e in ivs:
        if out and s <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


def gff_attrs(d: dict) -> str:
    return ";".join(f"{k}={v}" for k, v in d.items() if v is not None and v != "")


def copy_attrs(copy: dict) -> dict:
    cov = copy.get("coverage")
    return {"copy_number": copy.get("copy_number"),
            "model_coverage": f"{cov:.2f}" if cov is not None else None,
            "partial": "true" if copy.get("partial") else None}


def emit_rna_gene(lines: list[str], gene: str, suffix: str, chain: list[Locus], feature: str, copy: dict) -> None:
    """gene -> tRNA/rRNA, plus exon rows when the RNA is split by an intron
    (plastid trnK/trnL/trnV/trnI/trnA/trnG carry group II introns).
    tool_support records which predictors called the locus, so a tRNA kept
    only because it is the gene's sole locus (no tRNAscan-SE call) stays
    distinguishable from one both tools agree on."""
    gid, fid = f"gene-{gene}{suffix}", f"{feature.lower()}-{gene}{suffix}"
    tools = "+".join(sorted({h.source for l in chain for h in l.hits}))
    contig, strand = chain[0].contig, chain[0].strand
    coords = [str(min(l.start for l in chain) + 1), str(max(l.end for l in chain)), f"{copy['score']:.3f}", strand, "."]
    lines.append("\t".join([contig, "gff_export", "gene", *coords,
                             gff_attrs({"ID": gid, "Name": gene, "annotation_tier": "raw", **copy_attrs(copy),
                                        "tool_support": tools})]))
    lines.append("\t".join([contig, "gff_export", feature, *coords,
                             gff_attrs({"ID": fid, "Parent": gid, "Name": gene,
                                        "n_exons": len(chain) if len(chain) > 1 else None})]))
    if len(chain) > 1:
        for k, l in enumerate(chain, start=1):
            lines.append("\t".join([contig, "gff_export", "exon", str(l.start + 1), str(l.end), f"{l.score:.3f}",
                                     strand, ".", gff_attrs({"ID": f"exon-{gene}{suffix}-{k}", "Parent": fid,
                                                             "exon_number": k})]))


def emit_protein_coding_gene(lines: list[str], gene: str, suffix: str, tier: str, parts: list[dict],
                              edits: list[dict], gene_score: float, n_edits_summary: int, copy: dict,
                              phase_from_lengths: bool) -> None:
    """parts: list of {contig, start, end, strand, score, exon_number} in exon order.
    edits: list of {genomic_pos, resulting_aa, [slot]} to attach as sequence_alteration
    children - `slot` (trans_splicing) names the exon, and so the contig, an edit sits on.

    A trans-spliced gene's exons can sit on different contigs; its gene/mRNA
    then get one row per contig (same ID - GFF3's discontinuous-feature
    mechanism), each spanning only that contig's exons.

    phase_from_lengths: cis-spliced raw models get each CDS's phase from the
    cumulative length of the exons before it - internally consistent, but
    only as exact as the HMM-envelope exon boundaries (flagged
    exon_boundaries=approximate). A part's own "phase" (trans_splicing's,
    from transsplice >= 0.2.0) wins; otherwise editing's single exons are
    codon-aligned by construction, so 0 is exact there."""
    strands = {p["strand"] for p in parts}
    gene_strand = strands.pop() if len(strands) == 1 else "."  # mixed-strand for real trans-spliced genes
    gid, mid, cid = f"gene-{gene}{suffix}", f"mRNA-{gene}{suffix}", f"cds-{gene}{suffix}"
    approximate = phase_from_lengths and len(parts) > 1

    by_contig: dict[str, list[dict]] = defaultdict(list)
    for p in parts:
        by_contig[p["contig"]].append(p)
    for contig, cparts in by_contig.items():
        span = [str(min(p["start"] for p in cparts) + 1), str(max(p["end"] for p in cparts))]
        lines.append("\t".join([contig, "gff_export", "gene", *span, f"{gene_score:.3f}", gene_strand, ".",
                                 gff_attrs({"ID": gid, "Name": gene, "annotation_tier": tier, **copy_attrs(copy)})]))
        lines.append("\t".join([contig, "gff_export", "mRNA", *span, f"{gene_score:.3f}", gene_strand, ".",
                                 gff_attrs({"ID": mid, "Parent": gid, "Name": gene,
                                            "n_exons": len(parts), "n_edits": n_edits_summary,
                                            "exon_boundaries": "approximate" if approximate else None})]))
    upstream = 0
    for p in parts:
        if p.get("phase") is not None:
            phase = p["phase"]
        else:
            phase = (3 - upstream % 3) % 3 if phase_from_lengths else 0
        upstream += p["end"] - p["start"]
        lines.append("\t".join([p["contig"], "gff_export", "CDS", str(p["start"] + 1), str(p["end"]),
                                 f"{p['score']:.3f}", p["strand"], str(phase),
                                 gff_attrs({"ID": cid, "Parent": mid, "exon_number": p["exon_number"]})]))
    by_exon = {p["exon_number"]: p for p in parts}
    for i, e in enumerate(edits, start=1):
        pos = int(e["genomic_pos"])
        exon = by_exon.get(int(e["slot"]), parts[0]) if pd.notna(e.get("slot")) else parts[0]
        lines.append("\t".join([exon["contig"], "gff_export", "sequence_alteration", str(pos + 1), str(pos + 1),
                                 ".", exon["strand"], ".",
                                 gff_attrs({"ID": f"sa-{gene}{suffix}-{i}", "Parent": cid,
                                            "edited_from": "C", "edited_to": "T",
                                            "resulting_aa": e.get("resulting_aa", "")})]))


def overlaps(a_contig, a_start, a_end, b_contig, b_start, b_end) -> bool:
    return a_contig == b_contig and a_start < b_end and b_start < a_end


def build_species_gff(hits: pd.DataFrame, organelle: str, editing_calls: pd.DataFrame,
                       editing_edits: pd.DataFrame, ts_genes: pd.DataFrame, ts_exons: pd.DataFrame,
                       ts_edits: pd.DataFrame, stats: Counter, expected: dict) -> list[str]:
    """Every table must already be restricted to one species AND one
    organelle - genes like atpA/rpl16/rps3 exist in both genomes, and
    matching on species+gene alone let one organelle's editing call or
    reconstruction leak into the other's GFF.

    Every hit is used, not just the best per gene: overlapping hits
    collapse to one locus, protein-coding loci chain into cis-spliced
    genes, and every resulting locus/chain is emitted as its own copy
    (`gene-X`, `gene-X-2`, ...; strongest first)."""
    lines: list[str] = []
    is_trna = hits.gene.str.startswith("trn")
    hits = pd.concat([hits[~is_trna], resolve_trna_identity(hits[is_trna])])

    for gene, gh in hits.groupby("gene"):
        loci = collapse_hits(gh)
        is_rna = gene.startswith(("trn", "rrn"))
        models = [[l] for l in loci] if gene.startswith("rrn") else chain_exons(loci, organelle)
        copies = [{"loci": m, "score": sum(l.score for l in m),
                   "coverage": model_coverage(m, expected.get((organelle, gene)))} for m in models]
        for c in copies:
            c["partial"] = c["coverage"] is not None and c["coverage"] < PARTIAL_COVERAGE
        if gene.startswith("trn"):
            # the strongest locus always stays (a split plastid tRNA tRNAscan
            # misses must not vanish); further copies need both tools' support
            copies.sort(key=lambda c: -c["score"])
            supported = copies[:1] + [c for c in copies[1:] if trna_supported(c)]
            stats["trna_copies_unsupported_dropped"] += len(copies) - len(supported)
            copies = supported

        # Corrections come first and claim the loci they were computed on.
        corrected: list[dict] = []
        if gene in RECONSTRUCTED_GENES:
            g = ts_genes[ts_genes.gene == gene]
            if not g.empty and int(g.iloc[0].n_filled) >= 1:
                exons = ts_exons[ts_exons.gene == gene].sort_values("slot")
                # transsplice >= 0.2.0 reports each exon's CDS phase (split codons)
                parts = [{"contig": e.contig, "start": int(e.start), "end": int(e.end), "strand": e.strand,
                          "score": e.score, "exon_number": int(e.slot),
                          "phase": int(e.phase) if "phase" in exons.columns and pd.notna(e.phase) else None}
                         for e in exons.itertuples()]
                claimed = [c for c in copies if any(
                    overlaps(l.contig, l.start, l.end, p["contig"], p["start"], p["end"])
                    for l in c["loci"] for p in parts)]
                raw_loci = [l for c in copies for l in c["loci"]]
                inconsistent = gene in CIS_SPLICED_MITO and not all(
                    any(overlaps(l.contig, l.start, l.end, p["contig"], p["start"], p["end"])
                        and p["start"] >= l.start - CIS_MAX_EXTENSION and p["end"] <= l.end + CIS_MAX_EXTENSION
                        for l in raw_loci)
                    for p in parts)
                if inconsistent:
                    stats["cis_reconstruction_rejected_inconsistent_with_hits"] += 1
                    parts = None
                elif any(len(c["loci"]) > len(parts) for c in claimed):
                    # e.g. rps3 is cis-spliced in most plants: a 2-exon cis
                    # model beats a reconstruction that filled only 1 slot
                    stats["reconstruction_superseded_by_cis_model"] += 1
                    parts = None
                else:
                    copies = [c for c in copies if c not in claimed]
            if not g.empty and int(g.iloc[0].n_filled) >= 1 and parts is not None:
                corrected.append({"tier": "reconstructed", "parts": parts,
                                  "edits": ts_edits[ts_edits.gene == gene].to_dict("records"),
                                  "score": float(g.iloc[0].whole_gene_score),
                                  "n_edits": int(g.iloc[0].whole_gene_n_edits), "coverage": None, "partial": False})
        e = editing_calls[editing_calls.gene == gene]
        if not corrected and not e.empty:
            erow = e.iloc[0]
            # editing corrected the single best-scoring hit; find the copy holding it
            target = next((c for c in copies if any(
                l.contig == h.contig_id and l.start <= h.start and h.end <= l.end
                for l in c["loci"] for h in [gh.loc[gh.score.idxmax()]])), None)
            if target is not None and len(target["loci"]) == 1:
                l = target["loci"][0]
                corrected.append({"tier": "edited", "score": float(erow.score), "n_edits": int(erow.n_edits),
                                  "parts": [{"contig": l.contig, "start": int(erow.corrected_start),
                                             "end": int(erow.corrected_end), "strand": l.strand,
                                             "score": erow.score, "exon_number": 1}],
                                  "edits": editing_edits[editing_edits.gene == gene].to_dict("records"),
                                  "coverage": target["coverage"], "partial": target["partial"]})
                copies.remove(target)
            elif target is not None:
                # an ORF search over one exon of a cis-spliced gene can't
                # correct the gene's boundaries - keep the multi-exon model raw
                stats["editing_skipped_multi_exon"] += 1

        ordered = corrected + sorted(copies, key=lambda c: -c["score"])
        n_full = sum(1 for c in ordered if not c["partial"])
        for i, c in enumerate(ordered, start=1):
            c["copy_number"] = n_full
            suffix = "" if i == 1 else f"-{i}"
            stats["copies"] += 1
            stats["partial_copies"] += c["partial"]
            if is_rna:
                emit_rna_gene(lines, gene, suffix, c["loci"], "tRNA" if gene.startswith("trn") else "rRNA", c)
                continue
            if "tier" in c:
                emit_protein_coding_gene(lines, gene, suffix, c["tier"], c["parts"], c["edits"], c["score"],
                                          c["n_edits"], c, phase_from_lengths=False)
                stats[c["tier"]] += 1
                continue
            parts = [{"contig": l.contig, "start": l.start, "end": l.end, "strand": l.strand, "score": l.score,
                      "exon_number": k} for k, l in enumerate(c["loci"], start=1)]
            emit_protein_coding_gene(lines, gene, suffix, "raw", parts, [], c["score"], 0, c,
                                      phase_from_lengths=True)
            stats["raw"] += 1
            stats["multi_exon_raw"] += len(parts) > 1

    def sort_key(line: str):
        f = line.split("\t")
        return (f[0], int(f[3]))
    lines.sort(key=sort_key)
    return lines


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organelle", default="both", choices=["mito", "pltd", "both"])
    ap.add_argument("--species-list", default=None)
    ap.add_argument("--gene-calls", default=None, help="raw gene_calls.tsv to use for tier-3 fallback + "
                     "tRNA/rRNA + contig/strand lookups (default: denovo_annotation/results/gene_calls.tsv)")
    args = ap.parse_args()

    species_filter = sd.load_species_list(Path(args.species_list)) if args.species_list else None

    raw_path = Path(args.gene_calls) if args.gene_calls else ANALYSIS_DIR / "denovo_annotation" / "results" / "gene_calls.tsv"
    if not raw_path.exists():
        print(f"[err] missing {raw_path} - run denovo_annotation/src/04_build_gene_calls.py first "
              f"(or pass --gene-calls)", file=sys.stderr)
        sys.exit(1)
    raw = pd.read_csv(raw_path, sep="\t")
    if "source" not in raw.columns:
        # a gene_calls.tsv from before the model columns existed: loci still
        # collapse, but nothing chains into cis-spliced genes without them
        print(f"[warn] {raw_path} has no source/hmm_from/hmm_to/model_len columns - "
              f"re-run denovo_annotation/src/04_build_gene_calls.py", file=sys.stderr)
        raw = raw.assign(source="oatkDB", hmm_from=np.nan, hmm_to=np.nan, model_len=np.nan)
    if species_filter is not None:
        raw = raw[raw.species.isin(species_filter)]

    def load(path: Path, columns: list[str]) -> pd.DataFrame:
        # Always returns a DataFrame with `columns` present (empty if the
        # upstream file doesn't exist yet, e.g. plastid before editing/
        # trans_splicing ever cover it) so every `.species`/`.gene` filter
        # downstream is safe regardless of which modules have been run.
        if not path.exists():
            return pd.DataFrame(columns=columns)
        df = pd.read_csv(path, sep="\t")
        return df[df.species.isin(species_filter)] if species_filter is not None and not df.empty else df

    editing_calls = load(ANALYSIS_DIR / "editing" / "results" / "edited_gene_calls.tsv",
                          ["species", "organelle", "gene", "corrected_start", "corrected_end", "score", "n_edits"])
    editing_edits = load(ANALYSIS_DIR / "editing" / "results" / "predicted_edits.tsv",
                          ["species", "organelle", "gene", "genomic_pos", "resulting_aa"])
    ts_genes = load(ANALYSIS_DIR / "trans_splicing" / "results" / "reconstructed_genes.tsv",
                     ["species", "organelle", "gene", "n_filled", "whole_gene_score", "whole_gene_n_edits"])
    ts_exons = load(ANALYSIS_DIR / "trans_splicing" / "results" / "reconstructed_exons.tsv",
                     ["species", "organelle", "gene", "slot", "contig", "start", "end", "strand", "score"])
    ts_edits = load(ANALYSIS_DIR / "trans_splicing" / "results" / "reconstructed_edits.tsv",
                     ["species", "organelle", "gene", "slot", "genomic_pos", "resulting_aa"])

    organelles = ["mito", "pltd"] if args.organelle == "both" else [args.organelle]
    out_dir = ANALYSIS_DIR / "gff_export" / "results"
    unitig_out_dir = out_dir / "unitig"
    unitig_out_dir.mkdir(parents=True, exist_ok=True)

    tables = [editing_calls, editing_edits, ts_genes, ts_exons, ts_edits]
    empty = [t.iloc[0:0] for t in tables]
    split = [{k: g for k, g in t.groupby(["species", "organelle"])} if not t.empty else {} for t in tables]
    raw_split = {k: g for k, g in raw.groupby(["species", "organelle"])}
    # from the full table, so a --species-list run gets the same coverages
    expected = expected_model_spans(pd.read_csv(raw_path, sep="\t") if species_filter is not None else raw) \
        if "hmm_from" in raw.columns and raw.hmm_from.notna().any() else {}
    stats: Counter = Counter()

    n_written = n_unitig = n_no_fasta = n_noted = 0
    written: set[Path] = set()
    linearisation_notes = load_linearisation_notes()
    for organelle in organelles:
        species_here = sorted(raw[raw.organelle == organelle].species.unique())
        ctg_fastas = {r.species: r.ctg_fasta for r in
                      sd.discover_all(sd.repo_data_root(ROOT_DIR, organelle), organelle,
                                      species_filter=set(species_here))
                      if r.status == "ok"}
        for species in species_here:
            key = (species, organelle)
            lines = build_species_gff(raw_split[key], organelle, *[s.get(key, e) for s, e in zip(split, empty)],
                                      stats=stats, expected=expected)
            if not lines:
                continue
            out_path = out_dir / f"{species}.{organelle}.gff"
            if species not in ctg_fastas:
                n_no_fasta += 1
                out_path.write_text("##gff-version 3\n" + "\n".join(lines) + "\n")
                written.add(out_path)
                n_written += 1
                continue

            ctx = load_seq_context(Path(ctg_fastas[species]), UNITIG_MAP_DIR / f"{species}.{organelle}.tsv")
            ctg_lines, unitig_lines = add_unitig_layer(lines, ctx)
            low_frac = linearisation_notes.get(species) if organelle == "mito" else None
            n_noted += low_frac is not None
            out_path.write_text(assemble(*ctg_header(ctx, list(ctx.lengths), low_frac), ctg_lines))
            written.add(out_path)
            n_written += 1
            if unitig_lines:
                unitigs = sorted({l.split("\t", 1)[0] for l in unitig_lines}, key=lambda u: (len(u), u))
                unitig_path = unitig_out_dir / f"{species}.{organelle}.unitig.gff"
                unitig_path.write_text(assemble(*unitig_header(ctx, unitigs, low_frac), unitig_lines))
                written.add(unitig_path)
                n_unitig += 1
    # A full run defines the result set: a GFF from an earlier run for a
    # species no longer annotated (e.g. dropped from QC-pass, or its contigs
    # replaced) describes contigs that may not exist any more - remove it.
    n_stale = 0
    if species_filter is None:
        for f in list(out_dir.glob("*.gff")) + list(unitig_out_dir.glob("*.unitig.gff")):
            if f not in written:
                f.unlink()
                n_stale += 1
    print(f"[info] build_gff: removed {n_stale} stale GFFs from earlier runs", file=sys.stderr)
    print(f"[info] build_gff: {n_written} species -> {out_dir} ({n_unitig} with unitig-level GFF -> "
          f"{unitig_out_dir}; {n_no_fasta} without a resolvable ctg.fasta, no region/unitig layer; "
          f"{n_noted} marked linearisation=uncertain)",
          file=sys.stderr)
    print("[info] build_gff: " + " ".join(f"{k}={v}" for k, v in sorted(stats.items())), file=sys.stderr)


if __name__ == "__main__":
    main()
