# topology_plots

**The question:** which flagged/failed assemblies are worth a human look, and
what does the actual assembly graph look like? `qc_basic_stats` tells you
*that* something looks off (dead ends, size outlier, elevated cross-organelle
alignment, ...); this module lets you *see* it — a gene-labeled picture of
the assembly graph, generated automatically, no manual editing required.

Reproduces the approach in the sibling `mito_structural_variation` repo's
`figures/fig1b_genome_topology_examples/` (its own README documents the same
method) against this repo's data layout, fully automated end to end.

## Method

1. **`01_make_gene_queries.py`** — for each species×organelle you ask for
   (default: everything with QC status `fail` or `flag` — the point of this
   tool is triaging what's *already* flagged, not re-looking at everything),
   builds a small BLAST query FASTA: one sequence per gene name (the longest
   interval if a gene has multiple hits, reverse-complemented for `-`
   strand), skipping tRNAs (too small/numerous to label legibly). This only
   needs **one representative sequence per gene** — Bandage's own internal
   BLAST search finds every location that sequence's gene occurs in the
   graph, including duplicated copies on other nodes/contigs.
2. **`02_render_topology.py`** — runs `Bandage image` once per genome:

   ```bash
   Bandage image <species>.gfa <species>.genes.png --height 900 \
     --query <species>.genes.fasta --colour blastsolid --blasthits --fontsize 12
   ```

   This is what does the BLAST — Bandage runs `makeblastdb`/`blastn`
   internally against the query FASTA, colors each node by its best hit,
   and labels it with the gene name. (An earlier version also rendered a
   plain, no-gene-labels variant, matching the reference figures'
   `<topology>` + `<topology>_genes` pair — dropped as pure duplication:
   the gene-labeled plot shows the same topology plus gene positions, so
   the plain one added nothing.)

   Output is **PNG by default** (`--format png`, viewable in any image viewer,
   about half the file size of the equivalent SVG here) — Bandage's `image`
   command writes PNG/JPG/SVG natively, so no separate SVG→PNG conversion
   step (e.g. `resvg`) is needed. Pass `--format svg` for vector output if
   you want it (matches the reference figures' own `_genes.svg` naming).

## Setup

Needs `Bandage` (the original Qt-based tool — **not** BandageNG, which lacks
gene-query BLAST coloring) with a bundled BLAST+. There's no bioconda access
from this cluster's compute nodes, so this was installed as a self-contained
conda env at `/nfs/users/nfs_m/mb39/miniconda3/envs/bandage` (copied from a
working install rather than freshly `conda create`d) — see
`analysis/common/tool_paths.sh` for the path and the `QT_QPA_PLATFORM=offscreen`
env var it needs for headless rendering (no Xvfb required). There's also a
Bandage singularity image on this cluster
(`/software/tola/images/bandage-0.8.1.sif`) used elsewhere for quick
unlabeled PNGs — **don't use it here**, its headless SVG text rendering is
broken (glyphs silently fail to shape).

## Known limitation: unseeded layout

Bandage's force-directed graph layout has no seed/reuse option — two renders
of the same GFA can land on visibly different (but equally valid)
arrangements. This only matters for a genuinely tangled multi-node
single-component graph; single-node and clearly-separate multipartite
topologies have nothing to tangle. If a specific plot looks messy, just
re-render that one species with `--force` — you may land on a cleaner
layout. This isn't hidden or hand-fixed here; if you want to compare against
this repo's plots, that's the one thing they won't reproduce pixel-for-pixel.

## Reading the output

`results/plots/{mito,pltd}/<species>.genes.png` — open in any image viewer.
Node color = best BLAST hit; label = gene name. A normal, complete genome should show sensible
core-gene coverage; dead ends, tangled multi-node components, or an
unexpectedly large fraction of a plastid's genes' sequences turning up when
BLASTed against the *mito* graph (worth checking manually, this module
doesn't cross-check automatically) are the kinds of things worth a look for
a `fail`/`flag` species.

Validated examples already worth knowing about:
- `Juncus_inflexus` (flagged for `high_cross_organelle_alignment`) has a
  completely normal, gene-complete, textbook plastid graph — nothing broken
  about *that* assembly. The flag more likely reflects a genuinely large
  MTPT event landing in its mitochondrial assembly, not a plastid assembly
  problem — exactly the kind of distinction a QC number alone can't tell
  you but a picture can.
- `Clematis_viticella` visibly shows what its `n_subgraphs=10` QC value
  means: lots of small disconnected linear segments, not a coherent genome
  — this is exactly what led to adding the fragmentation check in
  `qc_basic_stats` (nothing checked `n_subgraphs` directly before it was
  spotted here).
- `Cerastium_semidecandrum` (27 dead-end nodes) visibly shows a tangled,
  unresolved graph, consistent with its `fail` status.
