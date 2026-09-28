# phylogeny

**The question:** given hundreds of assembled organelle genomes, can we
build a species tree from them? Organelle genomes are small and mostly
uniparentally inherited, so a handful of conserved marker genes is enough
— no need for whole-genome alignment.

## Method

1. **`01_select_core_genes.py`** — reads
   [`qc_basic_stats`](../qc_basic_stats/README.md)'s `core_genes_*.tsv`
   and finalizes the marker-gene set: present in ≥90% of species **and**
   single-copy in ≥80% of those. (Single-copy threshold is 80%, not 90% —
   verified at full-dataset scale that the best mito genes top out around
   89% single-copy, e.g. `ccmB` at 89.2%, `nad5` at 88.8%; plastid genes
   comfortably clear 93%+. A 90% bar excluded every mito gene entirely,
   which is a real biological signal (mitochondrial multipartite structure
   duplicates genes far more than plastids do) hitting an arbitrary
   cutoff, not evidence there's no usable mito marker gene.) Also writes
   `excluded_genes_*.txt`
   listing genes that were present but excluded for being multi-copy —
   worth reading, since it's basically a data-driven list of plant
   organellar genes with multiple exons/trans-splicing/duplication, e.g.
   `nad1/2/4/5/7`, `cox2`, `ccmFc/Fn`, `rps3` typically show up here for
   mito. This is real biology, not a limitation to apologize for — naively
   concatenating exon fragments as if they were one gene would produce a
   wrong tree, so those genes are left out rather than mishandled.
2. **`02_align_genes.py`** — `mafft --auto` per marker gene (cheap: a
   couple dozen short genes, no LSF needed).
3. **`03_concat_alignment.py`** — concatenates all per-gene alignments
   into one supermatrix with `AMAS`, writing a RAxML-style partition file.
   Partitioning (rather than naive concatenation) lets IQ-TREE fit a
   separate substitution model per gene — the correct way to combine loci
   evolving at different rates, and a second, independent safeguard around
   the multi-exon caveat above.
4. **`04_build_tree.py`** — one `iqtree2 -m MFP -bb 1000` run **per
   organelle** (not per species) — cheap even at ~1200 taxa, since the
   alignment is a handful of short genes, not whole genomes. `-bb 1000`
   gives ultrafast bootstrap support values on every internal branch for
   free. No dating step (that needs external fossil calibration and
   `treePL`, which is out of scope for a lightweight teaching resource) —
   branch lengths are in substitutions/site, not time.

5. **`05_tree_discordance.py`** (run after both trees exist) — a second,
   independent contamination signal on top of `qc_basic_stats`'s assembly-
   graph checks: for each species in both trees, finds its nearest-neighbor
   *set* (all species tied at the minimum pairwise ML distance — both
   `.mldist` matrices have many exact ties, especially mito with only
   ~9 marker genes vs plastid's ~60, so a single arbitrary nearest neighbor
   would be noisy) and checks whether the two sets share a genus. No
   overlap at all → `discordant=True`: the mito and plastid trees disagree
   about this species' closest relative, which a clean-looking single-
   subgraph assembly can still hide (unlike the graph-topology checks,
   this doesn't care whether the assembly *looks* structurally fine).
   Reuses the `.mldist` files IQ-TREE already wrote - no new tree-building.
   **Not merged into `qc_summary.tsv`**: it needs both trees, which
   themselves depend on a first QC pass having already run, so folding it
   into that same gate would be circular - treat `tree_discordance.tsv` as
   a second, independent lens, not a replacement for the main gate.
   Validated on real data: only 3/612 species flagged, one of them
   (`Ophrys_apifera`) independently also has 17 mito subgraphs (a hard
   `fail` under the fragmentation check) - two unrelated methods agreeing
   this specific assembly isn't trustworthy is a good sign both are doing
   their job rather than just generating noise.

## Reading the tree

`results/trees/{mito,pltd}.treefile` is a standard Newick tree; the
numbers on internal branches are ultrafast bootstrap support (0-100, not a
percentage of a fixed total — values ≥95 are conventionally treated as
well-supported). Load it in R (`ape::read.tree`), Python (`ete3`/`dendropy`),
or a GUI viewer like FigTree.

## If the species set changes

Re-run the whole module (extraction → alignment → concatenation → tree)
after any change to `annotation/results/genes/` — there's no fine-grained
incremental realignment here (not worth the complexity at this scale), so
a stale alignment missing newly-added taxa would silently produce a
tree that's out of date rather than wrong, but it's easy to forget.
`analysis/run_all.sh` always runs the full chain in order, so this is only
a concern if you invoke the phylogeny scripts directly and skip a step.

## Try this

Compare the mito and plastid trees for the same species set — plant
mitochondrial and plastid genomes are both maternally inherited in most
angiosperms, so you'd expect broadly congruent topologies; where they
disagree is often either weak signal (check bootstrap support) or a
genuinely interesting biological signal (e.g. hybridisation,
introgression, or an assembly-quality problem — check that species'
`qc_summary.tsv` row first).
