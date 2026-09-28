# repeats

**The question:** plant mitochondrial genomes are often "multipartite" —
not one simple circle but several chromosomes linked by recombination
across repeated sequence. Where are the repeats that could be mediating
this, and does repeat content predict the structural complexity we
already see in the assembly graph?

Reimplemented fresh (not ported) from the cleanest idea in the sibling
`mito_structural_variation` repo's repeat-finding code (see Wynn &
Christensen 2019 for the biology): self-alignment, filter, cluster. Plant
organelle contigs (tens of kb to ~1.5 Mb) make a full self-`blastn` cheap
— no genome-scale masking or indexing needed.

## Method

1. **`01_find_repeats.py`** — self-`blastn` each species' resolved
   contigs against themselves (short-seed parameters: `-word_size 50
   -reward 1 -penalty -4`, no masking), drop the trivial whole-sequence
   self-hit and mirrored duplicate pairs, then **merge overlapping/nested
   hit intervals into their maximal block** (blastn reports many
   near-duplicate overlapping HSPs for what is really one repeat locus —
   without merging, a single repeat region shows up as dozens of
   near-identical rows a few bp apart). The resulting maximal blocks are
   clustered into repeat families with `cd-hit-est` (90% identity).
2. **`02_flag_recomb_repeats.py`** — flags copies within 2kb (configurable
   via `--end-distance`) of a contig end as `putative_recomb_repeat`. A
   repeat near the end of a **linear** contig is a real candidate
   recombination point; near the end of a **circular** contig, "the end"
   is just an arbitrary linearisation choice, so it's flagged with a
   weaker note rather than treated the same way.
3. **`03_repeats_summary.py`** → `recomb_repeats_summary.tsv`: one row per
   species×organelle (family/copy counts, % of genome repetitive,
   putative recombination-repeat count).

## Reading the output

`per_species/<Species>.<organelle>.repeats.tsv`: one row per repeat
**block** (not per pairwise hit — a family with N copies produces N rows,
one per copy, all sharing a `repeat_family` ID). `strand` is always `+`
here: it describes the interval on the contig's own coordinates, not
orientation relative to a partner copy — a copy can be direct or inverted
relative to different partners, which is exactly what makes flip-flop
repeats interesting, but isn't a property of one copy in isolation.

## Try this

Join `recomb_repeats_summary.tsv` against `qc_basic_stats/results/gfa_stats.tsv`'s
`dead_end_nodes`/`n_subgraphs` (both keyed on `species`+`organelle`) — does
repeat content predict assembly graph complexity? Plastid genomes should
show a very consistent `pct_genome_repetitive` around 35-40% across almost
every species — that's the well-known large inverted repeat (IR); look for
outliers.
