# qc_basic_stats

**The question:** is this assembly real and complete enough to trust? With
~1250 species assembled semi-automatically, some fraction will be empty,
fragmented, or otherwise suspect — before building an exercise (or a
phylogeny!) on top of one, it's worth knowing that up front rather than
discovering it later as an unexplained outlier.

This module computes per-species metrics and combines them into a single
QC gate that every other module filters against by default.

## What it computes

| script | output | what |
|---|---|---|
| `01_gfa_stats.py` | `gfa_stats.tsv` | `gfatk stats` per subgraph: node/edge counts, dead-ends, length, GC%, coverage. A genome can legitimately have >1 subgraph (multipartite mitochondria are real biology). |
| `02_contig_stats.py` | `contig_stats.tsv` | Per-contig length/circularity/GC% from `.ctg.fasta` headers + `seqkit`. Falls back to `gfatools gfa2fa` (unitig sequences, not a resolved path) for the ~5% of plastid assemblies missing a `.ctg.fasta`, labelled `contig_source=gfa_fallback`. |
| `03_gene_matrix.py` | `gene_matrix_{mito,pltd}.tsv` | Long-format `(species, gene) → hits/score/aligned length/n_contigs`, straight from oatk's own `.ctg.bed` gene calls. |
| `04_core_gene_list.py` | `core_genes_{mito,pltd}.tsv` | Per-gene prevalence across non-empty assemblies — **computed from the data, never a hardcoded gene list.** |
| `05_qc_summary.py` | `qc_summary.tsv` | Combines everything into one `pass`/`flag`/`fail` per species×organelle, with a `status_reasons` column explaining why. |

## Reading `qc_summary.tsv`

Each check is either a **flag** (worth a look, not disqualifying) or a
**fail** (don't trust this assembly by default):

- **missing/empty assembly** → fail
- **dead-end nodes** in the graph → flag if any, fail if >2
- **fragmentation** (`n_subgraphs`, disconnected components in the
  assembly graph) → flag if ≥5, fail if ≥10. The dataset shows a clean
  natural break: 1154/1250 mito species have 1-4 subgraphs (consistent
  with real single-chromosome or modest multipartite biology — the
  literature describes real multipartite plant mitogenomes as essentially
  always 2-4 sub-genomic circles), then a long tail of 46 species with
  5-20. Confirmed against `topology_plots` output: e.g. `Clematis_viticella`
  (10 subgraphs) visibly shows "lots of small linear segments" — a failed
  assembly, not previously caught since nothing checked `n_subgraphs`
  before this.
- **non-circular** → flag only (many real plant mitochondria are
  legitimately multipartite/linear — this is informational, not a defect).
  Note: this reads the *resolved contig's* `circular=` flag from
  `contig_stats.tsv`, not `gfatk`'s raw-graph-level flag — verified these
  can disagree (a plastid's raw unitig graph can be "non-circular" while
  its resolved path still closes a genuine circle via a repeat, which is
  exactly the normal single-copy-repeat plastid structure).
- **genome size outlier** → flag if more than 4 MADs from the
  genus-level median (falling back to the dataset-wide median when a
  species has fewer than 3 sequenced congeners). Genus-relative because
  organelle genome size varies enormously across plant lineages — a fixed
  absolute cutoff would be meaningless.
- **core-gene completeness** below 80%/50% → flag/fail. The core-gene set
  is whichever genes are present in ≥90% of non-empty assemblies for that
  organelle (see `core_genes_*.tsv`) — not a textbook list, so it
  naturally adapts if oatk's gene family database changes.
- **no resolved contig fasta** (only the `gfa_fallback`, or an
  `unjoined`-tagged variant — see below) → flag.
- **cross-organelle alignment fraction** above 15%/35% → flag/fail. See
  below.

### `no_resolved_ctg_fasta(unjoined)`: a non-empty `.ctg.fasta` isn't proof anything resolved

Found by checking a specific user question ("ferns definitely wouldn't
linearize properly, right?") against `Azolla_filiculoides`: its
`.mito.ctg.fasta` has 175 contigs, but **every one of them is `nv=1`** —
one contig per raw graph segment, meaning Pathfinder never joined a single
pair of nodes into a path. The file being non-empty made `has_ctg_fasta`
true, so this was invisible to the QC gate entirely.

Scanning the whole dataset for the same signature (all contigs `nv=1`,
≥5 contigs — the ≥5 floor matches the fragmentation check's own precedent
that real multipartite plant mitogenomes are 2-4 subgenomic circles, so
5+ never-joined pieces isn't plausible biology) found **78 species**, heavily
concentrated in ferns/pteridophytes (`Vandenboschia_speciosa`,
`Asplenium_marinum`, `Polystichum_lonchitis`, `Dryopteris_filix_mas`,
`Pteridium_aquilinum`, and most other ferns in the dataset) — consistent
with their reputation for large, repeat-rich mitogenomes that trip up
graph-topology heuristics. `has_ctg_fasta` now treats this case as *not*
resolved, same as a genuinely empty file, with `n_contigs`/
`pct_unjoined_contigs` columns added to `qc_summary.tsv` so the two
subtypes stay distinguishable (`no_resolved_ctg_fasta` vs.
`no_resolved_ctg_fasta(unjoined)`).

This had a second, non-obvious effect: the cross-organelle contamination
check was previously running its mito-vs-plastid `minimap2` alignment
against these 78 species' fragmented mito *junk* too (since it only
gated on the old, wrong `has_ctg_fasta`), producing spurious high
alignment fractions against noise — 20 species' **plastid** status
corrected from `flag`/`fail` back to `pass` once their mito side was
correctly excluded from that comparison. Their plastid genomes were never
actually contaminated; they were being compared against garbage.

See [`../linearize/README.md`](../linearize/README.md) for how these 78
species (plus the pre-existing `no_resolved_ctg_fasta` ones) get a second
attempt at linearization.

### Why the contamination check is a minimap2 alignment, not a gene-name check

An earlier version of this check looked for plastid-only gene names
turning up in a species' mito annotation (or vice versa) as a
contamination signal. That was verified empirically to never fire: oatk
only ever searches a bin against its own family HMM database
(`acrogymnospermae_mito.fam` / `angiosperm_pltd.fam`), so a plastid gene
name can never appear in a `.mito.ctg.bed` regardless of what's actually
in the sequence.

Instead, for every species with both organelles assembled, this runs one
`minimap2 -x asm5` alignment between its mito and plastid contigs and
checks what fraction of one organelle's assembly aligns to the other.
Genuine MTPT (mitochondrial-plastid DNA transfer) is a small, isolated
insert — a few percent at most. A large aligned fraction instead suggests
a mis-binned contig. This is a coarse pre-screen; see
[`../synteny/README.md`](../synteny/README.md) for the full per-species
MTPT picture (PAFs, dotplots, a proper summary table).

## Try this

Run with `--qc-status pass,flag,fail` instead of the default `pass` and
compare the annotation/repeats/synteny/phylogeny output — what changes,
and does it match what the QC reasons predicted? This is a genuinely good
exercise in "what does bad input data actually do to downstream analysis,"
not just an abstract warning.
