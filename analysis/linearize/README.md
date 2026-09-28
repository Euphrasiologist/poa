# linearize

**The question:** every assembly in this dataset started as a graph
(`data/{mito,plastid}/<Species>/*.gfa`), not a genome — oatk's own
Pathfinder tries to resolve one linear (or circular) sequence out of that
graph, and for most species it succeeds (that's `*.ctg.fasta`). For the
species where it doesn't (`no_resolved_ctg_fasta` in `qc_basic_stats`),
there was no linearized genome at all. This module adds two independent,
increasingly capable ways to linearize *any* assembly graph, on top of
(not instead of) oatk's own resolution.

## Two tiers

### Tier 1 — `gfatk linear` (universal, cheap, graph-only)

`01_gfatk_linear.py` runs `gfatk linear -e` on every species x organelle
that has a GFA at all — seconds each, the whole ~1250-species dataset in
a few minutes, no new dependencies, no raw read access. It greedily picks
the highest-cumulative-coverage path through each subgraph. This is a
genuinely independent second opinion on path resolution, not a rerun of
oatk's own logic — verified on `Acaena_ovalifolia` (QC status `pass`,
`any_circular=True`): oatk's own gene-annotation-aware Pathfinder calls it
circular, but `gfatk linear` (no annotation weighting) picks a *linear*
2-node path and leaves a 112kb second segment as an orphan. Treat tier 1
as a second data point to compare against `data/*/*.ctg.fasta`, not a
strictly-better replacement for it.

Output:
- `results/linear/<organelle>/<species>.linear.fasta` — one sequence per
  subgraph (the chosen best path; see parsing note below for what counts
  as "the" path in a tangled subgraph).
- `results/linear_summary.tsv` — one row per species x organelle x
  subgraph: `path`, `n_segments`, `coverage`, `is_circular`, `length_bp`,
  and `n_orphan_nodes`/`orphan_bp`/`pct_subgraph_explained` — how much of
  that subgraph's sequence the chosen path actually accounts for. A low
  `pct_subgraph_explained` is itself a QC signal: `Cerastium_semidecandrum`
  (27 dead-end nodes, `fail` status) explains only 2% of its main
  subgraph this way — the graph is genuinely too tangled for a
  coverage-greedy heuristic to make sense of, consistent with its status.

Run it: `python3 analysis/linearize/src/01_gfatk_linear.py --organelle both`
(add `--species-list`/`--force` as usual).

### Tier 2 — `gfatk resolve` fallback (targeted, uses real read evidence)

For species where tier 1 *also* isn't good enough — scoped, by design, to
exactly `no_resolved_ctg_fasta` species (oatk's Pathfinder produced
nothing there at all) — `gfatk resolve` uses actual PacBio HiFi read-path
evidence (aligned to the graph with `GraphAligner`) to pick a circuit,
which can succeed where a pure coverage/topology heuristic can't.

**Why this is scoped narrowly, and why it needs a recruitment step first:**
whole-genome raw HiFi read sets for this dataset total **~24TB compressed**
(verified directly: `du` across all paths in `meta/file.txt`, ~3250 files,
mean 7.6GB each). Running GraphAligner against a species' full read set
just to align the handful of MB that actually belong to one small
organelle graph would be a huge, unjustified compute cost for what's meant
to be a lightweight teaching resource. So before GraphAligner ever runs,
`_resolve_one.sh` recruits only the organelle-relevant reads:

1. `gfatk fasta <gfa>` → segment-level reference (works even with no
   `ctg.fasta` — exactly this tier's target population).
2. `minimap2 -a -x map-hifi` against that reference, piped straight into
   `samtools fasta -F 0x904` (drop unmapped/secondary/supplementary) —
   recruitment and sequence extraction in **one pass** over the raw reads.
   (An earlier version ran `minimap2 -x map-hifi` for a PAF, then a
   second full `zcat | seqkit grep` pass to pull the matching reads -
   measured directly on `Carlina_vulgaris` (~30GB compressed reads):
   33.5min + 31.7min ≈ 65min, almost entirely spent re-decompressing the
   same reads twice - minimap2 itself only got ~2.7x parallelism out of
   16 threads, consistent with gzip decompression being the bottleneck,
   not alignment. Piping cuts that to one pass.)
3. `GraphAligner -x vg` aligns the recruited reads to the graph (GAF
   output) - measured at ~20min for ~412k recruited reads, <1GB RAM.
4. `gfatk resolve --gaf ...` uses that read-path evidence to resolve a
   circuit - effectively instant (<1s) once the GAF exists.

Recruiting a surprisingly large fraction of a species' total reads is
*expected*, not a bug: organelle genomes exist in hundreds of copies per
cell, so a large share of total whole-genome HiFi yield is genuinely
organelle-derived even against a reference that's only ~100-200kb
(verified: 412,256 reads recruited for `Carlina_vulgaris`, 85% at mapq=60,
i.e. high-confidence unique alignments, not noise).

**Validating result** (`Carlina_vulgaris`, plastid, `no_resolved_ctg_fasta`):
oatk's Pathfinder produced nothing; `gfatk linear` found *a* path (130kb,
3 segments, not circular); `gfatk resolve`, using real read-path evidence,
resolved a genuinely **circular** 153kb genome where one segment (`u16`)
correctly appears twice in the path (`u17+,u16-,u15+,u16+`) - the classic
plastid inverted-repeat architecture, which a topology-only heuristic has
no way to identify but read evidence can confirm directly. This is
exactly the case this tier exists for.

Scope grew from 46 to **124 species** partway through this tier's own
rollout: testing tier 2 on `Azolla_filiculoides` (prompted by a direct
question about whether ferns would linearize properly) surfaced a
`qc_basic_stats` blind spot - 78 species where oatk's `.ctg.fasta` is
non-empty but every single contig is an unjoined singleton node, which
`has_ctg_fasta` was wrongly treating as resolved (see
[`../qc_basic_stats/README.md`](../qc_basic_stats/README.md)'s
"`no_resolved_ctg_fasta(unjoined)`" section for the full story). Fixing
that check automatically folded those 78 - mostly fern mitogenomes,
exactly the repeat-rich case read-path evidence is best suited to - into
this tier's target list with no code changes needed here.

Steps:
```bash
python3 analysis/linearize/src/02_select_resolve_targets.py   # -> work/resolve_targets.tsv
bash analysis/linearize/src/03_recruit_and_resolve.sh --queue basement --threads 16
python3 analysis/linearize/src/04_summarize_resolve.py        # -> results/resolve_summary.tsv
```

`03_recruit_and_resolve.sh` names each bsub job `linearize_resolve.<species>.<organelle>`
and checks `bjobs -J` for that name before submitting, not just the output
file - an in-flight job hasn't written its output yet, so checking the
file alone let a rerun (e.g. after `--targets` grows, as it did here)
double-submit 35 genuine duplicate jobs the first time this happened.

`02_select_resolve_targets.py --qc-reason <other reason>` can target a
different QC failure category, but `no_resolved_ctg_fasta` is the
justified default — for other flag reasons (dead-ends, fragmentation,
etc.) oatk *did* produce a resolved sequence, so the read-recruitment cost
is much harder to justify without first checking whether a cheaper fix
(coverage-adjusted reassembly, see the top-level README's "Revising the
outputs" section) already works.

Output:
- `results/resolve/<organelle>/<species>.resolve.fasta` — one sequence
  per resolved circuit, gfatk resolve's own descriptive headers kept as-is.
- `results/resolve_summary.tsv` — one row per resolved circuit
  (`n_circuits`, `circuit_id`, `n_segments`, `length_bp`); a species with
  `n_circuits=0` means it was attempted but nothing resolved (distinct
  from a species not yet attempted at all, which won't appear in the
  table - cross-reference against `work/resolve_targets.tsv`).

`work/` holds the large per-species intermediates (segment references,
recruited reads, GAFs) - not committed (see `.gitignore`), safe to delete
once `results/` looks right; a rerun regenerates them.

## Resource footprint (per species x organelle, tier 2 only)

Measured directly, not estimated: ~30-65min wall time dominated by
recruitment (proportional to that species' total raw HiFi read volume,
not organelle size), <2GB peak RAM either step. Submitted as independent
parallel LSF jobs (`03_recruit_and_resolve.sh`, one bsub per target), so
aggregate wall-clock is bounded by the slowest single job, not the sum -
consistent with how this repo's other per-species batch work
(reassembly/promotion, see top-level README) is run.
