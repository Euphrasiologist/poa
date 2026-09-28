# plant_organelle_annotator

A gold-standard, `oatk`-independent annotation pipeline for plant
mitochondrial and plastid genomes: oatk-independent core gene calling
(direct `nhmmscan` against oatkDB's own gene-family HMMs, plus tRNA/rRNA
calling), RNA-editing-aware ORF finding, trans-spliced gene reconstruction
(`nad1`/`nad2`/`nad5`/`rps3`), and a merged GFF3 export - plus a wider
suite of QC, repeat, synteny, phylogeny and non-core ORF (TE/mitovirus)
analyses.

This repo is code only - no data, no per-dataset results. It was split out
of a specific 1250-species dataset repo so it can be installed and run
against *any* set of oatk-assembled organelle genomes, not just that one.

## Install

```bash
mamba env create -f environment.yml
conda activate plant_organelle_annotator
./install.sh   # cargo-installs this project's own Rust tools
```

See `environment.yml` for exactly what's pinned (minimap2, samtools,
seqkit, blast, cd-hit, mafft, iqtree, gfatools, hmmer, tRNAscan-SE,
barrnap, Bandage, GraphAligner, AMAS, rust) and `install.sh` for this
project's own Rust tools (`gfatk`, `orfedit`, `transsplice`, `hmm_to_gff`,
`filter_tblout`).

Reference data (oatkDB's mito/plastid gene-family `.fam` databases, and
the 5 curated GenBank genomes `trans_splicing` bootstraps its
trans-splice-junction profiles from) isn't bundled here yet - see
`analysis/denovo_annotation/README.md` and `analysis/trans_splicing/README.md`
for how to build them yourself in the meantime.

## Pointing this at your data

Every script expects a **data root**: a directory containing `data/`
(your oatk-assembled `<species>/*.ctg.fasta`+`.gfa` per species, laid out
as `data/mito/<species>/` and `data/plastid/<species>/`) plus an
`analysis/` tree that each module reads/writes its own `results/`/`work/`
into. Point at yours with:

```bash
export PLANT_ORGANELLE_DATA_ROOT=/path/to/your/data/repo
```

Left unset, scripts fall back to treating *this* repo as the data root
too (useful only for a quick local smoke test with a handful of files
alongside the code - not how a real dataset should be organised).

## Modules

See each module's own README under `analysis/<module>/` for full detail
- this is deliberately terse:

| module | what it does |
|---|---|
| `qc_basic_stats` | per-assembly `pass`/`flag`/`fail` - every other module filters against this by default |
| `annotation` | oatk's own gene calls, reshaped into one tidy table |
| `denovo_annotation` | oatk-*independent* core gene calling (nhmmscan + tRNAscan-SE + barrnap) |
| `editing` | RNA-editing-aware ORF finding/correction (`orfedit`) |
| `trans_splicing` | reconstructs `nad1`/`nad2`/`nad5`/`rps3` from scattered exons (`transsplice`) |
| `gff_export` | merges the above into one standardised GFF3 per species |
| `orf_scan` | non-core ORFs vs Pfam - TE/mitovirus content (opt-in, heaviest module) |
| `repeats` | repeat families and putative recombination-mediating repeats |
| `synteny` | MTPT (plastid-to-mito transfer) detection, dotplots, pairwise comparison |
| `phylogeny` | partitioned ML tree per organelle from single-copy marker genes |
| `linearize` | `gfatk`-based assembly-graph resolution to a single linear sequence |
| `topology_plots` | supporting visualisation for `synteny`/`linearize` |

Every module works on a full dataset by default, and takes a
`--species-list` for a subset. None of the LSF/`bsub` examples in these
READMEs are required - they're this cluster's own convention for running
the heavier stages at scale; anything here runs fine as a plain foreground
command on a small dataset.
