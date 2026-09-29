# plant_organelle_annotator

A gold-standard annotation pipeline for plant
mitochondrial and plastid genomes: core gene calling
(direct `nhmmscan` against oatkDB's own gene-family HMMs, plus tRNA/rRNA
calling), RNA-editing-aware ORF finding, spliced gene reconstruction
(trans-spliced `nad1`/`nad2`/`nad5`/`rps3` and cis-spliced `nad4`/`nad7`/
`ccmFc`/`cox2`/`rpl2`/`rps10`), unitig-level coordinates, and a merged GFF3
export - plus a wider
suite of QC, repeat, synteny, phylogeny and non-core ORF (TE/mitovirus)
analyses.

This repo is code and reference data only - no per-dataset results. It was
split out of a specific 1250-species dataset repo (`plant_organellar_database`)
so it can be installed and run against *any* set of oatk-assembled organelle
genomes, not just that one. `dev/SYNCED_FROM` records the dataset-repo
commit the code was last synced from (see "Re-syncing" below).

## Install

```bash
mamba env create -f environment.yml
conda activate plant_organelle_annotator
./install.sh   # cargo-installs this project's own Rust tools
```

See `environment.yml` for exactly what's pinned (minimap2, samtools,
seqkit, blast, cd-hit, mafft, iqtree, gfatools, hmmer, tRNAscan-SE,
barrnap, Bandage, GraphAligner, AMAS, rust) and `install.sh` for this
project's own Rust tools, pinned to the versions validated here: `gfatk`
0.4.0, `orfedit` 0.1.0 and `transsplice` 0.2.1.

## Reference data

Bundled:

- `analysis/trans_splicing/reference/` - the 11 angiosperm RefSeq
  mitogenomes (`genbank/`), their verified exon sets (`raw/`, with
  `junctions.tsv`) and the exon/whole-gene templates `transsplice` runs
  with (`profiles/`). `src/00a_fetch_references.py` regenerates `genbank/`
  and the cis-gene `raw/` files from NCBI; `src/00_build_exon_profiles.py`
  rebuilds templates from `raw/`. The trans-spliced genes' `whole_gene.pssm`
  were refined (`02_bootstrap_refine.py`) on the 1250-species dataset;
  rebuilding them from `raw/` alone resets that.
- `analysis/editing/reference/profiles/` - `orfedit` profiles built on the
  same dataset. `editing/src/01_scan_editing.py` uses them whenever the data
  root has none of its own; a dataset with 5+ species per gene can build
  its own with `00_build_reference_profiles.py`.

Not bundled - oatkDB's gene-family `.fam` databases (hmmpress'd,
`denovo_annotation/src/00_run_nhmmscan.py`). Point at them with:

```bash
export OATKDB_MITO_FAM=/path/to/viridiplantae_mito.fam
export OATKDB_PLTD_FAM=/path/to/viridiplantae_pltd.fam
```

Either use the prebuilt databases from
[OatkDB](https://github.com/c-zhou/OatkDB) (not the versions validated
here), or build current ones with its `oatkdb` script, as this pipeline's
were (`viridiplantae_mito_v20250217`, `viridiplantae_pltd_v20260928`;
taxid 33090 = Viridiplantae). The plastid one, about 3 hours on a cluster
node:

```bash
oatkdb -j 5 -t 10 -c 11 -T ./TEMP_PLTD -o viridiplantae_pltd 33090 chloroplast
```

The mito database is the same call with `mitochondrion` (check `oatkdb`'s
usage for the genetic-code option `-c`; 11 is the plastid code).

`oatkdb` downloads through NCBI EDirect; the plastid build here needed a
patched `nquire` on PATH. Other cluster-only paths (tRNAscan-SE config,
Pfam for `orf_scan`, mafft's `MAFFT_BINARIES`) are env-var overridable
(`TRNASCAN_CONF`, `PFAM_DB`, `MAFFT_BINARIES`) and ignored when they don't
exist, so a conda install needs none of them.

## Pointing this at your data

Every script expects a **data root**: a directory containing `data/`
(your e.g. oatk-assembled `<species>/*.ctg.fasta`+`.gfa` per species, laid out
as `data/mito/<species>/` and `data/plastid/<species>/`) plus an
`analysis/` tree that each module reads/writes its own `results/`/`work/`
into. Point at yours with:

```bash
export PLANT_ORGANELLE_DATA_ROOT=/path/to/your/data/repo
```

Results and work files go under the data root; shared modules and
bundled reference data are always read from this checkout, so one install
can serve several data roots. Left unset, scripts treat *this* repo as the
data root too (the original same-tree layout).

The annotation pipeline, in order (see `tests/smoke/run_smoke.sh` for the
exact commands):

```
qc_basic_stats 01-05 -> denovo_annotation 00-04 -> editing 01 ->
trans_splicing 01 -> unitig_coords 00-01 -> qc_basic_stats 06 -> gff_export 01
```

## Smoke test

```bash
tests/smoke/run_smoke.sh [DATA_ROOT] [THREADS]
```

Builds a fresh data root from `tests/data/` (the dataset's Arabidopsis
thaliana mito and plastid oatk assemblies, DToL ddAraThal4), runs the
pipeline above on it, and checks the mito GFF against RefSeq NC_037304.1
with `tests/smoke/check_arabidopsis.py`: gene content, exon structure and
editing-aware protein identity, with thresholds set just below the
dataset repo's own results so regressions fail. About 7 minutes on 8 cores.

## Modules

See each module's own README under `analysis/<module>/` for full detail
- this is deliberately terse:

| module | what it does |
|---|---|
| `qc_basic_stats` | per-assembly `pass`/`flag`/`fail` - every other module filters against this by default |
| `annotation` | oatk's own gene calls, reshaped into one tidy table |
| `denovo_annotation` | oatk-*independent* core gene calling (nhmmscan + tRNAscan-SE + barrnap) |
| `editing` | RNA-editing-aware ORF finding/correction (`orfedit`) |
| `trans_splicing` | reconstructs spliced genes (trans: `nad1`/`nad2`/`nad5`/`rps3`; cis: `nad4`/`nad7`/`ccmFc`/`cox2`/`rpl2`/`rps10`) from their exons (`transsplice`) |
| `unitig_coords` | maps every GFA unitig onto the linearised contigs, so annotations can be given in unitig coordinates |
| `gff_export` | merges the above into one standardised GFF3 per species, contig- and unitig-level |
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

## Re-syncing from the dataset repo

The code under `analysis/` is generated from the dataset repo, not edited
here:

```bash
python3 dev/sync_from_dataset.py ../plant_organellar_database [<commit>]
```

It exports a committed revision (`git archive`), drops per-dataset
material (results, work, the dataset's own `analysis/README.md` and
`run_all.sh`), applies the code/data-root split and the path overrides
above, and records the source commit in `dev/SYNCED_FROM`. Fix code in the
dataset repo and re-sync; an edit made here is overwritten on the next
sync. Then re-run the smoke test.
