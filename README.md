# poa - plant organelle annotator

A gold-standard annotation pipeline for plant
mitochondrial and plastid genomes: core gene calling
(direct `nhmmscan` against oatkDB's gene-family HMMs, plus tRNA/rRNA
calling), RNA-editing-aware ORF finding, spliced gene reconstruction
(trans-spliced `nad1`/`nad2`/`nad5`/`rps3` and cis-spliced `nad4`/`nad7`/
`ccmFc`/`cox2`/`rpl2`/`rps10`), unitig-level coordinates, and a merged GFF3
export at contig and unitig level.

It grew out of a 1250-species dataset repo (`plant_organellar_database`),
which now uses poa rather than the other way round.

## Install

```bash
mamba env create -f environment.yml
conda activate poa
./install.sh        # gfatk, orfedit, transsplice (Rust) via cargo
pip install .
poa check           # every tool found, oatkDB databases ready
```

`environment.yml` pins the third-party tools (hmmer, tRNAscan-SE, barrnap,
mafft, minimap2, samtools, seqkit, gfatools, GraphAligner); `install.sh`
pins this project's own Rust tools to the versions validated here: `gfatk`
0.4.0, `orfedit` 0.1.0 and `transsplice` 0.2.1. A bioconda package (and
with it Docker/Singularity images) is in progress.

## Reference data

All bundled - nothing to download, no oatk install needed:

- `src/poa/data/oatkdb/` - oatkDB's Viridiplantae mito and plastid
  gene-family HMM databases (`viridiplantae_mito_v20250217`,
  `viridiplantae_pltd_v20260928`), built with
  [OatkDB](https://github.com/c-zhou/OatkDB)'s `oatkdb` (MIT; see `NOTICE`
  there). Decompressed and `hmmpress`'d into `~/.cache/poa` (or
  `$POA_CACHE`) on first use. To use your own build instead, set
  `OATKDB_MITO_FAM` / `OATKDB_PLTD_FAM` to a pressed `.fam`.
- `src/poa/pipeline/trans_splicing/reference/` - the 11 angiosperm RefSeq
  mitogenomes (`genbank/`), their verified exon sets (`raw/`, with
  `junctions.tsv`) and the exon/whole-gene templates `transsplice` runs
  with (`profiles/`). `src/00a_fetch_references.py` regenerates `genbank/`
  and the cis-gene `raw/` files from NCBI; `src/00_build_exon_profiles.py`
  rebuilds templates from `raw/`. The trans-spliced genes' `whole_gene.pssm`
  were refined (`02_bootstrap_refine.py`) on the 1250-species dataset;
  rebuilding them from `raw/` alone resets that.
- `src/poa/pipeline/editing/reference/profiles/` - `orfedit` profiles built
  on the same dataset. `editing/src/01_scan_editing.py` uses them whenever
  the data root has none of its own.

Cluster-specific paths (tRNAscan-SE config, mafft's `MAFFT_BINARIES`) are
env-var overridable (`TRNASCAN_CONF`, `MAFFT_BINARIES`) and ignored when
they don't exist, so a conda install needs none of them.

## Running it

A single `poa run` command (FASTA and/or GFA in, GFF3 out) is being built.
Until then the stages run as scripts from the installed package
(`python3 -c 'import poa; print(poa.PIPELINE_DIR)'`) against a **data
root**: a directory containing `data/mito/<species>/` and
`data/plastid/<species>/` with each species' e.g. oatk-assembled
`*.ctg.fasta` + `.gfa`. Results and work files go under the data root's
`analysis/`. Point at it with `PLANT_ORGANELLE_DATA_ROOT`, or run from
inside it.

The annotation pipeline, in order (`tests/smoke/run_smoke.sh` has the exact
commands):

```
qc_basic_stats 01-05 -> denovo_annotation 00-04 -> editing 01 ->
trans_splicing 01 -> unitig_coords 00-01 -> qc_basic_stats 06 -> gff_export 01
```

## Tests

```bash
python3 -m unittest discover tests            # unit tests
tests/smoke/run_smoke.sh [DATA_ROOT] [THREADS] # end to end, installed poa
```

The smoke test builds a fresh data root from `tests/data/` (Arabidopsis
thaliana mito and plastid oatk assemblies, DToL ddAraThal4), runs the
pipeline above on it, and checks the mito GFF against RefSeq NC_037304.1
with `tests/smoke/check_arabidopsis.py`: gene content, exon structure and
editing-aware protein identity, with thresholds set just below the
dataset repo's own results so regressions fail. About 7 minutes on 8 cores.

## Modules

Each has its own README under `src/poa/pipeline/<module>/`:

| module | what it does |
|---|---|
| `qc_basic_stats` | per-assembly `pass`/`flag`/`fail` - the annotation stages filter against this by default |
| `denovo_annotation` | core gene calling (nhmmscan + tRNAscan-SE + barrnap) |
| `editing` | RNA-editing-aware ORF finding/correction (`orfedit`) |
| `trans_splicing` | reconstructs spliced genes (trans: `nad1`/`nad2`/`nad5`/`rps3`; cis: `nad4`/`nad7`/`ccmFc`/`cox2`/`rpl2`/`rps10`) from their exons (`transsplice`) |
| `linearize` | `gfatk`-based assembly-graph resolution to a linear sequence |
| `unitig_coords` | maps every GFA unitig onto the linearised contigs, so annotations can be given in unitig coordinates |
| `gff_export` | merges the above into one standardised GFF3 per species, contig- and unitig-level |

Every module works on a whole data root by default and takes a
`--species-list` for a subset. The LSF/`bsub` examples in these READMEs
are one cluster's convention for running at scale, not a requirement.
