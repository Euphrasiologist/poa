# poa - plant organelle annotator

A gold-standard annotation pipeline for plant
mitochondrial and plastid genomes: core gene calling
(direct `nhmmscan` against oatkDB's gene-family HMMs, plus tRNA/rRNA
calling), RNA-editing-aware ORF finding, spliced gene reconstruction
(trans-spliced `nad1`/`nad2`/`nad5`/`rps3` and cis-spliced `nad4`/`nad7`/
`ccmFc`/`cox2`/`rpl2`/`rps10`), unitig-level coordinates, and a merged GFF3
export at contig and unitig level.

It grew out of a need to annotate over 300 mitogenomes for a recent
<a href="https://www.biorxiv.org/content/10.64898/2026.10.26.761860v1.abstract">paper</a>.

## Install

```bash
mamba env create -f environment.yml
conda activate poa
pip install .
poa check           # every tool found, oatkDB databases ready
```

`environment.yml` pins every external tool from bioconda/conda-forge,
including this project's own `gfatk` 0.6.1, `orfedit` 0.1.0 and
`transsplice` 0.2.1, built for linux-64, linux-aarch64 and macOS (Intel and
Apple silicon). To build those three from source instead, remove them from
`environment.yml` and run `./install.sh` (cargo). A bioconda recipe for
poa itself is in `recipe/` and awaiting submission; once merged, poa (and
with it Docker/Singularity images) will install with
`mamba install -c conda-forge -c bioconda poa`.

## Reference data

All bundled so nothing to download, no oatk install needed:

- `src/poa/data/oatkdb/` - oatkDB's Viridiplantae mito and plastid
  gene-family HMM databases (`viridiplantae_mito_v20250217`,
  `viridiplantae_pltd_v20260928`), built with
  [OatkDB](https://github.com/c-zhou/OatkDB)'s `oatkdb` (MIT; see `NOTICE`
  there). Decompressed and `hmmpress`'d into `~/.cache/poa` (or
  `$POA_CACHE`) on first use. To use your own build instead, set
  `OATKDB_MITO_FAM` / `OATKDB_PLTD_FAM` to a pressed `.fam`.
- `src/poa/data/rrna/plant.hmm` - the rRNA models (barrnap's bacterial
  16S/23S/5S plus Rfam's RF00001 5S; see `NOTICE` there), run with
  `nhmmer` by `denovo_annotation/rrna_nhmmer.py`, a port of the
  `barrnap --kingdom plant` set-up this pipeline was validated with.
  barrnap itself isn't needed.
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
- `src/poa/pipeline/editing/reference/profiles_refseq/` - `orfedit`
  profiles for the other single-exon mito genes, built from the 11
  angiosperm RefSeq mitogenomes, and for the two pieces of split ccmFn
  (ccmFn1/ccmFn2), from 10 Brassicaceae and Allium cepa; used for any
  gene the dataset profiles don't cover.
- Moss references, under `lineages/moss/` in both of the above: exon
  templates and editing profiles built from 20 moss RefSeq mitogenomes
  (cached in the same `genbank/`). See "Reference lineage" below.

Cluster-specific paths (tRNAscan-SE config, mafft's `MAFFT_BINARIES`) are
env-var overridable (`TRNASCAN_CONF`, `MAFFT_BINARIES`) and ignored when
they don't exist, so a conda install needs none of them.

## Running it

```bash
poa run --mito sample.mito.ctg.fasta --mito-gfa sample.mito.gfa \
        --pltd sample.pltd.ctg.fasta --pltd-gfa sample.pltd.gfa \
        -o out/ -t 8
```

Give the contig FASTA for one or both organelles, from oatk or any other
assembler (`.gz` is fine). The assembly graph is optional: with it you
also get unitig-level coordinates (`unitig_loc=`/`unitig_span=`
attributes and a `.unitig.gff`). Without it, the contig GFF is the same
apart from those attributes. A graph alone isn't enough, because poa
annotates contigs and doesn't choose a path through the graph.

Output in `out/`:

- `<name>.<mito|pltd>.gff` - contig-level GFF3
- `<name>.<mito|pltd>.unitig.gff` - unitig-level GFF3 (with a GFA)
- `<name>.qc.tsv` - `pass`/`flag`/`fail` per organelle, with the reasons
  and the measurements behind them (see "QC" below)
- `tables/` - gene calls, predicted RNA edits, frameshifts, reconstructed spliced
  genes, linearisation QC
- `logs/` - one log per stage; `work/` - intermediate files

`<name>` defaults to the first input's file name up to its first `.`
(`-n` to set it). About 5 minutes for Arabidopsis mito + plastid on 4
cores, most of it tRNAscan-SE.

### QC

`poa run` ends with a `pass`/`flag`/`fail` line per organelle, also
written to `<name>.qc.tsv`:

```
[poa] QC mito: pass -> Arabidopsis_thaliana.qc.tsv
[poa] QC pltd: flag (core_gene_pct<0.80) -> Arabidopsis_thaliana.qc.tsv
```

and without a GFA or the other organelle, e.g.
`[poa] QC mito: pass; not assessed: graph, cross_organelle`.

The rules, thresholds and reason strings are those of the dataset QC gate
(`qc_basic_stats`, the gate behind the DToL dataset's pass/flag/fail
counts): dead-end nodes (flag >0, fail >2), fragmented graphs (flag >=5,
fail >=10 disconnected pieces, not counting closed circles at genome
depth), no circular contig (flag), a contig file of unjoined graph
segments (flag), core genes present (flag <80%, fail <50%), genome size
more than 4 MADs from its genus (flag), and, with both organelles, the
share of one assembly aligning to the other (`minimap2 -x asm5`; flag
>15%, fail >35%). The parts that need a dataset - which genes are core
(present in 90% of DToL assemblies) and each genus's genome sizes - ship
with poa as reference data (`qc_basic_stats/reference/`), so a single
sample is judged the way the dataset judged its own assemblies.

The genus comes from the first word of `--name`
(`Arabidopsis_thaliana` -> Arabidopsis); `--genus` sets it. A genus with
fewer than 3 DToL assemblies is compared with the whole dataset, as the
dataset gate does. The graph checks need the GFA (needing gfatk), and
the cross-organelle check both organelles (needing minimap2); whatever
can't be checked is listed in the `not_assessed` column rather than
counted against the assembly.

With a mito GFA, `poa run` also runs the low-depth side-path check
(`qc_basic_stats 06`, in `tables/low_depth_paths.tsv`). A flagged graph
gets a "linearisation uncertain" note in the contig GFF.

### Frameshifts

The editing correction (orfedit >= 0.2.0) aligns each single-exon
protein-coding gene to its reference profile allowing a change of reading
frame inside the gene, and reports where it happens: an assembly indel
(e.g. a miscounted homopolymer) or a real pseudogene, which poa can't
tell apart - check it against the reads. Each one is a row of
`tables/frameshifts.tsv` (the skipped bases' 0-based `start`/`end`,
`n_bases`: 1 = one base too many for the frame, e.g. a 1-bp insertion; 2 =
one too few, e.g. a 1-bp deletion), a `frameshift` feature under the
gene's CDS in the GFF, and a count in the QC table (`n_frameshifts`,
`frameshifted_genes`; not part of the status). `poa run` prints them:

```
[poa] frameshifts mito: 2 in atp1, cox1 -> tables/frameshifts.tsv
```

Multi-exon (trans- and cis-spliced) genes aren't checked. With orfedit
0.1.0 there are no frameshift calls.

### Circular contigs

A circular genome given as linear sequence starts at an arbitrary base,
and a gene through that base would be called as two partial genes, one
at each end. When `poa run` finds the same gene at both ends of a contig
(one not marked `circular=false` by oatk), it moves that contig's origin
into its longest gene-free stretch, annotates it again, and maps
everything back onto the sequence you gave it: the gene is one feature
with its end past the contig's length (GFF3's convention for features
across the origin of an `Is_circular=true` sequence). With a GFA it
leaves the contig as it is and warns, since the unitig coordinates
depend on the assembler's origin.

### Reference lineage

Moss mitochondria have their own intron set (nad1 has 3 exons, not 5;
nad2 2; nad4 2) and barely edit, so on the angiosperm references only 12%
of moss spliced genes came out complete and 98% of moss mito calls were
"corrected" towards the angiosperm protein. poa therefore picks a
reference lineage per species (`--lineage auto`, the default): the editing
stage aligns five single-exon mito genes (atp4 ccmB ccmC mttB rps14)
against both the angiosperm and the moss profiles, and a species goes to
moss only when at least 4 of them are scored and moss wins by 50 or more
per gene. The call is written to `tables/lineage.tsv`
(`analysis/editing/results/lineage.tsv` stage by stage), and the editing
and trans_splicing stages use it. `--lineage angiosperm|moss` forces a
lineage.

On the 956-species dataset, all 174 mosses were called moss and no other
species was. Moss spliced genes went from 12% to 99% complete, and the
share of moss calls changed by editing fell from 98% to 44%. Output for
every other species is unchanged. Liverworts, hornworts, ferns and
lycophytes still use the angiosperm references. See
`trans_splicing/README.md` and `editing/README.md`.

### Stage by stage, on a whole dataset

The stages are scripts in the installed package
(`python3 -c 'import poa; print(poa.PIPELINE_DIR)'`) that run against a
**data root**: a directory containing `data/mito/<species>/` and
`data/plastid/<species>/` with each species' `*.ctg.fasta` (and,
optionally, `.gfa`). Results and work files go under the data root's
`analysis/`. Point at it with `PLANT_ORGANELLE_DATA_ROOT`, or run from
inside it. In order (`tests/smoke/run_smoke.sh` has the exact commands):

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
dataset repo's own results so regressions fail. It then runs `poa run` on
the same inputs and requires identical GFFs. About 12 minutes on 8 cores.

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
