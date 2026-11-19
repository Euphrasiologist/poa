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
including this project's own `gfatk` 0.6.2, `orfedit` 0.1.0 and
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

Every protein-coding gene poa corrects is aligned to its reference
profile allowing a change of reading frame inside it, and each change is
reported: an assembly indel (e.g. a miscounted homopolymer) or a real
pseudogene, which poa can't tell apart - check it against the reads.
Single-exon genes go through orfedit (>= 0.2.0); spliced genes - mito
nad1/2/4/5/7, ccmFc, cox2, rpl2, rps3, rps10 and plastid atpF, clpP,
ndhA, ndhB, rpl2, rpoC1, rps12, rps16, ycf3 - through transsplice
(>= 0.3.0), which aligns the joined exons and names the exon.

Each frameshift is a row of `tables/frameshifts.tsv` (single-exon) or
`tables/reconstructed_frameshifts.tsv` (spliced, with the exon `slot`):
the skipped bases' 0-based `start`/`end`; `n_bases` (1 = one base too
many for the frame, e.g. a 1-bp insertion; 2 = one too few, e.g. a 1-bp
deletion); `window_start`/`window_end`, the stretch the change could
equally well sit in (a homopolymer run, say); and `dtol_assemblies`, how
many DToL assemblies have one at the same place in that gene. It is also
a `frameshift` feature under the gene's CDS in the GFF (`window=` when
wider than the skipped bases), and counted in the QC table
(`n_frameshifts`, `n_new_frameshifts` - seen in at most 2 DToL
assemblies - and `frameshifted_genes`; not part of the status). `poa run`
prints them:

```
[poa] frameshifts mito: 2, 0 not seen in DToL - rps14 (in 12 DToL), rps19 (in 4 DToL) -> tables/frameshifts.tsv
```

(Arabidopsis's mitochondrial rps14 and rps19 are pseudogenes.) A
frameshift seen in many DToL assemblies is a feature of the gene; one
seen nowhere else in a gene that otherwise aligns well is the first thing
to check.

How well it works, on 1-bp indels planted at random in the coding
sequence of real DToL genes (scripts in `build/fs_val`, `build/ts_val`):
single-exon genes 92% found (96% plastid, 88% mito; 94% in genes' best,
near-complete copies), 99% with the right type, the indel at the call or
inside its window. Spliced genes: mito 84% found (652 of 778; rpl2, rps10
and cox2 are the weak ones), plastid 92% when the indel is in the copy
reconstructed (74% overall, as the second inverted-repeat copy of ndhB,
rpl2 and rps12 isn't reconstructed), 99-100% with the right type and all
in the right exon (the spliced-gene figures were measured before the
exon-placement changes below and not repeated since). The misses are mostly short genes (rps19, atp9, rps14)
and indels near a gene end, with too little gene on one side to tell.
Across the DToL assemblies there's a mean of 2 single-exon calls per
assembly, nearly all in partial copies or recurring in other assemblies;
79 are one-off calls in a gene's best copy. With older
orfedit/transsplice there are no frameshift calls (a warning says so).

Read frameshifts in spliced genes with care. In the DToL mitochondria
12% of reconstructed genes carry one (ccmFc 49%, rpl2 27%, nad5 16%, nad7
6%), most at a handful of positions shared by hundreds of assemblies
(ccmFc profile column 470 in 193, nad1 column 132 in 184, nad5 column 248 in 143)
- so a high `dtol_assemblies` there says the gene differs from the
template at that point in that lineage, not that the assembly is wrong.
In plastids about 1% of reconstructed genes have one (rps16 3%,
clpP 3%).

### Exon placement

Where transsplice puts exons was checked on the 24 RefSeq genomes the
templates come from (11 mito, 13 plastid), each run with templates built
without it (`build/exon_val`, `build/exon_impl`): 89% of 477 exons have
both ends within 3 bp of GenBank's (mito 86%, plastid 92%) and 68% are
base-exact; exons of 30 bp or less 15 of 16; last exons end exactly at
the stop codon in 142 of 169 and first exons start exactly at the start
codon in 166 of 178. Before these changes the same test gave 65% (no
exon of 30 bp or less right). What changed: the exon-slot profiles
(`00_build_exon_profiles.py`, see CHANGELOG), and in transsplice a
whole-contig and a next-to-the-neighbours search for every short exon,
decided by the whole-gene alignment and intron ends, and gene ends held
to the start and stop codon. Single-exon genes corrected by orfedit now
end at the stop codon too when the profile stops a few residues short
of it: on the same genomes 89% of plastid and 80% of mito genes end
exactly where GenBank's CDS does (17% and 67% before).

What is still wrong, mostly: genes whose exon count differs from the
template's (cox2 with one intron or two, an intron lost in a lineage),
and gene ends the references disagree on (ccmFc, rpl2). Liverworts,
hornworts, ferns and lycophytes are reconstructed on the angiosperm
templates.

The recurrence table leaves out the 54 DToL assemblies resolved with
gfatk before 0.6.3, whose contigs repeat each graph overlap: they are
2.5% of the assemblies and had 8% of the single-exon calls.

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
