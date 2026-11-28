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
including this project's own `gfatk` 0.7.0, `orfedit` 0.2.1 and
`transsplice` 0.3.3, built for linux-64, linux-aarch64 and macOS (Intel and
Apple silicon). To build those three from source instead, remove them from
`environment.yml` and run `./install.sh` (cargo); do that too while
bioconda does not have them (on 2026-11-27 it had gfatk 0.6.3, orfedit
0.1.0 and transsplice 0.2.1, with updates awaiting review). gfatk 0.7.0's
GitHub release also has a Linux binary. On bioconda poa is called
`plant-organelle-annotator`, because `poa` there is the POA v2 aligner. Its
recipe is in `recipe/`; once it is merged, poa (and with it
Docker/Singularity images) will install with
`mamba install -c conda-forge -c bioconda plant-organelle-annotator`, and
the command is still `poa`.

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
- `src/poa/pipeline/trans_splicing/reference/` - angiosperm RefSeq
  mitogenomes (`genbank/`), their verified exon sets (`raw/`, with
  `junctions.tsv`), spliced genes of DToL angiosperms with the editing
  their RNA-seq shows applied (`raw_rna/`), and the exon/whole-gene
  templates `transsplice` runs with (`profiles/`), built from the loci
  whose proteins are the edited ones. `src/00a_fetch_references.py`
  regenerates `genbank/` and the cis-gene `raw/` files from NCBI;
  `src/00_build_exon_profiles.py` rebuilds the templates. `plastid/` is
  the same for the plastid genes.
- `src/poa/pipeline/editing/reference/profiles_refseq/` - `orfedit`
  profiles for the single-exon mito genes, built from the RefSeq
  angiosperm mitogenomes annotated with edited proteins and from
  RNA-edited genes of DToL angiosperms
  (`reference/lineages/angiosperm/`, which also holds the plastid
  profiles built from such genes), and for the two pieces of split ccmFn
  (ccmFn1/ccmFn2), from 10 Brassicaceae and Allium cepa. These come
  first.
- `src/poa/pipeline/editing/reference/profiles/` - `orfedit` profiles built
  from the DToL dataset's own gene calls, for a gene the others lack.
- Moss and liverwort references, under `lineages/moss/` and
  `lineages/liverwort/` in both of the above: exon templates and editing
  profiles built from 20 moss and 27 liverwort RefSeq mitogenomes (cached
  in the same `genbank/`), and plastid ones from 26 moss and 30 liverwort
  RefSeq plastomes (`plastid/`, `pltd/`, `refseq_proteins.pltd.fasta`).
  See "Reference lineage" below.
- Fern, horsetail and clubmoss references, under `lineages/fern/`,
  `lineages/horsetail/` and `lineages/lycophyte/`: genes with the editing
  RNA-seq shows applied (`raw_rna/`, `rna_proteins.fasta`, each with a
  README) - of Darwin Tree of Life assemblies for fern (mostly) and
  horsetail, of a RefSeq genome and Darwin Tree of Life assemblies for
  clubmoss. Each has plastid references as well (`lineages/<x>/plastid/`,
  `lineages/<x>/pltd/`). Picked by `--lineage auto` or forced.

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
  genes, linearisation QC, and with `--cross-search` / `--cross-align` the
  other organelle's DNA (see "DNA from the other organelle" below)
- `logs/` - one log per stage; `work/` - intermediate files

`<name>` defaults to the first input's file name up to its first `.`
(`-n` to set it). About a minute for Arabidopsis mito + plastid on 2
cores, most of it nhmmscan and the spliced-gene reconstruction; Vitis'
773 kb mitogenome with its plastome takes under 2 minutes. tRNAscan-SE
runs behind its `--mid` filter, which takes seconds instead of minutes
and misses only very weak hits (25 bits or less); `--trnascan-full`
runs the full search (about 5 more minutes for Arabidopsis).

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
79 are one-off calls in a gene's best copy. poa needs orfedit 0.2.0 and
transsplice 0.3.0: older ones cannot read the profiles that record
editing sites, and poa stops with a message saying so.

Read frameshifts in spliced genes with care. In the DToL mitochondria
9% of reconstructed genes carry one (ccmFc 29%, rpl2 26%, nad5 16%, nad7
7%), most at a handful of positions shared by hundreds of assemblies
(ccmFc profile column 470 in 193, nad1 column 132 in 184, nad5 column 248 in 143)
- so a high `dtol_assemblies` there says the gene differs from the
template at that point in that lineage, not that the assembly is wrong.
In plastids fewer than 1% of reconstructed genes have one (rps16 3%,
clpP 1%).

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
to the start and stop codon.

Where a single-exon gene starts. The alignment starts a gene wherever
the profile's first columns match, which is not always a start codon. A
gene whose first aligned codon is not ATG (or ACG, which an edit makes
AUG) now starts at the nearest ATG within 10 codons: upstream, in an
open frame, or downstream at a column whose references can start with
Met; failing that, the nearest ACG. On 19 RefSeq angiosperm mitogenomes
and 13 RefSeq plastomes, starts within 3 bp of the annotation went from
384 to 401 of 457 (mito) and from 564 to 613 of 671 (plastid), with one
made wrong.

Where a single-exon gene ends. A reading frame ends at its first stop
codon, and orfedit no longer lets an alignment step over one to reach
the profile's last columns, as it did wherever a gene ended a few
residues before its profile. With no stop right after the aligned
protein, the first one within 50 codons ends the gene. A stop the
protein goes on matching beyond (a pseudogene, a stop the lineage edits
away, an error) is kept inside the gene and reported: `*` in the
protein and `n_internal_stops` in `tables/edited_gene_calls.tsv`
(`reconstructed_genes.tsv` for spliced genes), `n_internal_stops` on
the gene's mRNA in the GFF. And a gene is not read through an internal
stop or a frame change with fewer than 25 profile columns before or
after it: it starts after it, or ends before it in its own frame, since
one that buys only the profile's last few columns cannot be told from a
gene that ends differently there. Such a frame change is still reported
- a row of `tables/frameshifts.tsv` with `near_end` true, a `frameshift`
feature with `near_end=true` in the GFF - because an indel there is as
worth checking against the reads as any other (leaving them out, 83% of
planted indels were found in place of 95%). A gene that runs out of the window it is
looked for in (the hit plus 300 bp each side) is looked for again in a
longer one: rbcL, whose HMM is the gene's first half, used to stop at
codon 341.

On 19 RefSeq angiosperm mitogenomes and 12 RefSeq plastomes, with the
same profiles before and after (`build/rnaedit/ends`), single-exon
genes ending within 3 bp of the annotated stop went from 397 to 443 of
457 (mito) and from 558 to 621 of 626 (plastid), and the frame changes
called in those mito genes from 23 to 3; one end that was right is now
wrong (Liriodendron sdh4, an edit-made stop a codon past the profile's
end; with the references rebuilt since, its profile ends at that stop
and the end is right again). On 153 DToL
assemblies, single-exon genes ending at a stop went from 10,728 to
11,479 of 11,608 and the frame changes their genes are read through
from 234 to 184 (105 more, near a gene's end, are reported and not
followed); spliced genes complete went from 3,179 to 3,187 of 3,252;
81 single-exon and 25 spliced genes have an internal stop. Of 1-bp
indels planted in real genes as many are found as before (single-exon
94.9% of 1,500 against 94.7%; spliced mito 90.2% of 768, plastid 77.3%
of 894, against 90.2% and 77.0%). A
gene that ends before its profile's tail no longer counts as aligned
end to end - it used to by the step over its stop - so that count falls
where a profile's tail is longer than a species' gene (ccmC in every
moss; sdh4, rps19, rps14 and ccmC in angiosperms).

What is still wrong, mostly: genes whose exon count differs from the
template's (cox2 with one intron or two, an intron lost in a lineage),
and gene ends the references disagree on (ccmFc, rpl2). Liverworts and
hornworts are reconstructed on the angiosperm templates, as are ferns,
horsetails and clubmosses unless their lineage is given (below).

The recurrence table leaves out the 54 DToL assemblies resolved with
gfatk before 0.6.3, whose contigs repeat each graph overlap: they are
2.5% of the assemblies and had 8% of the single-exon calls.

### Circular contigs

A circular genome given as linear sequence starts at an arbitrary base,
and a gene through that base would be called as two partial genes, one
at each end. When `poa run` finds the same gene at both ends of a contig
(one not marked `circular=false` by oatk), it moves that contig's origin
into its longest gene-free stretch, annotates it again, and maps
everything back onto the sequence you gave it. It looks for such genes in
the nhmmscan and rRNA calls, before tRNAscan-SE runs (oatkDB has tRNA
models too), so only that organelle's nhmmscan is run again; a tRNA only
tRNAscan-SE calls across an origin moves it a second time. In the output the gene is one feature
with its end past the contig's length (GFF3's convention for features
across the origin of an `Is_circular=true` sequence). With a GFA it
leaves the contig as it is and warns, since the unitig coordinates
depend on the assembler's origin.

### DNA from the other organelle

Seed-plant mitogenomes carry plastid DNA, from a few hundred bases to tens
of kilobases, and a plastome now and then carries mitochondrial DNA.
Gene calling searches each organelle only with its own gene models, so
`poa run` finds none of this unless asked. Two options, both off by
default, look for it, and their findings go into the GFFs as
`biological_region` features and into `tables/organelle_transfers.tsv`.

`--cross-search` searches each organelle's contigs with the other
organelle's oatkDB models as well (one more nhmmscan). A hit counts only
when it beats every one of the organelle's own models on the same stretch
by at least 10 bits. A plastid psbA in a mitogenome is reported, but the
mitogenome's own rrn18 hit by the plastid rrn16 model, or its rpl2 by the
plastid rpl2 model, is not, because its own models fit it better. Each
feature found this way is the other organelle's gene, or a piece of it:
`plastid_derived=true` in a mitogenome and `mito_derived=true` in a
plastome, with the model's name, its share of the model
(`model_coverage`), the E-value and the best score the organelle's own
models reach there (`own_best_score`). It needs only one organelle.

`--cross-align` needs both organelles and BLAST+. It aligns the plastid
contigs to the mito contigs (`blastn -task dc-megablast`, E-value below
1e-10) and reports every alignment of at least 100 bp at 80% identity or
more, in both GFFs, as a stretch `shared_with` the other organelle, with
where it lies there (`other_region`, `orientation`, `identity`). Below 80%
identity the alignments are almost all between the two organelles' rRNA
genes, which are related from long before any transfer. An alignment does
not say which way the DNA moved, so these stretches are not marked
derived. A `--cross-search` gene inside one says which way that part went.
On five test samples the alignments held 132 of the 143 plastid genes
(150 bp or longer) that `--cross-search` found in the mitogenomes;
minimap2 held 81 of them, which is why this step uses BLAST.

| sample | plastid genes or pieces in the mitogenome | shared stretches (mito bp) | mito genes in the plastome |
| --- | --- | --- | --- |
| *Arabidopsis thaliana* | 15 | 13 (5,178) | none |
| *Daucus carota* | 39 | 21 (47,858) | a cox1 piece in both inverted-repeat copies |
| *Vitis vinifera* | 82 | 25 (60,549) | none |

The carrot cox1 piece is consistent with the mitochondrial DNA known to
sit in carrot's plastid inverted repeat. The searches took 15 to 31 seconds
of these runs' 6 to 19 minutes (2 cores). The genes both organelles
have (rpl2, rpl16, rps3, rps4, rps12, rps14, rps19) are checked whatever
the options: a call whose protein fits the other organelle's references at
least twice as well as its own has `plastid_derived=true` or
`mito_derived=true` on its gene.

### Reference lineage

poa's first references were angiosperm ones, and other land plants differ
in both things they encode: where the introns are, and how much the
transcripts are edited. On the angiosperm references only 12% of moss
spliced genes came out complete, and 98% of moss mito calls were
"corrected" towards the angiosperm protein. poa therefore picks a
reference lineage per species (`--lineage auto`, the default; any lineage
can be forced):

| lineage | for | mito references | plastid references |
|---|---|---|---|
| angiosperm | seed plants, and anything not recognised | 14 RefSeq mitogenomes annotated with edited proteins, and RNA-edited genes of DToL angiosperms | RNA-edited genes of DToL angiosperms |
| moss | mosses | 20 RefSeq mitogenomes | 26 RefSeq plastomes |
| liverwort | liverworts | 27 RefSeq mitogenomes | 30 RefSeq plastomes |
| fern | ferns | RNA-edited genes of 11 DToL ferns, on top of angiosperm | the same 11 ferns |
| horsetail | *Equisetum* | RNA-edited genes of 8 DToL horsetails, on top of angiosperm | the same horsetails |
| lycophyte | clubmosses | RNA-edited genes of two RefSeq and three DToL clubmosses, on top of angiosperm | two DToL clubmosses |

The editing stage aligns the species' mito genes to the angiosperm
profiles and to each other lineage's (plain C-to-U editing), and the
species takes the lineage whose profiles beat the angiosperm ones by the
most per gene, provided that reaches the lineage's margin over at least
4 genes (215 for moss, 275 for liverwort, 60 for fern, 250 for horsetail
and lycophyte; `common/lineage.py` has the calibration). Anything less
clear stays angiosperm. Moss and liverwort replace the angiosperm
references outright; fern, horsetail and lycophyte keep them and add
exon templates and editing profiles of their own, and are looked for
among the species still on the angiosperm references. The call and the
margins are written to `tables/lineage.tsv`
(`analysis/editing/results/lineage.tsv` stage by stage), and the editing
and trans_splicing stages use it. The lineage is read from the
mitochondrial genes: a plastid run on its own uses the angiosperm
references unless `--lineage` names one.

On fresh runs of 153 Darwin Tree of Life assemblies every call was
right: 64 of 64 liverworts, 30 of 30 ferns, 28 of 28 mosses, 9 of 9
horsetails and 4 of 4 clubmosses, with 12 angiosperms, 3 gymnosperms, a
hornwort, *Isoetes* and *Selaginella* left on the angiosperm references.
The ferns furthest from the fern references, *Ophioglossum* and
*Cryptogramma*, have margins of 87 and 140. On the dataset's stored gene
calls the highest of 716 angiosperms reaches 50 on the liverwort
profiles, 28 on the moss ones and 12 on the fern ones. Output for a species is the same whether its lineage
is detected or forced (15 species compared).

Angiosperm. An editing profile corrects a gene towards its reference
proteins, so those have to be the edited ones. The first angiosperm
references were RefSeq proteins, most of them the genome's own
translation, and poa found one annotated edit in twenty. They are now
the proteins of the 14 RefSeq mitogenomes that are annotated edited, and
genes of DToL angiosperms with the editing their own public RNA-seq
shows applied; the plastid ones are RNA-derived throughout. Those
proteins also say where editing happens, and each profile records the
columns edited in at least two genera. A C-to-U edit then costs 1.5 at
such a column and 10 anywhere else, where a difference from the profile
is far more often a residue of the species' own; without that, most
edits called in fast-evolving plastid genes (rpoC2, matK) were wrong.

The two costs suit a plant that edits where angiosperms do. Gymnosperms,
hornworts, quillworts and spikemosses use the same references and edit
heavily elsewhere, so each species is first corrected at one cost (2.5)
everywhere, and takes the two only if at least 65% of the C-to-U edits
then called in its mitochondrial genes are at known columns (10% of the
plastid ones, when there is no mitochondrion to judge by). Of 89
angiosperms 85 do, with a median of 95% there; four sedges, whose
mitochondrial genes have drifted furthest, have 58 to 62% and keep the
one cost. Ginkgo has 42%, a yew 38%, a hornwort 32%, a spikemoss 25%
and a quillwort 25%, and these too are annotated at the one cost. `tables/lineage.tsv` gives the counts and the outcome
(`edit_sites`).

With each species' genus left out of the references, against the RNA-seq
of 109 DToL angiosperms and the annotation of 19 RefSeq mitogenomes:

| | RefSeq edited proteins alone | with the RNA-derived genes | and editing sites |
|---|---|---|---|
| mito single-exon genes: C-to-U sites in the RNA found | 56% | 71% | 77% |
| ... edit calls that are in the RNA | 94% | 93% | 97% |
| mito spliced genes: sites found | 69% | 78% | 80% |
| ... edit calls that are in the RNA | 95% | 96% | 98% |
| plastid single-exon genes aligned end to end (of 4,561) | 3,426 | 4,044 | 4,044 |
| ... sites found, 32 species whose reads are the assembly's own | 14% | 70% | 75% |
| ... edit calls that are in the RNA | 8% (259 right) | 40% (1,474) | 93% (1,520) |
| plastid spliced genes: edit calls that are in the RNA | 12% (64 right) | 59% (729) | 94% (835) |
| RefSeq coding sequences with every exon placed (of 655) | 505 | 535 | 536 |
| RefSeq annotated edits found | 68% | 84% | 89% |
| ... edit calls that are annotated | 85% | 86% | 87% |

Before any of this, with references of unedited proteins and genes not
yet ended at their stops, poa placed 424 of the 655 RefSeq coding
sequences and found 6% of their annotated edits.
A silent edit changes no residue and cannot be called, and the RNA-seq is
mostly poly-A selected, so "found" is out of a list that holds
unprocessed transcripts' partial edits. It is also often another
plant's than the assembly, and where the two differ by a C/T the reads
show an edit that is not one: hence the plastid row for the 32 species
without such differences (over all 109, the same calls are 42% of the
"sites"). Checked site by site against the reads, 93% of the mito calls
with at least 10 reads are edited there and another 3% in under a fifth
of the reads; of the plastid ones, 90% and 4%.

atp6. Its annotated starts lie 18 to 146 codons before the conserved
part, which its profile begins at: a presequence too variable to align,
annotated from the first ATG of the open reading frame. With orfedit
0.2.1 poa starts angiosperm-reference atp6 there (`--upstream-start
atp6`, scan.UPSTREAM_START_GENES), unless the gene starts at an ACG an
edit makes into AUG. On 19 RefSeq angiosperm mitogenomes atp6 starts
within 3 bp of the annotation went from 8 to 14 of 19, single-exon mito
starts from 401 to 407 of 457, none made wrong; four of the five atp6
still wrong now start further upstream than RefSeq does. Over 1,785 atp6
copies in 836 DToL angiosperm mitogenomes 861 start further upstream
(median 87 codons); their alignment, edits, frame changes and ends are
as before.

Still open: starts RefSeq's own records disagree on. Of the single-exon
genes of ten RefSeq mitogenomes, 200 of 245 had both ends within 3 bp of
the annotation (164 before the start and stop rules above; the atp6 rule
adds 3 in the gene-end harness); most of the rest are rps4, mttB, rpl16
and matR, whose annotated starts differ from record to record. A
plastid-derived rpl2 in a mitogenome is reported as an incomplete mito
rpl2 where it used to pass for a complete one.

Moss. Moss mitochondria have their own intron set (nad1 has 3 exons, not
5; nad2 2; nad4 2) and edit little. On the 956-species dataset, all 174
mosses were called moss and no other species was, and moss spliced genes
went from 12% to 99% complete. Moss genes are aligned at an edit cost of
6, not 2.5 (`lineage.EDIT_PENALTY`): against RNA-seq of 16 DToL mosses,
poa at 2.5 called 579 edits in 626 plastid genes and 127 of them were in
the RNA; at 6, on the moss plastid references, it calls 31 and 16 are.
That gives up the real edits of the few mosses that do edit their
plastids (*Brachythecium rutabulum* and *Fissidens dubius*, about 80
sites each, of which 10 to 15 are still found), for proteins that are
otherwise the genome's own.

Liverwort. Liverwort mitochondrial genes are split differently again
(cob 4 exons, cox1 up to 10, cox2 and cox3 3, atp1 3 or 1; none
trans-spliced) and the complex thalloids do not edit at all. With each
genome's genus left out of the references, on 27 RefSeq liverwort
mitogenomes against their own annotation:

| | angiosperm references | liverwort references |
|---|---|---|
| annotated exons placed within 3 bp | 176 of 1,696 | 1,542 of 1,696 |
| intron-containing genes with every exon placed | 0 of 337 | 305 of 337 |
| single-exon genes placed | 172 of 697 | 579 of 697 |
| annotated genes with a frameshift call | 144 | 17 |
| annotated genes not called | 126 | 0 |
| spliced genes complete, 64 DToL liverworts | 372 of 635 | 1,072 of 1,082 |
| frameshift calls, same | 590 | 20 |
| mito edits called in 9 complex thalloid liverworts | 942 | 12 |

Most single-exon genes still off are ccmB, ccmFn and nad6, whose
annotated start differs between the RefSeq records themselves. nad7, a
pseudogene in nearly all liverworts, is not reconstructed, and rpl6 and
rps8, which liverworts keep, have no gene model and are not called. The
references are the RefSeq proteins, which are unedited, so in the leafy
and simple thalloid liverworts, which do edit, poa finds part of the
editing only: against RNA-seq of 14 DToL species, 226 edits called in
132 mito genes and 142 of them right, of 581 in the RNA (on the
angiosperm references, 594 called in the 122 genes placed and 124
right, of 406).

Plastids. Each lineage but angiosperm has plastid exon templates and
editing profiles of its own, used for the species that take the lineage.
Against the angiosperm-built profiles poa "edited" these plastids
towards the angiosperm protein - 713 edits called in 431 plastid genes
of 8 horsetails, whose RNA shows a handful.

| plastid single-exon genes, against RNA-seq | aligned end to end | edits called | of them in the RNA |
|---|---|---|---|
| horsetail, 8 species (each left out in turn) | 263 of 431 -> 439 of 445 | 713 -> 59 | 6 -> 9 |
| clubmoss, 2 species (each left out in turn) | 57 of 89 -> 76 of 79 | 401 -> 282 | 214 -> 186 |
| moss, 16 species | 405 of 626 -> 597 of 624 | 579 -> 31 | 127 -> 16 |
| leafy liverwort, 15 species | 230 of 358 -> 333 of 358 | 541 -> 286 | 198 -> 153 |
| complex thalloid liverwort, 4 species | 119 of 178 -> 167 of 178 | 112 -> 36 | 1 -> 1 |

The moss and liverwort rows are not held out (their references are
RefSeq plastomes of other individuals, often of the same genus). See
`trans_splicing/README.md`, `editing/README.md` and the READMEs under
`reference/lineages/`.

`--lineage fern` is for ferns, and, forced, for quillworts (*Isoetes*)
and hornworts, which are not detected. These lineages edit U to C as well as C to
U, and their organelle genes have genomic in-frame stop codons that
editing turns into Gln or Arg - about nine per mitochondrial gene in
Polypodiales. Without the lineage set, treat a fern's proteins and
frameshift calls as unreliable: each unread stop is a skipped column, an
apparent frame change or an early end. Horsetails and clubmosses have
lineages of their own (further down); *Selaginella* has none and belongs
on the default.

With it, orfedit and transsplice read such stops through at the cost of a
C-to-U edit, and mitochondrial genes are corrected against fern
references: exon templates for the spliced genes (nad1, nad4, nad5 and
nad7 have four exons in ferns, rps3 three, rps10 one, and there is no
ccmFc; atp1, cox1 and rps14, which have no intron in angiosperms, have
one in ferns and are reconstructed from templates too) and editing
profiles for 22 single-exon genes. On those fern references other U-to-C
edits are called too, at a cost of 5; plastid genes, which have fern
references of their own, and the mitochondrial genes still on angiosperm
profiles (matR and nad9 among them) get stop read-through only. It needs
builds of orfedit
and transsplice with `--stop-edit-penalty` and `--reverse-edit-penalty`.
U-to-C edits are marked `U>C` in the `edit` column of the edit tables and
`edited_from=T;edited_to=C` in the GFF.

The fern references are built from edited proteins only. Few ferns have a
deposited mitochondrial editome (*Haplopteris ensiformis* and four
Salviniales; Kwok van der Giezen et al. 2025), so most come from eleven
Darwin Tree of Life ferns - nine Polypodiales in six genera, *Azolla
filiculoides* - whose own RNA-seq was mapped to their assemblies and the
editing it shows applied to their genes (`lineages/fern/raw_rna/` and
`editing/reference/lineages/fern/rna_proteins.fasta`; each has a README
with the species, runs and method). The six GenBank ferns whose proteins
only have their stop codons edited out are not used for profiles: with
them in the exon templates, transsplice found 9% of the C-to-U edits the
RNA shows.

Scored against that RNA, with each genus held out of the references in
turn (found = share of the RNA's edits called, right = share of calls
the RNA confirms; an RNA edit is one in at least 20% of reads):

| | C-to-U found / right | U-to-C found / right |
|---|---|---|
| single-exon genes, angiosperm profiles | 46% / 70% | 45% / 89% |
| single-exon genes, fern profiles | 80% / 96% | 89% / 97% |
| spliced genes, fern templates | 82% / 97% | 92% / 96% |

On the 18 DToL ferns that supplied no RNA, spliced mitochondrial genes
aligned over 99% of their template go from 74 to 137 of 159 against the
previous fern templates, and genes with a frameshift call from 29 to 7;
they get as many edit calls per gene as the eleven that did. Planted 1-bp
indels are found in 82% of spliced genes and in all 60 planted single-exon
ones. Exon junctions were checked against the RNA too: of 144 with
intron-spanning reads, 135 give the transcript the reads show, with
transsplice allowed to slide a junction a few bases (`--junction-slide`,
which the fern lineage sets where transsplice has it; 126 without). The
other nine are off by a base or two - nad2's first junction in eight
species, where the alternatives give nearly the same protein and neither
has a recognisable intron end. Two of 89 spliced genes run past a stop
codon that is never edited. Nine of the 29 ferns had their assemblies
replaced after their genes were first called; all figures here are from
calls on the current ones.

atp1, cox1 and rps14 were single-exon calls across an intron until they
got exon templates, from gene models built on the same RNA (the introns
are the junctions intron-spanning reads show). Held out by genus on the
eleven ferns, atp1 is complete in 11 of 11 (none of 19 single-exon calls
aligned end to end, and 6 drew a frameshift call), cox1 in 9 of 11 (0 of
12), and rps14 no longer draws the frameshift call it drew in all 12
copies; of the RNA edits in the three genes 84%, 77% and 91% are found,
with 95% or more of the calls right, against 49%, 39% and 34%.

Fern plastid genes are corrected against fern plastid references as
well, from the same eleven ferns: 613 edited proteins of 56 single-exon
genes and 79 loci of the 9 spliced ones
(`trans_splicing/reference/lineages/fern/plastid/`). Held out by genus,
against the land-plant plastid references every other lineage uses (as
first built, before rbcL and matK were read whole):

| | land-plant references | fern plastid references |
|---|---|---|
| single-exon genes aligned end to end | 339 of 614 | 596 of 614 |
| RNA edits found / calls right, single-exon | 74% / 59% | 82% / 83% |
| RNA edits found / calls right, spliced | 65% / 50% | 81% / 83% |
| genes running past a stop codon the RNA never edits | 81 | 36 |

Measured again the same way with poa 0.3.2's tools (orfedit 0.2.1,
transsplice 0.3.3) and the psbF profiles trimmed to their Met
(`build/rnaedit/pilot/remeasure033nb.sh`): single-exon genes aligned to
99% or more 589 of 614, RNA edits found / calls right 81% / 81%
(single-exon) and 77% / 83% (spliced), and 49 genes running past a stop
codon the RNA never edits. Genes ending at their first stop more often
read through one by the U-to-C edit ferns make at stops (ndhG, atpF,
matK, rpl32, ndhK, 2 to 28 codons). Held out the same way, horsetail
plastid genes align to 99% or more in 455 of 500 (9 species; the RNA
shows almost no plastid editing) and clubmoss ones in 80 of 87 (2
species, RNA edits found / calls right 79% / 60%).

Plastid calls are less often right than mitochondrial ones: about one in
seven is a C that other fern genera edit, or have as T, and this species
leaves alone. Other U-to-C edits are not called in plastid genes: they
found 91% of the plastid U-to-C sites against 77%, but only 42% of those
calls were in the RNA.

`--lineage horsetail` is for *Equisetum*, the one living genus of
horsetails. No horsetail mitogenome is annotated in GenBank, so its exon
templates (11 genes) and editing profiles (23) come from gene models
built on RNA-seq of eight *Equisetum* species in the Darwin Tree of Life
dataset: the introns are the junctions intron-spanning reads show, and
the proteins carry the C-to-U editing the reads show. Horsetails do not
edit U to C, so the editing model is the default one. *Equisetum* genes
are split differently from angiosperm ones - cox1 has three introns,
nad7 four exons, rps3 three, rps10 no intron, and nad1 five exons that
are not the angiosperm five. Poly-A RNA-seq holds few mitochondrial
reads, and two species (*E. sylvaticum*, *E. telmateia*) supply most of
the references; a profile is built from as little as one species
(`trans_splicing/reference/lineages/horsetail/raw_rna/README.md`).

With each species' own loci left out in turn, on nine *Equisetum*
assemblies (the default-reference figures are for the eight that have
them):

| | default references | horsetail references |
|---|---|---|
| single-exon genes aligned end to end | 119 of 218 | 172 of 189 |
| spliced genes complete | 64 of 80 | 86 of 86 |
| genes with a frameshift call | 45 | 4 |
| RNA edits found / calls right, single-exon | 44% / 50% | 80% / 96% |
| RNA edits found / calls right, spliced | 47% / 58% | 79% / 99% |

All 110 exon junctions with intron-spanning reads agree with the reads.
The four frameshift calls left are in matR, whose *Equisetum* reading
frame is short and differs between species. ccmFn aligns over 89% of its
profile: its HMM hit is in two pieces and the correction window stops
short of the start.

`--lineage lycophyte` is for clubmosses (Lycopodiaceae). Their genes
have an intron set of their own: atp6, atp9, cob, cox1, cox3, nad3 and
rps14 have introns that angiosperms lack, nad2 has three exons, nad4 two,
cox2 four, rps3 three, and the ccm genes and nad7 are gone. A few
in-frame stop codons are edited away, so the lineage reads stops through
as the fern one does. Four clubmoss mitogenomes are in RefSeq, but with
unedited proteins: templates built from them placed exons and found a
fifth of the edits RNA-seq shows. The references are therefore the
RefSeq annotations of *Phylloglossum drummondii* and *Lycopodium
japonicum* with the editing their own RNA-seq shows applied (all 35
genes, 289 C-to-U and 5 U-to-C edits; all 34 genes, 266 C-to-U edits),
plus the genes that poly-A RNA-seq covers in three Darwin Tree
of Life clubmosses
(`trans_splicing/reference/lineages/lycophyte/README.md`).

| | default references | lycophyte references |
|---|---|---|
| annotated exons placed within 3 bp, 3 RefSeq genomes, genus held out | 51 of 198 | 176 of 198 |
| spliced genes with every exon placed | 4 of 44 | 41 of 44 |
| genes with a frameshift call | 14 | 3 |
| spliced genes complete, 3 DToL clubmosses, own loci held out | 15 of 27 | 45 of 45 |
| RNA edits found / calls right, same | 65 of 119 / 65 of 147 | 88 of 98 / 88 of 100 |

All 21 exon junctions that intron-spanning reads cross in the three
dataset clubmosses are as the reads show. The edit figures rest on the
18 genes with enough reads, and the references on two deeply sequenced
species, so treat this lineage as less tested than the fern one. atp9's
third exon, 7 bp, is rarely placed; extra copies of atp4, atp8 and rps2
draw frameshift calls that may well be real; and nad6 draws one near its
end in most clubmosses, where the species' proteins differ in length,
which probably is not. It is not for
*Selaginella* or *Isoetes*, whose mitogenomes are unlike these
(*Isoetes* edits like a fern).

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

`denovo_annotation 05` (`05_organelle_transfers.py`, what `poa run
--cross-search` and `--cross-align` run) is optional, after
`denovo_annotation 00` and before `gff_export 01`.

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
