#import "@local/report:0.1.0": *
// figures quoted in the prose, from the tables in data/ (set when the tables are regenerated)
#let n-species = [111]       // Tree of Life angiosperms in the references
#let n-rna = [109]           // ...with RNA-seq, scored held out
#let pct-now = [82%]         // CDS with every exon placed, ten RefSeq genomes, now
#let single-now = [200]      // single-exon CDS placed, now
#let ss-now = [295]          // CDS placed allowing start/stop choice, now
#let found-gb = [88%]        // annotated edits found, five RefSeq genomes, now
#let right-gb = [89%]
#let found-rna = [80%]       // RNA edits found, five species, now
#let right-rna = [98%]
#let right-pltd = [93%]      // plastid single-exon edit calls in the RNA, now
#show: report.with(
  title: "poa against PMGA",
  subtitle: "A head-to-head on angiosperm mitogenomes, November 2026",
  abstract: [
    PMGA (Li et al. 2025, _Plant Communications_ 6:101191) annotates angiosperm mitogenomes
    with MAKER against a curated reference set and predicts RNA editing with Deepred-mt. I ran
    its published container and poa on the same 22 angiosperm mitogenomes and scored both the
    same way. As poa stood on 23 November, PMGA was ahead on everything measured, and above all
    on editing, where poa found one edit in twenty. poa's angiosperm references have since been
    rebuilt from RNA-edited proteins, and its profiles now record where those proteins are
    edited. With each test genome's genus kept out of the references, poa now finds nearly as
    many edits as PMGA with Deepred-mt (88% against 93% of those annotated in RefSeq, 80%
    against 84% of those in RNA-seq) and, against the reads, calls about a quarter as many
    wrong ones. Its genes now start at a start codon and end at their stop codon, and on gene
    models it is level with PMGA (every exon placed in 82% of coding sequences against PMGA's
    81%), though PMGA still has more of them exact to the base. Outside angiosperm mitochondria
    there is no contest: PMGA does not go there.
  ],
)

= What was run

PMGA's source is distributed as a Singularity sandbox on figshare (3.2 GB). It ran unpatched
on all 22 genomes, one core and 2.6 GB each, in 4 to 26 minutes (one run hung on its node and
finished in nine minutes when run again). Nothing was uploaded to PMGA's web server, since
twelve of the genomes are unpublished. The 22 genomes are ten RefSeq mitogenomes
(_Arabidopsis_, _Liriodendron_, rice, maize, tobacco, beet, grape, _Brassica napus_, soybean,
_Silene latifolia_) and twelve Darwin Tree of Life assemblies picked at regular intervals from
the dataset's angiosperms.

The container does not run its editing step as shipped: the driver script never calls it, and
the Deepred-mt it bundles has no TensorFlow. With your go-ahead I fetched the TensorFlow 2.4.1
image and three Python packages, and ran PMGA's own editing script on its own intermediate
files with one path line changed, at its own probability cut-off of 0.9. The container's two
reference sets are not the 29-species and 319-species sets of the paper: both hold the same
mitochondrial references, 5 to 15 sequences a gene, and the second adds plastid genes.

poa is shown at three stages. "Before" is poa as committed on 23 November, with angiosperm
references built from 11 RefSeq mitogenomes, most of them annotated with the genome's own
translation. "Edited RefSeq references" is the first step of the rebuild: the 14 RefSeq
mitogenomes annotated with RNA-edited proteins. "Now" adds genes of #n-species Darwin Tree of
Life angiosperms with the editing their own public RNA-seq shows applied, for the plastid as
well as the mitochondrion; profiles that record which columns are edited in the references of
at least two genera; a cost for an edit that depends on whether it falls at such a column
(1.5) or not (10), in place of 2.5 everywhere; and genes started at a start codon and ended
at their stop codon (below). The last of these is in the working tree of poa, orfedit and
transsplice, not yet committed.

Every "now" figure for a genome or species is from references built without its genus, since
most of the ten RefSeq genomes are among poa's references and the Tree of Life species are
now among them too. I could not do the same for PMGA or for Deepred-mt, and would expect most
of the ten RefSeq genomes to be among their training genomes.

= Gene structure against the RefSeq annotation

Each annotated coding sequence of the canonical mitochondrial genes (333 across the ten
genomes) was matched to the tool's best model of that gene. An exon counts as placed when
both its ends are within 3 bp of the annotated ones.

#figure(
  tsv-table(read("data/pmga_structure.tsv")),
  caption: [Gene models of ten RefSeq angiosperm mitogenomes against their annotation.],
) <tab-structure>

PMGA places every exon of 81% of the coding sequences (@tab-structure); poa placed 65% before
and places #pct-now now. On spliced genes poa is level or ahead, trans-spliced nad1, nad2 and
nad5 included. On single-exon genes it has #single-now of 245, as PMGA does. It had 164 until
two changes. Its genes now end at their stop codon: the aligner had been reading on over a
gene's stop wherever the reference profile ran a few residues further (189). And they start
at a start codon: the alignment starts a gene wherever the profile's first columns match,
and a gene that began elsewhere now starts at the nearest ATG within ten codons (200).
Allowing for the choice of start codon and for the stop, poa places #ss-now of 333 and PMGA
292; exact to the base, PMGA has 258 and poa 247. What is left is mostly genes whose annotated
start lies far from the conserved part (atp6) or differs from record to record (rps4, mttB,
rpl16, matR); PMGA follows the annotation's own conventions there, because its references
are curated from such annotations.

= Editing against curated edited proteins

Five of the ten RefSeq records (_Arabidopsis_, rice, beet, _Brassica napus_, _Silene_) are
annotated with RNA-edited proteins. A codon whose annotated residue differs from the genome's
own translation is an edit; a predicted edit is right when it falls in such a codon and gives
that residue.

#figure(
  tsv-table(read("data/pmga_edits.tsv")),
  caption: [Predicted edits against the edits annotated in five RefSeq mitogenomes, codon by
  codon. Calls outside the annotated coding sequences, and silent ones, cannot be judged and
  are left out.],
) <tab-edits>

Deepred-mt, as PMGA runs it, finds 93% of the annotated edits and 88% of its calls are right
(@tab-edits). poa found 6% before: at a conserved edit site its profile held the unedited
residue, and there was nothing to correct towards. It now finds #found-gb with #right-gb of
its calls right, with these five genera out of its references.

= Against RNA-seq

For five of the twelve dataset assemblies there is public RNA-seq with mitochondrial reads
(_Filago_, _Levisticum_, _Populus alba_, _Bellardia_, _Nardus_; the libraries are mostly
poly-A selected, so about 55 genes are covered end to end). Neither PMGA nor Deepred-mt has
seen these assemblies, and poa's figures are again with each genus held out.

#figure(
  tsv-table(read("data/pmga_rna.tsv")),
  caption: [Gene models and predicted edits of five Darwin Tree of Life angiosperms against
  their RNA-seq. A coding sequence reads cleanly when, with the edits the RNA shows applied,
  it translates without an internal stop. The RNA's edits include silent ones, which poa
  cannot call.],
) <tab-rna>

Against the reads (@tab-rna) PMGA with Deepred-mt finds 84% of the edits, with 93% of its
calls right. poa found 2% before and now finds #found-rna, with #right-rna of its calls right:
nearly as many edits, and about a quarter as many wrong calls (9 against 33). On gene models
the two are level.

= What changed in poa

Four things, the first three in the order they mattered for editing.

The references. An editing profile corrects a gene towards its reference proteins, so they
have to be the edited ones. RefSeq has 14 angiosperm mitogenomes annotated that way; the rest
of the new references are read from RNA-seq, which the European Nucleotide Archive holds for
several hundred of the Tree of Life angiosperms. On 19 RefSeq genomes, each with its genus
held out, that took poa from 6% of annotated edits to 84%.

A fault in how profiles were scored. A residue that a protein hardly uses, such as cysteine
in a small ribosomal protein, scored as a match in every column, and poa called an edit to it
wherever one C-to-U would make it: arginine to cysteine at all three arginine codons of
plastid rps19, in 44 of 46 angiosperms. This had been there from the start, in every lineage's
profiles, and all of them are rebuilt.

Editing sites. RNA-derived proteins say not only what the edited protein is but where the
edits are. A profile now records the columns edited in at least two genera, and an edit away
from them costs four times what it did. In plastids, where a plant has about thirty edits and
fast-evolving genes differ from the profile at many more places, this took the share of
poa's calls that are in the RNA from 40% to #right-pltd without losing a correct one
(@tab-poa); the plastid was never PMGA's subject, so there is no comparison to make there. In
the mitochondrion it also let an edit at a known site cost less, which is where the last five
points of annotated edits came from.

The two costs suit a plant that edits where angiosperms do. Gymnosperms, hornworts and
quillworts use the same references and edit heavily elsewhere, so poa decides for each
species from its own genes whether the sites apply: of 89 angiosperms annotated with the
references as built, 85 take them, four sedges do not, and neither do Ginkgo, a yew, a
hornwort, a spikemoss or a quillwort, which keep the single cost and their edit calls.

Gene starts and ends. A reading frame ends at its first stop codon, but the aligner treated a stop as
one more codon it could leave out, so a gene that ended a few residues before its reference
profile was read on over its stop. Genes now end at their first stop; a stop the protein
goes on matching beyond (a pseudogene, mostly) is kept in the gene and reported; and a frame
change in the last few codons, which cannot be told from a gene that simply ends differently,
is reported but no longer followed. On 19 RefSeq mitogenomes and 12 plastomes, with the same
references, genes ending within 3 bp of the annotated stop went from 397 to 443 of 457 and
from 558 to 621 of 626, and rbcL, which a short search window had been cutting at codon 341,
is whole. A gene whose alignment began at a codon other than ATG or ACG (an edit away from
AUG) now starts at the nearest ATG within ten codons: starts within 3 bp of the annotation
went from 384 to 401 of 457 and, on 13 plastomes, from 564 to 613 of 671, with one made wrong. Neither change
costs frameshift detection: of 1-bp indels planted in real genes as many are found as
before.

#figure(
  tsv-table(read("data/poa_heldout.tsv")),
  caption: [poa on #n-rna Darwin Tree of Life angiosperms with RNA-seq and 19 RefSeq
  mitogenomes, each annotated with references built without its genus, all three columns
  with genes started at a start codon and ended at their stop. "Edited RefSeq references" has no Tree of Life species in
  it. "RNA-derived references" is at one edit cost everywhere; "with editing sites" is poa
  now.],
) <tab-poa>

= What each tool does that the other does not

PMGA annotates angiosperm mitochondria, and its authors say so. Beyond the protein genes it
writes tRNAs, rRNAs, repeats, plastid-derived genes, a GenBank file and a genome map. For an
input in several contigs it joins them end to end and writes one circular record, and a
trans-spliced gene's own row can have its start after its end. Its editing prediction needs a
runtime the published container does not have.

poa annotates plastids as well, and has references for mosses, liverworts, ferns, horsetails
and clubmosses, picked from the genes themselves. On liverworts the angiosperm references
place 176 of 1,696 annotated exons and the liverwort ones 1,542, so a tool with angiosperm
references alone is not usable there. poa also reports frameshifts, takes an assembly graph,
gives a quality table per sample, and has been run over 1,454 assemblies. It applies the edits
it calls and writes the edited protein, which PMGA does not.

= Where that leaves the claim

On angiosperm mitochondria poa is now within a few points of PMGA on editing, by a different
route: Deepred-mt predicts an edit from the sequence around a cytidine, poa from what the edit
does to the protein and from where edits have been seen. Deepred-mt finds a few more edits and
can call a silent one, which poa cannot; against RNA-seq poa's calls are more often right, and
it gives the edited protein. On gene models the two are now level, with PMGA ahead on how
many ends are exact to the base. Beyond angiosperm mitochondria, in the plastid
and in the other lineages, poa is the only one of the two that runs.

#note(title: "Caveats")[
  Ten RefSeq genomes and five assemblies with RNA are a small test, though poa's own held-out
  figures on #n-rna species and 19 genomes say the same (@tab-poa). The RefSeq annotations are themselves
  imperfect; the PMGA paper found 9% of published splice sites wrong, and some of the "edited"
  annotations are predictions. Public RNA-seq is often from another plant than the one
  assembled, and where the two differ the reads show a difference that is not an edit; poa
  leaves such species out of its evidence for editing sites, but they are in the counts here.
  PMGA's supplement was not read, its web server may behave differently from the container,
  and the editing step was run outside the container's own driver.
]
