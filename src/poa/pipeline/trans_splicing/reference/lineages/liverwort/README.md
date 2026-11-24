# Liverwort (Marchantiophyta) exon templates

Liverwort mitochondrial genes are split differently from angiosperm ones,
and none of their introns is trans-spliced. Across the 27 RefSeq
mitogenomes used here (`reference_genes.LIVERWORT_MITOGENOMES`, one per
genus: Haplomitriopsida, complex thalloids, simple thalloids and leafy
liverworts) cob has four exons, cox2 and cox3 three, nad2, nad3, nad4,
nad5, atp9, rpl2 and rps14 two, and nad1, rps3, rps10 and ccmFc one,
where the angiosperm templates expect 1, 2, 1, 5, 1, 4, 5, 1, 2, 1, 5,
2, 2 and 2. Three genes vary within liverworts by loss of introns: cox1
has ten exons in 17 genomes and eight, seven or six in the rest, atp1
three or none (12 genomes each, two with one), nad4L three or two. A
template has the full exon set; a species without an intron has the two
exons end to end, and transsplice fills both slots from the one stretch
(0.3.0 on; before that it filled two slots at most from one genomic
exon).

nad7 is a pseudogene in every liverwort but *Haplomitrium* and *Treubia*
and gets no template (`reference_genes.LINEAGE_NOT_RECONSTRUCTED`): one
built from those two drew 187 frameshift calls in 55 of 64 Darwin Tree of
Life liverworts. rpl6 and rps8, which liverworts keep, are not called at
all: poa's gene models are the angiosperm set.

`raw/` holds the loci (`src/00a_fetch_references.py --lineage
liverwort`): 458 verified against the RefSeq protein, 1 not. The
annotated proteins are the genome's own translation, which is right for
the complex thalloids (no RNA editing) and a simplification for the
leafy and simple thalloid liverworts, which edit C to U.

Held out (2026-11). Each of the 27 genomes annotated with its genus left
out of the references, against its RefSeq annotation (genes poa has a
model for): 1,542 of 1,696 annotated exons placed within 3 bp (176 on
the angiosperm references), 305 of 337 intron-containing genes with every
exon placed (0), 17 annotated genes with a frameshift call (144), none
not called (126). Most single-exon genes still off are ccmB, ccmFn and
nad6, whose annotated start differs between the RefSeq records
themselves (ccmB runs from 62 to 277 residues). On 64 Darwin Tree of
Life liverworts, all called liverwort by `--lineage auto`: spliced genes
complete 1,072 of 1,082 (372 of 635 on the angiosperm references),
frameshift calls 20 (590).

## Plastid

`plastid/raw/` holds the spliced plastid genes of 30 RefSeq liverwort
plastomes (`reference_genes.LIVERWORT_PLASTOMES`, one per genus, not
*Haplomitrium*; `src/00a_fetch_references.py --lineage liverwort
--organelle pltd`, 236 loci verified), and `plastid/profiles/` the exon
templates built from them. The GenBank records are cached in
`plastid/genbank/` when fetched and are not shipped.
