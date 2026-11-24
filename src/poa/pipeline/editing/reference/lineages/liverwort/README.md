# Liverwort (Marchantiophyta) editing profiles

`mito/*.pssm` are built by `editing/src/00b_build_lineage_profiles.py
--lineage liverwort` from the single-exon genes of 27 RefSeq liverwort
mitogenomes (`reference_genes.LIVERWORT_MITOGENOMES`): 22 genes. The
genes with introns in liverworts (atp1, atp9, cob, cox1, cox3, nad3,
nad4L, rps14) are exon templates instead
(`trans_splicing/reference/lineages/liverwort/README.md`, which has the
held-out figures).

`pltd/*.pssm` are built with `--organelle pltd` from
`refseq_proteins.pltd.fasta` (headers `accession|gene`): 1,532 proteins
of 53 genes from 30 RefSeq liverwort plastomes
(`reference_genes.LIVERWORT_PLASTOMES`), written by `--organelle pltd
--fetch`.

The RefSeq proteins are the genome's own translation. Complex thalloid
liverworts do not edit their organelle transcripts, and on these
profiles poa calls next to no edits in them (12 mito edits in nine
Darwin Tree of Life species, against 942 on the angiosperm profiles).
Leafy and simple thalloid liverworts do edit, and unedited references
find only part of it: in *Nowellia curvifolia*, whose RNA-seq shows 137
mito edits in 39 genes, 63 are called and 33 of those are right.
