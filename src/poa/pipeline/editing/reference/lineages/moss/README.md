# Moss (Bryophyta) editing profiles

`mito/*.pssm` are built by `editing/src/00b_build_lineage_profiles.py
--lineage moss` from 20 RefSeq moss mitogenomes
(`reference_genes.MOSS_MITOGENOMES`).

`pltd/*.pssm` are built with `--organelle pltd` from
`refseq_proteins.pltd.fasta` (headers `accession|gene`): 1,339 proteins
of 53 genes from 26 RefSeq moss plastomes
(`reference_genes.MOSS_PLASTOMES`, one per genus, not *Takakia*), written
by `--organelle pltd --fetch`. The spliced plastid genes of the same
plastomes are exon templates
(`trans_splicing/reference/lineages/moss/plastid/`).

Moss genes are aligned to these, and to the moss exon templates, at an
edit cost of 6, not 2.5 (`lineage.EDIT_PENALTY`): mosses edit little, and
at 2.5 poa called 579 edits in the plastid genes of 16 Darwin Tree of
Life mosses, 127 of them in their RNA-seq.
