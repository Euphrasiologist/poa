# Horsetail (Equisetum) editing references

`mito/*.pssm` are built by `editing/src/00b_build_lineage_profiles.py
--lineage horsetail` from `rna_proteins.fasta` (headers `species|gene`):
89 proteins of 23 single-exon mitochondrial genes from eight *Equisetum*
species of the Darwin Tree of Life dataset, each the gene's reading frame
with the C-to-U editing the species' own RNA-seq shows applied. The
species, runs and method are in
`trans_splicing/reference/lineages/horsetail/raw_rna/README.md`, with the
spliced genes.

cox1 has introns in *Equisetum* and is an exon template there. matR is
included although the *Equisetum* reading frame is two thirds the length
of the angiosperm protein and shows no editing: against the angiosperm
profile it drew 24 edit calls a copy, none of them in the RNA.
