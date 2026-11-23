# Clubmoss (Lycopodiaceae) editing profiles

`mito/*.pssm` are built by `editing/src/00b_build_lineage_profiles.py
--lineage lycophyte` from `rna_proteins.fasta` (headers `species|gene`):
50 RNA-edited proteins of 19 single-exon mitochondrial genes, 18 and 19
from the RefSeq annotations of *Phylloglossum drummondii* and
*Lycopodium japonicum* with the editing their own RNA-seq shows applied,
the rest from three Darwin Tree of Life
clubmosses. Sources, runs, method and how the lineage does are in
`trans_splicing/reference/lineages/lycophyte/README.md`, with the
spliced genes (atp6, atp9, cob, cox1, cox3, nad3 and rps14 have introns
in clubmosses and are exon templates there).

The four RefSeq mitogenomes' own proteins are not used: they are
unedited.
