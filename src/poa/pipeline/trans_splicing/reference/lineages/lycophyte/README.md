# Clubmoss (Lycopodiaceae) exon templates

Clubmoss genes are split differently from angiosperm ones. In all four
RefSeq Lycopodiaceae mitogenomes (`reference_genes.LYCOPHYTE_MITOGENOMES`:
*Huperzia squarrosa* NC_017755.1, *Huperzia crispata* NC_071971.1,
*Phylloglossum drummondii* NC_086557.1, *Lycopodium japonicum*
NC_080981.1) atp6 and cob have three exons, atp9 and cox1 four, cox3,
nad3 and rps14 two, where angiosperms have one; nad2 has three exons,
nad4 two, cox2 four, rps3 three, where angiosperm templates expect five,
four, two and two. The ccm genes and nad7 are missing.

`raw/` holds those genomes' spliced genes as annotated
(`src/00a_fetch_references.py --lineage lycophyte`), as a record of the
exon structures. Their proteins are the genome's own translation except
at four in-frame stop codons the records mark as edited (cob codon 46,
nad1 177, nad5 544, rpl2 36), so they are left out of the templates
(`LINEAGE_UNEDITED_MITOGENOMES`): templates built from them placed exons
but found a fifth of the edits RNA-seq shows, fewer than the angiosperm
references do.

`raw_rna/` holds what the templates are built from, loci with the RNA
editing applied:

| source | RNA-seq | spliced-gene loci | single-exon proteins |
|---|---|---|---|
| *Phylloglossum drummondii* NC_086557.1, exons as annotated | SRR12832746 (total RNA) | 15 | 18 |
| *Huperzia selago*, Darwin Tree of Life assembly | 44268_8#11, SRR1994288, ERR2040865, ERR2040868 | 2 | 8 |
| *Diphasiastrum alpinum*, Darwin Tree of Life assembly | 49772_1#22 | 1 | 4 |
| *Spinulum annotinum*, Darwin Tree of Life assembly | ERR2040871 | 1 | 1 |

*Phylloglossum* is the one clubmoss with deep RNA-seq of its
mitochondrion: a total-RNA library that covers the RefSeq genome about
2,000 times. Its reads were mapped to NC_086557.1 and every annotated
gene edited by the rules in `../fern/raw_rna/README.md` (an edit is a C
read as U, or a U read as C, in at least 20% of 10 or more reads); all 35
came out clean, with 289 C-to-U and 5 U-to-C edits (`dev/rna_editome/
gb_edit.py`); sdh3 and matR, which have introns in *Phylloglossum* alone,
are not used. The three Darwin Tree of Life clubmosses have poly-A
libraries with few mitochondrial reads and give the genes that happen to
be covered, as gene models built from the reads
(`dev/rna_editome/boot.py`, see `../horsetail/raw_rna/README.md`). A
template is built from as little as one locus
(`exon_profiles.LINEAGE_MIN_LOCI`), and most are *Phylloglossum* alone.

Held out (2026-11). Three RefSeq genomes annotated with their genus left
out of the references (*Phylloglossum* cannot be, being the reference):
88% of 198 annotated exons placed within 3 bp (26% on the angiosperm
templates), 41 of 44 spliced genes with every exon placed (4), 5 genes
with a frameshift call (14). Three Darwin Tree of Life clubmosses, each
with its own loci left out, against their RNA-seq: in the 18 genes with
enough reads, 87 of 100 calls right and 87 of 98 edits found (65 of 147
and 65 of 119 on the angiosperm references); all 21 exon junctions that
intron-spanning reads cross as the reads show; spliced genes complete 45
of 45 (15 of 27). What stays unplaced is mostly atp9's third exon, which
is 7 bp.
