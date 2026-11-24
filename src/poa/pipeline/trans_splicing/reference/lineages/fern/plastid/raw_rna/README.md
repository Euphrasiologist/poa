# RNA-edited fern plastid loci

Plastid genes of the same eleven Darwin Tree of Life ferns as
`../../raw_rna/`, from the same RNA-seq runs and by the same rules (see
the README there): each species' reads mapped to its own mitochondrial
and plastid contigs together, the plastid contigs piled up, and the
editing the reads show applied to the genes poa had called. Plastid
depth is 700 to 9,000 reads, and of the 715 plastid genes poa called in
the eleven species 708 gave a clean protein.

Here: the spliced genes (`reference_genes.PLASTID_SPLICED`) transsplice
had reconstructed completely on the land-plant templates, 79 loci of 9
genes, exons cut to the protein's span. The 613 single-exon proteins (56
genes) are in `editing/reference/lineages/fern/rna_proteins.pltd.fasta`.
(Read again on 2026-11-24 from gene models that end at their stop codons
and are not cut off by the search window: the first set, 94 loci and 614
proteins, had rbcL stopping at codon 341 and matK at 243. The held-out
figures below are the first set's.)
`src/00_build_exon_profiles.py --lineage fern --organelle pltd` builds
`../profiles/` from this directory, and `editing/src/00b_build_lineage_profiles.py
--lineage fern --organelle pltd` the editing profiles.

They are used for a fern's plastid genes with `--lineage fern`, with stop
codons read through but without the other U-to-C edits the
mitochondrial genes get. Held out by genus on the eleven ferns
(2026-11), against the land-plant references:

| | land-plant references | fern plastid references |
|---|---|---|
| single-exon genes aligned end to end | 339 of 614 | 596 of 614 |
| RNA edits found / calls right, single-exon | 74% / 59% | 82% / 83% |
| RNA edits found / calls right, spliced | 65% / 50% | 81% / 83% |
| genes running past a stop codon the RNA never edits | 81 | 36 |

With other U-to-C edits allowed as well, 91% of the plastid U-to-C sites
were found against 77%, but 42% of the U-to-C calls were in the RNA
against 89%, and calls right overall fell to 77%. About one call in seven is a C the
other genera edit, or have as T, and this species leaves alone.
