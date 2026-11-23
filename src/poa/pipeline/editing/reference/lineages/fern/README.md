# Fern editing references

`mito/*.pssm` are built by `editing/src/00b_build_lineage_profiles.py --lineage fern` from the GenBank ferns with RNA-seq editomes (`reference_genes.FERN_EDITED_MITOGENOMES`) and from `rna_proteins.fasta` (headers `species|gene`).

`pltd/*.pssm` are the plastid profiles, built with `--organelle pltd` from `rna_proteins.pltd.fasta`: 614 proteins of 56 single-exon plastid genes from the same eleven ferns and runs, by the same rules. See `trans_splicing/reference/lineages/fern/plastid/raw_rna/README.md` for how they do.

These references are not GenBank records. They are mitochondrial genes of
eleven Darwin Tree of Life fern assemblies with the RNA editing their own
RNA-seq shows applied, made because no leptosporangiate fern outside
Salviniales and *Haplopteris* has a deposited mitochondrial editome, and
profiles built from unedited proteins give orfedit and transsplice
nothing to correct towards.

| species | RNA-seq run (Sanger, DToL poly-A library) | spliced-gene loci | single-exon proteins |
|---|---|---|---|
| Asplenium scolopendrium | 51987_7#38 | 8 | 17 |
| Asplenium trichomanes | 50409_7#54 | 8 | 18 |
| Asplenium viride | 47846_3#40 | 8 | 1 |
| Azolla filiculoides | 53189_7#64 | 8 | 4 |
| Dryopteris borreri | 49913_2#23 | 8 | 17 |
| Dryopteris dilatata | 49913_2#24 | 8 | 12 |
| Dryopteris filix mas | 50726_2#46 | 7 | 19 |
| Onoclea sensibilis | 49913_2#29 | 8 | 18 |
| Phegopteris connectilis | 47503_1#12 | 7 | 16 |
| Polystichum setiferum | 48017_1#27 | 8 | 19 |
| Woodsia ilvensis | 48017_1#35 | 1 | 12 |

Method (2026-11): each species' reads were mapped with bwa-mem2 to its
own mitochondrial and plastid contigs together, and the mitochondrial
contigs piled up (base quality >= 20). In a gene poa had called:

- an edit is a C read as U, or a U read as C, on the gene's strand in at
  least 20% of reads (at least 2). Poly-A libraries carry unprocessed
  transcripts, so sites are often edited in well under half the reads;
- a stop codon left after that is read through if a U-to-C site in it is
  edited in at least 5% of reads (at least 2);
- a stop codon with 20 or more reads at each base and no such evidence is
  a real stop: the first in the second half of the gene ends the protein
  (poa's call ran past it), one in the first 15% moves the start to after
  it;
- the gene is used if no other stop is left, at least 90% of it has 10 or
  more reads, and its length is a multiple of three.

For spliced genes the exons are transsplice's, from a run on templates
that did not include the species' own genus, with each exon boundary
moved to the junction the intron-spanning reads show where at least 5 of
them, and three times as many as agree with the call, show another.
