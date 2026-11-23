# RNA-edited fern loci

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

Files: `RNA_<species>.fasta` in the format of `../raw/` (one record per
exon, genomic sequence, and the edited protein as
`translation_REFSEQ_ANNOTATED`), and `junctions.tsv`. `exon_profiles.py`
reads this directory next to `../raw/`; `fetch_references.py` does not
touch it.

## atp1, cox1 and rps14

These three genes have one intron each in ferns and none in angiosperms,
so poa's call for them was a single-exon one that ran across the intron
and there were no transsplice exons to start from. Their loci here are
gene models built from the reads (2026-11, `dev/rna_editome/boot.py`):

- the introns are the junctions intron-spanning reads show (bwa-mem2
  split alignments), with at least 5 reads that are at least a tenth of
  the reads at the donor;
- the transcript is followed from the middle of poa's call in both
  directions, across every such intron, and the edits above applied;
- the gene is the open reading frame holding most of the called bases,
  from the ATG that a GenBank fern's protein starts at (*Haplopteris*
  rps14 and cox1, *Azolla rubra* atp1) to its stop;
- a model is used if at least 90% of it has 10 or more reads, its exon
  count is the one most species show, and its protein is within 10% of
  the median length.

That gave atp1 in 8 species, cox1 in 7 and rps14 in 6, two exons each,
which is also what the GenBank ferns annotate (*Haplopteris* atp1 alone
has no intron). Held out by genus, on the eleven ferns: atp1 complete in
11 of 11 (none of 19 single-exon calls aligned end to end before, and 6
drew a frameshift call), cox1 in 9 of 11 (0 of 12), and rps14 without the
frameshift call it drew in 12 of 12 copies; of the RNA edits in them 84%,
77% and 91% found, with 99%, 99.6% and 95% of calls right (before: 49%,
39% and 34% found). Of the 22 junctions with intron-spanning reads, 21
agree with them.
