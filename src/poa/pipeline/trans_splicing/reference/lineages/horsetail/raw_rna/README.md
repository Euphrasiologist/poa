# Horsetail (Equisetum) gene models from RNA-seq

There is no annotated horsetail mitogenome in GenBank, so these
references are all made here: mitochondrial genes of nine *Equisetum*
assemblies from the Darwin Tree of Life dataset, with their exons placed
by, and the RNA editing taken from, RNA-seq of the same species.
*Equisetum* is the only living genus of horsetails.

| species | RNA-seq | spliced-gene loci | single-exon proteins |
|---|---|---|---|
| Equisetum sylvaticum | 50726_2#64 | 11 | 23 |
| Equisetum telmateia | 50812_5#40, SRR4061752, SRR4061753 | 8 | 21 |
| Equisetum ramosissimum | SRR5499399 | 2 | 9 |
| Equisetum palustre | 45657_1#3 | 2 | 8 |
| Equisetum pratense | 49913_2#21 | 2 | 8 |
| Equisetum variegatum | 47503_1#17 | 0 | 8 |
| Equisetum hyemale | 50726_2#63, SRR4061756, SRR4061757 | 1 | 7 |
| Equisetum arvense | 46781_2#93, SRR4061754, SRR4061755 | 1 | 5 |
| Equisetum fluviatile | 51142_1#43 | 0 | 0 |

Runs named like 50726_2#64 are Sanger poly-A libraries of the Darwin Tree
of Life specimen; SRR accessions are public runs of the same species from
other plants (PRJNA340020, PRJNA384992). Poly-A libraries hold few
mitochondrial reads: *E. sylvaticum* and *E. telmateia* reach 150 to 280
reads over the mitochondrial contigs and give nearly every gene, the rest
a handful each.

Method (2026-11, `dev/rna_editome/boot.py`, `mk_lineage.py`). Each run
was mapped with bwa-mem2 to the species' own mitochondrial and plastid
contigs together and the runs of a species pooled. poa, on whatever
references it had, said roughly where each gene is; then

- the introns are the junctions intron-spanning reads show (split
  alignments), with at least 5 reads that are at least a tenth of the
  reads at the donor (3 next to a junction poa had called);
- the transcript is followed from the middle of poa's call in both
  directions across every such intron, 1.5 kb past the called ends;
  poa's own junction is used only where the reads show none nearby (a
  22 bp exon no read piece aligns to), and where an exon lies in a repeat
  the reads place in another copy, the copy poa called is kept;
- a C read as U in at least 20% of 10 or more reads (at least 2) is an
  edit. Horsetails show no U-to-C editing: over the mitochondrial contigs
  of *E. sylvaticum* 700 sites read C-to-U against 15 U-to-C, about
  what any other kind of mismatch comes to;
- the gene is the open reading frame holding most of the called bases,
  from the ATG nearest where a homologue's protein starts (*Marchantia
  polymorpha*, NC_001660) to its stop;
- a model is used if at least 90% of it has 10 or more reads, its exon
  count is the one most species show and its protein is within 10% of
  their median length.

What that shows of *Equisetum* genes: cox1 has three introns
(angiosperms have none); nad1 has five exons but not the angiosperm five
(angiosperm exons 2 and 3 are one exon, and there is an intron
angiosperms lack); nad7 has four exons (five), rps3 three (two), rps10
no intron; nad4's first exon lies on the other strand from the rest in
these assemblies, so is trans-spliced; nad2, nad5, ccmFc, cox2 and rpl2
are split as in angiosperms.

Files: `RNA_<species>.fasta` in the format of a lineage's `raw/` (one
record per exon, genomic sequence, and the edited protein as
`translation_REFSEQ_ANNOTATED`) and `junctions.tsv`; the single-exon
proteins are in `editing/reference/lineages/horsetail/rna_proteins.fasta`.
A template or profile is built from as little as one species
(`exon_profiles.LINEAGE_MIN_LOCI`), with fewer pseudocounts than the
other lineages' (`LINEAGE_PSEUDOCOUNT`): the species of the one genus
differ by a few percent.
