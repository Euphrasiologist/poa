# RNA-seq editomes for reference building

The scripts that made the fern references in
`src/poa/pipeline/trans_splicing/reference/lineages/fern/raw_rna/` and
`src/poa/pipeline/editing/reference/lineages/fern/rna_proteins.fasta`
(2026-11). They are kept as a record of the method and are not part of
the package: paths are those of the Sanger farm and of this dataset, and
they read tables from earlier validation runs under `build/`.

- `map.sh SPECIES TAG irods:PATH.cram` - fetch an RNA-seq CRAM, map it
  with bwa-mem2 to the species' mito + plastid contigs, pile up the mito
  contigs (`TAG.pile.tsv.gz`) and keep the mito alignments
  (`TAG.mito.bam`).
- `rna.py` - an RNA edit, a real stop codon, and the edited protein of a
  called gene, from a pileup.
- `rebuild.py` - single-exon genes: edited proteins, and editing profiles
  cross-validated by genus against the RNA.
- `tpl_rebuild.py` - spliced genes: edited loci, exon boundaries moved to
  the junctions intron-spanning reads show, exon templates
  cross-validated by genus.
- `junctions.py` - exon junctions against intron-spanning reads.
- `splice_score.py`, `ends.py` - edit calls in exons, and gene ends,
  against the RNA.
- `wire.sh` - rebuild the bundled profiles and templates with poa's own
  builders and validate.
