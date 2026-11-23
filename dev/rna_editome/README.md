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

## Gene models, more lineages, plastid (2026-11, second round)

For horsetails there was no annotated genome to start from, and fern
atp1, cox1 and rps14 had no exon calls, so these build the gene from the
reads instead of correcting a call:

- `map4.sh`, `map5.sh` - a public paired run (ENA fastq) mapped to a
  dataset species' contigs, or to a GenBank genome in `extra_genomes/`.
- `boot.py` - gene models from RNA-seq, seeded by a poa run: introns
  from intron-spanning reads, edits applied, the reading frame holding
  the call, start by homology. Pools a species' runs.
- `mk_lineage.py` - a lineage's `raw_rna/` and `rna_proteins.fasta` from
  the models that pass (coverage, modal exon count, length).
- `pltd_refs.py` - RNA-edited plastid proteins and spliced loci from a
  poa run and the plastid pileup (the rules of `rna.py`).
- `cv.py` - a lineage's references held out by species or genus:
  templates and profiles rebuilt without the group, its species
  annotated with them (mito, and plastid for a lineage with plastid
  references), tables in poa's layout.
- `assess.py`, `assess_pltd.py`, `calls_support.py`, `summ.py` - a run's
  calls against the RNA: edits found and right, proteins, exon junctions
  against intron-spanning reads; calls by how much of the RNA shows them;
  gene-level counts without RNA.
- `score_gb.py` - calls on a GenBank genome against its own annotation
  (the clubmoss RefSeq genomes).
- `gb_edit.py` - a GenBank genome's annotated genes with the editing its
  own RNA-seq shows applied (*Phylloglossum* and *Lycopodium japonicum*, the clubmoss references).

`assess.py`, `assess_pltd.py`, `score_gb.py`, `calls_support.py` and
`summ.py` sat in a `pilot/` directory under the others and import them
from the directory above.
