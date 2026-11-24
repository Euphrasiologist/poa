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
- `map4.sh` - a public paired RNA-seq run (ENA) mapped to a dataset
  species' own contigs, pileups for both organelles. ENA holds the poly-A
  libraries of most Darwin Tree of Life bryophytes and many angiosperms
  (portal API, `tax_tree(<taxid>)`, `library_source="TRANSCRIPTOMIC"`).
- `detect_cal.py`, `detect_cal_sum.py` - every dataset species' mito
  genes scored against each lineage's profiles and the default ones, and
  the per-gene margins by taxonomic group: the calibration of
  `lineage.MIN_MARGIN`.
- `scope_sum.py`, `rna_sum.py` - runs summed by taxonomic group: without
  RNA (genes aligned, complete, edits and frameshift calls per species)
  and from the RNA-scored tables of `assess.py`/`assess_pltd.py`.
- `gb2fa.py` - GenBank records to FASTA, for annotating RefSeq genomes as
  tests (`score_gb.py`; the liverwort held-out test).
- `score_models.py`, `prot_score.py`, `pmga2poa.py` - poa against PMGA (Li
  et al. 2025) on the same genomes: gene models against a RefSeq
  annotation, proteins against RefSeq's RNA-edited ones, and PMGA's GFF
  as poa-shaped tables so `assess.py` scores it against RNA-seq.

## Angiosperm references and editing sites (2026-11, third round)

- `rna_refs.py` - RNA-edited references for one organelle from poa runs
  and the species' pileups: single-exon proteins with the residues an
  edit made in lower case, and spliced loci. `RNA_REVERSE=0` for a
  lineage that edits C to U only: the U-to-C differences the reads then
  show measure how far the reads are from the assembly's own plant, and a
  species with too many gives proteins but no editing sites, or is left
  out (`species.tsv`).
- `map4.sh` now starts a download again when the file arrives whole but
  with the wrong checksum (resuming such a file never repairs it).
- `cv.py` takes the edit costs and the site settings from the
  environment (`CV_KNOWN`, `CV_NOVEL`, `CV_SITE_GROUPS`,
  `CV_PSEUDOCOUNT_MASS`, plastid ones apart).
- `call_frac.py` - every predicted edit against the reads at its site:
  edited, edited in under a fifth of reads, not edited, or not covered.
  Fairer than counting against a list of RNA edits, which for public
  RNA-seq holds the differences between two plants as well.
- `pair_sum.py` - RNA-scored tables of the same species under different
  references, summed over the genes every one of them has.
- `edit_score.py` - predicted edits against the edits a GenBank record's
  annotation applies, codon by codon.
- `cv_sites.sh`, `site_score.sh` - the held-out runs behind the figures in
  the README (angiosperms, genus held out, at given edit costs) and their
  scoring; `heldout_table.py` makes the table of
  `reports/pmga_comparison` from them.
- `cmp_runs.py` - poa runs of the same species before and after a
  change, by taxonomic group: lineage calls, genes aligned and complete,
  edits and frameshift calls.
- `detect_cal_sum2.py` - `detect_cal.py`'s scores as per-gene margins by
  taxonomic group (the 2026-11-24 calibration of `lineage.MIN_MARGIN`).
- `angio_readme.py` - writes
  `editing/reference/lineages/angiosperm/README.md` (species, runs, and
  whose reads count as evidence of editing sites) from `rna_refs.py`'s
  output.

`assess.py`, `assess_pltd.py`, `score_gb.py`, `calls_support.py`,
`call_frac.py`, `pair_sum.py` and `summ.py` sat in a `pilot/` directory
under the others and import them from the directory above.
