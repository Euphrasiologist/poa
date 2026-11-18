# Changelog

## Unreleased

- Frameshifts: with orfedit >= 0.2.0, the editing correction allows a
  change of reading frame inside a single-exon gene and reports where it
  is - `tables/frameshifts.tsv` (the skipped bases' coordinates, 1- or
  2-base type, profile column), a `frameshift` child of the gene's CDS in
  the GFF, and `n_frameshifts`/`n_new_frameshifts`/`frameshifted_genes` in
  the QC table, printed at the end of `poa run`. Each frameshift carries
  `dtol_assemblies`, how many DToL assemblies have one at the same place
  in the same gene (`editing/reference/frameshifts_dtol.tsv`): mito rpl16
  has one in 239 of 979 mito assemblies, so that's the gene, while one
  seen nowhere else is the first thing to check against the reads. On the
  DToL dataset, 1-bp indels planted in 900 genes are found in 83% (94% of
  plastid genes, 72% of mito), 99.9% with the right type, and every gene
  without a frameshift is corrected exactly as before. With orfedit 0.1.0
  there are no frameshift calls (a warning says so); `scan_editing
  --frameshift-penalty inf` turns them off.
- Circular contigs: a gene running through the start of a circular contig
  (an oatk contig not marked `circular=false`, or a finished, linearised
  organelle FASTA) was called as two partial genes, one at each end.
  `poa run` now moves such a contig's origin into its longest gene-free
  stretch, calls genes again, and maps the GFF and tables back onto the
  input, the gene as one feature ending past the contig's length (GFF3's
  convention for `Is_circular=true` sequences). With a GFA the contig is
  left as it is, with a warning, as unitig coordinates need the
  assembler's origin.
- `poa run` ends with single-sample QC: a `pass`/`flag`/`fail` line per
  organelle, with reasons, in `<name>.qc.tsv`. It applies the dataset QC
  gate's rules (dead ends, fragmentation, circularity, unjoined contigs,
  core genes, genome size against the genus, mito/plastid
  cross-alignment), with the dataset-derived parts - the core-gene set and
  per-genus genome sizes from the DToL assemblies - shipped as reference
  data (`qc_basic_stats/reference/`). Re-run over all 2,906 DToL
  assemblies, it gives the dataset gate's status and reason strings for
  every assembly that has files, apart from two where the dataset gate is
  wrong: a stale graph-stats row (Silene vulgaris mito) and an assembly
  with no core genes recorded as unknown rather than 0% (Neottia ovata
  mito). `--genus` sets the genus for the size check (default: the first
  word of `--name`). Checks that need a GFA (gfatk) or both organelles
  (minimap2) are listed as not assessed when they can't run.
- RNA editing correction now covers 26 single-exon mito genes, up from 7.
  Genes the dataset-built profiles don't cover (they need a gene to be
  single-copy in 80% of species) use bundled profiles built from the 11
  angiosperm RefSeq mitogenomes (`editing/reference/profiles_refseq/`).
  On Arabidopsis, with Arabidopsis held out of the profiles, the median
  protein identity of the 19 newly corrected genes against RefSeq rises
  from 0.837 to 1.000 (rps12 0.57 -> 1.00, atp9 0.51 -> 1.00). The 7
  existing profiles, the plastid annotation and the moss lineage call are
  unchanged.
- `00b_build_lineage_profiles.py` takes `--lineage angiosperm`, and
  `--exclude ACC` / `--out-dir` for held-out validation.
- Editing correction covers every substantial mito copy of a gene, not just
  the best-scoring hit: any other hit covering at least half the HMM that
  isn't another gene's stronger locus, kept if orfedit aligns at least 80%
  of the profile. Arabidopsis's second atp6 (a different protein, 349 aa)
  goes from raw (0.785 identity to RefSeq) to edited (0.975). Across the
  DToL gene calls this is up to 1,366 extra mito copies (atp1 163, atp9
  107). Plastids still correct the best hit only.
  `edited_gene_calls.tsv` and `predicted_edits.tsv` gain `contig`,
  `hit_start` and `hit_end` columns naming the copy.
- Split ccmFn: where mitochondrial ccmFn is two genes (Brassicaceae,
  Allium, Galeopsis - 37 species in the DToL dataset), its two nhmmscan
  hits are named ccmFn1 and ccmFn2, as in RefSeq, and each is corrected
  against its own profile, built from 11 RefSeq Brassicaceae and Allium
  mitogenomes. Before, gff_export could chain the pieces into one
  cis-spliced ccmFn when they were close, and editing corrected the first
  piece against the full-length gene (Arabidopsis ccmFN1 1380 bp vs
  RefSeq's 1149). With Arabidopsis held out, ccmFN1 goes from 0.826 to
  0.995 identity and ccmFN2 from 0.784 (raw) to 1.000.
- The Arabidopsis smoke check keeps every RefSeq copy of a gene (atp6 x2,
  ccmFN1/ccmFN2) and pairs copies with GFF models one-to-one; before, the
  second RefSeq atp6 overwrote the first.

## 0.2.0

Annotations are unchanged from 0.1.0. On the full dataset each reworked
stage gave the same output as before (trans_splicing's tables up to row
order and ties between equal-scoring repeat copies), and `poa run`'s GFFs
are byte-identical. The changes are in how poa is built and called.

- Every stage `poa run` uses is now an importable module, and `poa run`
  calls them as functions in a single process instead of running each stage
  script as a subprocess. Each function takes its paths as arguments and
  raises `StageError` on failure. The `src/NN_*.py` stage scripts are now
  thin wrappers with the same flags as before.
- `poa.run.run()` takes keyword arguments
  (`run(mito="x.fasta", outdir="out", lineage="moss")`), so a sample can be
  annotated from Python.
- The command line uses click. Bad options (a missing input file, an
  unknown `--lineage`, `-t 0`, no `-o`) are usage errors with exit status 2,
  reported before anything runs, and `-h` works as well as `--help`.
- Pipeline modules report through `logging`, in the same `[info]`,
  `[warn]` and `[err]` format as before.
- New dependency: click >= 8.
- Code formatted with black (line length 120); tests clean up the
  temporary directories they create.

## 0.1.0

First release.
