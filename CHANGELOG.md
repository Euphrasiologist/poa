# Changelog

## Unreleased

- Fern references rebuilt from RNA-seq editomes of eleven Darwin Tree of
  Life ferns (nine Polypodiales, Azolla filiculoides): each species'
  RNA-seq mapped to its own assembly, the editing it shows applied to the
  genes poa calls, bundled as `trans_splicing/reference/lineages/fern/
  raw_rna/` (79 spliced-gene loci) and `editing/reference/lineages/fern/
  rna_proteins.fasta` (148 single-exon proteins), each with a README of
  species, runs and method. With `--lineage fern`:
  - the six GenBank ferns whose proteins only have stops edited out are
    left out of the exon templates (`LINEAGE_UNEDITED_MITOGENOMES`);
  - editing profiles use 2 pseudo-observations per column shared out by
    composition, as the exon slots do (`LINEAGE_PSEUDOCOUNT`), and
    template columns under half the references have are dropped
    (`exon_profiles.LINEAGE_GAP_THRESHOLD`);
  - other U-to-C edits are called at a cost of 5 on fern profiles and
    templates (`lineage.REVERSE_EDIT_PENALTY`; orfedit and transsplice
    `--reverse-edit-penalty`), not on plastid genes or genes on borrowed
    profiles.
  Held out by genus against the RNA: single-exon genes 80% of C-to-U
  edits found with 96% of calls right (42% and 92% on the previous fern
  profiles, 46% and 70% on the angiosperm ones) and 88% / 97% for U-to-C;
  spliced genes 82% / 97% and 92% / 96% (9% / 74% and 50% / 93% before).
  On 18 DToL ferns that supplied no RNA: spliced genes aligned over 99%
  74 -> 137 of 159, genes with a frameshift call 29 -> 7. Spliced genes
  running past a never-edited stop codon 52 -> 13 of 89. The whole-gene
  column mask is now empty. Exon junctions agree with intron-spanning
  reads at 116 of 144, as before.
  The single-exon figures in the two entries below counted nine ferns
  whose assemblies had been replaced after their gene calls were made;
  without them the previous fern profiles align 329 of 397 copies over
  99% (the new ones 351) and the angiosperm ones 202.
- `--lineage fern`, for ferns, quillworts and hornworts (never picked by
  `auto`): the angiosperm references, with orfedit and transsplice
  allowed to read a genomic in-frame stop codon as Gln or Arg by a
  U-to-C RNA edit, at the cost of a C-to-U edit (`common/lineage.py`,
  `STOP_EDIT_PENALTY`; needs orfedit and transsplice builds with
  `--stop-edit-penalty`, and stops with an error otherwise). These
  lineages edit both ways - a Polypodiales mito gene has about nine such
  stops - and without it each one cost a skipped column, a frameshift
  call or an early gene end. On 30 DToL fern assemblies: single-exon
  mito copies aligned over 99% of their profile 152 -> 292 of 702,
  frameshift calls 277 -> 113; spliced mito genes aligned over 99%
  0 -> 80 of 249, plastid 74 -> 117 of 233; U-to-C comes to 32% of mito
  and 6% of plastid edits, against 33-35% and 1-10% measured in fern
  transcriptomes. Edit tables gain an `edit` column (`C>U`/`U>C`) and the
  GFF's `sequence_alteration` features say `edited_from=T;edited_to=C`
  for the reverse ones. In a dataset run, species recorded as fern in
  `lineage.tsv` go through transsplice as a batch of their own.
- Fern editing profiles (`editing/reference/lineages/fern/mito/`, 22
  single-exon mito genes), used with `--lineage fern` in place of the
  angiosperm ones: from the five ferns whose deposited proteins carry
  their RNA-seq editome (`FERN_EDITED_MITOGENOMES`: Haplopteris
  ensiformis and the per-gene records of Azolla rubra, A. pinnata,
  Marsilea mutica and Salvinia molesta), one protein per species and
  gene (`lineage_profiles.one_per_species`). On 538 single-exon mito
  copies in 29 DToL ferns, against the angiosperm profiles with stop
  read-through: aligned over 99% of the profile 231 -> 412, frameshift
  calls 46 -> 44, planted 1-bp indels found 43% -> 61%, C-to-U calls per
  copy 11.0 -> 7.3. Cross-validated on the five species' own editomes
  (each held out): 29-41% of real C-to-U edits found with 60-94% of
  calls right, against 47-60% and 49-71% on the angiosperm profiles;
  stop read-throughs found 85-100% on both. matR and nad9 are left on
  the angiosperm profiles. `column_mask.py --editing` can mask editing
  profile columns the dataset's ferns lack; none reached the 80%
  threshold, so `profile_mask.tsv` is empty.
- Fern exon templates (`trans_splicing/reference/lineages/fern/`), used
  for spliced mito genes with `--lineage fern`: from seven annotated
  leptosporangiate fern mitogenomes (`FERN_MITOGENOMES`: three
  Woodwardia, Blechnopsis, Dryopteris, Osmunda, and Haplopteris's nine
  chromosomes; none RefSeq, and Dryopteris drops out because its
  annotation does not match its sequence). nad1, nad4, nad5 and nad7 have
  four exons, nad2 five, rps3 three, cox2 and rpl2 two, rps10 one; no
  ccmFc; every intron cis. Reference proteins are checked against their
  genomes allowing U-to-C as well as C-to-U edits (78 of 91 loci pass),
  and slot profiles are cut from the annotated, edited proteins, since a
  fern exon's own translation has stop codons in it. Whole-gene profile
  columns that at least 80% of the DToL fern reconstructions leave
  unaligned are dropped (`whole_gene_mask.tsv`, `column_mask.py`;
  `00_build_exon_profiles.py --no-mask` builds without): 94 columns over
  seven genes, most in nad5 and nad2. On 29 DToL fern assemblies, against
  the angiosperm templates with stop read-through on: complete spliced
  mito genes 149 -> 231 of 248 (rpl2 0 -> 28, rps10 1 -> 29, rps3
  0 -> 26), genes aligned over 99% of their profile 80 -> 109 (62 before
  masking); planted 1-bp indels found 74% (77% on the angiosperm
  templates, different planted sets). nad7 still aligns end to end in
  only 8 of 29 (27 on the angiosperm template): its first 22 columns come
  from the Woodwardia and Blechnopsis annotations and fall just short of
  the masking rule. Horsetails and lycophytes do worse on these templates
  and stay on the default. Single-exon genes still use the angiosperm
  editing profiles.
- Planted-indel recall, re-measured after the exon-placement and
  reverse-editing work (`build/recall`, 54 gfatk assemblies left out):
  mito spliced genes 88.2% -> 88.0% of 768 with other calls in the same
  genes 143 -> 75; plastid spliced 74.5% -> 77.2% of 894; single-exon
  93.3% of 1,500, unchanged. With `--lineage fern` (stop read-through off
  -> on, angiosperm templates): fern mito spliced 45% -> 77%, fern mito
  single-exon 39% -> 46%.
- Exon placement in spliced genes, and gene ends. Checked on the 24
  RefSeq genomes the exon templates come from, each with its own genome
  held out of the templates, 89% of 477 exons now have both ends within
  3 bp of GenBank's (65% before; base-exact 68% against 45%), 15 of 16
  exons of 30 bp or less (none before), 142 of 169 last exons end exactly
  at the stop codon (73) and 166 of 178 first exons start exactly at the
  start codon (120). Three changes:
  - Exon-slot profiles (`00_build_exon_profiles.py`) score against the
    whole protein's composition with 2 pseudo-observations per column
    (`--slot-pseudocount`, now the default; `none` for the old way), and
    lose end columns on their intron side that most references lack
    (one reference annotating an intron boundary 7 codons out put every
    other species' ccmFc exon 2 that far into its intron). The old slot
    profiles were inverted for short exons - a residue the exon lacks
    scored as a match in every column. All shipped slot profiles and
    floors (mito, moss, plastid) are rebuilt; `--slots-only` rebuilds
    them without touching `whole_gene.pssm`.
  - transsplice (0.3.0, unreleased): every short exon slot gets a
    whole-contig search and a search beside its neighbouring exons, and
    the whole-gene alignment plus intron-end motifs choose between what
    they find and what filled the slot; the first exon is cut back to
    the start codon and the last ends at the stop codon.
  - orfedit (0.2.0, unreleased; `--stop-search`, default 10 codons): a
    corrected single-exon gene whose profile stops short of the protein's
    end is extended to the stop codon. On the same 24 genomes 89% of
    plastid and 80% of mito edited genes end exactly at GenBank's stop
    (17% and 67% before); 511 ends fixed, none made wrong.
  The whole-gene profiles of nad1, nad2, nad5 and rps3 are back to the
  RefSeq-only ones: the versions refined from DToL reconstructions
  (`02_bootstrap_refine.py`) had learned nad5's misplaced exon 3, and
  with them the whole-gene alignment preferred the wrong exon (like a
  reference exon in 569 of 972 DToL assemblies, against 880 with the
  RefSeq profile; none before). Their `whole_gene_score`s are lower for
  it (nad5 median 2,607 -> 1,663). Refine again only from reconstructions
  made with this version.
  Slot floors are 0.075 of the lowest reference self-score
  (`FLOOR_FRACTION`; half before): the new profiles are sharper, and at
  half, ferns and lycophytes, which run on the angiosperm templates, lost
  most of their exons (117 complete mito genes in 29 assemblies fell to
  27). Lower was better at every value tried down to 0.075, with held-out
  placement unchanged; the cost is a few more short exons placed that
  resemble no reference exon.
  On the DToL assemblies (transsplice alone, angiosperm templates,
  the 54 assemblies re-resolved with gfatk 0.6.3 left out):
  complete mito genes 5,957 -> 7,423 of 8,566 (cox2 +401, nad4 +302,
  ccmFc +275, rps3 +238, nad1 +205, rpl2 +113; rps10 -94, in angiosperms
  placements unlike any reference exon, in ferns because their rps10 has
  one exon), plastid 9,373 -> 9,868 of 9,954 (rps12 +303, clpP +105,
  rps16 +41); nad5 exon 3 like a reference exon in 861 of 927 (0 of 775),
  nad1 exon 4 928 of 928 (79 of 834), plastid rps12 exon 3 973 of 1,139
  (0 of 837), rps16 exon 1 812 of 852 (2 of 777). Genes with a frameshift
  call: mito 1,521 -> 762 (ccmFc 568 -> 268, nad1 215 -> 12, rpl2 292 ->
  171, cox2 118 -> 5; nad7 37 -> 65), plastid 362 -> 74 (rps12 136 -> 0,
  rps16 106 -> 23). On the 174 moss assemblies with the moss templates:
  1,555 of 1,566 genes complete (1,550), 24 with a frameshift call (165).
  `frameshifts_dtol.tsv`'s reconstruction
  rows are rebuilt from these runs. On Arabidopsis (smoke test): nad1
  0.936 -> 0.997 identical to RefSeq, nad5 0.990 -> 0.997, cox2 and rpl2
  reconstructed (were raw), plastid rps16 0.85 -> 0.99, rps12 0.94 ->
  1.00.
- Frameshifts: every protein-coding gene poa corrects is aligned allowing
  a change of reading frame, and each change is reported - single-exon
  genes through orfedit >= 0.2.0 (`tables/frameshifts.tsv`), spliced ones
  through transsplice >= 0.3.0 (`tables/reconstructed_frameshifts.tsv`,
  naming the exon) - with the skipped bases' coordinates, 1- or 2-base
  type, a `window` the change could equally well sit in (a homopolymer
  run, say), and `dtol_assemblies`, how many DToL assemblies have one at
  the same place in the same gene (`editing/reference/frameshifts_dtol.tsv`,
  built without the 54 assemblies resolved by gfatk before 0.6.3). Also a
  `frameshift` child of the gene's CDS in the GFF, and
  `n_frameshifts`/`n_new_frameshifts`/`frameshifted_genes` in the QC
  table, printed at the end of `poa run`. Frame-change penalty 20, no
  filters. On 1-bp indels planted in the coding sequence of real DToL
  genes: single-exon genes 92% found (96% plastid, 88% mito), 99% with the
  right type, and genes without a frameshift corrected exactly as with
  frame changes off; mito spliced genes 84% (652 of 778); plastid spliced
  genes 92% when the indel is in the copy reconstructed. Known weak
  point, in the README: calls in spliced genes are common and mostly
  recurrent (9% of DToL mito reconstructions; ccmFc 29%). With older orfedit/transsplice
  there are no frameshift calls (a warning says so);
  `--frameshift-penalty inf` turns them off.
- transsplice 0.3.0 aligns the joined exons in frame 0 only, which changes
  results with frame changes off too: on the 8,993 DToL mito spliced
  genes 528 differ from 0.2.1, nearly all nad5, 174 of which become
  complete and none stop being complete.
- `00_build_exon_profiles.py` takes `--exclude ACC`, `--out-dir` and
  `--slot-pseudocount` (see the exon-placement entry above).
- Plastid spliced genes (atpF, clpP, ndhA, ndhB, rpl2, rpoC1, rps12, rps16,
  ycf3) are reconstructed by transsplice like the mito ones, from exon
  templates built from 13 RefSeq plastomes across land plants
  (`trans_splicing/reference/plastid/`; `00a_fetch_references.py
  --organelle pltd`, `00_build_exon_profiles.py --organelle pltd`), and
  replace the raw multi-exon models in the plastid GFF when complete (73%
  of DToL plastid rps12, 87-100% for the rest). On Arabidopsis, all nine come
  out identical or nearly identical (0.99-1.00) to RefSeq.
  One copy per gene is reconstructed, so the second inverted-repeat copy
  of ndhB, rpl2 and rps12 isn't checked for frameshifts separately (in a
  graph assembly both copies come from one unitig).
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
- gfatk 0.6.2 (environment.yml, install.sh, recipe). `poa run` only uses
  `gfatk stats`, which is unchanged; 0.6.2 fixes `gfatk resolve` walking
  junctions the graph has no link for. `linearize/06_reannotate_promoted.py`
  reads each circuit's walk from its contig header (`path=`, gfatk >= 0.6),
  falling back to the resolve log for older outputs.
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
