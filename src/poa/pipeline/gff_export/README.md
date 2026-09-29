# gff_export

**The question:** `denovo_annotation` produces one combined GFF per
species, but it's flat - isolated single-line features, no `gene`/`mRNA`/
`CDS` hierarchy, and it doesn't reflect anything `editing` or
`trans_splicing` found (their results live in separate TSVs). There was
no single file you could open and see "here is `nad2`, here are its 5
exons, here is where it's trans-spliced, here is where it's RNA-edited".
This module builds that file: one standardized, GFF3-spec-correct,
hierarchical annotation per species, merging all three modules by
precedence.

## Gene models: every hit, every copy, every exon

Until 2026-09-28 this kept only the best-scoring hit per gene, which
silently dropped 190,077 of 348,663 hits: every second copy (both
inverted-repeat copies of plastid rrn16/rrn23/ndhB/..., mito repeat
copies) and every exon but one of every cis-spliced gene (mito nad4/nad7/
ccmFc/cox2/rpl2, plastid ndhA/ndhB/ycf3/clpP/petB/rpoC1/...). Models are
now built from all hits, using where each hit lies along its gene's HMM
(`hmm_from`/`hmm_to`, added to `gene_calls.tsv` from nhmmscan's tblout):

- **Collapse:** a gene's overlapping hits become one locus; opposite-strand
  members (weak antisense rRNA/tRNA model hits) lose to the strongest
  hit's strand. tRNA loci are first resolved *across* names - one real
  trnK is also hit by the trnI/trnM/trnN/trnQ/trnR/trnT/trnV/trnW models -
  keeping tRNAscan-SE's anticodon-based identity (oatkDB's for CAU, which
  tRNAscan-SE always calls Met: trnfM and trnI-CAU are distinct genes).
- **Chain:** protein-coding and tRNA loci on one strand, 250 bp - 8 kb
  (mito) / 4 kb (plastid) apart, each picking up the gene model where the
  previous one left off, are exons of one gene. A second copy restarts the
  model, so it never extends a chain. A gap under 250 bp is not an intron
  (organellar group I/II introns are longer) but one exon whose HMM hit
  broke over a divergent stretch - the pieces merge (this was ~800 bogus
  "introns" in ycf1/ycf2/rpoC2 before 2026-09-28's fix). Chains in
  bryophyte atp9/atp1/cox1/cox3/cob/nad9/sdh3 are expected: moss and
  liverwort mitogenomes carry introns there that angiosperms lack.
- **Copies:** everything left is a copy - `gene-X`, `gene-X-2`, ... -
  strongest first. `model_coverage` is measured against the model span
  the gene's hits usually cover across the dataset (oatkDB models often
  run well past the CDS - rps12's is 1119 nt for a ~380 bp gene); under
  0.8 is `partial=true`, and `copy_number` counts the rest. An extra tRNA
  copy is kept only if tRNAscan-SE (score >= 35) and oatkDB agree, or it
  is a near-complete intron-split chain (tRNAscan-SE misses those - the
  second IR copies of plastid trnA-UGC/trnI-GAU look exactly like that).
  The strongest locus of every gene is always kept; `tool_support` says
  which predictors called each tRNA/rRNA.

Checked against references: Arabidopsis mito (GenBank NC_037304) -
every protein-coding gene's copy count and exon count match (atp6 x2,
nad4 4 exons, nad7 5, ccmFc/cox2/rpl2/rps3 2) apart from nad1 (the
trans_splicing result fills 4 of 5 slots), rpl2 (tagged partial - in
Brassicaceae the mito gene encodes only the N-terminal part) and loci the
reference leaves unannotated (rps14/rps19/sdh4 pseudogenes, atpA/atp1
naming); all 22 reference tRNAs recovered with the right copy numbers.
Arabidopsis plastid - every IR gene x2, ndhA/ndhB/rpl2/rps16/rpoC1/atpF/
trnK/trnL/trnV/trnA/trnI 2 exons, ycf3/clpP 3; petB/petD/rpl16/trnG-UCC
come out one exon short because their first exon (6-25 bp) is too short
for any HMM hit.

Full dataset: 201,364 gene models (34,641 partial), 16,815 raw cis-spliced
multi-exon models; 2,638 unsupported extra tRNA loci dropped.

## Mito genomes whose linearisation is uncertain

Some mito assembly graphs carry minor configurations: low-depth side
paths (well under the ~2-fold abundance range of real mitochondrial
chromosomes - Wu et al. 2015, PNAS 112:10185) embedded in the graph.
They may be low-abundance recombination products, nuclear copies of
mitochondrial DNA, or artefacts - not yet told apart (that needs
read-level evidence: junction-spanning reads, what flanks the side-path
reads). Any single contig sequence is then one reading of the graph, and
a gfatk circuit can walk the real genome several times through them.

For the species `qc_basic_stats/src/06_low_depth_paths.py` flags
(`results/low_depth_paths.tsv`), both GFFs start with a
`# linearisation uncertain` comment and the contig `region` rows carry
`linearisation=uncertain;linearisation_reason=low_depth_side_paths;
graph_low_depth_frac=...`. Prefer the unitig-level GFF for these: it
annotates the graph's own sequences and doesn't depend on a path choice.
Gene-free sequence alone is never grounds for the flag - Silene-type
multichromosomal mitogenomes carry whole gene-free chromosomes at normal
abundance.

## Known issue: trans_splicing exon junctions (being fixed)

`trans_splicing`'s reconstructions are not codon-exact at junctions where
a codon is split between two exons (junction phase 1 or 2 - e.g. nad1
exons 2/5, nad2 exons 3-5, nad5 exons 2/4, rps3 exon 2). Its template
builder translates every exon from frame 0, and exons are joined on codon
boundaries, so the split codon's 1-2 bases are dropped and the joined CDS
is out of frame downstream. Measured on Arabidopsis (a training genome):
nad2's five exons are each placed within 2 bp, yet the protein matches
RefSeq at only 50% of residues (nad5 80%, nad1 61%); tiny exons (nad5
exon 3, 22 bp; nad1 exon 4, 59 bp) can land on the wrong locus. So
`annotation_tier=reconstructed` means right exons, approximately right
boundaries - not exact ones. Phase-aware templates and junction
refinement in `transsplice` are the fix in progress.

## Precedence: which module's correction wins

Corrections claim the loci they were computed on:

1. **`trans_splicing`** - the genes in `trans_splicing/src/reference_genes.py`
   (trans-spliced nad1/nad2/nad5/rps3; since 2026-09-29 also cis-spliced
   nad4/nad7/ccmFc/cox2/rpl2/rps10) with >=1 filled slot, unless a raw
   chain of the same loci has more exons (the template under-segments
   that species, e.g. an extra intron). For the cis-spliced genes a
   reconstruction must also agree with the HMM evidence - every exon
   overlaps a raw hit for the gene and extends at most 100 bp past it -
   or the raw chain is kept (Brassicaceae rpl2 is truncated, and a
   full-length template aligned 573 bp past its end). Either fallback
   leaves the gene `raw`: no edit sites, approximate boundaries.
2. **`editing`** - on the single-exon copy holding the hit it corrected.
   An ORF search over one exon of a cis-spliced gene can't correct that
   gene's boundaries, so multi-exon copies stay raw (716 cases).
3. **`denovo_annotation`** - everything else (and all tRNA/rRNA).

Recorded per gene as `annotation_tier=reconstructed|edited|raw` - an
explicit confidence signal, not implicit in which script you happened to
run last. Raw multi-exon models take CDS phase from the exon lengths
before each exon and say `exon_boundaries=approximate` on the mRNA - HMM
envelope bounds aren't splice sites.

## How GFF3 actually encodes this (confirmed, not assumed)

**Multi-exon / trans-spliced genes**: the spec's real mechanism is
multiple `CDS` rows sharing one `ID`, each `Parent=`-linked to one `mRNA`.
Nothing in the spec requires those rows to share a strand or even a
seqid - which is exactly what a trans-spliced gene needs, since a real
exon can be on the opposite strand from its neighbours (confirmed
directly in the 5 GenBank reference genomes `trans_splicing` fetched:
their `CDS` lines use nested `join(complement(a..b),c..d,...)` location
strings for precisely this reason - GFF3's equivalent is just multiple
same-`ID` rows instead of one nested string). One real consequence worth
knowing, not a bug: the parent `gene`/`mRNA` feature's own start/end must
span every exon's extremes, so a trans-spliced gene's `gene` row can
cover a huge region (e.g. `nad2` below spans 233kb) that's mostly
unrelated intervening sequence - inherent to representing trans-splicing
in GFF3 at all, and its strand is written as `.` (mixed/not applicable)
rather than a misleading single value.

**RNA editing sites**: no single certain, widely-adopted SO term was
available with full confidence, so this uses `sequence_alteration` (a
real, general term - "this position differs from a reference") as a
pragmatic, documented choice, not a claim that it's the one definitive
standard.

**Phase column** (GFF3 field 8, required for `CDS`): `0` for
`editing`/`trans_splicing`-derived exons. For `editing`'s single-exon
calls that is a verified guarantee (orfedit's DP only emits codon-aligned
boundaries); for `trans_splicing` it is wrong wherever a codon is split
across a junction - see "Known issue" above. Raw single-exon calls also
get `0`, and raw cis-spliced models get each exon's phase from the
lengths of the exons before it; both are best-effort, not confirmed -
the coordinates are HMMER envelope bounds with no codon-boundary
enforcement (multi-exon ones say `exon_boundaries=approximate`).

## A real coordinate bug this caught

Building this surfaced an off-by-one between the two `gene_calls.tsv`-
shaped tables: `annotation/`'s (oatk's own `.ctg.bed`) is 0-based
half-open (BED-style), but `denovo_annotation/`'s was copying coordinates
straight from `nhmmscan`/`tRNAscan-SE`/`barrnap`'s native GFF3 output
(1-based inclusive) with no conversion - a silent 1bp inconsistency
between two tables documented as interchangeable. Fixed at the source
(`denovo_annotation/src/04_build_gene_calls.py`, `start - 1` at the one
place raw GFF coordinates enter that table) and regenerated. Given the
~60bp flank used everywhere downstream, this almost certainly didn't
change any of `denovo_annotation`'s already-reported comparison numbers,
but it needed fixing before coordinates here could be trusted at bp
resolution.

## Worked example: `Acer_campestre`, real output

```
gene   382417 383028 -  annotation_tier=edited        # ccmB
mRNA   382417 383028 -  n_exons=1;n_edits=11
CDS    382417 383028 -  ID=cds-ccmB
sequence_alteration  382457 ... (x11, one per real predicted edit site)

gene   442932 676371 .  annotation_tier=reconstructed # nad2
mRNA   442932 676371 .  n_exons=5;n_edits=8
CDS  442932-443090 +  exon_number=3
CDS  670221-670610 -  exon_number=2
CDS  671813-671965 -  exon_number=1
CDS  674100-674672 +  exon_number=4
CDS  676123-676371 +  exon_number=5

gene  rrn26  216102-219285 -  annotation_tier=raw
rRNA  216102-219285 -
```

`ccmB`'s 11 `sequence_alteration` children match exactly the 11 real
edits already documented for this species in `editing/README.md` - a
useful cross-check that the merge is reading the right rows, not just
that it runs without error. `nad2`'s exon strands (1,2 on `-`; 3,4,5 on
`+`) match `trans_splicing`'s own `reconstructed_junctions.tsv` for this
species exactly (cis 1->2, trans 2->3, trans 3->4, cis 4->5).

Verified structurally sound, not just eyeballed: every `Parent=` in every
one of the 14 smoke-test files resolves to a real `ID=` in the same file
(0 orphan references).

## Assembly provenance and unitig coordinates

Every ctg gets a `##sequence-region` pragma and a `region` feature
recording which linearisation produced it - `resolver=oatk_pathfinder`
or `resolver=gfatk_resolve` (the `linearize` fallback promoted into
`data/` by `05_promote_resolve.py`; read straight from the ctg.fasta
header) - plus GFF3's own `Is_circular=true` where applicable.

Using `unitig_coords`' placement map, every feature also records where it
sits on the raw assembly-graph unitigs:

```
ctg000001c ... CDS 442932 443090 + ... exon_number=3;unitig_loc=u72:1107-1265:-
ctg000001c ... gene 105511 107090 + ... unitig_span=u309,u1593
```

`unitig_loc` (1-based inclusive, strand relative to the unitig) when the
feature lies inside one unitig placement; `unitig_span` listing the
unitigs it crosses when it spans a junction between placements - flagged,
not guessed (`unitig_span=unplaced` if it sits on ctg sequence no
verified placement covers). A position inside a link overlap belongs to
two placements; the one where it's core is preferred.

`results/unitig/<species>.<organelle>.unitig.gff` is the same annotation
on the unitig sequences themselves (`unitig_coords/results/unitig_fasta/`).
A junction-crossing feature becomes several rows sharing one `ID` (GFF3's
discontinuous-feature mechanism, as for multi-exon CDS), each tagged
`unitig_split=<ctg>:<start>-<end>` with the ctg stretch it came from, and
CDS phase recomputed per piece. This view collapses repeat copies onto
one sequence - which is why the ctg-level GFF stays primary - so each
unitig's `region` row carries `n_ctg_copies` to keep that visible.

For `gfatk_resolve` contigs there is no recorded unitig walk, so
placements come from exact sequence matching; a match found only after
trimming the unitig's ends (typically one of several near-identical
repeat variants) claims only its verified middle.

Verified on a 7-species test set spanning both resolvers, linear and
circular contigs, cross-contig trans-splicing and an origin-wrapping
single-unitig genome: every `unitig_loc` and every split piece is
byte-identical between ctg and unitig sequence (bar SNP/1bp-indel
differences confined to link overlaps, where the ctg carries one
neighbour's copy), every split feature's pieces sum to its full length,
and every `Parent=` resolves in both files.

## Two merge bugs fixed (2026-09-28)

- **Cross-organelle leakage.** Editing calls, edits and
  `trans_splicing` reconstructions were matched on species + gene only.
  Genes present in both genomes (`atpA`, `rpl16`, `rps19`, `rps4`,
  `rps14`, ...) could pick up the *other* organelle's call - 3911 edited
  gene models in 673 species had the other genome's coordinates and edit
  sites, and plastid `rps3` could take the mito trans-splicing
  reconstruction. Every table is now split by species AND organelle
  before lookup (with organelle matched, all 54278 editing calls
  correspond exactly to their raw best hit).
- **Cross-contig trans-spliced genes.** 258 of 3558 reconstructions have
  exons on more than one contig; the `gene`/`mRNA` row used the first
  exon's contig but spanned every exon's coordinates, running off the
  end of the contig. These now get one `gene`/`mRNA` row per contig
  (shared `ID`), and each edit site sits on its own exon's contig.

## Running

Stage 4 of `denovo_annotation/work/run_full_suite.sh`; not in
`run_all.sh`. Depends on `denovo_annotation`'s `gene_calls.tsv` and
benefits from `editing`/`trans_splicing` having been run first, but
degrades gracefully to tier-3-only output if they haven't (every
upstream table is loaded with its expected columns even when the file
doesn't exist yet). The unitig layer needs `unitig_coords` run first and
is skipped per species when its map is absent.
