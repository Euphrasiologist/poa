# trans_splicing

**The question:** `nad1`, `nad2`, `nad5`, and (occasionally) `rps3` are
trans-spliced in angiosperm mitochondria - their exons are transcribed as
separate RNA molecules from genomic loci that can be hundreds of kb apart,
sometimes on opposite strands, then spliced together at the RNA level.
Neither oatk nor anything else in this pipeline reconstructs that: oatk's
`.ctg.bed` calls surface these genes as multiple independent per-exon hits
under the same gene name, and `phylogeny` explicitly excludes them from its
marker-gene set as "multi-copy/fragmented" rather than mishandle them (see
`annotation/README.md`'s "Multi-copy / multi-exon genes" section and
`phylogeny/results/excluded_genes_mito.txt`). This module joins them, using
a purpose-built tool, [`transsplice`](https://github.com/tolkit/transsplice)
(a standalone Rust CLI that depends on
[`orfedit`](https://github.com/tolkit/orfedit) - see `../editing/README.md`
- as a library, reusing its edit-tolerant alignment engine rather than
duplicating it).

## What it computes

| script | output | what |
|---|---|---|
| `00_build_exon_profiles.py` | `reference/profiles/<gene>/exon{N}.pssm` + `.threshold`, `whole_gene.pssm` | per-exon-slot and whole-gene reference profiles, from 5 real GenBank mitochondrial genomes |
| `01_reconstruct.py` | `results/reconstructed_genes.tsv`, `results/reconstructed_exons.tsv`, `results/reconstructed_junctions.tsv` | runs `transsplice scan-batch` over every species' raw `.ctg.bed` fragments for these genes |
| `02_bootstrap_refine.py` | overwrites `whole_gene.pssm` | stage 2 - rebuilds the whole-gene profile from this dataset's own high-confidence reconstructions instead of just the 5 GenBank references |

## Reference data: real GenBank genomes, independently verified

`reference/raw/` holds exon coordinates/sequences from 5 complete, curated
RefSeq mitochondrial genomes (Arabidopsis NC_037304.1, Beta vulgaris
NC_002511.2, maize NC_007982.1, rice NC_011033.1, Vitis NC_012119.1),
fetched directly (NCBI efetch) and parsed from their `CDS join()` features -
not literature paraphrase. **Independently re-verified**: every exon was
re-translated from its own raw coordinates and checked against the
genome's own curated RefSeq protein; all 23 files matched except one
flagged, confirmed-real exception (rice `rps3`, ~16% identity past residue
34 - likely an assembly artifact in that specific older RefSeq record, not
a transcription error on our end) - see `reference/raw/junctions.tsv`'s
`notes` column, marked `UNRESOLVED` and excluded from training.

**What this data changed about the design** (worth remembering when
reading the algorithm below): `nad2` and `nad5` have a *perfectly
conserved* 5-exon topology across all 5 genomes (`nad2`: intron 2 trans,
1/3/4 cis; `nad5`: introns 2 **and** 3 trans, reconstituting a **22-nt**
central exon - Wissinger et al. 1991, EMBO J). `nad1`'s exon *count* is
always 5, but *which* junctions are cis vs trans genuinely varies by
lineage (confirmed via hard strand-switches, not just large gaps - e.g.
Beta vulgaris trans at junctions 3/4 where Arabidopsis/Vitis are trans at
1/3). `rps3` is mostly a plain 2-exon cis-spliced gene - Beta vulgaris has
lost the intron entirely (single exon); only maize confirms a genuinely
derived trans-spliced copy (`ZeamMp001`, exons ~97kb apart) sitting
alongside an ordinary cis paralog (`ZeamMp024`, ~1.8kb intron) in the same
genome. **Consequence**: the algorithm never assumes a fixed cis/trans
*pattern* per junction - only fixed exon *identity/order*, which is
structural. Cis vs trans is a per-species, per-junction *output label*
(distance + strand check), never an input assumption.

## Method

1. **`00_build_exon_profiles.py`** - per gene, the *canonical* exon count K
   is the modal `n_exons` among included reference loci (5 for
   `nad1`/`nad2`/`nad5`, 2 for `rps3`). For each slot 1..K, only loci whose
   *own* exon count matches K contribute training sequences - this excludes
   Beta vulgaris' intron-less `rps3` from `exon1`/`exon2` (it has neither -
   its single CDS is stored as `exon1of1`, a different topology, not a
   35 vs. real ~24aa leader exon) while still using its full protein for
   the whole-gene profile. Each `.threshold` (`min_self_score`) is not an
   invented cutoff - it's computed by running the slot's own training
   sequences back through a real `orfedit scan-batch` call against the
   profile just built from them, and taking the minimum genuine score.
2. **`01_reconstruct.py`** builds one manifest row per raw `.ctg.bed`
   fragment (not best-hit-only, unlike `editing`) and calls
   `transsplice scan-batch`, which per (species, gene):
   - **Merges** same-strand fragments within 5kb into an additional
     candidate spanning the whole run (oatk's HMM sometimes splits one
     true exon into adjacent hits - confirmed directly in real data, e.g.
     `nad2` pairs ~1.2-1.4kb apart in `Acer_campestre`) - both the pieces
     and the merged span go into the pool; scoring decides which is real.
   - **Scores** every candidate against every exon-slot profile via
     `orfedit::scan_gene` reused as-is (so scoring *is* a boundary-corrected,
     edit-tolerant alignment, not a raw heuristic).
   - **De novo fallback**: if a slot has no candidate above its floor
     (the routine case for `nad5` - see caveat below), scans the same
     contig(s) the other exons were found on, both strands, for that one
     slot's profile - bounded to gene-known contigs, not a genome-wide scan.
   - **Assigns** candidates to slots via best-scoring injective partial
     matching (brute force - K≤5, pool capped, trivially fast).
   - **Concatenates** assigned exons in slot order, each independently
     strand-corrected, classifies each junction cis/trans for *this*
     species, and runs one final whole-gene alignment pass for an
     aggregate score/protein.
3. **`02_bootstrap_refine.py`** - once `01_reconstruct.py` has produced
   real output, rebuilds each gene's whole-gene profile from this
   dataset's own `complete=True, whole_gene_score>0` reconstructions
   (gap characters stripped) instead of just the 5 GenBank references, and
   re-run `01_reconstruct.py` to benefit. Same two-stage bootstrap pattern
   as `editing/README.md`'s.

## A real bug this caught, worth knowing about

Validating on a 14-species smoke test surfaced a genuine correctness bug
in `orfedit`'s core (shared by both tools, fixed there, not worked around
here): the alignment DP could "align" a stop-codon translation to a real
profile column at a small finite penalty when every other reading scored
even worse, producing literal `*` characters inside output protein
sequences instead of treating that position as an unexplained deletion.
Fixed by scoring stop/unknown as `NEG_INFINITY` in `orfedit`'s
`Profile::score` (a real coding column can never truly be a stop), with a
regression test (`orfedit/src/align.rs`). Confirmed this materially
improved reconstruction quality, not just cosmetics: e.g.
`Potamogeton_compressus nad5` went from a 4/5-filled, negative-scoring
partial reconstruction to a fully complete (5/5), strongly positive one
after the fix.

## Known limitation: `nad5` recovery rate

`nad5`'s modal fragment count in this dataset is **1** - oatk's own gene
model mostly only ever calls the single largest exon; the other 4 (one
only 22nt!) usually aren't independently called *at all*, not just
mis-joined. The de novo fallback exists specifically for this, but its
floor for the 22nt micro-exon (`exon3`) is unusually weak
(`min_self_score` came out slightly *negative* - a real, small-profile
consequence of that exon being only ~7 residues, not a bug) - treat any
`nad5` reconstruction that only completed via that one slot's de novo hit
as lower-confidence than one built from real `.ctg.bed` calls. On the
smoke test, `nad2`/`nad1` completion looked more constrained by genuinely
missing raw fragments than by assignment errors; `nad5` completion should
be watched at full scale for the opposite failure mode - the fallback
firing too permissively - rather than assumed to be a clean win if it
looks unexpectedly high.

## Ground-truth validation found a real, precisely-localized bug

`Arabidopsis_thaliana` is in this dataset as a real oatk-assembled species
(not the NCBI reference - a different physical PacBio-sequenced genome of
the same species), which makes it a genuine ground-truth check: its
reconstructed proteins can be compared directly against the curated
RefSeq protein sequences already fetched for `nad1`/`nad2`/`nad5`/`rps3`
(a proper global alignment, BLOSUM62 - a naive positional/`difflib`
comparison is *not* valid here and will make correct reconstructions look
wrong, since a single real indel shifts every downstream position).

**Result, gene identification**: 54/55 real `/gene` features in the
curated GenBank record were found (the one apparent miss, `trnF`, is a
parsing artifact - GenBank names it without an anticodon suffix; we
correctly call it `trnF-GAA`). `editing`'s single-exon corrected proteins
(`atp4`/`ccmB`/`ccmC`/`mttB`/`nad9`) are essentially perfect against the
real curated proteins: **99.5-100% identity** - `orfedit`'s core
edit-tolerant alignment engine is validated, not just plausible.

**Result, trans-spliced gene reconstruction - a real, precisely-localized
problem, not previously known at this resolution**: identity is excellent
specifically for *large, well-powered* exon slots and poor for *small*
ones, in every gene checked - `nad1`'s big first exon (~128aa, ~40% of the
protein) aligns near-perfectly, then everything after it (three much
shorter exons) falls apart; `nad2`'s first two exons (~180aa combined)
align near-perfectly, exons 3-5 (all short) don't; `rps3` (a small 24aa
leader + one long main exon) is poor throughout, consistent with its
already-documented weak whole-gene profile, not a new issue. This is not
the same thing as the `nad5` micro-exon caveat above (that's about
*finding* a slot at all); this is about slots that *are* found and
assigned a plausible-looking region, but the wrong one - almost certainly
because a short exon-slot profile has little power to reject a
locally-plausible wrong frame/region in favour of the true one when
scored via `scan_gene` against a whole candidate window.

**Fixed** (`transsplice/src/assign.rs`'s `best_assignments_beam` +
`lib.rs`'s `reconstruct_gene`): per-slot-independent scoring alone can't
tell a permutation error from a correct assignment when several real
candidates score comparably across slots - confirmed directly (`nad2`,
above: all 5 raw fragments had decent independent nhmmscan scores, so
this wasn't a weak-evidence problem, it was a genuine "right fragments,
wrong slots" error). Fixed by generating a small beam of plausible
assignments (each slot's top-2 candidates, all valid injective
combinations - bounded, ~32 at most for K<=5) and picking whichever one's
*whole-gene* alignment score is best, not just the locally-greedy
per-slot sum - a new regression test
(`whole_gene_score_resolves_a_real_permutation_ambiguity`) mirrors this
exact scenario synthetically.

**Result on the same `Arabidopsis_thaliana` ground truth, honestly mixed,
not uniformly fixed**: `nad5` improved sharply (46%->82% identity) -
the fix working as designed. `nad1` is unchanged (~80%) - expected,
since its problem (slot 4's raw evidence was genuinely weak, bitscore
27-48) was never a permutation error to begin with, just correctly
diagnosed as a different, harder failure mode. `nad2` - the specific
case this fix was built around - is **unchanged** (still 60.3%, identical
exon coordinates to before the fix). Root cause: the beam-search fix only
helps when its arbiter (the whole-gene profile) has enough power to
prefer the correct permutation over the wrong one; `nad2`'s whole-gene
profile, like the per-slot ones, is still built from only 5-7 GenBank
sequences and apparently isn't discriminating in the exons 3-5 region
either. This was predicted to be fixable by `02_bootstrap_refine.py`
(stage 2 of the profile bootstrap) once enough dataset-derived
reconstructions existed - **checked directly, and it wasn't fixed, for a
real and important reason.**

**`02_bootstrap_refine.py` run on this 15-species set**: `nad2` cleared
its >=5-complete-reconstructions gate for the first time (exactly 5) and
got a refined `whole_gene.pssm`; `nad5` too (6). Re-running
`01_reconstruct.py` afterward: `nad2`'s `whole_gene_score` jumped sharply
(278 -> 946), but **identity against the real curated protein did not
improve** (57.5% vs 60.3% before - if anything slightly worse), and the
exon coordinates for `Arabidopsis_thaliana` are byte-identical to the
pre-bootstrap run. `nad5` did improve further (82% -> 85%).

**Why `nad2` didn't improve - a real methodological gap in
`02_bootstrap_refine.py`, not a dead end**: `Arabidopsis_thaliana`'s own
(partially wrong) reconstruction was one of the only 5 sequences used to
build `nad2`'s refined profile. The profile then learned to *expect*
exactly that same wrong pattern at those positions, so re-scoring
`Arabidopsis_thaliana` against it just confirms its own prior error
rather than being corrected by independent evidence - the much higher
score is an artifact of that circularity, not evidence of correctness.
This is a leave-one-out problem: `02_bootstrap_refine.py` should exclude
each species from its own refined profile (or, more practically at full
scale, simply require a much larger and more diverse training pool than
5 sequences before trusting the refinement - with hundreds of species
each contributing a small fraction of the training set, any single
species' self-influence shrinks to negligible, which is exactly why this
is expected to behave better once run on the full ~1250-species dataset
rather than this 15-species smoke test). Re-check `nad2` again once a
full-scale run gives the bootstrap a properly large, diverse training
pool - this 5-sequence result doesn't yet test the mechanism this was
designed for.

## Reading the output

`reconstructed_genes.tsv`: one row per species x gene - `n_filled`/`n_slots`,
`complete`, aggregate `whole_gene_score`/`whole_gene_n_edits`, and the
joined protein (`-` = an unfilled/deleted position). `reconstructed_exons.tsv`:
one row per assigned exon (slot, corrected coordinates, strand, its own
score, `via_denovo` flag, raw-fragment provenance - `merged:...` for a
run of merged fragments). `reconstructed_junctions.tsv`: one row per
junction between consecutive assigned exons (`cis`/`trans`, gap in bp) -
this is the actual answer to "is this really trans-spliced in this
species", derived fresh per species, not inherited from the reference
topology.

## Not yet run at full-dataset scale / not in `run_all.sh`

New, unvalidated at scale - same opt-in precedent as `orf_scan`/`editing`.
Run `00_build_exon_profiles.py` once, then `01_reconstruct.py
--species-list <subset>` and sanity-check against this README's smoke-test
numbers before a full run; `02_bootstrap_refine.py` is only useful after a
first full-ish pass has produced enough `complete`/positive-score
reconstructions per gene (5+) to be worth rebuilding from.

## Try this

Feed successfully reconstructed (`complete=True`, decent `whole_gene_score`)
genes back into `phylogeny`'s marker-gene set - that module currently
excludes these genes entirely rather than mishandle them; this is what
closes that gap. Not done in this phase.

## Cis-spliced genes (stage 2a, 2026-09-29)

`transsplice` never assumed a junction is trans - it classifies each one
per species after assignment - so cis-spliced mito genes are just more
templates. `src/reference_genes.py` is the single gene list for the
fetch, template and reconstruction scripts and `gff_export`: nad1/nad2/
nad5/rps3 plus nad4, nad7, ccmFc, cox2, rpl2, rps10.

References: `src/00a_fetch_references.py` fetches 11 angiosperm RefSeq
mitogenomes (the original 5 plus Nicotiana, Glycine, Daucus,
Liriodendron, Phoenix, Silene latifolia; GenBank cached in
`reference/genbank/`) and writes each cis gene's exons in the existing
`reference/raw/` format. Each locus is verified: the genomic translation
must equal RefSeq's curated protein except at C-to-U-editable residues
(an edit-created stop counts). 52/54 verify; rice nad7 and rice rpl2 are
UNRESOLVED and excluded. Build only these genes' templates with
`00_build_exon_profiles.py --genes ...` - rebuilding a trans gene resets
the whole_gene.pssm 02_bootstrap_refine.py refined from this dataset.

Exon count varies between references for some genes (cox2: 1/2/3 exons;
nad4 3 in Beta and Silene, 4 elsewhere; nad7 4 in Nicotiana); templates
use the modal count, so a species with an extra intron is under-segmented
by its template - gff_export then keeps the raw chain (more exons).

Arabidopsis vs RefSeq with transsplice 0.2.1, end to end in the GFF:
nad4 4/4 and nad7 5/5 exons base-exact (94% protein identity = its
editing sites), nad2 5/5, rps3 2/2; cox2 falls back to the raw chain
(its 83 bp exon 2 is below what a profile places reliably) and rpl2 is
rejected by gff_export's consistency guard (Brassicaceae rpl2 is
truncated; a full-length template aligned 573 bp past its end).

Full dataset (job 709710, transsplice 0.2.1), complete reconstructions:
nad1 76%, nad2 75%, nad4 66%, nad5 55%, nad7 68%, ccmFc 69%, cox2 57%,
rpl2 59%, rps10 81%, rps3 71% of species with the gene. gff_export
rejected 868 cis reconstructions as inconsistent with the HMM hits:
mostly slots placed de novo where no hit exists (rpl2, rps10 - often
lost or nuclear-transferred) and cox2 exon 1 running ~340 bp 5' of the
hit in ~220 mosses (Sphagnum, Orthotrichaceae, Hypnales), whose cox2
diverges from the angiosperm template - stage 2b's lineage templates
would be needed there. In every case the raw chain is kept.
