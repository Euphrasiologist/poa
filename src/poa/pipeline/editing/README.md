# editing

**The question:** plant mitochondrial (and, more sparingly, plastid) mRNAs
undergo extensive post-transcriptional C-to-U RNA editing - including,
routinely, at start and stop codons. Neither `annotation` (a reshape of
oatk's own DNA-level HMM gene calls) nor `orf_scan` (plain `ORFfinder`,
DNA-level, standard code) has any way to see this: a genomic sense codon
that only becomes a stop after editing (`CGA`/`CAA`/`CAG` ->
`UGA`/`UAA`/`UAG`) reads straight through at the DNA level, so a DNA-only
caller over-extends past the true C-terminus; conversely a start created
only by editing (`ACG` -> `AUG`) can leave the true N-terminus uncalled.
This module corrects for that, using a purpose-built tool,
[`orfedit`](https://github.com/tolkit/orfedit) (a standalone Rust CLI,
`cargo install`ed like `gfatk` - see `analysis/common/tool_paths.sh`), not
Python.

If you're coming from `orf_scan`: that module finds *extra*, non-core
content (TE/mitovirus-like ORFs) that oatk's targeted gene search doesn't
report at all. This module does the opposite - it doesn't find new genes,
it corrects the boundaries of ones oatk already found.

## What it computes

| script | output | what |
|---|---|---|
| `00_build_reference_profiles.py` | `results/profiles/<gene>.pssm` | per-gene amino acid profile (PSSM), built from this dataset's own already-annotated orthologs |
| `01_scan_editing.py` | `results/edited_gene_calls.tsv`, `results/predicted_edits.tsv`, `results/lineage.tsv` | runs `orfedit scan-batch` over every (species, gene) anchor that has a profile |
| `00b_build_lineage_profiles.py` | `reference/lineages/<lineage>/mito/<gene>.pssm`, `reference/profiles_refseq/mito/<gene>.pssm` | bundled profiles from RefSeq: per non-angiosperm lineage (see "Moss profiles" below), and angiosperm profiles for the mito genes the dataset doesn't profile (see "RefSeq profiles") |

The scripts in `src/` are command-line wrappers. The code is in importable
modules - `poa.pipeline.editing.scan`, `.profiles` and `.lineage_profiles` -
whose functions take the data root and every option as arguments and raise
`StageError` (`poa.pipeline.common.errors`) on failure:

```python
from poa.pipeline.editing.scan import scan

result = scan(root, ["mito"], threads=8)  # -> ScanResult (output paths, species -> lineage)
```

`poa run` calls `scan()` in-process. The profile maths (mafft, log-odds
columns) is `common/pssm.py`'s, shared with `trans_splicing`.

## Method

1. **`00_build_reference_profiles.py`** - for each marker gene that already
   has a FASTA in `annotation/results/genes/{mito,pltd}/` (restricted to
   protein-coding genes - the marker-gene set also includes `trn*`/`rrn*`
   non-coding genes, which are skipped here, not translated), translates
   every species' sequence, drops anything short or with an internal stop
   (a very low bar - see the script; it's a training-set sanity filter, not
   a correction step), aligns the survivors with `mafft --auto`, and writes
   a per-column amino-acid frequency profile (log-odds vs. this alignment's
   own background frequency, columns gappier than 90% dropped).

   90%, not a more obvious 50%: verified empirically on a 14-species smoke
   test that a 50% cutoff systematically trimmed almost every profile by a
   few N-/C-terminal columns and made ~98% of *zero-edit* calls look
   "changed" by a few bp - not real editing, just the training alignment
   itself containing plenty of boundary-truncated raw calls (the exact
   problem this tool exists to fix), so a naive per-column gap-majority
   filter was throwing away genuine terminal signal, not just rare
   lineage-specific insertions.

   **Caveat worth remembering** (documented in the script, repeating here
   because it matters for interpreting results): this is trained on the
   dataset's own existing calls, so if the *majority* of species share the
   same boundary truncation, the profile bakes that in as "correct" too.
   This approach finds and fixes genes truncated in the *minority* of
   species relative to the dataset's own majority pattern - it is not a
   substitute for an externally-validated editing-site atlas (PREP-Mt-
   style), and results should be spot-checked against literature-known
   cases before being trusted wholesale, especially early on.

   One specific instance of this worth flagging: `nad5` is canonically
   trans-spliced in angiosperm mitochondria, yet it clears this dataset's
   80% single-copy marker-gene threshold (`phylogeny/README.md`:
   88.8% single-copy) and so gets a profile here. That almost certainly
   means oatk's `nad5` HMM hit only ever targets one (the largest) exon,
   silently, not the full trans-spliced gene - so this module's "corrected"
   `nad5` calls are edit-corrected boundaries of that one exon, not a
   full-length reconstruction. Full-length `nad5`/`nad1`/`nad2`/`rps3`
   reconstruction is a separate, not-yet-built module (see below).

2. **`01_scan_editing.py`** - joins `annotation/results/gene_calls.tsv`
   (already QC-gated) against the profiles above (one manifest row per
   species x gene, best-scoring `.ctg.bed` hit only), and makes one call to
   `orfedit scan-batch`. `orfedit` extracts a flanked genomic window around
   each existing call, and for each of the 3 reading frames runs a
   semi-global affine-gap alignment of the profile against the window's
   codons - at every codon, scored as the *best* of a literal translation or
   any subset of its C positions hypothesized edited to T (a small penalty
   per hypothesized edit, so it isn't free to invent edits). It also checks
   the codon immediately after the aligned protein for a genuine (possibly
   edit-created) stop and folds it into the corrected span if found -
   without this, `changed` was true on almost every call regardless of
   real editing, since the profile (built from stop-stripped protein
   translations) is one codon short of `.ctg.bed`'s own stop-inclusive
   convention by construction, not because of anything biological. See
   `orfedit`'s own README/`src/align.rs`/`src/lib.rs` for the algorithm and
   its unit tests (both of these issues are covered by regression tests
   there, not just fixed ad hoc).

## Reading the output

`edited_gene_calls.tsv`: one row per species x gene, with both the raw
(`raw_start`/`raw_end`) and corrected (`corrected_start`/`corrected_end`)
coordinates, `n_edits`, `includes_stop`, alignment `score`, and the
edit-corrected protein sequence. `predicted_edits.tsv`: one row per
predicted edit site (genomic position, which profile/protein column it
falls in - `n_columns_total` itself, i.e. one past the last real column,
for a predicted edit-created stop - and the resulting amino acid).

**`changed` means "editing was invoked" (`n_edits > 0`), not "coordinates
moved at all".** Coordinate-equality was tried first and dropped: a
profile-alignment boundary almost never lands on the exact same bp as an
independently-called boundary even when nothing biological changed, so
`changed` needs to mean something an alignment-noise nudge can't trigger -
`n_edits > 0` is that signal, and it's what this module exists to produce.
If you want raw boundary drift regardless of editing, compare
`raw_start`/`raw_end` to `corrected_start`/`corrected_end` yourself.

## Validated on a 14-species smoke test (2026-09-24)

Real numbers, not a specification - useful as a sanity baseline before a
full run, not a target to reproduce exactly. `--organelle both
--species-list <14 species incl. Acer_campestre>`: 55 mito + 504 pltd
gene calls, `changed` (`n_edits>0`) in 58% of mito calls vs. 37% of
plastid calls, edit density 3.98/kb (mito) vs. 1.24/kb (plastid) over the
corrected coding span - directionally right (mitochondrial editing denser
than plastid) and the plastid rate is the right order of magnitude versus
the literature (~30-40 total sites genome-wide in well-studied species,
concentrated in ~50-60kb of CDS -> roughly 0.6-0.7/kb; this run's 1.24/kb
is about 2x that, plausible for an uncalibrated MVP on real predictions
rather than literature-verified sites, not the 10-100x-off signature of a
real bug).

### RefSeq profiles for the other mito genes

`00_build_reference_profiles.py` only profiles genes single-copy in at
least 80% of species, which in angiosperm mitochondria leaves 7: atp4,
ccmB, ccmC, mttB, nad5, nad9 and rps14. atp6, atp9, cox1, rps12 and the
rest are too often duplicated or missing, so their calls kept the raw
HMM boundaries, which ran well past the gene (Arabidopsis rps12: 657 bp
against RefSeq's 378).

`src/00b_build_lineage_profiles.py --lineage angiosperm` builds profiles
for every single-exon mito protein gene from the 11 angiosperm RefSeq
mitogenomes (`trans_splicing/reference_genes.py`'s
`REFERENCE_MITOGENOMES`), into `reference/profiles_refseq/mito/`: 26
genes; rpl10, rps2, rps11, sdh3 and sdh4 have fewer than 5 verified
RefSeq proteins and get none. RefSeq's `/translation` is the edited
protein, so these profiles point straight at the target. `scan` uses one
for a gene only when the dataset-built (or bundled dataset) profiles have
none, so the 7 above are unchanged, and the lineage call, which compares
only genes both lineages have profiles for, is unchanged.

These profiles drop alignment columns that are gaps in more than half
the proteins, not 90% as for the dataset profiles: with ~10 curated
proteins, 90% keeps a column only one species fills, and orfedit then
over-extends towards that species' extension (atp1 1845 bp vs RefSeq's
1524 at 0.9).

Validated on the Arabidopsis smoke assembly against RefSeq with
Arabidopsis (NC_037304.1) held out of the profiles
(`--exclude NC_037304.1 --out-dir ...`): the 19 mito genes that were
previously raw have median protein identity 0.837 before and 1.000 after,
and 13 match RefSeq's length exactly. ccmFn is the one that got worse
(0.875 -> 0.826): in Brassicaceae it is split into ccmFN1 and ccmFN2,
and a full-length profile can only partly fit either piece. atp6 stays
raw in the check because Arabidopsis has two copies and `scan` corrects
the best-scoring one, not the one RefSeq is compared against.

### Moss profiles and the lineage call

The dataset-built profiles are angiosperm profiles in practice, because
the dataset is mostly angiosperms. Aligned against them, mosses were
"corrected" towards the angiosperm protein: 98% of moss mito calls
changed (dicots 37%), at a median score of 391 (dicots 555), with 15% of
calls aligning under 80% of the profile (dicots 0.9%).

`src/00b_build_lineage_profiles.py --lineage moss` builds bundled moss
profiles (`reference/lineages/moss/mito/`) from RefSeq's curated proteins
in 20 moss mitogenomes (`trans_splicing/src/reference_genes.py`'s
`MOSS_MITOGENOMES`). It builds one for each bundled mito profile whose
gene is single-exon in mosses: atp4, ccmB, ccmC, mttB and rps14. nad5 and
nad9 are multi-exon in mosses and get none.

`01_scan_editing.py --lineage auto` (the default) runs orfedit on the
angiosperm profiles and then, for mito, on the moss ones, and calls each
species' lineage with `common/lineage.py`: moss when at least 4 of the 5
genes are scored and moss wins by 50 or more per gene on average. The call
goes to `results/lineage.tsv`. A moss-called species keeps its moss mito
calls only, so its nad5/nad9 are left uncorrected rather than corrected
against an angiosperm profile. Its plastid calls are unchanged.

Calibration, 956-species mito dataset with ToL-ID clades as the truth:
mosses win by 88-150 per gene, while dicots and monocots lose by 57-430
and liverworts by -25 to +6. Ferns and lycophytes have lost ccmB/ccmC,
and 16 of them preferred moss by 55-93 on the 3 genes left, which is why
the call needs 4 genes. With that rule, 174/174 mosses are called moss
and nothing else is. Moss results (auto vs `--lineage angiosperm`):

| | angiosperm profiles | moss profiles |
|---|---|---|
| calls changed | 98% | 44% |
| edits per call | 2.4 | 0.7 |
| median score | 391 | 592 |
| aligned < 80% | 15% | 0.5% |

Every non-moss row is byte-identical between the two runs. Liverworts
(87% changed), ferns and lycophytes still align against angiosperm
profiles and have the same low-score symptoms. Treat their corrections
as lower-confidence.

## Not yet run at full-dataset scale / not in `run_all.sh`

Like `orf_scan`, this is new and not yet validated broadly enough to be
part of the default chain - run it explicitly on a `--species-list` subset
first and compare against the smoke-test numbers above before trusting a
full run.

## Try this

Compare `n_edits` and `changed` rates between `nad5`/`nad9`/`ccmB` (mito
genes already close to the 80%/90% marker-gene thresholds, i.e. genes with
more real biological irregularity to begin with) and comfortably-clean
plastid genes (`atpA`, `atpB`, ...) - plastid editing is much sparser than
mitochondrial, so a similarly high correction rate there would be a red
flag for the profile/scoring, not a finding.
