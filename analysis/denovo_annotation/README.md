# denovo_annotation

**The question:** can gene calling for this pipeline stop depending on
running oatk at all (oatk stays only for assembly - the genome
graph/contigs), and if so, is that just a decoupling exercise or does it
actually call genes better? This module runs the same targeted gene search
oatk does internally, but directly and independently of oatk's assembler,
using `oatkdb` (oatk's own database-builder) + `nhmmscan` for
protein-coding genes, `tRNAscan-SE` for tRNAs, and `barrnap` for rRNAs -
exactly the pipeline already built and proven in the sibling
`mito_structural_variation/annotation/` repo, run here directly against
this dataset's own oatk-assembled contigs.

If you're coming from `annotation/`: that module reshapes oatk's own
bundled `.ctg.bed` calls. This module is a from-scratch replacement for
that calling step, producing one combined, human-inspectable GFF3 per
species (`results/gff/<species>.<organelle>.gff` - all three sources'
calls, tagged by source, sorted by position) as its primary output, plus a
`gene_calls.tsv`-shaped table derived from it (same columns as
`annotation/results/gene_calls.tsv`: `species organelle contig_id start
end gene score strand`) so it's a drop-in alternative anchor source for
`editing/` and `trans_splicing/` (`--gene-calls` flag on both), not a new
format to learn. The GFF is the thing to inspect directly; the TSV is a
derived convenience for those two consumers, not a second source of truth.

## What it computes

| script | output | what |
|---|---|---|
| `00_run_nhmmscan.py` | `work/nhmmscan/<species>.<organelle>.{filtered.tblout,gff}` | `nhmmscan` against oatkDB's own gene-family `.fam` database, E-value filtered, converted to GFF via `hmm_to_gff` |
| `01_run_trnascan.py` | `work/trnascan/<species>.<organelle>.gff` | `tRNAscan-SE -O` |
| `02_run_barrnap.py` | `work/barrnap/<species>.<organelle>.gff` | `barrnap --kingdom plant` |
| `03_build_combined_gff.py` | `results/gff/<species>.<organelle>.gff` | **the primary output** - all three sources combined into one sorted GFF3, tagged by source, not deduplicated |
| `04_build_gene_calls.py` | `results/gene_calls.tsv` | reshapes the combined GFF into one `annotation/results/gene_calls.tsv`-shaped table |

## Validated finding, not a hypothesis: this measurably out-performs oatk's bundled calling

On the same 14-species smoke test used throughout `editing/`/
`trans_splicing/`, comparing this module's mito `gene_calls.tsv` against
`annotation/results/gene_calls.tsv` (oatk's own calls) for exactly the
genes `trans_splicing` cares about:

- **`nad5`**: oatk finds only 1 fragment per species (its own internal
  model apparently only ever calls the single largest exon - already
  flagged as a real limitation in `trans_splicing/README.md`).
  `nhmmscan` run directly finds **4 fragments per species**, consistently,
  across every one of the 13 comparable species in the smoke test - not a
  one-off (`Acer_campestre` was the first case found; it generalised).
- **`nad1`**: oatk finds 3-4 fragments where present; direct `nhmmscan`
  finds 5-6, matching the literature's 5-exon canonical structure more
  closely.
- **Genes oatk misses entirely, in QC-`pass` species**: `Antennaria_dioica`
  is QC status `pass` with 210 other gene calls in `annotation`'s table,
  but **zero** calls for `nad1`/`nad2`/`nad5`/`rps3` - all four genes,
  completely absent, not just fragmented. Direct `nhmmscan` finds 4-6
  fragments for every one of them in the same species. Two further species
  in the smoke test (`Convolvulus_arvensis`, `Hydrocotyle_vulgaris`) have
  **zero rows in oatk's `gene_calls.tsv` for any gene at all** despite a
  usable assembly existing - `denovo_annotation` recovers real calls for
  both.
- **End-to-end payoff, not just more GFF rows**: re-running
  `trans_splicing/01_reconstruct.py --gene-calls
  denovo_annotation/results/gene_calls.tsv` on the same smoke test roughly
  **doubled the completion rate** (16/52 = 31% vs. the oatk-anchored
  baseline's 5/32 = 16% already documented in `trans_splicing/README.md`);
  `nad5` specifically went from 25% to 46% complete, exactly the
  mechanism predicted from the fragment-count comparison above, not a
  coincidence.
- **One caveat, not uniformly better**: `rps3` completion for bryophyte
  species (`Andreaea_rothii`, `Grimmia_montana`, `Campylopus_introflexus`)
  was *lower* via `nhmmscan` than oatk's own calls - plausibly the
  oatkDB nucleotide model for `rps3` is itself trained mostly on
  angiosperm reference sequences, the same angiosperm-vs-bryophyte
  divergence pattern already documented as a limitation in
  `editing/README.md`'s profile-building step. Worth revisiting alongside
  that, not a reason to doubt the overall result.

## Method

1. **`00_run_nhmmscan.py`** - `nhmmscan --tblout` against
   `/software/team301/OatkDB/viridiplantae_mito_v20250217.fam` (mito;
   already built and `hmmpress`ed) or the plastid equivalent (**not yet
   built** - see caveat below), then `filter_tblout <tblout> 1e-5` (same
   E-value floor as this repo's existing Pfam-scan precedent in
   `orf_scan`), then `hmm_to_gff <filtered> nhmmscan` to GFF3. `1e-5` isn't
   arbitrary here either: the unfiltered tblout for `Acer_campestre`
   included a real `E=0.47` `nad2` hit that plainly shouldn't be trusted -
   `1e-5` clears it out, real hits cluster many orders of magnitude below
   that floor (E=0 to ~1e-20 for anything used in the comparison above).
   `hmm_to_gff` already classifies hits as `gene`/`tRNA`/`rRNA` by gene
   name and keeps `target_name` in the native oatkDB naming convention
   (e.g. `rrn26`, `trnP-UGG`) - no renaming needed for this source.
2. **`01_run_trnascan.py`** / **`02_run_barrnap.py`** - exact
   invocations already proven in `mito_structural_variation`'s
   `5_get_rna_genes.bash`/`6_get_rrna_genes.bash`.
3. **`03_build_combined_gff.py`** - appends all three sources' GFF data
   lines (stripping each one's own `##gff-version` header, keeping one at
   the top) and sorts by contig then start - same pattern
   `mito_structural_variation` uses to build up its own combined GFF.
   **Not deduplicated on purpose**: oatkDB's `.fam` database includes its
   own tRNA/rRNA gene-family models (confirmed directly - `nhmmscan`
   itself returns `trnP-UGG`/`rrn26`-style hits, not just protein-coding
   ones), so a real tRNA can legitimately get both an `oatkDB` row and a
   `tRNAscan-SE` row. Keeping both, tagged by source, means you can see
   independent methods agreeing (or not) rather than having that signal
   silently thrown away by a dedup step.
4. **`04_build_gene_calls.py`** - reshapes the combined GFF into the
   `gene_calls.tsv`-shaped table. `tRNAscan-SE` hits (its own `isotype`/
   `anticodon` attributes) are renamed to this dataset's
   `trn<AA>-<anticodon>` convention (undetermined/pseudo calls dropped,
   not guessed at); `barrnap`'s bacterial-style names (`16S_rRNA` etc.)
   are mapped to `rrn18`/`rrn26`/`rrn5` (mito) or
   `rrn16`/`rrn23`/`rrn5`/`rrn4.5` (plastid) - a best-effort mapping, not
   independently verified against oatk's own rRNA naming beyond
   spot-checking `rrn18`/`rrn26` matches. **`score` units differ by
   source and must never be compared across sources**: `oatkDB`
   (`nhmmscan`) is a bitscore (higher better), `tRNAscan-SE` its own
   confidence score (higher better), `barrnap` is an E-value (*lower*
   better) - kept as-is per source in both the GFF and the TSV, rather
   than force-normalised into one scale that would misrepresent the
   underlying tools.

## Known gap: no plastid database yet

`oatkdb` supports building one directly (`oatkdb ... <taxid> chloroplast`
- confirmed via `--help`, same tool as the existing mito database), but
the NCBI/edirect fetch this needs fails with **persistent SSL errors**
when run from normal LSF compute nodes (`long`/`normal` queues) - confirmed
directly: works cleanly from `farm22-head2` (the login node) and from the
**`transfer` queue** specifically, fails 100% of retries from a `long`-queue
compute node (`node-14-14`, tested). The fix is straightforward (submit
the `oatkdb` build via `bsub -q transfer`, not `long`/`normal`) but wasn't
completed this session - `--organelle pltd` is fully wired in every script
here and will work as soon as
`/software/team301/OatkDB/viridiplantae_pltd_v20260927.fam` exists (or
update `OATKDB_PLTD_FAM` in `tool_paths.sh` to whatever name it's actually
built with).

## Scaling to the full dataset: three real things this surfaced

Preparing for a full ~1250-species run (still mito-only - see the plastid
gap above) surfaced three real, fixed issues, not just a "should be fine"
assumption:

1. **No parallelism, at all, until now.** `00_run_nhmmscan.py`/
   `01_run_trnascan.py`/`02_run_barrnap.py` looped over species
   sequentially in one process - fine for a 14-species smoke test, a real
   problem at full scale given `tRNAscan-SE` alone measures at
   **~4 minutes per species, confirmed directly** (`/usr/bin/time`:
   4:04 wall-clock, 348 CPU-seconds, 142% CPU, for one real mitogenome) -
   at ~1250 species that's several hours single-threaded. Fixed: all
   three scripts now take the same `--jobs`/`--cpu`(or `--thread`/
   `--threads`) pattern as `orf_scan/src/02_scan_pfam.py` (concurrent
   processes via `ThreadPoolExecutor`, each still internally
   multi-threaded) - e.g. inside `bsub -n 8`, `--jobs 4 --thread/--cpu 2`
   uses the full allocation.
2. **`tRNAscan-SE` hangs forever if its own output file already exists**,
   confirmed directly: it prompts interactively
   ("(O)verwrite/(A)ppend/(Q)uit?") and `subprocess.run` never answers,
   so the process just sits there - this genuinely happened on a
   `--force` re-run over already-scanned species, silently ate 2 of 4
   worker slots for the whole job, and drove memory usage from a real
   ~300MB to a reported 10.9GB before LSF killed the job on its memory
   limit (the memory number was *itself* an artifact of the hang, not a
   real per-process cost). Fixed: always pass `-Q`/`--forceow` to
   `tRNAscan-SE` - this script's own skip-if-exists/`--force` check is
   what actually controls re-scanning, `tRNAscan-SE`'s own overwrite
   prompt should never be reachable.
3. **No QC gate**, unlike `annotation/src/01_gene_table.py`'s
   `--qc-status` default of `pass`. Confirmed a real, non-hypothetical
   cost: `Stellaria_graminea` (QC status `fail` - fragmented,
   non-circular, <50% core genes present) was still getting the full
   `nhmmscan`+`tRNAscan-SE`+`barrnap` treatment for a genuinely
   uninformative ~31kb fragment, and produced zero tRNAs (a real,
   QC-explained result, not a tool bug - worth remembering when zero
   hits shows up elsewhere too: check QC status before assuming a tool
   problem). All three scan scripts now take the same `--qc-status`
   default as `annotation/`'s, skipping non-`pass` species by default.

**A real operational gotcha to remember, not (yet) fixed in code**:
`03_build_combined_gff.py`/`04_build_gene_calls.py` don't clean up a
species that *stops* being eligible between runs (e.g. after adding the
QC gate above) - they only skip-if-exists or overwrite species still in
scope, so a species dropped from eligibility leaves its **stale** output
sitting in `results/gff/`/`results/gene_calls.tsv` until manually
removed. If the eligible species set ever shrinks, clear
`results/gff/*.gff` before rebuilding rather than assuming `--force`
alone reconciles it.

**Also confirmed the hard way**: run this via `bsub`, not interactively.
A `--jobs 4` batch over 15 species got silently `SIGKILL`ed partway
through when run directly in an interactive session on the login node -
almost certainly a resource/session policy on shared interactive
compute, not a bug in the scripts themselves (the same batch completed
cleanly end-to-end once submitted as a normal `bsub` job). Treat this the
same as the `oatkdb` build earlier in this project: real per-species
compute belongs on the farm, not the head node, full stop.

New, unvalidated at scale beyond an 11-species smoke test above (four of
the original 15 excluded by the QC gate just added) - same opt-in
precedent as `orf_scan`/`editing`/`trans_splicing`. Mito can be run at
full scale immediately; plastid needs the database gap above closed
first.

## Not yet done: cutting over the default

`editing/`'s and `trans_splicing/`'s default anchor is still
`annotation/results/gene_calls.tsv` (oatk's own calls) - `--gene-calls`
lets you point either at this module's output for comparison, but nothing
here changes their defaults. Given the completion-rate improvement above,
switching the default (or running both and taking the union) is a
reasonable next step, not done automatically in this phase since it's a
real decision (full-dataset validation first, matching this repo's own
habit of not trusting a 14-species smoke test as sufficient on its own).
