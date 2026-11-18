# Changelog

## Unreleased

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
