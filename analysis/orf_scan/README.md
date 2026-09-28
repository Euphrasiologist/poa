# orf_scan

**The question:** oatk's core-gene annotation only looks for the genes it's
told to look for (a targeted HMM family database of known mitochondrial/
plastid genes). Plant organelle genomes also routinely carry other
content — transposable-element remnants, and (mitochondria specifically)
mitovirus-derived sequences — that a targeted core-gene search will never
report, because it isn't looking for it. This module finds ORFs
independently with NCBI's `ORFfinder`, keeps only the ones **not** already
explained by oatk's own gene calls, and screens those against the full
Pfam-A library to see what's actually there.

If you're coming from the `mito_structural_variation` repo: this
deliberately reuses its `ORFfinder → hmmscan` idea (rather than its fuller
from-scratch reannotation pipeline, which is redundant for core genes —
see [`../annotation/README.md`](../annotation/README.md)), because this
specific question — what's the *non-core* content — genuinely isn't
answered anywhere else in this dataset.

## This is the heaviest module in the suite

Everything else here runs on the full ~1250-species dataset in minutes.
This one doesn't: `hmmscan` against Pfam-A's ~20,000 profiles measures at
roughly 15-90 seconds *per species* (mito genomes, with more ORFs, cost
more than plastid genomes) even against the pressed/binary database
already set up on this cluster. A full run is a multi-hour job. Because of
that, **it is not part of the default `analysis/run_all.sh` chain** —
run it explicitly:

```bash
# quick look at a handful of species
analysis/orf_scan/src/01_find_orfs.py --species-list my_species.txt
analysis/orf_scan/src/02_scan_pfam.py --species-list my_species.txt
analysis/orf_scan/src/03_orf_scan_summary.py

# full dataset - wrap in an LSF job with real parallelism
bsub -n 16 -q normal -o orf_scan.out -e orf_scan.err \
  bash -c '
    analysis/orf_scan/src/01_find_orfs.py --organelle both
    analysis/orf_scan/src/02_scan_pfam.py --jobs 4 --cpu 4
    analysis/orf_scan/src/03_orf_scan_summary.py
  '
```

`--jobs` (in `02_scan_pfam.py`) runs several `hmmscan` processes
concurrently, each still multi-threaded via `--cpu` — size these to your
LSF allocation (`--jobs 4 --cpu 4` for a 16-core job).

## Method

1. **`01_find_orfs.py`** — `ORFfinder -g 1` (standard genetic code — plant
   mito/plastid, unlike animal mitochondria, use the standard code)
   `-ml 300` (100 aa minimum) `-n true` (drop nested ORFs) on each
   species' resolved contigs. Verified empirically: the default 75nt
   cutoff finds ~4000 ORFs on a single mitogenome, almost all spurious
   short background ORFs with no biological meaning; 300nt cuts that to
   ~200 while still comfortably keeping real domains (which are typically
   hundreds of amino acids). Any ORF overlapping an oatk core-gene call
   (`.ctg.bed`) by ≥50% of its own length is dropped — verified
   empirically that `ORFfinder` independently re-finds real annotated
   genes (NADH dehydrogenase, COX3, ribosomal proteins...), so without
   this filter the scan would mostly just rediscover known genes rather
   than surface anything new.
2. **`02_scan_pfam.py`** — one `hmmscan --cut_ga` call per species (not
   per ORF, to avoid reloading the ~1.7GB database thousands of times)
   against the full Pfam-A library, keeping the best domain hit per ORF.
3. **`03_orf_scan_summary.py`** — categorizes each hit as `te_like`,
   `mitovirus_like`, or `other_pfam_hit` against curated keyword lists
   (see the script — things like `transposase`, `reverse transcriptase`,
   `zf-CCHC` for TEs; `Mitovir_RNA_pol`, `RdRp`, `mitovirus` for
   mitoviruses). **These lists are curated, not data-driven** — unlike
   the core-gene thresholds elsewhere in this suite, this needs real
   biological knowledge of which Pfam families mean what, and they are
   not exhaustive. Treat `category` as a lead worth checking, not a final
   verdict — the raw Pfam hit (family name + description) is always kept
   regardless of category, so nothing is hidden.

## Reading the output

`per_species/<Species>.<organelle>.orfs.tsv`: one row per non-core ORF,
with its best Pfam hit (if any) and category. `orf_scan_summary.tsv`: one
row per species×organelle with counts. A validated real example already
found in this dataset: *Solanum nigrum*'s mitochondrial genome carries
several `Mitovir_RNA_pol` (PF05919, "Mitovirus RNA-dependent RNA
polymerase") and `RVT_2`/`zf-RVT` (reverse transcriptase) hits — genuine,
recognizable mitovirus and retroelement signal, not noise.

## Try this

Join `orf_scan_summary.tsv` against `repeats/results/recomb_repeats_summary.tsv`
— TE-derived sequences are often themselves repetitive, so do
`n_te_like`-heavy species also show higher repeat content? And against
`qc_basic_stats/results/gfa_stats.tsv`'s `dead_end_nodes` — could TE
insertions be contributing to assembly graph complexity?
