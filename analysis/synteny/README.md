# synteny

**The question:** the top-level README already suggested a classic
exercise — plastid DNA inserted into mitochondria (MTPT, "mitochondrial
plastid DNA transfer"), detected with:

```bash
minimap2 -x asm5 mito.fasta plastid.fasta > mito_vs_plastid.paf
```

This module automates exactly that across the whole dataset, plus an
opt-in mode for comparing contigs between different species.

## Two modes

1. **MTPT mode (default)** — for every QC-eligible species with both a
   mito and plastid assembly, aligns mito (reference) vs plastid (query)
   with `minimap2 -x asm5` and saves the PAF to `results/mtpt/`.
   `03_mtpt_summary.py` turns every PAF into one row
   (`total_aligned_length_bp`, `n_alignment_blocks`, `longest_block_bp`,
   `pct_plastid_genome_in_mito`) in `mtpt_summary.tsv` — directly
   answering the exercises the top-level README already proposed
   ("quantify total aligned length... compare across species... test
   whether MTPT burden correlates with genome size") without writing a
   one-off script per question.
2. **Pairwise mode** — for a student-chosen `--species-list` (e.g. one
   genus) and a single `--organelle`, aligns every pair within that list.
   **Deliberately guarded**: it refuses more than 30 species unless you
   pass `--allow-all-pairs`, because an all-vs-all alignment across the
   full ~1250-species dataset would no longer be lightweight. Run it like:

   ```bash
   analysis/synteny/src/01_pairwise_paf.py --mode pairwise \
     --organelle pltd --species-list my_genus.txt
   ```

## Dotplots

`02_dotplot.py` renders every saved PAF as a simple matplotlib dotplot —
each alignment block drawn as a line from `(target_start, query_start)` to
`(target_end, query_end)` (mirrored for `-` strand blocks), colored by
strand. If you haven't read a PAF-based dotplot before: a straight
diagonal line means the two sequences are collinear over that block; a
line running the "wrong way" (anti-diagonal) means an inversion; a gap
means no alignment there (either genuinely absent, or too diverged for
`asm5`'s divergence tolerance).

## Note on the QC contamination pre-screen

[`qc_basic_stats`](../qc_basic_stats/README.md) runs a similar
mito-vs-plastid alignment internally as a quick contamination check
(flagging species where the aligned fraction looks like mis-binning rather
than genuine MTPT). This module's MTPT mode is the fuller, persisted
version of that same idea — expect the two `pct_*` numbers to closely
agree, since they're the same alignment normalized by the same resolved
contig length.
