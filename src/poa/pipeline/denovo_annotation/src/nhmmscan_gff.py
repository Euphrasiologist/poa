"""nhmmscan --tblout (against an oatkDB .fam) -> E-value-filtered tblout + GFF3.

Python port of the two small Rust tools this step used to shell out to,
`filter_tblout` 0.1.1 and `hmm_to_gff` 0.1.1 (nhmmscan mode), so the
pipeline has no compiled glue of its own. Checked against both tools on
all 2096 tblouts of the 1250-species dataset: identical GFF (their float32
arithmetic and Rust number formatting included) and identical filtered
rows; tests/test_nhmmscan_gff.py keeps the Arabidopsis pair as a fixture.

  - a hit is kept if its E-value < `evalue` (filter_tblout) and <= 1e-5
    (hmm_to_gff's own floor), both compared as float32;
  - start/end are the envelope (envfrom/envto), smaller first;
  - feature type: `trn` in the model name -> tRNA, else `rrn` -> rRNA,
    else gene;
  - coverage = envelope length (end - start, not +1) / model length,
    `partial:<pct>` below 0.8, else `full`. Model lengths are read from the
    .fam's own LENG lines (hmm_to_gff embedded the same numbers for the two
    database versions it knew; reading the .fam works for any version).

The filtered tblout keeps the input's lines verbatim (filter_tblout
re-rendered them; 04_build_gene_calls.py reads it whitespace-split, so the
columns are what matter).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

COVERAGE = np.float32(0.8)
HMM_TO_GFF_MAX_EVALUE = np.float32(1e-05)
# nhmmscan --tblout columns (0-based, whitespace-split)
TARGET, QUERY, ENV_FROM, ENV_TO, STRAND, EVALUE, SCORE = 0, 2, 8, 9, 11, 12, 13


def fam_model_lengths(fam: Path) -> dict[str, int]:
    """Model name -> length (LENG) for every HMM in an HMMER3 .fam file."""
    lengths, name = {}, None
    with open(fam) as fh:
        for line in fh:
            if line.startswith("NAME "):
                name = line.split()[1]
            elif line.startswith("LENG ") and name is not None:
                lengths[name] = int(line.split()[1])
    return lengths


def rust_f32(x: np.float32) -> str:
    """Rust's `{}` for f32: shortest round-trip digits, no exponent, no `.0`."""
    return np.format_float_positional(x, unique=True, trim="-")


def rust_exp2(x: np.float32) -> str:
    """Rust's `{:.2e}` for f32: `1.23e-5`, `0.00e0` - no `+`, no zero padding."""
    mantissa, exp = f"{float(x):.2e}".split("e")
    return f"{mantissa}e{int(exp)}"


def filter_and_convert(tblout: Path, fam: Path, evalue: str) -> tuple[str, str]:
    """(filtered tblout text, GFF3 text) for one nhmmscan --tblout file."""
    threshold = np.float32(evalue)
    lengths = fam_model_lengths(fam)
    kept_lines, records = [], []
    for line in Path(tblout).read_text().splitlines(keepends=True):
        if line.startswith("#") or not line.strip():
            kept_lines.append(line)
            continue
        f = line.split()
        e = np.float32(f[EVALUE])
        if not e < threshold:
            continue
        kept_lines.append(line)
        if e > HMM_TO_GFF_MAX_EVALUE:
            continue
        records.append((f, e))

    gff_rows = []
    for n, (f, e) in enumerate(records, start=1):
        target, seqid = f[TARGET], f[QUERY]
        a, b = int(f[ENV_FROM]), int(f[ENV_TO])
        start, end = min(a, b), max(a, b)
        feature = "tRNA" if "trn" in target else "rRNA" if "rrn" in target else "gene"
        ratio = np.float32(end - start) / np.float32(lengths[target])
        coverage = f"partial:{float(ratio * np.float32(100.0)):.2f}" if ratio < COVERAGE else "full"
        attrs = f"target_name={target};ID={seqid}_{n};coverage={coverage};e_value={rust_exp2(e)}"
        gff_rows.append((seqid, start, "\t".join([
            seqid, "oatkDB", feature, str(start), str(end), rust_f32(np.float32(f[SCORE])), f[STRAND][0], ".", attrs])))
    gff_rows.sort(key=lambda r: (r[0], r[1]))  # stable, like Rust's sort_by
    gff = "##gff-version 3\n" + "".join(row + "\n" for _, _, row in gff_rows)
    return "".join(kept_lines), gff
