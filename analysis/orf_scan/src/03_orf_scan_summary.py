#!/usr/bin/env python3
"""Categorize each non-core ORF's best Pfam hit as TE-like, mitovirus-like,
or other, and aggregate to one summary row per species x organelle.

The keyword lists below are deliberately curated/documented rather than
data-driven (unlike the core-gene thresholds elsewhere in this suite) -
this one check specifically needs real biological knowledge of which Pfam
families indicate TE or mitovirus content, which can't be inferred from
prevalence statistics the way core-gene completeness can. They are not
exhaustive; treat `category` as a starting point for manual inspection,
not a final verdict - and the raw per-ORF table always keeps the full
Pfam hit (name + description) regardless of category, so nothing is hidden
by the categorization.

Output: analysis/orf_scan/results/orf_scan_summary.tsv
        (+ recategorizes the `category` column into every per-species .orfs.tsv, in place)
"""
from __future__ import annotations

import os
import re
import sys
from pathlib import Path

import pandas as pd

# CODE_DIR is this checkout's analysis/ (shared modules, bundled reference
# data); data/ and every module's results/ and work/ live under the data
# root - PLANT_ORGANELLE_DATA_ROOT, or this checkout if unset.
CODE_DIR = Path(__file__).resolve().parents[2]
ROOT_DIR = Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT", str(CODE_DIR.parent)))
ANALYSIS_DIR = ROOT_DIR / "analysis"
PER_SPECIES_DIR = ANALYSIS_DIR / "orf_scan" / "results" / "per_species"
FILENAME_RE = re.compile(r"^(.+)\.(mito|pltd)\.orfs\.tsv$")

# Matched case-insensitively against the Pfam family name + free-text description.
TE_KEYWORDS = [
    "transposase", "transposon", "retrotransposon", "reverse transcriptase",
    "integrase", "gag-pol", "gag polyprotein", "zf-CCHC", "DDE", " rve ",
    "mariner", "hAT family", "CACTA", "PIF/harbinger", "helitron", "mutator",
    "LTR retrotransposon",
]
MITOVIRUS_KEYWORDS = [
    "mitovirus", "narnavirus", "rna-dependent rna polymerase", "rna dependent rna polymerase",
    "viral rna-directed rna polymerase", "rdrp",
]


def categorize(name: str, description: str) -> str:
    text = f"{name} {description}".lower()
    if any(k.lower() in text for k in MITOVIRUS_KEYWORDS):
        return "mitovirus_like"
    if any(k.lower() in text for k in TE_KEYWORDS):
        return "te_like"
    return "other_pfam_hit"


def main():
    rows = []
    for path in sorted(PER_SPECIES_DIR.glob("*.orfs.tsv")):
        m = FILENAME_RE.match(path.name)
        if not m:
            continue
        species, organelle = m.groups()
        # dtype=str for evalue/score: 02_scan_pfam.py already formats these as
        # clean strings (e.g. "1.5e-29"); reading them back as float64 here and
        # writing again would silently reintroduce float round-trip artifacts
        # (e.g. "1.4999999999999998e-29") since we don't need numeric evalue here.
        df = pd.read_csv(path, sep="\t", dtype={"evalue": str, "score": str})

        if "pfam_name" in df.columns:
            has_hit = df["pfam_name"].notna()
            df["category"] = "no_hit"
            df.loc[has_hit, "category"] = [
                categorize(n, d) for n, d in zip(df.loc[has_hit, "pfam_name"], df.loc[has_hit, "description"].fillna(""))
            ]
            df.to_csv(path, sep="\t", index=False)
        else:
            has_hit = pd.Series([False] * len(df))
            df["category"] = "not_scanned"

        rows.append({
            "species": species, "organelle": organelle,
            "n_noncore_orfs": len(df),
            "n_pfam_hits": int(has_hit.sum()),
            "n_te_like": int((df["category"] == "te_like").sum()),
            "n_mitovirus_like": int((df["category"] == "mitovirus_like").sum()),
        })

    out_path = ANALYSIS_DIR / "orf_scan" / "results" / "orf_scan_summary.tsv"
    columns = ["species", "organelle", "n_noncore_orfs", "n_pfam_hits", "n_te_like", "n_mitovirus_like"]
    pd.DataFrame(rows, columns=columns).to_csv(out_path, sep="\t", index=False)
    print(f"[info] orf_scan_summary: {len(rows)} species x organelle rows -> {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
