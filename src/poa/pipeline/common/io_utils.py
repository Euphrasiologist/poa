"""Shared incremental/atomic TSV helpers for analysis/ result tables.

Philosophy: every module builds one or more flat TSVs with one row per
species (or per species x organelle, or per species x gene). Re-running a
script should be cheap: skip species already present in the output unless
their input file has changed since the output was last written, or --force
is set. Writes are atomic (write to .tmp, then os.replace) so a job killed
mid-write never corrupts a shared table other scripts/students may be
reading concurrently.
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import Iterable

import pandas as pd


def load_existing(path: Path, key_cols: list[str]) -> tuple[pd.DataFrame, set]:
    """Return (existing_dataframe, set_of_key_tuples). Empty if file absent."""
    path = Path(path)
    if not path.exists():
        return pd.DataFrame(), set()
    df = pd.read_csv(path, sep="\t", dtype=str)
    if df.empty:
        return df, set()
    keys = set(df[key_cols].itertuples(index=False, name=None))
    return df, keys


def needs_recompute(key: tuple, existing_keys: set, output_path: Path, input_path: Path, force: bool) -> bool:
    """Should this key be (re)computed?

    - Always True if --force.
    - Always True if the key isn't in the existing output yet.
    - Otherwise True only if input_path is newer than output_path (the
      species' assembly was regenerated/fixed since we last computed it).
    """
    if force:
        return True
    if key not in existing_keys:
        return True
    output_path = Path(output_path)
    input_path = Path(input_path)
    if not output_path.exists() or not input_path.exists():
        return True
    return input_path.stat().st_mtime > output_path.stat().st_mtime


def _tidy_numeric_dtypes(df: pd.DataFrame) -> pd.DataFrame:
    """Columns that are all-integral-or-missing get pandas' nullable Int64 dtype,
    so e.g. a subgraph_id column with some NaN doesn't render as '1.0' in the TSV."""
    for col in df.columns:
        if df[col].dtype != object and not pd.api.types.is_float_dtype(df[col]):
            continue
        numeric = pd.to_numeric(df[col], errors="coerce")
        if numeric.isna().equals(df[col].isna()) and numeric.dropna().astype(float).apply(float.is_integer).all():
            df[col] = numeric.astype("Int64")
    return df


def atomic_write_tsv(path: Path, rows: Iterable[dict], columns: list[str]) -> None:
    """Write rows (list of dicts) to path atomically, always with `columns` as header order."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    df = pd.DataFrame(list(rows), columns=columns)
    df = _tidy_numeric_dtypes(df)
    df.to_csv(tmp_path, sep="\t", index=False)
    os.replace(tmp_path, path)


def merge_rows(existing_df: pd.DataFrame, new_rows: list[dict], key_cols: list[str], columns: list[str]) -> pd.DataFrame:
    """Replace/append new_rows into existing_df keyed on key_cols, return the merged frame."""
    new_df = pd.DataFrame(new_rows, columns=columns) if new_rows else pd.DataFrame(columns=columns)
    if existing_df.empty:
        return new_df
    if new_df.empty:
        return existing_df
    existing_df = existing_df.set_index(key_cols)
    new_df_idx = new_df.set_index(key_cols)
    existing_df = existing_df.drop(index=new_df_idx.index, errors="ignore")
    merged = pd.concat([existing_df, new_df_idx]).reset_index()
    return merged[columns]


def merge_species_scoped_tsv(path: Path, new_path: Path, species_filter: set | None) -> None:
    """Merge a (possibly --species-list-scoped) freshly-computed result table
    at `new_path` into the existing one at `path`, atomically.

    Found the hard way: a script whose output covers every species that
    belongs in the table (gene_calls.tsv and friends) used to just write
    `new_path`'s rows straight to `path` unconditionally. That's correct for
    an unscoped run (`species_filter` is None - `new_path` already covers
    everyone), but a `--species-list`-scoped run only ever computes `new_path`
    for the listed species, so writing it straight to `path` silently
    discarded every OTHER species already in the table - confirmed directly
    against a real ~1250-species dataset (a 2-species-scoped run collapsed a
    418216-row table to 811 rows).

    So: rows in `path` for a species NOT in `species_filter` are kept
    untouched; rows for a species IN `species_filter` are dropped and
    replaced by whatever `new_path` has for it - including replaced with
    nothing, if that species genuinely no longer produces any row this run.
    """
    path = Path(path)
    new_path = Path(new_path)
    new_df = pd.read_csv(new_path, sep="\t", dtype=str) if new_path.exists() else pd.DataFrame()

    if species_filter is None or not path.exists():
        merged = new_df
    else:
        existing = pd.read_csv(path, sep="\t", dtype=str)
        kept = existing[~existing["species"].isin(species_filter)] if "species" in existing.columns else existing
        merged = pd.concat([kept, new_df], ignore_index=True) if not new_df.empty else kept

    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    merged.to_csv(tmp_path, sep="\t", index=False)
    os.replace(tmp_path, path)
