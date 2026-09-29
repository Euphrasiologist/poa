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
