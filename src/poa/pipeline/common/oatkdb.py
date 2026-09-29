"""Locate the oatkDB gene-family HMM database (.fam) for an organelle.

$OATKDB_MITO_FAM / $OATKDB_PLTD_FAM win if set (an already-hmmpress'd
.fam of your own). Otherwise the database bundled with poa
(poa/data/oatkdb/*.fam.gz, see NOTICE there) is decompressed and
hmmpress'd once into the cache - $POA_CACHE, else $XDG_CACHE_HOME/poa,
else ~/.cache/poa - and reused from then on. Preparing it is atomic
(built in a temp dir, then renamed), so parallel jobs can't see a
half-pressed database.
"""
from __future__ import annotations

import gzip
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

BUNDLED_DIR = Path(__file__).resolve().parents[2] / "data" / "oatkdb"
BUNDLED = {
    "mito": "viridiplantae_mito_v20250217",
    "pltd": "viridiplantae_pltd_v20260928",
}
ENV = {"mito": "OATKDB_MITO_FAM", "pltd": "OATKDB_PLTD_FAM"}
PRESSED_SUFFIXES = (".h3f", ".h3i", ".h3m", ".h3p")


def cache_dir() -> Path:
    if os.environ.get("POA_CACHE"):
        return Path(os.environ["POA_CACHE"])
    base = os.environ.get("XDG_CACHE_HOME") or str(Path.home() / ".cache")
    return Path(base) / "poa"


def fam_path(organelle: str) -> Path:
    """Path to a ready-to-use (hmmpress'd) .fam for 'mito' or 'pltd'."""
    if os.environ.get(ENV[organelle]):
        return Path(os.environ[ENV[organelle]])
    name = BUNDLED[organelle]
    dest = cache_dir() / "oatkdb" / name
    fam = dest / f"{name}.fam"
    if all(Path(f"{fam}{s}").exists() for s in PRESSED_SUFFIXES):
        return fam
    hmmpress = shutil.which("hmmpress")
    if hmmpress is None:
        raise FileNotFoundError("hmmpress (HMMER) not found on PATH - needed once to prepare the bundled oatkDB")
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = Path(tempfile.mkdtemp(prefix=f".{name}.", dir=dest.parent))
    try:
        with gzip.open(BUNDLED_DIR / f"{name}.fam.gz") as src, open(tmp / fam.name, "wb") as out:
            shutil.copyfileobj(src, out)
        subprocess.run([hmmpress, str(tmp / fam.name)], check=True, capture_output=True)
        try:
            tmp.rename(dest)
        except OSError:  # another job got there first
            if not all(Path(f"{fam}{s}").exists() for s in PRESSED_SUFFIXES):
                raise
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return fam
