"""Locating the Rust tools poa drives (orfedit, transsplice)."""

from __future__ import annotations

import shutil
from pathlib import Path

from poa.pipeline.common.errors import StageError


def cargo_tool(name: str) -> str:
    """`name` on PATH, else ~/.cargo/bin/<name> (a `cargo install`)."""
    exe = shutil.which(name) or str(Path.home() / ".cargo" / "bin" / name)
    if not Path(exe).exists():
        raise StageError(
            f"{name} binary not found at {exe} - install it from bioconda, or `cargo install --path .` in its repo"
        )
    return exe
