#!/usr/bin/env python3
"""Command-line wrapper for poa.pipeline.unitig_coords.unitig_map (see its docstring)."""

import sys
from pathlib import Path

# run by path: make the poa package this file belongs to importable
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))

from poa.pipeline.unitig_coords.unitig_map import main  # noqa: E402

if __name__ == "__main__":
    sys.exit(main())
