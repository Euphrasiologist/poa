"""poa - plant organelle annotator."""
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

try:
    __version__ = version("poa")
except PackageNotFoundError:  # running from a source checkout
    __version__ = "0.1.0.dev0"

PIPELINE_DIR = Path(__file__).resolve().parent / "pipeline"
