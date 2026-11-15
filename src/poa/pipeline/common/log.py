"""Logging for poa's stages.

Library code logs through `logging.getLogger(__name__)` and never
configures anything; entry points (`poa`, each stage's main()) call
setup_logging() once. Records go to stderr as "[info] ...", "[warn] ...",
"[err] ...": the format the stages printed before they used logging,
which `poa run` (and tests/smoke/run_smoke.sh) still scan stage logs for.
"""

from __future__ import annotations

import logging
import sys

TAGS = {logging.DEBUG: "debug", logging.INFO: "info", logging.WARNING: "warn", logging.ERROR: "err"}


class TagFormatter(logging.Formatter):
    def format(self, record: logging.LogRecord) -> str:
        record.tag = TAGS.get(record.levelno, "err" if record.levelno > logging.ERROR else "info")
        return super().format(record)


class StderrHandler(logging.StreamHandler):
    """Writes to whatever sys.stderr is at the time of each record, not when
    the handler was made, so redirecting or capturing stderr (`poa run`'s
    per-stage logs, test runners) also catches log output."""

    def __init__(self) -> None:
        super().__init__()

    @property
    def stream(self):
        return sys.stderr

    @stream.setter
    def stream(self, _value) -> None:
        pass


def setup_logging(level: int = logging.INFO) -> None:
    """Send log records to stderr in the [tag] format; safe to call twice."""
    root = logging.getLogger()
    root.setLevel(level)
    if not any(isinstance(h, StderrHandler) for h in root.handlers):
        handler = StderrHandler()
        handler.setFormatter(TagFormatter("[%(tag)s] %(message)s"))
        root.addHandler(handler)
