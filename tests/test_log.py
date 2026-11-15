"""Logging setup: the [info]/[warn]/[err] format stage logs are checked
for, and that log output follows stderr wherever it is redirected.

    python3 -m unittest tests.test_log
"""

from __future__ import annotations

import contextlib
import io
import logging
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from poa import run  # noqa: E402
from poa.pipeline.common.log import StderrHandler, setup_logging  # noqa: E402

log = logging.getLogger("poa.test")


class Logging(unittest.TestCase):
    def setUp(self):
        root = logging.getLogger()
        saved = root.handlers[:], root.level
        self.addCleanup(lambda: (setattr(root, "handlers", saved[0]), root.setLevel(saved[1])))
        root.handlers = []

    def test_tags(self):
        setup_logging()
        buf = io.StringIO()
        with contextlib.redirect_stderr(buf):
            log.debug("hidden")
            log.info("a: scanned=1 failed=0")
            log.warning("b")
            log.error("c")
        self.assertEqual(buf.getvalue(), "[info] a: scanned=1 failed=0\n[warn] b\n[err] c\n")

    def test_setup_twice_adds_one_handler(self):
        setup_logging()
        setup_logging()
        self.assertEqual(sum(isinstance(h, StderrHandler) for h in logging.getLogger().handlers), 1)

    def test_stage_log_captures_logging(self):
        """poa run's per-stage log (fd 2 redirected to a file) gets log
        records too, and a logged [err] fails the stage."""
        setup_logging()
        logs = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, logs)
        run.run_stage("ok", lambda: log.info("fine"), logs)
        self.assertEqual((logs / "ok.log").read_text(), "[info] fine\n")
        with self.assertRaises(run.StageFailed):
            run.run_stage("bad", lambda: log.error("broke"), logs)


if __name__ == "__main__":
    unittest.main()
