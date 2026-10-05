"""poa run's input handling: what it accepts, what it refuses before running
any stage, and how it lays out the data root. The full run is covered by
tests/smoke/run_smoke.sh.

    python3 -m unittest tests.test_run
"""
from __future__ import annotations

import argparse
import gzip
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from poa import run  # noqa: E402

FASTA = b">ctg1\nACGT\n"


def parse(*argv: str) -> argparse.Namespace:
    p = argparse.ArgumentParser()
    run.add_arguments(p)
    return p.parse_args(argv)


class DefaultName(unittest.TestCase):
    def test_oatk_names(self):
        self.assertEqual(run.default_name("x/Arabidopsis_thaliana.mito.ctg.fasta"), "Arabidopsis_thaliana")
        self.assertEqual(run.default_name("Arabidopsis_thaliana.k1001.s31.c100.pltd.ctg.fasta.gz"),
                         "Arabidopsis_thaliana")

    def test_unsafe_characters(self):
        self.assertEqual(run.default_name("my sample+1.fa"), "my_sample_1")


class Refusals(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.fa = self.tmp / "s.fasta"
        self.fa.write_bytes(FASTA)

    def refuses(self, *argv: str) -> str:
        with self.assertRaises(SystemExit) as cm:
            run.run(parse(*argv, "-o", str(self.tmp / "out")))
        return str(cm.exception.code)

    def test_no_fasta(self):
        self.assertIn("--mito and/or --pltd", self.refuses())

    def test_gfa_without_fasta(self):
        self.assertIn("--mito-gfa needs --mito", self.refuses("--pltd", str(self.fa), "--mito-gfa", str(self.fa)))

    def test_bad_name(self):
        self.assertIn("--name", self.refuses("--mito", str(self.fa), "-n", "a.b"))

    def test_missing_and_empty_inputs(self):
        self.assertIn("no such file", self.refuses("--mito", str(self.tmp / "nope.fasta")))
        empty = self.tmp / "empty.fasta"
        empty.write_bytes(b"")
        self.assertIn("empty input", self.refuses("--mito", str(empty)))
        self.assertFalse((self.tmp / "out" / "work").exists())

    def test_empty_after_decompressing(self):
        with gzip.open(self.tmp / "e.fasta.gz", "wb"):
            pass
        self.assertIn("empty input", self.refuses("--mito", str(self.tmp / "e.fasta.gz")))
        self.assertFalse((self.tmp / "out" / "work").exists())

    def test_existing_work_needs_force(self):
        (self.tmp / "out" / "work").mkdir(parents=True)
        self.assertIn("--force", self.refuses("--mito", str(self.fa)))


class StageInput(unittest.TestCase):
    def test_gz_is_decompressed(self):
        tmp = Path(tempfile.mkdtemp())
        with gzip.open(tmp / "s.fasta.gz", "wb") as fh:
            fh.write(FASTA)
        run.stage_input(str(tmp / "s.fasta.gz"), tmp / "out.fasta")
        self.assertEqual((tmp / "out.fasta").read_bytes(), FASTA)


if __name__ == "__main__":
    unittest.main()
