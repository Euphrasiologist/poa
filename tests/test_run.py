"""poa run's input handling: what it accepts, what it refuses before running
any stage, and how it lays out the data root. The full run is covered by
tests/smoke/run_smoke.sh.

    python3 -m unittest tests.test_run
"""

from __future__ import annotations

import gzip
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from click.testing import CliRunner  # noqa: E402

from poa import cli, run  # noqa: E402

FASTA = b">ctg1\nACGT\n"


class DefaultName(unittest.TestCase):
    def test_oatk_names(self):
        self.assertEqual(run.default_name("x/Arabidopsis_thaliana.mito.ctg.fasta"), "Arabidopsis_thaliana")
        self.assertEqual(
            run.default_name("Arabidopsis_thaliana.k1001.s31.c100.pltd.ctg.fasta.gz"), "Arabidopsis_thaliana"
        )

    def test_unsafe_characters(self):
        self.assertEqual(run.default_name("my sample+1.fa"), "my_sample_1")


class Refusals(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.tmp)
        self.fa = self.tmp / "s.fasta"
        self.fa.write_bytes(FASTA)

    def refuses(self, **kwargs) -> str:
        with self.assertRaises(SystemExit) as cm:
            run.run(outdir=self.tmp / "out", **kwargs)
        return str(cm.exception.code)

    def test_no_fasta(self):
        self.assertIn("--mito and/or --pltd", self.refuses())

    def test_gfa_without_fasta(self):
        self.assertIn("--mito-gfa needs --mito", self.refuses(pltd=str(self.fa), mito_gfa=str(self.fa)))

    def test_bad_name(self):
        self.assertIn("--name", self.refuses(mito=str(self.fa), name="a.b"))

    def test_missing_and_empty_inputs(self):
        self.assertIn("no such file", self.refuses(mito=str(self.tmp / "nope.fasta")))
        empty = self.tmp / "empty.fasta"
        empty.write_bytes(b"")
        self.assertIn("empty input", self.refuses(mito=str(empty)))
        self.assertFalse((self.tmp / "out" / "work").exists())

    def test_empty_after_decompressing(self):
        with gzip.open(self.tmp / "e.fasta.gz", "wb"):
            pass
        self.assertIn("empty input", self.refuses(mito=str(self.tmp / "e.fasta.gz")))
        self.assertFalse((self.tmp / "out" / "work").exists())

    def test_existing_work_needs_force(self):
        (self.tmp / "out" / "work").mkdir(parents=True)
        self.assertIn("--force", self.refuses(mito=str(self.fa)))


class Cli(unittest.TestCase):
    """The click layer: options reach run.run as keywords, and bad
    options are usage errors (exit 2) before anything runs."""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.tmp)
        self.fa = self.tmp / "s.fasta"
        self.fa.write_bytes(FASTA)
        self.calls = []
        real = run.run
        run.run = lambda **kw: self.calls.append(kw) or 0
        self.addCleanup(setattr, run, "run", real)

    def invoke(self, *argv: str):
        return CliRunner().invoke(cli.main, list(argv))

    def test_options_reach_run(self):
        res = self.invoke("run", "--mito", str(self.fa), "-o", str(self.tmp / "out"), "-t", "8", "--lineage", "moss")
        self.assertEqual(res.exit_code, 0, res.output)
        self.assertEqual(
            self.calls,
            [
                dict(
                    mito=str(self.fa),
                    mito_gfa=None,
                    pltd=None,
                    pltd_gfa=None,
                    outdir=str(self.tmp / "out"),
                    name=None,
                    lineage="moss",
                    threads=8,
                    force=False,
                )
            ],
        )

    def test_usage_errors(self):
        out = str(self.tmp / "out")
        for argv in (
            ("run", "--mito", str(self.fa)),  # no -o
            ("run", "--mito", str(self.tmp / "nope.fasta"), "-o", out),
            ("run", "--mito", str(self.fa), "-o", out, "--lineage", "fern"),
            ("run", "--mito", str(self.fa), "-o", out, "-t", "0"),
        ):
            self.assertEqual(self.invoke(*argv).exit_code, 2, argv)
        self.assertEqual(self.calls, [])

    def test_help_keeps_layout(self):
        res = self.invoke("run", "--help")
        self.assertEqual(res.exit_code, 0)
        self.assertIn("  denovo_annotation 00-04 -> editing 01", res.output)
        self.assertIn("  OUTDIR/logs/<stage>.log ", res.output)

    def test_version(self):
        self.assertIn(f"poa {cli.__version__}", self.invoke("--version").output)


class StageInput(unittest.TestCase):
    def test_gz_is_decompressed(self):
        tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, tmp)
        with gzip.open(tmp / "s.fasta.gz", "wb") as fh:
            fh.write(FASTA)
        run.stage_input(str(tmp / "s.fasta.gz"), tmp / "out.fasta")
        self.assertEqual((tmp / "out.fasta").read_bytes(), FASTA)


class InProcessStage(unittest.TestCase):
    """A stage's output - its own and its child processes' - goes to its
    log, and a failure becomes StageFailed."""

    def setUp(self):
        self.logs = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.logs)

    def log(self, label):
        return (self.logs / f"{label}.log").read_text()

    def test_output_captured(self):
        def stage():
            print("from python", file=sys.stderr)
            subprocess.run(["echo", "from a child"], check=True)

        run.run_stage("ok", stage, self.logs)
        self.assertEqual(self.log("ok").split("\n")[:2], ["from python", "from a child"])

    def test_stage_error(self):
        def stage():
            raise run.StageError("orfedit scan-batch failed (exit 2)")

        with self.assertRaises(run.StageFailed) as cm:
            run.run_stage("bad", stage, self.logs)
        self.assertIn("[err] orfedit scan-batch failed", str(cm.exception))

    def test_bug_gives_traceback(self):
        def stage():
            return {}["missing"]

        with self.assertRaises(run.StageFailed):
            run.run_stage("bug", stage, self.logs)
        self.assertIn("KeyError: 'missing'", self.log("bug"))


if __name__ == "__main__":
    unittest.main()
