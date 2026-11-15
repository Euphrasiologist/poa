"""trans_splicing.reconstruct's manifest: which exon templates each species'
fragments are scanned against.

    python3 -m unittest tests.test_reconstruct
"""

from __future__ import annotations

import shutil
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from poa import PIPELINE_DIR  # noqa: E402
from poa.pipeline.common.errors import StageError  # noqa: E402
from poa.pipeline.trans_splicing import reconstruct as rc  # noqa: E402

REFERENCE = PIPELINE_DIR / "trans_splicing" / "reference"


class Manifest(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.root)
        for sp in ("Moss", "Rose"):
            d = self.root / "data" / "mito" / sp
            d.mkdir(parents=True)
            (d / f"{sp}.mito.ctg.fasta").write_text(">c1\nACGT\n")
        rows = [(sp, "mito", g, "c1", 1, 100, "+") for sp in ("Moss", "Rose") for g in ("nad1", "rps10", "cox1")]
        rows.append(("Ghost", "mito", "nad1", "c1", 1, 100, "+"))  # no assembly: skipped
        self.calls = pd.DataFrame(rows, columns=["species", "organelle", "gene", "contig_id", "start", "end", "strand"])

    def test_templates_follow_lineage(self):
        rows = rc.build_manifest(self.root, self.calls, None, {"Moss": "moss"}, "angiosperm")
        got = {(r["species"], r["gene"]): r["template_dir"] for r in rows}
        # cox1 isn't a reconstructed gene; moss has no rps10 template
        self.assertEqual(set(got), {("Moss", "nad1"), ("Rose", "nad1"), ("Rose", "rps10")})
        self.assertEqual(got["Moss", "nad1"], str(REFERENCE / "lineages" / "moss" / "profiles" / "nad1"))
        self.assertEqual(got["Rose", "nad1"], str(REFERENCE / "profiles" / "nad1"))

    def test_species_filter_and_forced_lineage(self):
        rows = rc.build_manifest(self.root, self.calls, {"Rose"}, {}, "moss")
        self.assertEqual({(r["species"], r["gene"]) for r in rows}, {("Rose", "nad1")})

    def test_missing_gene_calls(self):
        with self.assertRaises(StageError):
            rc.reconstruct(self.root, gene_calls=self.root / "nope.tsv")


if __name__ == "__main__":
    unittest.main()
