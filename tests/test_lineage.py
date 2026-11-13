"""The per-species reference-lineage call (common/lineage.py) and how
editing/01_scan_editing.py splits its orfedit runs by it.

    python3 -m unittest tests.test_lineage
"""

from __future__ import annotations

import shutil
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
PIPELINE = REPO / "src" / "poa" / "pipeline"
sys.path.insert(0, str(REPO / "src"))
from poa.pipeline.common import lineage as lin  # noqa: E402
from poa.pipeline.editing import scan as scan_editing  # noqa: E402

GENES = ["atp4", "ccmB", "ccmC", "mttB", "rps14"]


def scores(angio: float, moss: float, genes=GENES) -> dict:
    return {"angiosperm": {g: angio for g in genes}, "moss": {g: moss for g in genes}}


class Choose(unittest.TestCase):
    def test_clear_moss(self):
        lineage, n, totals = lin.choose(scores(300, 600))
        self.assertEqual((lineage, n), ("moss", 5))
        self.assertEqual(totals, {"angiosperm": 1500, "moss": 3000})

    def test_clear_angiosperm(self):
        self.assertEqual(lin.choose(scores(600, 300))[0], "angiosperm")

    def test_narrow_win_stays_angiosperm(self):
        self.assertEqual(lin.choose(scores(500, 500 + lin.MIN_MARGIN - 1))[0], "angiosperm")
        self.assertEqual(lin.choose(scores(500, 500 + lin.MIN_MARGIN))[0], "moss")

    def test_too_few_shared_genes(self):
        few = GENES[: lin.MIN_SHARED_GENES - 1]
        self.assertEqual(
            lin.choose(scores(100, 900, few)),
            ("angiosperm", len(few), {"angiosperm": 100 * len(few), "moss": 900 * len(few)}),
        )

    def test_only_shared_genes_count(self):
        s = scores(300, 600)
        s["angiosperm"]["nad5"] = 10_000  # no moss profile for nad5: ignored
        self.assertEqual(lin.choose(s)[:2], ("moss", 5))

    def test_no_moss_run(self):
        self.assertEqual(lin.choose({"angiosperm": {g: 1 for g in GENES}})[0], "angiosperm")


class SplitRuns(unittest.TestCase):
    def setUp(self):
        def calls(rows):
            return pd.DataFrame(rows, columns=["species", "organelle", "gene", "score"]).astype(str)

        self.angio = calls(
            [("Moss", "mito", g, 300) for g in GENES + ["nad5"]]
            + [("Moss", "pltd", "psbA", 900)]
            + [("Rose", "mito", g, 600) for g in GENES + ["nad5"]]
        )
        self.moss = calls([("Moss", "mito", g, 600) for g in GENES] + [("Rose", "mito", g, 200) for g in GENES])

    def test_assign_and_keep(self):
        table = scan_editing.assign_lineages({"angiosperm": self.angio, "moss": self.moss}, None)
        lineages = dict(zip(table["species"], table["lineage"]))
        self.assertEqual(lineages, {"Moss": "moss", "Rose": "angiosperm"})
        self.assertEqual(list(table.columns), lin.COLUMNS)
        kept_a = scan_editing.keep_lineage(self.angio, lineages, "angiosperm")
        kept_m = scan_editing.keep_lineage(self.moss, lineages, "moss")
        # the moss keeps its plastid call from the default run, and loses nad5
        # (no moss profile) rather than keeping the angiosperm correction
        self.assertEqual(sorted(kept_a[kept_a.species == "Moss"].gene), ["psbA"])
        self.assertEqual(sorted(kept_m.gene), sorted(GENES))
        self.assertEqual(set(kept_a[kept_a.species == "Rose"].gene), set(GENES + ["nad5"]))

    def test_forced(self):
        table = scan_editing.assign_lineages({"angiosperm": self.angio, "moss": self.moss}, "angiosperm")
        self.assertEqual(set(table["lineage"]), {"angiosperm"})
        self.assertEqual(set(table["decided_by"]), {"forced"})


class Load(unittest.TestCase):
    def test_missing_and_present(self):
        root = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, root)
        self.assertIsNone(lin.load(root))
        lin.table_path(root).parent.mkdir(parents=True)
        lin.table_path(root).write_text("species\tlineage\nA\tmoss\n")
        self.assertEqual(lin.load(root), {"A": "moss"})


class References(unittest.TestCase):
    """The bundled moss references exist where the stages look for them."""

    def test_bundled(self):
        ts = lin.trans_splicing_reference(PIPELINE, "moss") / "profiles"
        self.assertTrue((ts / "nad1" / "exon3.pssm").exists())
        self.assertFalse((ts / "nad1" / "exon4.pssm").exists())  # moss nad1 has 3 exons
        self.assertFalse((ts / "rps10").exists())  # lost from moss mitogenomes
        ed = lin.editing_profiles(PIPELINE, "moss") / "mito"
        self.assertEqual(sorted(p.stem for p in ed.glob("*.pssm")), sorted(GENES))
        self.assertEqual(
            lin.trans_splicing_reference(PIPELINE, "angiosperm"), PIPELINE / "trans_splicing" / "reference"
        )


if __name__ == "__main__":
    unittest.main()
