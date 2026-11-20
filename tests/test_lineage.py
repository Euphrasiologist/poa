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

    def test_fern_borrows_the_default_references(self):
        # fern: no references of its own, never chosen automatically, and a
        # stop read-through cost
        self.assertEqual(lin.reference("fern"), lin.DEFAULT_LINEAGE)
        self.assertEqual(lin.reference("moss"), "moss")
        self.assertIn("fern", lin.CHOICES)
        self.assertNotIn("fern", lin.LINEAGES)
        self.assertEqual(lin.stop_edit_penalty("fern"), 2.5)
        self.assertIsNone(lin.stop_edit_penalty("angiosperm"))
        self.assertIsNone(lin.stop_edit_penalty("moss"))
        self.assertEqual(lin.choose(scores(100, 900))[0], "moss")  # auto never says fern

    def test_forced_fern_keeps_the_default_run(self):
        # a species forced to fern takes its mito rows from the default
        # (angiosperm-profile) run, not from another lineage's
        df = pd.DataFrame({"species": ["F", "F", "A"], "organelle": ["mito", "pltd", "mito"], "gene": ["atp4", "rbcL", "atp4"]})
        lineages = {"F": "fern", "A": "angiosperm"}
        kept = scan_editing.keep_lineage(df, lineages, lin.DEFAULT_LINEAGE)
        self.assertEqual(len(kept), 3)
        self.assertTrue(scan_editing.keep_lineage(df, lineages, "moss").empty)

    def test_forced(self):
        table = scan_editing.assign_lineages({"angiosperm": self.angio, "moss": self.moss}, "angiosperm")
        self.assertEqual(set(table["lineage"]), {"angiosperm"})
        self.assertEqual(set(table["decided_by"]), {"forced"})


class FernReferences(unittest.TestCase):
    def test_verify_allows_reverse_edits_only_when_asked(self):
        from poa.pipeline.common.refseq_mito import editing_explains, verify

        # ATG AAA TAA GTT: the annotated protein reads the TAA as Gln
        dna, protein = "ATGAAATAAGTT", "MKQV"
        self.assertFalse(verify(dna, protein)[0])
        ok, note = verify(dna, protein, reverse=True)
        self.assertTrue(ok, note)
        self.assertTrue(editing_explains("TGA", "R", reverse=True))
        self.assertFalse(editing_explains("TGA", "R"))
        # C-to-U still works either way, and nothing else is explained
        self.assertTrue(editing_explains("CCA", "L"))
        self.assertFalse(editing_explains("AAA", "E", reverse=True))

    def test_column_mask(self):
        from poa.pipeline.trans_splicing import column_mask as cm

        # column 0 unaligned in all, column 3 in 4 of 5 (80%), column 5 in 3 of 5
        proteins = ["-KV-F-", "-KV-F-", "-KV-F-", "-KV-FA", "-KVLFA"]
        self.assertEqual(cm.unsupported_columns(proteins), [0, 3])
        self.assertEqual(cm.unsupported_columns([]), [])
        tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, tmp)
        pssm = tmp / "whole_gene.pssm"
        pssm.write_text("A\tR\n" + "".join(f"{i}.0\t0.0\n" for i in range(6)))
        self.assertEqual(cm.apply_mask(pssm, {0, 3}), 2)
        self.assertEqual([line.split("\t")[0] for line in pssm.read_text().splitlines()[1:]], ["1.0", "2.0", "4.0", "5.0"])

    def test_fern_templates_are_bundled_and_used(self):
        from poa.pipeline.trans_splicing import reconstruct
        from poa.pipeline.trans_splicing.reference_genes import LINEAGE_MITOGENOMES

        self.assertIn("fern", LINEAGE_MITOGENOMES)
        profiles = lin.trans_splicing_reference(PIPELINE, "fern") / "profiles"
        self.assertTrue((profiles / "rps10" / "exon1.pssm").exists())
        self.assertFalse((profiles / "rps10" / "exon2.pssm").exists())  # one exon in ferns
        self.assertFalse((profiles / "ccmFc").exists())  # lost from fern mitogenomes
        # fern: its own exon templates, the angiosperm editing profiles
        self.assertEqual(reconstruct.template_lineage("fern"), "fern")
        self.assertEqual(lin.reference("fern"), lin.DEFAULT_LINEAGE)


class FernEditingProfiles(unittest.TestCase):
    def test_bundled_and_gated(self):
        from poa.pipeline.editing import lineage_profiles as lp

        mito = lin.editing_profiles(PIPELINE, "fern") / "mito"
        have = {p.stem for p in mito.glob("*.pssm")}
        self.assertIn("rps4", have)
        self.assertIn("atp4", have)
        # left out: their fern profiles did worse than the angiosperm ones
        self.assertFalse({"matR", "nad9"} & have)
        self.assertFalse({"matR", "nad9"} & set(lp.LINEAGE_GENES["fern"]))
        # multi-exon in ferns, so no single-locus profile
        self.assertFalse({"atp1", "cox1", "rps14"} & have)
        self.assertEqual(lp.LINEAGE_MIN_SEQUENCES["fern"], 4)

    def test_one_per_species(self):
        from poa.pipeline.editing import lineage_profiles as lp

        proteins = {"nad3": {"A1": "MK", "A2": "MK", "A3": "ML", "B1": "MK"}}
        species = {"A1": "a", "A2": "a", "A3": "a", "B1": "b"}
        self.assertEqual(list(lp.one_per_species(proteins, species)["nad3"]), ["A1", "A3", "B1"])


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
