"""RefSeq profile fallback in editing.scan and the pinned lineage gene lists."""

import tempfile
import unittest
from pathlib import Path

from poa.pipeline.common import lineage as lin
from poa.pipeline.editing import lineage_profiles as lp
from poa.pipeline.editing.scan import profile_paths


class ProfilePaths(unittest.TestCase):
    def test_first_directory_wins_and_others_fill_in(self):
        with tempfile.TemporaryDirectory() as tmp:
            dataset, refseq = Path(tmp) / "dataset", Path(tmp) / "refseq"
            for d, genes in ((dataset, ["atp4", "nad9"]), (refseq, ["atp4", "rps12"])):
                (d / "mito").mkdir(parents=True)
                for g in genes:
                    (d / "mito" / f"{g}.pssm").write_text("")
            found = profile_paths([dataset, refseq], "mito")
            self.assertEqual(sorted(found), ["atp4", "nad9", "rps12"])
            self.assertEqual(found["atp4"], dataset / "mito" / "atp4.pssm")
            self.assertEqual(found["rps12"], refseq / "mito" / "rps12.pssm")
            self.assertEqual(profile_paths([dataset, refseq], "pltd"), {})


class BundledProfiles(unittest.TestCase):
    def test_moss_genes_pinned(self):
        # the lineage call was calibrated on these; adding genes changes it
        self.assertEqual(lp.MOSS_GENES, ["atp4", "ccmB", "ccmC", "mttB", "nad5", "nad9", "rps14"])
        moss = sorted(p.stem for p in (lin.editing_profiles(lp.PIPELINE_DIR, "moss") / "mito").glob("*.pssm"))
        self.assertEqual(moss, ["atp4", "ccmB", "ccmC", "mttB", "rps14"])

    def test_refseq_profiles_bundled_for_single_exon_genes(self):
        bundled = {p.stem for p in (lp.refseq_profiles_dir() / "mito").glob("*.pssm")}
        self.assertEqual(len(bundled), 26)
        self.assertTrue({"atp6", "atp9", "cox1", "rps12"} <= bundled)
        self.assertTrue(bundled <= set(lp.REFSEQ_GENES))
