"""Profile building: pseudocounts shared out by the background, and the
editing sites a profile takes from reference proteins that mark them."""

import tempfile
import unittest
from pathlib import Path

from poa.pipeline.common import pssm


class Pseudocounts(unittest.TestCase):
    def test_a_residue_the_proteins_never_use_does_not_score_as_a_match(self):
        # 30 copies of a protein without Cys or Trp: with a flat pseudocount per
        # residue their background was far below their pseudocount, and both
        # scored +5 bits or more in every column
        alignment = {f"s{i}": "MKLNTEKEIIVTSRASTIIPMIGHTIAIHNGREHLPVYITDRMVGHK" for i in range(30)}
        columns = pssm.build_profile(alignment)
        for col in columns:
            self.assertLess(col["C"], 0)
            self.assertLess(col["W"], 0)
        # and the residue a column does have is still its best
        self.assertEqual(max(pssm.AA_ORDER, key=lambda aa: columns[0][aa]), "M")

    def test_a_missing_residue_scores_the_same_whatever_its_background(self):
        alignment = {f"s{i}": "MKLLLLLLLLW" for i in range(10)}
        col = pssm.build_profile(alignment)[0]
        self.assertAlmostEqual(col["L"], col["C"])  # neither is in column 0; L is common, C absent


class EditingSites(unittest.TestCase):
    def test_mark_edited(self):
        # Ser to Leu, Pro to Leu and a Thr to Met start: C-to-U edits. The stop codon is not a residue
        self.assertEqual(pssm.mark_edited("ACGTCACCAAAATAA", "MLLK"), "mllK")
        # a GTG start read as Met is not an edit, nor is a residue no C-to-U explains
        self.assertEqual(pssm.mark_edited("GTGTCACGA", "MLW"), "MlW")
        # proteins that do not line up with the gene are left unmarked
        self.assertEqual(pssm.mark_edited("ATGTCA", "MLK"), "MLK")
        self.assertEqual(pssm.edited_residues("mLlK"), {0, 2})

    def test_a_site_needs_two_groups(self):
        alignment = {"RNA_A_x": "MLK-F", "RNA_A_z": "MLKSF", "RNA_B_y": "MLKSF", "NC_1": "MLK-F"}
        edited = {"RNA_A_x": {1}, "RNA_A_z": {1, 3}, "RNA_B_y": {1}}
        groups = {"RNA_A_x": "A", "RNA_A_z": "A", "RNA_B_y": "B", "NC_1": "N"}
        columns = pssm.build_profile(alignment, edited=edited, groups=groups)
        # column 1: genera A and B of the three present; column 3: genus A alone (twice)
        self.assertEqual([round(c[pssm.EDIT_COLUMN], 2) for c in columns], [0.0, 0.67, 0.0, 0.0, 0.0])
        # without groups every protein counts for itself
        columns = pssm.build_profile(alignment, edited=edited)
        self.assertEqual([c[pssm.EDIT_COLUMN] > 0 for c in columns], [False, True, False, False, False])

    def test_the_edit_column_is_written_only_for_profiles_with_sites(self):
        alignment = {"a": "MLKSF", "b": "MLKSF"}
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "x.pssm"
            pssm.write_profile(pssm.build_profile(alignment), path)
            self.assertEqual(path.read_text().splitlines()[0].split("\t"), list(pssm.AA_ORDER))
            pssm.write_profile(pssm.build_profile(alignment, edited={"a": {1}, "b": {1}}), path)
            lines = [line.split("\t") for line in path.read_text().splitlines()]
            self.assertEqual(lines[0], [*pssm.AA_ORDER, pssm.EDIT_COLUMN])
            self.assertEqual([row[20] for row in lines[1:]], ["0.000", "1.000", "0.000", "0.000", "0.000"])


if __name__ == "__main__":
    unittest.main()
