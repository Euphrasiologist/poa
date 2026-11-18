"""Frameshift recurrence across the DToL dataset."""

import tempfile
import unittest
from pathlib import Path
from unittest import mock

import pandas as pd

from poa.pipeline.editing import scan


class Recurrence(unittest.TestCase):
    def test_counts_assemblies_within_the_window(self):
        with tempfile.TemporaryDirectory() as tmp:
            ref = Path(tmp) / "ref.tsv"
            pd.DataFrame(
                {"organelle": ["mito"] * 3, "gene": ["rpl16"] * 3, "profile_column": [20, 24, 40],
                 "n_assemblies": [300, 17, 2]}  # fmt: skip
            ).to_csv(ref, sep="\t", index=False)
            shifts = pd.DataFrame(
                {"species": ["A_b", "A_b", "M_n"], "organelle": ["mito"] * 3, "gene": ["rpl16", "cox1", "rpl16"],
                 "profile_column": [22, 100, 22]}  # fmt: skip
            )
            with mock.patch.object(scan, "FRAMESHIFT_REFERENCE", ref):
                out = scan.add_recurrence(shifts, {"A_b": "angiosperm", "M_n": "moss"})
            # columns 20 and 24 are within 5 of 22; cox1 has none; a moss
            # gene's profile columns aren't the reference's
            self.assertEqual(list(out.dtol_assemblies), [317, 0, ""])


if __name__ == "__main__":
    unittest.main()
