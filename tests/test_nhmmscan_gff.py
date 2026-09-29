"""nhmmscan_gff.py must reproduce the Rust filter_tblout + hmm_to_gff output
it replaced. Fixtures: Arabidopsis nhmmscan tblouts and the GFFs the Rust
tools made from them (fixtures/nhmmscan/*.rust.gff).

    python3 -m unittest tests.test_nhmmscan_gff
"""
from __future__ import annotations

import os
import sys
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
FIXTURES = REPO / "tests" / "fixtures" / "nhmmscan"
sys.path.insert(0, str(REPO / "analysis" / "denovo_annotation" / "src"))
from nhmmscan_gff import filter_and_convert, rust_exp2, rust_f32  # noqa: E402

import numpy as np  # noqa: E402

FAM = {
    "mito": os.environ.get("OATKDB_MITO_FAM", "/software/team301/OatkDB/viridiplantae_mito_v20250217.fam"),
    "pltd": os.environ.get("OATKDB_PLTD_FAM", "/software/team301/OatkDB/viridiplantae_pltd_v20260928.fam"),
}


class RustFormatting(unittest.TestCase):
    def test_f32_display(self):
        self.assertEqual(rust_f32(np.float32("1705.0")), "1705")
        self.assertEqual(rust_f32(np.float32("44.5")), "44.5")
        self.assertEqual(rust_f32(np.float32("3843.1")), "3843.1")

    def test_exp2(self):
        self.assertEqual(rust_exp2(np.float32("0")), "0.00e0")
        self.assertEqual(rust_exp2(np.float32("3.4e-12")), "3.40e-12")
        self.assertEqual(rust_exp2(np.float32("1e-5")), "1.00e-5")


class MatchesRustTools(unittest.TestCase):
    def check(self, organelle: str):
        fam = Path(FAM[organelle])
        if not fam.exists():
            self.skipTest(f"no {organelle} .fam at {fam} (set OATKDB_{organelle.upper()}_FAM)")
        tblout = FIXTURES / f"Arabidopsis_thaliana.{organelle}.tblout"
        filtered, gff = filter_and_convert(tblout, fam, "1e-5")
        self.assertEqual(gff, (FIXTURES / f"Arabidopsis_thaliana.{organelle}.rust.gff").read_text())
        # every hit in the GFF is a filtered-tblout row
        n_rows = sum(1 for line in filtered.splitlines() if line.strip() and not line.startswith("#"))
        self.assertEqual(n_rows, gff.count("\n") - 1)

    def test_mito(self):
        self.check("mito")

    def test_pltd(self):
        self.check("pltd")


if __name__ == "__main__":
    unittest.main()
