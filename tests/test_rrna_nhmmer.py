"""rrna_nhmmer.py must reproduce `barrnap --kingdom plant` (barrnap 0.9 with
the plant.hmm database poa now bundles). Fixtures: what that barrnap wrote
for the Arabidopsis assemblies in tests/data (fixtures/rrna/*.barrnap.gff).

    python3 -m unittest tests.test_rrna_nhmmer
"""
from __future__ import annotations

import gzip
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src" / "poa" / "pipeline" / "denovo_annotation" / "src"))
from rrna_nhmmer import call_rrna, tblout_to_gff  # noqa: E402

ASSEMBLIES = {
    "mito": REPO / "tests/data/mito/Arabidopsis_thaliana/Arabidopsis_thaliana.mito.ctg.fasta.gz",
    "pltd": REPO / "tests/data/plastid/Arabidopsis_thaliana/Arabidopsis_thaliana.k1001.s31.c100.pltd.ctg.fasta.gz",
}


class Parsing(unittest.TestCase):
    def test_partial_and_reject(self):
        row = "ctg1 - {g} - 1 100 {a} {b} 1 100 9999 + {e} 50.0 0.1 -"
        tbl = "\n".join([
            row.format(g="16S_rRNA", a=1, b=1585, e="0"),            # full length
            row.format(g="16S_rRNA", a=3000, b=2001, e="1.2e-30"),   # 1000 bp, minus strand: partial
            row.format(g="5S_rRNA", a=10, b=30, e="0.001"),          # 21 bp < int(0.25 x 119): rejected
        ])
        gff = tblout_to_gff(tbl).splitlines()
        self.assertEqual(gff[0], "##gff-version 3")
        self.assertEqual(gff[1], "ctg1\tbarrnap:0.9\trRNA\t1\t1585\t0\t+\t.\tName=16S_rRNA;product=16S ribosomal RNA")
        self.assertEqual(gff[2], "ctg1\tbarrnap:0.9\trRNA\t2001\t3000\t1.2e-30\t-\t.\t"
                                 "Name=16S_rRNA;product=16S ribosomal RNA (partial);"
                                 "note=aligned only 63 percent of the 16S ribosomal RNA")
        self.assertEqual(len(gff), 3)


class MatchesBarrnap(unittest.TestCase):
    def check(self, organelle: str):
        if shutil.which("nhmmer") is None:
            self.skipTest("nhmmer (HMMER) not on PATH")
        with tempfile.NamedTemporaryFile(suffix=".fa") as fa:
            fa.write(gzip.decompress(ASSEMBLIES[organelle].read_bytes()))
            fa.flush()
            gff = call_rrna(fa.name, "1e-5", 1)
        expected = (REPO / "tests/fixtures/rrna" / f"Arabidopsis_thaliana.{organelle}.barrnap.gff").read_text()
        self.assertEqual(gff, expected)

    def test_mito(self):
        self.check("mito")

    def test_pltd(self):
        self.check("pltd")


if __name__ == "__main__":
    unittest.main()
