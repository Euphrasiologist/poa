"""06_low_depth_paths' gene-calls fallback (no oatk .annot_mito.txt): which
unitigs carry a gene, given poa's calls on the contigs and unitig_coords'
map. On the dataset it agrees with oatk's own unitig hits on every flag
(956 species).

    python3 -m unittest tests.test_low_depth_paths
"""
from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "low_depth_paths", REPO / "src/poa/pipeline/qc_basic_stats/src/06_low_depth_paths.py")
ldp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ldp)

HEADER = "unitig\tunitig_length\tctg_seqid\tctg_start\tctg_end\tstrand\tsource\tmatch_trim\n"


def write_map(rows: list[str]) -> Path:
    path = Path(tempfile.mkdtemp()) / "sp.mito.tsv"
    path.write_text(HEADER + "".join(r + "\n" for r in rows))
    return path


class GeneUnitigsFromCalls(unittest.TestCase):
    def test_linear_placements(self):
        m = write_map(["u1\t1000\tctg1\t0\t1000\t+\tpath\t0",
                       "u2\t1000\tctg1\t990\t1990\t+\tpath\t0"])
        # mostly on u2, only 10bp on u1
        self.assertEqual(ldp.gene_unitigs_from_calls([("ctg1", 980, 1500)], m, {"ctg1": 1990}), {"u2"})
        # straddling both evenly: both carry it
        self.assertEqual(ldp.gene_unitigs_from_calls([("ctg1", 700, 1300)], m, {"ctg1": 1990}), {"u1", "u2"})

    def test_single_unitig_circle_wraps_origin(self):
        # Achillea_maritima's shape: one circular unitig, placed past the contig end
        m = write_map(["u6\t202131\tctg000001c\t202059\t404190\t-\tpath\t0"])
        calls = [("ctg000001c", 90819, 92375)]
        self.assertEqual(ldp.gene_unitigs_from_calls(calls, m, {"ctg000001c": 202095}), {"u6"})

    def test_other_contig_and_no_map(self):
        m = write_map(["u1\t1000\tctg1\t0\t1000\t+\tpath\t0"])
        self.assertEqual(ldp.gene_unitigs_from_calls([("ctg2", 0, 500)], m, {}), set())
        self.assertEqual(ldp.gene_unitigs_from_calls([("ctg1", 0, 500)], m.with_name("none.tsv"), {}), set())


if __name__ == "__main__":
    unittest.main()
