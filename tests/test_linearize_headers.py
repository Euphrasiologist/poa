"""gfatk resolve's circuit headers through linearize (04 summarise, 05
promote, 06 re-annotate) and into gff_export's region feature, for gfatk
0.6 and 0.7 headers.

    python3 -m unittest tests.test_linearize_headers
"""

from __future__ import annotations

import importlib.util
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from poa.pipeline.gff_export import build_gff as bg  # noqa: E402

SRC = REPO / "src" / "poa" / "pipeline" / "linearize" / "src"
PATH = "u0+,u10+,u14+,u4-,u14-,u10-"
V063 = f"Sp.pltd_circuit1:n_segments=4:bp=120:circular=true:path={PATH}"
V070 = f"Sp.pltd_circuit1:n_segments=4:bp=120:circular=true:repeat_runs=2:unspanned_repeat_runs=1:unspanned_at=3:path={PATH}"
ARM = "bubble_arm:u11:between=u10-_to_u23-:support_reads=131:strand=+:edit_distance=8:pct_identity=99.6"


def script(name: str):
    spec = importlib.util.spec_from_file_location(f"linearize_{name[:2]}", SRC / name)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class CircuitHeaders(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.s04, cls.s05, cls.s06 = (script(n) for n in ("04_summarize_resolve.py", "05_promote_resolve.py", "06_reannotate_promoted.py"))

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.tmp)

    def resolve_fasta(self, header: str) -> Path:
        p = self.tmp / "Sp.resolve.fasta"
        p.write_text(f">{header}\n{'ACGT' * 30}\n>{ARM}\nACGT\n")
        return p

    def test_both_versions_parse_in_every_step(self):
        for header in (V063, V070):
            with self.subTest(header=header[:60]):
                f = self.resolve_fasta(header)
                rows, n_arms = self.s04.parse_resolve_fasta(f)
                self.assertEqual((len(rows), n_arms), (1, 1))
                circuits = self.s05.parse_circuits(f)
                self.assertEqual([(c["n_segments"], c["bp"]) for c in circuits], [(4, 120)])
                out = self.tmp / "promoted.ctg.fasta"
                self.s05.write_promoted_fasta(out, circuits, "0.7.0")
                promoted = out.read_text().splitlines()[0]
                self.assertEqual(self.s06.HEADER_PATH_RE.search(promoted).group(1), PATH)
                ctx = bg.load_seq_context(out, self.tmp / "no_map.tsv")
                self.assertEqual((ctx.resolver["ctg000001c"], ctx.resolver_version["ctg000001c"]), ("gfatk_resolve", "0.7.0"))
                self.assertTrue(ctx.circular["ctg000001c"])

    def test_the_version_comes_from_beside_the_resolve_fasta(self):
        f = self.resolve_fasta(V070)
        self.assertIsNone(self.s05.resolve_version(f))
        (self.tmp / "Sp.resolve.version").write_text("0.7.0\n")
        self.assertEqual(self.s05.resolve_version(f), "0.7.0")
        out = self.tmp / "old.ctg.fasta"
        self.s05.write_promoted_fasta(out, self.s05.parse_circuits(f))  # promoted before versions were recorded
        self.assertNotIn("gfatk_version", out.read_text())
        self.assertNotIn("ctg000001c", bg.load_seq_context(out, self.tmp / "no_map.tsv").resolver_version)


if __name__ == "__main__":
    unittest.main()
