"""The other organelle's DNA (denovo_annotation/transfers.py): which cross
hits are kept, what the alignment reports, how both reach the GFF, and how
they are mapped back after poa moves an origin.

    python3 -m unittest tests.test_transfers
"""

from __future__ import annotations

import random
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from poa.pipeline.common import origin  # noqa: E402
from poa.pipeline.denovo_annotation import transfers as T  # noqa: E402
from poa.pipeline.gff_export import build_gff as bg  # noqa: E402


def tbl_line(model, contig, env_from, env_to, strand, evalue, score):
    # nhmmscan --tblout: target acc query acc hmmfrom hmmto alifrom alito envfrom envto modlen strand E score bias desc
    return f"{model} - {contig} - 1 100 {env_from} {env_to} {env_from} {env_to} 100 {strand} {evalue} {score} 0.0 -"


def hits(*rows):
    return pd.DataFrame(
        [dict(model=m, contig=c, start=s, end=e, strand=st, e_value=1e-20, score=sc) for m, c, s, e, st, sc in rows],
        columns=["model", "contig", "start", "end", "strand", "e_value", "score"],
    )


class CrossSearch(unittest.TestCase):
    def test_tblout_is_read_0_based_and_filtered(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "x.tblout"
            p.write_text("# header\n" + tbl_line("psbA", "ctg1", 101, 1160, "+", "1e-50", "800.0") + "\n"
                         + tbl_line("psbA", "ctg1", 2000, 1901, "-", "0.01", "12.0") + "\n")
            t = T.read_tblout(p, 1e-5)
            self.assertEqual(t[["model", "contig", "start", "end", "strand", "score"]].values.tolist(), [["psbA", "ctg1", 100, 1160, "+", 800.0]])
            self.assertEqual(len(T.read_tblout(p)), 2)
            self.assertEqual(T.read_tblout(p).iloc[1][["start", "end", "strand"]].tolist(), [1900, 2000, "-"])

    def test_a_hit_only_the_other_organelles_models_explain_is_kept(self):
        cross = hits(("psbA", "c", 100, 1100, "+", 800.0))
        kept = T.foreign_hits(cross, hits(("cox1", "c", 5000, 6000, "+", 900.0)))
        self.assertEqual(kept[["model", "own_best_score"]].values.tolist(), [["psbA", 0.0]])

    def test_an_own_model_as_good_or_nearly_as_good_on_the_stretch_drops_it(self):
        cross = hits(("rrn16", "c", 100, 1500, "+", 300.0), ("trnS-GGA", "c", 3000, 3087, "+", 97.4))
        own = hits(("rrn18", "c", 150, 1600, "+", 900.0), ("trnS-GGA", "c", 3000, 3087, "+", 96.9))
        self.assertTrue(T.foreign_hits(cross, own).empty)

    def test_an_own_hit_on_the_other_strand_or_elsewhere_does_not_count(self):
        cross = hits(("psbA", "c", 100, 1100, "+", 800.0))
        own = hits(("orf", "c", 100, 1100, "-", 900.0), ("orf", "c", 1000, 2000, "+", 900.0))  # 100 bp of 1000 overlap
        self.assertEqual(len(T.foreign_hits(cross, own)), 1)

    def test_several_models_on_one_stretch_leave_the_best(self):
        cross = hits(("trnM-CAU", "c", 10, 83, "+", 60.0), ("trnfM-CAU", "c", 10, 85, "+", 82.0), ("psbA", "c", 500, 1500, "+", 700.0))
        self.assertEqual(sorted(T.foreign_hits(cross, hits()).model), ["psbA", "trnfM-CAU"])


def mutate(seq, rate, rng):
    return "".join(rng.choice("ACGT".replace(b, "")) if rng.random() < rate else b for b in seq)


def revcomp(seq):
    return seq[::-1].translate(str.maketrans("ACGT", "TGCA"))


@unittest.skipUnless(T.blast_exe("blastn") and T.blast_exe("makeblastdb"), "BLAST+ not installed")
class CrossAlign(unittest.TestCase):
    def test_a_plastid_stretch_in_the_mitogenome_is_reported_in_both(self):
        rng = random.Random(1)
        pltd = "".join(rng.choice("ACGT") for _ in range(4000))
        piece = revcomp(mutate(pltd[1000:1800], 0.05, rng))
        mito = "".join(rng.choice("ACGT") for _ in range(3000)) + piece + "".join(rng.choice("ACGT") for _ in range(3000))
        d = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, d)
        (d / "m.fa").write_text(f">m1\n{mito}\n")
        (d / "p.fa").write_text(f">p1\n{pltd}\n")
        rows = pd.DataFrame(T.cross_align_one("s", str(d / "m.fa"), str(d / "p.fa"), 1))
        self.assertEqual(sorted(rows.organelle), ["mito", "pltd"])
        m = rows[rows.organelle == "mito"].iloc[0]
        self.assertEqual((m.contig, m.strand, m.other_contig, m.other_organelle), ("m1", "-", "p1", "pltd"))
        self.assertLess(abs(m.start - 3000), 30)
        self.assertLess(abs(m.end - 3800), 30)
        self.assertLess(abs(m.other_start - 1000), 30)
        self.assertGreater(m.identity, 0.9)
        p = rows[rows.organelle == "pltd"].iloc[0]
        self.assertEqual((p.contig, p.start, p.end, p.other_start, p.other_end), ("p1", m.other_start, m.other_end, m.start, m.end))

    def test_unrelated_genomes_share_nothing(self):
        rng = random.Random(2)
        d = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, d)
        (d / "m.fa").write_text(">m1\n" + "".join(rng.choice("ACGT") for _ in range(5000)) + "\n")
        (d / "p.fa").write_text(">p1\n" + "".join(rng.choice("ACGT") for _ in range(5000)) + "\n")
        self.assertEqual(T.cross_align_one("s", str(d / "m.fa"), str(d / "p.fa"), 1), [])


ROWS = pd.DataFrame(
    [
        dict(species="s", organelle="mito", contig="m1", start=99, end=1160, strand="+", method="cross_search", name="psbA",
             score=812.3, e_value=1.5e-60, coverage=0.97, identity=None, aligned_bp=None, own_best_score=0.0,
             other_organelle="pltd", other_contig=None, other_start=None, other_end=None),
        dict(species="s", organelle="mito", contig="m1", start=50, end=2050, strand="-", method="cross_align",
             name="shared with the plastid", score=None, e_value=None, coverage=None, identity=0.981, aligned_bp=2004,
             own_best_score=None, other_organelle="pltd", other_contig="p1", other_start=4999, other_end=7000),
    ],
    columns=T.COLUMNS,
)


class Gff(unittest.TestCase):
    def test_both_kinds_are_biological_regions(self):
        lines = bg.transfer_lines(ROWS)
        f = [l.split("\t") for l in lines]
        self.assertEqual([x[2] for x in f], ["biological_region", "biological_region"])
        aln, hit = f  # sorted by start
        self.assertEqual(aln[3:7], ["51", "2050", ".", "."])
        self.assertIn("shared_with=pltd", aln[8])
        self.assertIn("other_region=p1:5000-7000", aln[8])
        self.assertIn("orientation=reverse", aln[8])
        self.assertNotIn("plastid_derived", aln[8])
        self.assertEqual(hit[3:7], ["100", "1160", "812.3", "+"])
        self.assertIn("Name=plastid psbA;plastid_derived=true", hit[8])
        self.assertIn("model_coverage=0.97", hit[8])
        self.assertEqual(len({x[8].split(";")[0] for x in f}), 2)  # IDs unique

    def test_a_mito_gene_in_a_plastome_is_mito_derived(self):
        row = ROWS.iloc[[0]].assign(organelle="pltd", other_organelle="mito", name="cox1")
        self.assertIn("Name=mito cox1;mito_derived=true", bg.transfer_lines(row)[0])

    def test_no_rows_no_lines(self):
        self.assertEqual(bg.transfer_lines(ROWS.iloc[0:0]), [])
        self.assertEqual(bg.transfer_lines(None), [])


class Unrotate(unittest.TestCase):
    def test_both_sides_go_back_onto_their_own_input(self):
        d = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, d)
        p = d / "organelle_transfers.tsv"
        t = pd.concat([ROWS, ROWS.iloc[[1]].assign(organelle="pltd", contig="p1", start=4999, end=7000, other_organelle="mito",
                                                   other_contig="m1", other_start=50, other_end=2050)])
        t.to_csv(p, sep="\t", index=False)
        origin.unrotate_table(p, "pltd", {"p1": 1000}, {"p1": 10000})
        u = pd.read_csv(p, sep="\t")
        self.assertEqual(u.iloc[1][["start", "end", "other_start", "other_end"]].tolist(), [50, 2050, 5999, 8000])  # mito row: other side moved
        self.assertEqual(u.iloc[2][["start", "end", "other_start", "other_end"]].tolist(), [5999, 8000, 50, 2050])  # pltd row: own side moved
        self.assertEqual(u.iloc[0][["start", "end"]].tolist(), [99, 1160])
        self.assertNotIn(".0", p.read_text().splitlines()[2])

    def test_other_region_attributes(self):
        d = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, d)
        g = d / "s.mito.gff"
        g.write_text("##gff-version 3\n" + "\n".join(bg.transfer_lines(ROWS)) + "\n")
        origin.unrotate_other_regions(g, {"p1": 9000}, {"p1": 10000})
        self.assertIn("other_region=p1:4000-6000", g.read_text())
        origin.unrotate_other_regions(g, {"x": 5}, {"x": 10})
        self.assertIn("other_region=p1:4000-6000", g.read_text())


if __name__ == "__main__":
    unittest.main()
