"""trans_splicing.reconstruct's manifest: which exon templates each species'
fragments are scanned against; and the templates' junctions.tsv, which says
how often the references join each pair of exons cis.

    python3 -m unittest tests.test_reconstruct
"""

from __future__ import annotations

import shutil
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from poa import PIPELINE_DIR  # noqa: E402
from poa.pipeline.common.errors import StageError  # noqa: E402
from poa.pipeline.trans_splicing import exon_profiles as ep  # noqa: E402
from poa.pipeline.trans_splicing import reconstruct as rc  # noqa: E402

REFERENCE = PIPELINE_DIR / "trans_splicing" / "reference"


class Manifest(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.root)
        for sp in ("Moss", "Rose"):
            d = self.root / "data" / "mito" / sp
            d.mkdir(parents=True)
            (d / f"{sp}.mito.ctg.fasta").write_text(">c1\nACGT\n")
        rows = [(sp, "mito", g, "c1", 1, 100, "+") for sp in ("Moss", "Rose") for g in ("nad1", "rps10", "cox1")]
        rows.append(("Ghost", "mito", "nad1", "c1", 1, 100, "+"))  # no assembly: skipped
        self.calls = pd.DataFrame(rows, columns=["species", "organelle", "gene", "contig_id", "start", "end", "strand"])

    def test_templates_follow_lineage(self):
        rows = rc.build_manifest(self.root, self.calls, None, {"Moss": "moss"}, "angiosperm")
        got = {(r["species"], r["gene"]): r["template_dir"] for r in rows}
        # cox1 isn't a reconstructed gene; moss has no rps10 template
        self.assertEqual(set(got), {("Moss", "nad1"), ("Rose", "nad1"), ("Rose", "rps10")})
        self.assertEqual(got["Moss", "nad1"], str(REFERENCE / "lineages" / "moss" / "profiles" / "nad1"))
        self.assertEqual(got["Rose", "nad1"], str(REFERENCE / "profiles" / "nad1"))

    def test_species_filter_and_forced_lineage(self):
        rows = rc.build_manifest(self.root, self.calls, {"Rose"}, {}, "moss")
        self.assertEqual({(r["species"], r["gene"]) for r in rows}, {("Rose", "nad1")})

    def test_missing_gene_calls(self):
        with self.assertRaises(StageError):
            rc.reconstruct(self.root, gene_calls=self.root / "nope.tsv")


class JunctionPriors(unittest.TestCase):
    def test_counted_over_references_with_the_templates_exon_count(self):
        tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, tmp)
        head = "accession\tgene\tlocus_tag\texon_num\tn_exons\tstart\tend\tstrand\tjunction_to_next\tgap_bp\tnotes\n"
        rows = [
            ("A", "clpP", "a1", 1, 3, "cis"), ("A", "clpP", "a1", 2, 3, "cis"), ("A", "clpP", "a1", 3, 3, "none"),
            ("B", "clpP", "b1", 1, 3, "trans"), ("B", "clpP", "b1", 2, 3, "cis"), ("B", "clpP", "b1", 3, 3, "none"),
            ("C", "clpP", "c1", 1, 1, "none"),  # no introns: numbered otherwise, left out
            ("D", "clpP", "d1", 1, 3, "trans"), ("D", "clpP", "d1", 2, 3, "trans"),  # held out
            ("A", "rps12", "a2", 1, 3, "trans"),  # another gene
        ]  # fmt: skip
        text = head + "".join(f"{a}\t{g}\t{t}\t{n}\t{k}\t1\t2\t+\t{j}\t0\tok\n" for a, g, t, n, k, j in rows)
        text += "E\tclpP\te1\t1\t3\t1\t2\t+\ttrans\t0\tUNRESOLVED\n"
        (tmp / "junctions.tsv").write_text(text)
        self.assertEqual(ep.junction_counts("clpP", 3, tmp / "junctions.tsv", ("D",)), {1: (1, 1), 2: (2, 0)})

    def test_bundled_templates_know_which_junctions_are_trans(self):
        def shares(d: Path) -> dict[int, float]:
            t = pd.read_csv(d / ep.JUNCTIONS_NAME, sep="\t")
            return {r.exon: r.n_cis / (r.n_cis + r.n_trans) for r in t.itertuples()}

        rps12 = shares(REFERENCE / "plastid" / "profiles" / "rps12")
        self.assertEqual((rps12[1], rps12[2]), (0.0, 1.0))
        self.assertEqual(set(shares(REFERENCE / "plastid" / "profiles" / "clpP").values()), {1.0})
        nad5 = shares(REFERENCE / "profiles" / "nad5")
        self.assertEqual((nad5[1], nad5[2], nad5[3], nad5[4]), (1.0, 0.0, 0.0, 1.0))
        # one row per junction, for every template that has the file
        for f in REFERENCE.rglob(ep.JUNCTIONS_NAME):
            if f.parent.parent.name != "profiles":
                continue  # raw/ tables
            n = max(int(p.stem[4:]) for p in f.parent.glob("exon*.pssm"))
            self.assertEqual(list(pd.read_csv(f, sep="\t").exon), list(range(1, n)), f)


if __name__ == "__main__":
    unittest.main()
