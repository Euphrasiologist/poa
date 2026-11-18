"""Genes across a circular contig's origin: find, rotate, map back."""

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from poa.pipeline.common import origin


def calls(rows):
    return pd.DataFrame(rows, columns=["contig_id", "start", "end", "gene", "strand"])


class Origin(unittest.TestCase):
    def test_split_gene_found_only_on_circular_contigs(self):
        c = calls([("c", 0, 800, "cox1", "-"), ("c", 9000, 10000, "cox1", "-"), ("c", 4000, 5000, "atp1", "+")])
        self.assertEqual(origin.split_genes(c, {"c": (10000, None)}), {"c": ["cox1"]})
        self.assertEqual(origin.split_genes(c, {"c": (10000, False)}), {})
        # different strands: two genes that happen to sit at the ends
        c.loc[0, "strand"] = "+"
        self.assertEqual(origin.split_genes(c, {"c": (10000, True)}), {})

    def test_new_origin_is_in_the_longest_gap_going_round(self):
        c = calls([("c", 0, 800, "a", "+"), ("c", 1000, 2000, "b", "+"), ("c", 9000, 10000, "a", "+")])
        self.assertEqual(origin.new_origin(c, 10000), 2000 + 7000 // 2)
        # longest gap wraps the origin: from 6000 round to 500
        c = calls([("c", 500, 1000, "a", "+"), ("c", 3000, 6000, "b", "+")])
        self.assertEqual(origin.new_origin(c, 10000), (6000 + 4500 // 2) % 10000)

    def test_round_trip(self):
        seq = "".join("ACGT"[(i * 7) % 4] for i in range(100))
        with tempfile.TemporaryDirectory() as tmp:
            fa = Path(tmp) / "x.fa"
            fa.write_text(f">c circular=true\n{seq}\n")
            origin.rotate_fasta(fa, {"c": 30})
            rotated = origin.read_fasta(fa)[0][1]
            self.assertEqual(rotated, seq[30:] + seq[:30])
            self.assertEqual(origin.contigs(fa), {"c": (100, True)})
            # a feature at rotated [60, 80) is input [90, 110): across the origin
            self.assertEqual(origin.interval(60, 80, 30, 100), (90, 110))
            self.assertEqual(rotated[60:80], (seq + seq)[90:110])
            self.assertEqual(origin.point(75, 30, 100), 5)

    def test_gff_children_stay_inside_their_gene(self):
        gff = "\n".join(
            [
                "##gff-version 3",
                "##sequence-region c 1 100",
                "c\tx\tregion\t1\t100\t.\t+\t.\tID=region-c",
                "c\tx\tgene\t61\t80\t1\t+\t.\tID=gene-a;Name=a",
                "c\tx\tmRNA\t61\t80\t1\t+\t.\tID=mRNA-a;Parent=gene-a",
                "c\tx\tCDS\t61\t80\t1\t+\t0\tID=cds-a;Parent=mRNA-a",
                "c\tx\tsequence_alteration\t76\t76\t.\t+\t.\tID=sa-a-1;Parent=cds-a",
                "c\tx\tgene\t11\t20\t1\t+\t.\tID=gene-b;Name=b",
            ]
        )
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "x.gff"
            p.write_text(gff + "\n")
            origin.unrotate_gff(p, {"c": 30}, {"c": 100})
            rows = [l.split("\t") for l in p.read_text().splitlines() if not l.startswith("#")]
            got = {r[8].split(";")[0]: (int(r[3]), int(r[4])) for r in rows}
            self.assertEqual(got["ID=gene-a"], (91, 110))
            self.assertEqual(got["ID=cds-a"], (91, 110))
            self.assertEqual(got["ID=sa-a-1"], (106, 106))  # 6 on the input, kept inside gene-a
            self.assertEqual(got["ID=gene-b"], (41, 50))
            self.assertIn("Is_circular=true", rows[0][8])

    def test_table_rows_of_the_other_organelle_untouched(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "gene_calls.tsv"
            pd.DataFrame(
                {"species": ["s", "s"], "organelle": ["mito", "pltd"], "contig_id": ["c", "c"],
                 "start": [60, 60], "end": [80, 80], "gene": ["a", "b"]}  # fmt: skip
            ).to_csv(p, sep="\t", index=False)
            origin.unrotate_table(p, "mito", {"c": 30}, {"c": 100})
            df = pd.read_csv(p, sep="\t")
            self.assertEqual(list(zip(df.start, df.end)), [(90, 110), (60, 80)])


if __name__ == "__main__":
    unittest.main()
