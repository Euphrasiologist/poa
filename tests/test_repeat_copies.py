"""A reconstructed gene written again on an exact repeat copy
(gff_export.build_gff.repeat_copies): the second inverted-repeat copy of a
plastid gene, with its exons, edits and frameshifts mapped across.

    python3 -m unittest tests.test_repeat_copies
"""

from __future__ import annotations

import random
import sys
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from poa.pipeline.gff_export import build_gff as bg  # noqa: E402


def revcomp(s: str) -> str:
    return s.translate(bg.COMPLEMENT)[::-1]


def locus(contig, start, end, strand="+"):
    return bg.Locus(contig, start, end, strand, 100.0, *([None] * (len(bg.Locus.__dataclass_fields__) - 5)))


class RepeatCopies(unittest.TestCase):
    def setUp(self):
        rng = random.Random(3)
        r = lambda n: "".join(rng.choice("ACGT") for _ in range(n))  # noqa: E731
        self.gene = r(300)  # exons at 0-100 and 200-300 of it
        # the gene at 100-400, its reverse complement (the other IR copy) at 1000-1300
        self.seq = r(100) + self.gene + r(600) + revcomp(self.gene) + r(700)
        self.rec = {
            "parts": [
                {"contig": "c", "start": 100, "end": 200, "strand": "+", "score": 1.0, "exon_number": 1, "phase": 0},
                {"contig": "c", "start": 300, "end": 400, "strand": "+", "score": 1.0, "exon_number": 2, "phase": 2},
            ],
            "junctions": ["cis"],
            "edits": [{"slot": 1, "genomic_pos": 150, "resulting_aa": "L"}, {"slot": 2, "genomic_pos": 310, "resulting_aa": "W"}],
            "frameshifts": [{"slot": 2, "contig": "c", "start": 320, "end": 321, "n_bases": 1, "window_start": 318, "window_end": 324}],
            "score": 900.0,
        }

    def test_the_other_inverted_repeat_copy_gets_the_reconstruction(self):
        raw = {"loci": [locus("c", 1010, 1290, "-")]}
        [(claimed, again)] = bg.repeat_copies(self.rec, [raw], {"c": self.seq})
        self.assertIs(claimed, raw)
        self.assertEqual(
            [(p["start"], p["end"], p["strand"], p["exon_number"], p["phase"]) for p in again["parts"]],
            [(1200, 1300, "-", 1, 0), (1000, 1100, "-", 2, 2)],
        )
        # each exon's sequence is the original's reverse complement
        for p, q in zip(self.rec["parts"], again["parts"]):
            self.assertEqual(self.seq[q["start"]:q["end"]], revcomp(self.seq[p["start"]:p["end"]]))
        self.assertEqual([e["genomic_pos"] for e in again["edits"]], [1249, 1089])
        self.assertEqual(self.seq[1249], revcomp(self.seq[150]))
        f = again["frameshifts"][0]
        self.assertEqual((f["start"], f["end"], f["window_start"], f["window_end"]), (1079, 1080, 1076, 1082))
        self.assertIs(again["repeat_copy_of"], self.rec)
        self.assertEqual(self.rec["parts"][0]["start"], 100)  # the original is untouched

    def test_a_trans_spliced_gene_shares_the_block_with_no_duplicate(self):
        # rps12: exon 1 single-copy elsewhere, exons 2-3 (cis) in each IR copy
        rec = dict(self.rec, junctions=["trans", "cis"])
        rec["parts"] = [
            {"contig": "c", "start": 1500, "end": 1560, "strand": "-", "score": 1.0, "exon_number": 1},
            {"contig": "c", "start": 100, "end": 200, "strand": "+", "score": 1.0, "exon_number": 2},
            {"contig": "c", "start": 300, "end": 400, "strand": "+", "score": 1.0, "exon_number": 3},
        ]
        rec["edits"], rec["frameshifts"] = [{"slot": 1, "genomic_pos": 1520}], []
        [(_, again)] = bg.repeat_copies(rec, [{"loci": [locus("c", 1000, 1300, "-")]}], {"c": self.seq})
        self.assertEqual([(p["start"], p["end"], p["strand"]) for p in again["parts"]], [(1500, 1560, "-"), (1200, 1300, "-"), (1000, 1100, "-")])
        self.assertEqual(again["edits"], rec["edits"])

    def test_in_a_mitogenome_the_whole_gene_must_be_duplicated(self):
        rec = dict(self.rec, junctions=["trans", "cis"])
        rec["parts"] = [
            {"contig": "c", "start": 1500, "end": 1560, "strand": "-", "score": 1.0, "exon_number": 1},
            {"contig": "c", "start": 100, "end": 200, "strand": "+", "score": 1.0, "exon_number": 2},
            {"contig": "c", "start": 300, "end": 400, "strand": "+", "score": 1.0, "exon_number": 3},
        ]
        rec["edits"], rec["frameshifts"] = [], []
        raw = [{"loci": [locus("c", 1000, 1300, "-")]}]
        self.assertEqual(bg.repeat_copies(rec, raw, {"c": self.seq}, share_blocks=False), [])
        self.assertEqual(len(bg.repeat_copies(self.rec, raw, {"c": self.seq}, share_blocks=False)), 1)

    def test_a_gene_across_the_origin_is_read_round_it(self):
        # the gene's second exon past the end of a circular contig: rotate so it wraps
        k = 350  # the contig now starts inside exon 2
        rot = self.seq[k:] + self.seq[:k]
        n = len(rot)
        rec = dict(self.rec, edits=[], frameshifts=[])
        rec["parts"] = [
            {"contig": "c", "start": 100 - k + n, "end": 200 - k + n, "strand": "+", "score": 1.0, "exon_number": 1},
            {"contig": "c", "start": 300 - k + n, "end": 400 - k + n, "strand": "+", "score": 1.0, "exon_number": 2},
        ]
        [(_, again)] = bg.repeat_copies(rec, [{"loci": [locus("c", 1010 - k, 1290 - k, "-")]}], {"c": rot})
        self.assertEqual([(p["start"], p["end"], p["strand"]) for p in again["parts"]], [(1200 - k, 1300 - k, "-"), (1000 - k, 1100 - k, "-")])

    def test_no_copy_without_an_exact_duplicate_at_the_raw_copy(self):
        elsewhere = {"loci": [locus("c", 1500, 1700)]}
        self.assertEqual(bg.repeat_copies(self.rec, [elsewhere], {"c": self.seq}), [])
        diverged = self.seq[:1150] + ("A" if self.seq[1150] != "A" else "C") + self.seq[1151:]
        self.assertEqual(bg.repeat_copies(self.rec, [{"loci": [locus("c", 1010, 1290)]}], {"c": diverged}), [])
        self.assertEqual(bg.repeat_copies(self.rec, [{"loci": [locus("c", 1010, 1290)]}], {}), [])


if __name__ == "__main__":
    unittest.main()
