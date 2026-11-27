"""Spliced genes in the GFF: junction labels from transsplice, cis genes
across a circular contig's origin, and their unitig coordinates."""

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from poa.pipeline.gff_export import build_gff as bg


def part(n, start, end, strand="+", contig="c"):
    return {"contig": contig, "start": start, "end": end, "strand": strand, "score": 1.0, "exon_number": n}


class Junctions(unittest.TestCase):
    def test_labels_in_exon_order_and_another_contig_said_so(self):
        j = pd.DataFrame(
            {"from_slot": [2, 1, 3], "to_slot": [3, 2, 4], "kind": ["trans", "cis", "trans"], "gap_bp": [80000, 700, -1]}
        )
        self.assertEqual(bg.junction_labels(j), ["cis", "trans", "other_contig"])
        self.assertEqual(bg.junction_labels(j.iloc[0:0]), [])

    def test_a_gene_not_all_cis_is_marked_trans_splicing(self):
        for junctions, exception in ((["cis"], None), (["cis", "trans"], "trans-splicing"), (["other_contig"], "trans-splicing")):
            lines: list[str] = []
            parts = [part(1, 100, 200), part(2, 900, 1000)]
            bg.emit_protein_coding_gene(
                lines, "ndhB", "", "reconstructed", parts, [], 1.0, 0, {"junctions": junctions}, phase_from_lengths=False
            )
            rows = {l.split("\t")[2]: l.split("\t")[8] for l in lines}
            self.assertIn(f"junctions={','.join(junctions)}", rows["mRNA"])
            for kind in ("gene", "mRNA"):
                self.assertEqual("exception=trans-splicing" in rows[kind], exception is not None, (junctions, kind))


class AcrossOrigin(unittest.TestCase):
    def test_exons_either_side_of_the_origin_are_written_past_the_length(self):
        # ndhB of a fern plastome: exon 1 just after the origin, exon 2 just before it
        parts = [part(1, 382, 901, "-"), part(2, 148421, 149192, "-")]
        moved, which = bg.across_origin(parts, 149609)
        self.assertEqual(which, {1})
        self.assertEqual([(p["start"], p["end"]) for p in moved], [(382 + 149609, 901 + 149609), (148421, 149192)])

    def test_exons_on_one_side_stay(self):
        parts = [part(1, 1000, 1500), part(2, 2200, 2600)]
        self.assertEqual(bg.across_origin(parts, 150000), (parts, set()))

    def test_placements_take_a_second_lap_of_a_circular_contig(self):
        with tempfile.TemporaryDirectory() as tmp:
            fa, tsv = Path(tmp) / "c.fa", Path(tmp) / "map.tsv"
            fa.write_text(">c circular=true\n" + "A" * 1000 + "\n")
            pd.DataFrame(
                {"ctg_seqid": ["c", "c"], "ctg_start": [0, 600], "ctg_end": [600, 1000], "unitig": ["u1", "u2"],
                 "unitig_length": [600, 400], "strand": ["+", "+"]}  # fmt: skip
            ).to_csv(tsv, sep="\t", index=False)
            ctx = bg.load_seq_context(fa, tsv)
            # a CDS at 1101-1200 is 101-200, on u1
            line = "c\tx\tCDS\t1101\t1200\t.\t+\t0\tID=cds-a"
            out, unitig = bg.add_unitig_layer([line], ctx)
            self.assertIn("unitig_loc=u1:101-200:+", out[0])
            # a gene from 901 to 1100 crosses from u2 into u1
            out, unitig = bg.add_unitig_layer(["c\tx\tgene\t901\t1100\t.\t+\t.\tID=gene-a"], ctx)
            self.assertIn("unitig_span=u2,u1", out[0])
            self.assertEqual(sorted(l.split("\t")[0] for l in unitig), ["u1", "u2"])


if __name__ == "__main__":
    unittest.main()



class OrganelleOrigin(unittest.TestCase):
    """common.organelle_origin on consensus proteins of the bundled profiles."""

    @staticmethod
    def consensus(organelle, gene):
        from poa.pipeline.common import organelle_origin as oo

        index, scores = oo.load_profile(organelle, gene)
        residues = {i: a for a, i in index.items()}
        return "".join(residues[int(k)] for k in scores.argmax(axis=1))

    def test_each_organelles_own_rpl2_stays_and_the_others_is_flagged(self):
        from poa.pipeline.common import organelle_origin as oo

        mito, pltd = self.consensus("mito", "rpl2"), self.consensus("pltd", "rpl2")
        calls = pd.DataFrame([
            dict(species="a", organelle="mito", gene="rpl2", protein=mito),
            dict(species="a", organelle="pltd", gene="rpl2", protein=pltd),
            dict(species="b", organelle="mito", gene="rpl2", protein=pltd),  # plastid DNA in the mitogenome
            dict(species="c", organelle="pltd", gene="rps14", protein=self.consensus("mito", "rps14")),
        ])
        self.assertEqual(oo.foreign_calls(calls), {("b", "mito", "rpl2"): "pltd", ("c", "pltd", "rps14"): "mito"})

    def test_edited_copies_sharing_a_name_are_told_apart_by_location(self):
        from poa.pipeline.common import organelle_origin as oo

        calls = pd.DataFrame([
            dict(species="a", organelle="mito", gene="rps14", contig="c1", corrected_start=100, corrected_end=400, protein=self.consensus("mito", "rps14")),
            dict(species="a", organelle="mito", gene="rps14", contig="c1", corrected_start=900, corrected_end=1200, protein=self.consensus("pltd", "rps14")),
        ])
        self.assertEqual(oo.foreign_calls(calls), {("a", "mito", "rps14@c1:900-1200"): "pltd"})

    def test_genes_one_organelle_has_and_short_or_failed_calls_are_left_alone(self):
        from poa.pipeline.common import organelle_origin as oo

        calls = pd.DataFrame([
            dict(species="a", organelle="mito", gene="cox1", protein=self.consensus("pltd", "rpl2")),
            dict(species="a", organelle="mito", gene="rps14", protein="MKRLAQ"),
        ])
        self.assertEqual(oo.foreign_calls(calls), {})
        self.assertEqual(oo.foreign_calls(pd.DataFrame()), {})

    def test_the_flag_is_a_gene_attribute_named_for_the_other_organelle(self):
        self.assertEqual(bg.copy_attrs({"derived": "plastid_derived"})["plastid_derived"], "true")
        self.assertEqual(bg.copy_attrs({"derived": "mito_derived"})["mito_derived"], "true")
        self.assertNotIn("plastid_derived", bg.copy_attrs({"derived": None}))
