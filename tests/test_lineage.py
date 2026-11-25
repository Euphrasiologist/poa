"""The per-species reference-lineage call (common/lineage.py) and how
editing/01_scan_editing.py splits its orfedit runs by it.

    python3 -m unittest tests.test_lineage
"""

from __future__ import annotations

import shutil
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
PIPELINE = REPO / "src" / "poa" / "pipeline"
sys.path.insert(0, str(REPO / "src"))
from poa.pipeline.common import lineage as lin  # noqa: E402
from poa.pipeline.editing import scan as scan_editing  # noqa: E402

GENES = ["atp4", "ccmB", "ccmC", "mttB", "rps14"]


def scores(angio: float, moss: float, genes=GENES) -> dict:
    return {"angiosperm": {g: angio for g in genes}, "moss": {g: moss for g in genes}}


class Choose(unittest.TestCase):
    def pick(self, s):
        return lin.choose(lin.margins(s), lin.LINEAGES)

    def test_clear_moss(self):
        m = lin.margins(scores(300, 600))
        self.assertEqual(m, {"moss": (300.0, 5)})
        self.assertEqual(lin.choose(m), "moss")

    def test_clear_angiosperm(self):
        self.assertEqual(self.pick(scores(600, 300)), "angiosperm")

    def test_narrow_win_stays_angiosperm(self):
        self.assertEqual(self.pick(scores(500, 500 + lin.MIN_MARGIN["moss"] - 1)), "angiosperm")
        self.assertEqual(self.pick(scores(500, 500 + lin.MIN_MARGIN["moss"])), "moss")

    def test_too_few_shared_genes(self):
        few = GENES[: lin.MIN_SHARED_GENES - 1]
        s = scores(100, 900, few)
        self.assertEqual(lin.margins(s), {"moss": (800.0, len(few))})
        self.assertEqual(self.pick(s), "angiosperm")

    def test_only_shared_genes_count(self):
        s = scores(300, 600)
        s["angiosperm"]["nad5"] = 10_000  # no moss profile for nad5: ignored
        self.assertEqual(lin.margins(s), {"moss": (300.0, 5)})
        self.assertEqual(self.pick(s), "moss")

    def test_no_moss_run(self):
        s = {"angiosperm": {g: 1 for g in GENES}}
        self.assertEqual(lin.margins(s), {})
        self.assertEqual(self.pick(s), "angiosperm")

    def test_no_default_run(self):
        self.assertEqual(lin.margins({"moss": {g: 1 for g in GENES}}), {})

    def test_each_lineage_over_its_own_genes(self):
        # liverwort profiles cover more genes than the moss ones; each is compared with the default over its own
        s = scores(300, 600)
        s["liverwort"] = {"atp4": 700, "atp6": 700, "atp8": 700, "nad6": 700}
        s["angiosperm"].update(atp6=300, atp8=300, nad6=300)
        self.assertEqual(lin.margins(s), {"moss": (300.0, 5), "liverwort": (400.0, 4)})


class ChooseAmongLineages(unittest.TestCase):
    """A margin per lineage: profiles from a few edited proteins are sharper, and bryophytes fit each other's."""

    def test_takes_the_best_lineage_over_its_margin(self):
        # a horsetail also fits the fern profiles, less well
        self.assertEqual(lin.choose({"fern": (85.0, 17), "horsetail": (499.0, 23), "lycophyte": (52.0, 14)}), "horsetail")
        self.assertEqual(lin.choose({"fern": (271.0, 16), "horsetail": (11.0, 18), "lycophyte": (25.0, 12)}), "fern")
        self.assertEqual(lin.choose({"fern": (8.0, 17), "horsetail": (-47.0, 19), "lycophyte": (429.0, 14)}), "lycophyte")

    def test_moss_and_liverwort(self):
        # a moss fits the liverwort profiles too, and a liverwort the moss ones: each short of the other's margin
        self.assertEqual(lin.choose({"moss": (310.0, 5), "liverwort": (210.0, 20)}), "moss")
        self.assertEqual(lin.choose({"moss": (169.0, 5), "liverwort": (361.0, 20)}), "liverwort")
        # a clubmoss on the liverwort profiles, a hornwort
        self.assertEqual(lin.choose({"liverwort": (151.0, 18)}), "angiosperm")
        self.assertEqual(lin.choose({"liverwort": (157.0, 18)}), "angiosperm")

    def test_leaves_the_rest_on_the_default(self):
        # a liverwort on the lycophyte profiles, an angiosperm, a hornwort
        self.assertEqual(lin.choose({"fern": (18.0, 18), "horsetail": (-43.0, 23), "lycophyte": (144.0, 14)}), "angiosperm")
        self.assertEqual(lin.choose({"fern": (-275.0, 18), "horsetail": (-417.0, 23), "lycophyte": (-380.0, 14)}), "angiosperm")
        self.assertEqual(lin.choose({"fern": (13.0, 15), "horsetail": (-24.0, 15), "lycophyte": (100.0, 15)}), "angiosperm")
        self.assertEqual(lin.choose({}), "angiosperm")
        # the ferns furthest from the fern references (Ophioglossum, Cryptogramma) are still ferns
        self.assertEqual(lin.choose({"fern": (85.0, 20), "horsetail": (41.0, 20)}), "fern")
        self.assertEqual(lin.choose({"fern": (139.0, 11)}), "fern")

    def test_needs_enough_genes_and_the_margin(self):
        few = lin.MIN_SHARED_GENES - 1
        self.assertEqual(lin.choose({"horsetail": (600.0, few)}), "angiosperm")
        for name, margin in lin.MIN_MARGIN.items():
            self.assertEqual(lin.choose({name: (margin - 1, 10)}), "angiosperm")
            self.assertEqual(lin.choose({name: (margin, 10)}), name)

    def test_among(self):
        m = {"moss": (900.0, 10), "fern": (300.0, 10)}
        self.assertEqual(lin.choose(m, lin.DETECTED), "fern")
        self.assertEqual(lin.choose(m, lin.LINEAGES), "moss")
        self.assertEqual(set(lin.DETECTED), set(lin.BORROWS))
        self.assertEqual(lin.MARGIN_COLUMNS, [f"{x}_margin" for x in ("moss", "liverwort", *lin.DETECTED)])
        self.assertEqual(lin.COLUMNS[4:], [*lin.MARGIN_COLUMNS, *lin.SITE_COLUMNS])
        self.assertEqual(set(lin.MIN_MARGIN), set(lin.LINEAGES[1:]) | set(lin.DETECTED))


class EditingSites(unittest.TestCase):
    def test_follows_sites(self):
        # judged on the mitochondrial calls when there are enough of them
        self.assertTrue(lin.follows_sites({"mito": (190, 200), "pltd": (5, 60)}))
        self.assertFalse(lin.follows_sites({"mito": (300, 700), "pltd": (25, 300)}))  # Ginkgo
        self.assertEqual(lin.follows_sites({"mito": (13, 20)}), 13 / 20 >= lin.SITE_MIN_SHARE["mito"])
        # else on the plastid ones, which have many more calls away from the known columns
        self.assertTrue(lin.follows_sites({"mito": (2, 5), "pltd": (15, 40)}))
        self.assertFalse(lin.follows_sites({"pltd": (22, 600)}))  # a hornwort's plastid
        # and with too few of either the species keeps the ordinary cost
        self.assertFalse(lin.follows_sites({"mito": (5, 5), "pltd": (3, 3)}))
        self.assertFalse(lin.follows_sites({}))

    def test_known_site_calls(self):
        with tempfile.TemporaryDirectory() as tmp:
            header = "\t".join("ARNDCQEGHILKMFPSTWYV")
            row = "\t".join(["0.0"] * 20)
            sites, plain = Path(tmp) / "mito" / "nad3.pssm", Path(tmp) / "mito" / "rps2.pssm"
            sites.parent.mkdir()
            sites.write_text(f"{header}\tedit\n" + "".join(f"{row}\t{e}\n" for e in ("0", "0.5", "0", "1")))
            plain.write_text(f"{header}\n{row}\n{row}\n")
            self.assertEqual(scan_editing.site_columns(sites), [False, True, False, True])
            self.assertIsNone(scan_editing.site_columns(plain))
            rows = [dict(organelle="mito", profile=str(p)) for p in (sites, plain)]
            edits = pd.DataFrame(
                [
                    ("Rose", "mito", "nad3", "1", "C>U"),
                    ("Rose", "mito", "nad3", "2", "C>U"),
                    ("Rose", "mito", "nad3", "3", "U>C"),  # not a C-to-U edit
                    ("Rose", "mito", "nad3", "4", "C>U"),  # an edit-made stop, past the last column
                    ("Rose", "mito", "rps2", "1", "C>U"),  # a profile without sites
                ],
                columns=["species", "organelle", "gene", "profile_column", "edit"],
            )
            self.assertEqual(scan_editing.known_site_calls(edits, rows), {"Rose": {"mito": (1, 2)}})

    def test_load_site_species(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(lin.load_site_species(Path(tmp)), set())
            path = lin.table_path(Path(tmp))
            path.parent.mkdir(parents=True)
            path.write_text("species\tlineage\tedit_sites\nRose\tangiosperm\tyes\nGinkgo\tangiosperm\tno\nMoss\tmoss\t\n")
            self.assertEqual(lin.load_site_species(Path(tmp)), {"Rose"})


class SplitRuns(unittest.TestCase):
    def setUp(self):
        def calls(rows):
            return pd.DataFrame(rows, columns=["species", "organelle", "gene", "score"]).astype(str)

        self.angio = calls(
            [("Moss", "mito", g, 300) for g in GENES + ["nad5"]]
            + [("Moss", "pltd", "psbA", 900)]
            + [("Rose", "mito", g, 600) for g in GENES + ["nad5"]]
        )
        self.moss = calls([("Moss", "mito", g, 600) for g in GENES] + [("Rose", "mito", g, 200) for g in GENES])

    def test_assign_and_keep(self):
        table = scan_editing.assign_lineages({"angiosperm": self.angio, "moss": self.moss}, None)
        lineages = dict(zip(table["species"], table["lineage"]))
        self.assertEqual(lineages, {"Moss": "moss", "Rose": "angiosperm"})
        self.assertEqual(list(table.columns), lin.COLUMNS)
        kept_a = scan_editing.keep_lineage(self.angio, lineages, "angiosperm")
        kept_m = scan_editing.keep_lineage(self.moss, lineages, "moss")
        # the moss keeps its plastid call from the default run, and loses nad5
        # (no moss profile) rather than keeping the angiosperm correction
        self.assertEqual(sorted(kept_a[kept_a.species == "Moss"].gene), ["psbA"])
        self.assertEqual(sorted(kept_m.gene), sorted(GENES))
        self.assertEqual(set(kept_a[kept_a.species == "Rose"].gene), set(GENES + ["nad5"]))

    def test_fern_borrows_the_default_references(self):
        # fern: no references of its own, never chosen automatically, and a
        # stop read-through cost
        self.assertEqual(lin.reference("fern"), lin.DEFAULT_LINEAGE)
        self.assertEqual(lin.reference("moss"), "moss")
        self.assertIn("fern", lin.CHOICES)
        self.assertNotIn("fern", lin.LINEAGES)
        self.assertEqual(lin.stop_edit_penalty("fern"), 2.5)
        # other U-to-C edits: fern only, and only on its own references
        self.assertEqual(lin.reverse_edit_penalty("fern"), 5.0)
        self.assertIsNone(lin.reverse_edit_penalty("angiosperm"))
        self.assertIsNone(lin.reverse_edit_penalty("moss"))
        self.assertEqual(lin.junction_slide("fern"), 6)
        self.assertIsNone(lin.junction_slide("angiosperm"))
        self.assertIsNone(lin.stop_edit_penalty("angiosperm"))
        self.assertIsNone(lin.stop_edit_penalty("moss"))
        self.assertEqual(lin.choose(lin.margins(scores(100, 900)), lin.LINEAGES), "moss")  # fern is not looked for here

    def test_forced_fern_keeps_the_default_run(self):
        # a species forced to fern takes its mito rows from the default
        # (angiosperm-profile) run, not from another lineage's
        df = pd.DataFrame({"species": ["F", "F", "A"], "organelle": ["mito", "pltd", "mito"], "gene": ["atp4", "rbcL", "atp4"]})
        lineages = {"F": "fern", "A": "angiosperm"}
        kept = scan_editing.keep_lineage(df, lineages, lin.DEFAULT_LINEAGE)
        self.assertEqual(len(kept), 3)
        self.assertTrue(scan_editing.keep_lineage(df, lineages, "moss").empty)

    def test_forced(self):
        table = scan_editing.assign_lineages({"angiosperm": self.angio, "moss": self.moss}, "angiosperm")
        self.assertEqual(set(table["lineage"]), {"angiosperm"})
        self.assertEqual(set(table["decided_by"]), {"forced"})


class FernReferences(unittest.TestCase):
    def test_verify_allows_reverse_edits_only_when_asked(self):
        from poa.pipeline.common.refseq_mito import editing_explains, verify

        # ATG AAA TAA GTT: the annotated protein reads the TAA as Gln
        dna, protein = "ATGAAATAAGTT", "MKQV"
        self.assertFalse(verify(dna, protein)[0])
        ok, note = verify(dna, protein, reverse=True)
        self.assertTrue(ok, note)
        self.assertTrue(editing_explains("TGA", "R", reverse=True))
        self.assertFalse(editing_explains("TGA", "R"))
        # C-to-U still works either way, and nothing else is explained
        self.assertTrue(editing_explains("CCA", "L"))
        self.assertFalse(editing_explains("AAA", "E", reverse=True))

    def test_column_mask(self):
        from poa.pipeline.trans_splicing import column_mask as cm

        # column 0 unaligned in all, column 3 in 4 of 5 (80%), column 5 in 3 of 5
        proteins = ["-KV-F-", "-KV-F-", "-KV-F-", "-KV-FA", "-KVLFA"]
        self.assertEqual(cm.unsupported_columns(proteins), [0, 3])
        self.assertEqual(cm.unsupported_columns([]), [])
        tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, tmp)
        pssm = tmp / "whole_gene.pssm"
        pssm.write_text("A\tR\n" + "".join(f"{i}.0\t0.0\n" for i in range(6)))
        self.assertEqual(cm.apply_mask(pssm, {0, 3}), 2)
        self.assertEqual([line.split("\t")[0] for line in pssm.read_text().splitlines()[1:]], ["1.0", "2.0", "4.0", "5.0"])

    def test_fern_templates_are_bundled_and_used(self):
        from poa.pipeline.trans_splicing import reconstruct
        from poa.pipeline.trans_splicing.reference_genes import LINEAGE_MITOGENOMES

        self.assertIn("fern", LINEAGE_MITOGENOMES)
        profiles = lin.trans_splicing_reference(PIPELINE, "fern") / "profiles"
        self.assertTrue((profiles / "rps10" / "exon1.pssm").exists())
        self.assertFalse((profiles / "rps10" / "exon2.pssm").exists())  # one exon in ferns
        self.assertFalse((profiles / "ccmFc").exists())  # lost from fern mitogenomes
        # fern: its own exon templates, the angiosperm editing profiles
        self.assertEqual(reconstruct.template_lineage("fern"), "fern")
        self.assertEqual(lin.reference("fern"), lin.DEFAULT_LINEAGE)


class FernEditingProfiles(unittest.TestCase):
    def test_bundled_and_gated(self):
        from poa.pipeline.editing import lineage_profiles as lp

        mito = lin.editing_profiles(PIPELINE, "fern") / "mito"
        have = {p.stem for p in mito.glob("*.pssm")}
        self.assertIn("rps4", have)
        self.assertIn("atp4", have)
        # left out: their fern profiles did worse than the angiosperm ones
        self.assertFalse({"matR", "nad9"} & have)
        self.assertFalse({"matR", "nad9"} & set(lp.LINEAGE_GENES["fern"]))
        # multi-exon in ferns, so no single-locus profile
        self.assertFalse({"atp1", "cox1", "rps14"} & have)
        self.assertEqual(lp.LINEAGE_MIN_SEQUENCES["fern"], 4)

    def test_rna_references(self):
        """the RNA-edited DToL loci and proteins are bundled and read; the
        GenBank ferns with unedited proteins are left out of the templates"""
        from poa.pipeline.editing import lineage_profiles as lp
        from poa.pipeline.trans_splicing import exon_profiles as ep
        from poa.pipeline.trans_splicing import reference_genes as rg

        proteins = lp.rna_proteins("fern")
        self.assertIn("RNA_Onoclea_sensibilis", proteins["nad3"])
        self.assertNotIn("RNA_Onoclea_sensibilis", lp.rna_proteins("fern", ("RNA_Onoclea_sensibilis",))["nad3"])
        self.assertEqual(lp.rna_proteins("moss"), {})
        self.assertEqual(lp.LINEAGE_PSEUDOCOUNT["fern"], 2.0)

        ref = lin.trans_splicing_reference(PIPELINE, "fern")
        dirs = ep.raw_dirs(ref)
        self.assertEqual([d.name for d in dirs], ["raw", ep.RNA_RAW_NAME])
        raw, junctions = ep.load_all_raw_data(*dirs), ep.load_junctions(*dirs)
        key = ("RNA_Onoclea_sensibilis", "rps3", "rps3_1")
        self.assertEqual(junctions[key]["n_exons"], 3)
        self.assertEqual(sorted(k for k in raw[key] if k != "refseq"), [1, 2, 3])
        self.assertNotIn("*", raw[key]["refseq"])
        self.assertIn("MW732172.1", rg.LINEAGE_UNEDITED_MITOGENOMES["fern"])
        self.assertFalse(any(a.startswith("OM8675") for a in rg.LINEAGE_UNEDITED_MITOGENOMES["fern"]))
        self.assertEqual(ep.LINEAGE_GAP_THRESHOLD["fern"], 0.5)

    def test_one_per_species(self):
        from poa.pipeline.editing import lineage_profiles as lp

        proteins = {"nad3": {"A1": "MK", "A2": "MK", "A3": "ML", "B1": "MK"}}
        species = {"A1": "a", "A2": "a", "A3": "a", "B1": "b"}
        self.assertEqual(list(lp.one_per_species(proteins, species)["nad3"]), ["A1", "A3", "B1"])


class Load(unittest.TestCase):
    def test_missing_and_present(self):
        root = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, root)
        self.assertIsNone(lin.load(root))
        lin.table_path(root).parent.mkdir(parents=True)
        lin.table_path(root).write_text("species\tlineage\nA\tmoss\n")
        self.assertEqual(lin.load(root), {"A": "moss"})


class LineageSplicedGenes(unittest.TestCase):
    """Genes with introns in one lineage only: reconstructed there, and
    left out of that lineage's single-exon editing."""

    def test_gene_lists(self):
        from poa.pipeline.trans_splicing import reference_genes as rg

        self.assertEqual(rg.mito_spliced("angiosperm"), rg.GENES)
        self.assertEqual(rg.mito_spliced("moss"), rg.GENES)
        self.assertEqual(rg.mito_spliced("fern")[len(rg.GENES) :], ["atp1", "cox1", "rps14"])
        self.assertIn("atp9", rg.mito_spliced("lycophyte"))
        every = rg.mito_spliced()
        self.assertEqual(len(every), len(set(every)))
        self.assertTrue({"atp1", "cox1", "atp6", "nad3"} <= set(every))

    def test_editing_leaves_them_out(self):
        self.assertEqual(scan_editing.lineage_spliced("angiosperm"), set())
        self.assertEqual(scan_editing.lineage_spliced("auto"), set())
        self.assertEqual(scan_editing.lineage_spliced("fern"), {"atp1", "cox1", "rps14"})
        self.assertEqual(scan_editing.lineage_spliced("horsetail"), {"cox1"})
        rows = [
            {"organelle": "mito", "gene": "cox1|ctg1:5-900"},
            {"organelle": "mito", "gene": "cox3|ctg1:1000-1800"},
            {"organelle": "pltd", "gene": "cox1|ctg1:5-900"},
        ]
        self.assertEqual(scan_editing.without_genes(rows, set()), rows)
        self.assertEqual(scan_editing.without_genes(rows, {"cox1"}), rows[1:])

    def test_templates(self):
        from poa.pipeline.editing import lineage_profiles as lp
        from poa.pipeline.trans_splicing import reconstruct

        fern = lin.trans_splicing_reference(PIPELINE, "fern") / "profiles"
        for gene in ("atp1", "cox1", "rps14"):
            self.assertTrue((fern / gene / "exon2.pssm").exists(), gene)
            self.assertFalse((fern / gene / "exon3.pssm").exists(), gene)
            self.assertNotIn(gene, lp.LINEAGE_GENES["fern"])
        # no template for them in a lineage that does not splice them
        self.assertFalse((lin.trans_splicing_reference(PIPELINE, "moss") / "profiles" / "cox1").exists())
        self.assertFalse((lin.trans_splicing_reference(PIPELINE, "angiosperm") / "profiles" / "cox1").exists())
        self.assertEqual(reconstruct.template_lineage("fern"), "fern")


class HorsetailAndLycophyte(unittest.TestCase):
    """Two more lineages chosen only with --lineage, each with its own
    exon templates and editing profiles on top of the angiosperm ones."""

    def test_registered(self):
        from poa.pipeline.trans_splicing import reconstruct

        for name in ("horsetail", "lycophyte"):
            self.assertIn(name, lin.CHOICES)
            self.assertNotIn(name, lin.LINEAGES)  # never picked automatically
            self.assertEqual(lin.reference(name), lin.DEFAULT_LINEAGE)
            self.assertEqual(reconstruct.template_lineage(name), name)
            self.assertIsNone(lin.reverse_edit_penalty(name))
            self.assertIsNone(lin.junction_slide(name))
        # horsetails edit C to U only; clubmosses have a few stop codons edited away
        self.assertIsNone(lin.stop_edit_penalty("horsetail"))
        self.assertEqual(lin.stop_edit_penalty("lycophyte"), 2.5)

    def test_horsetail_references(self):
        from poa.pipeline.editing import lineage_profiles as lp
        from poa.pipeline.trans_splicing import exon_profiles as ep
        from poa.pipeline.trans_splicing import reference_genes as rg

        self.assertEqual(rg.LINEAGE_MITOGENOMES["horsetail"], {})  # RNA-derived only
        ref = lin.trans_splicing_reference(PIPELINE, "horsetail")
        self.assertEqual([d.name for d in ep.raw_dirs(ref)], [ep.RNA_RAW_NAME])
        self.assertTrue((ref / "profiles" / "cox1" / "exon4.pssm").exists())
        self.assertTrue((ref / "profiles" / "nad1" / "whole_gene.pssm").exists())
        self.assertIn("RNA_Equisetum_sylvaticum", lp.rna_proteins("horsetail")["atp9"])
        have = {p.stem for p in (lin.editing_profiles(PIPELINE, "horsetail") / "mito").glob("*.pssm")}
        self.assertTrue({"atp9", "matR"} <= have)
        self.assertNotIn("cox1", have)  # three introns in Equisetum
        # one genus: a single locus is enough, with fewer pseudocounts
        self.assertEqual(ep.LINEAGE_MIN_LOCI["horsetail"], 1)
        self.assertLess(ep.LINEAGE_PSEUDOCOUNT["horsetail"], 1.0)

    def test_lycophyte_references(self):
        from poa.pipeline.editing import lineage_profiles as lp
        from poa.pipeline.trans_splicing import exon_profiles as ep
        from poa.pipeline.trans_splicing import reference_genes as rg

        self.assertEqual(len(rg.LINEAGE_MITOGENOMES["lycophyte"]), 4)
        # the RefSeq proteins are unedited: kept as a record, not built from
        self.assertEqual(set(rg.LINEAGE_UNEDITED_MITOGENOMES["lycophyte"]), set(rg.LINEAGE_MITOGENOMES["lycophyte"]))
        self.assertEqual(rg.LINEAGE_EDITED_MITOGENOMES["lycophyte"], {})
        ref = lin.trans_splicing_reference(PIPELINE, "lycophyte")
        self.assertEqual([d.name for d in ep.raw_dirs(ref)], ["raw", ep.RNA_RAW_NAME])
        raw = ep.load_all_raw_data(*ep.raw_dirs(ref))
        self.assertIn(("RNA_Phylloglossum_drummondii", "cox1", "cox1_1"), raw)
        self.assertIn("RNA_Phylloglossum_drummondii", lp.rna_proteins("lycophyte")["atp1"])
        ts = lin.trans_splicing_reference(PIPELINE, "lycophyte") / "profiles"
        self.assertTrue((ts / "atp9" / "exon4.pssm").exists())
        self.assertTrue((ts / "nad4" / "exon2.pssm").exists())
        self.assertFalse((ts / "nad4" / "exon3.pssm").exists())  # two exons in clubmosses
        self.assertFalse((ts / "nad7").exists())  # lost from their mitogenomes
        have = {p.stem for p in (lin.editing_profiles(PIPELINE, "lycophyte") / "mito").glob("*.pssm")}
        self.assertIn("atp1", have)
        self.assertFalse(set(rg.LINEAGE_EXTRA_SPLICED["lycophyte"]) & have)

    def test_single_sequence_profile(self):
        from poa.pipeline.common import pssm

        tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, tmp)
        columns = pssm.align_and_build({"only": "MKLV"}, tmp, "one", 0.5, pssm.composition(["MKLV"]), 0.25)
        self.assertEqual(len(columns), 4)
        self.assertGreater(columns[0]["M"], 0)
        self.assertLess(columns[0]["K"], 0)


class Liverwort(unittest.TestCase):
    """Liverwort: a lineage with references of its own, like moss, from RefSeq mitogenomes."""

    def test_is_a_full_lineage(self):
        self.assertIn("liverwort", lin.LINEAGES)
        self.assertNotIn("liverwort", lin.BORROWS)
        self.assertEqual(lin.reference("liverwort"), "liverwort")
        self.assertIn("liverwort", lin.CHOICES)
        self.assertIsNone(lin.stop_edit_penalty("liverwort"))
        self.assertIsNone(lin.edit_penalty("liverwort"))

    def test_spliced_genes(self):
        from poa.pipeline.trans_splicing import reference_genes as rg

        genes = rg.mito_spliced("liverwort")
        self.assertTrue({"atp1", "atp9", "cob", "cox1", "cox3", "nad3", "nad4L", "rps14"} <= set(genes))
        # a pseudogene in nearly every liverwort: no template
        self.assertNotIn("nad7", genes)
        self.assertIn("nad7", rg.mito_spliced("angiosperm"))
        self.assertIn("nad7", rg.mito_spliced())
        templates = lin.trans_splicing_reference(PIPELINE, "liverwort") / "profiles"
        for gene in genes:
            self.assertTrue((templates / gene / "whole_gene.pssm").exists(), gene)
        self.assertFalse((templates / "nad7").exists())
        # cox1 has ten exons in most liverworts; the others have lost introns and are not in the slots
        self.assertTrue((templates / "cox1" / "exon10.pssm").exists())
        self.assertEqual(scan_editing.lineage_spliced("liverwort"), set(rg.LINEAGE_EXTRA_SPLICED["liverwort"]))

    def test_profiles(self):
        from poa.pipeline.editing import lineage_profiles as lp
        from poa.pipeline.trans_splicing import reference_genes as rg

        have = {p.stem for p in (lin.editing_profiles(PIPELINE, "liverwort") / "mito").glob("*.pssm")}
        self.assertTrue({"atp4", "atp6", "ccmB", "nad6", "rps2", "sdh3", "sdh4"} <= have)
        # the genes liverworts splice are templates, not profiles; matR is not a liverwort gene
        self.assertFalse(have & set(rg.LINEAGE_EXTRA_SPLICED["liverwort"]))
        self.assertNotIn("matR", have)
        self.assertEqual(set(have) - set(lp.LIVERWORT_GENES), set())
        self.assertEqual(len(rg.LIVERWORT_MITOGENOMES), 27)


class EditCost(unittest.TestCase):
    """A lineage's own cost of a C-to-U edit (moss)."""

    def test_moss_only(self):
        self.assertEqual(lin.edit_penalty("moss"), 6.0)
        for name in ("angiosperm", "liverwort", "fern", "horsetail", "lycophyte", "auto"):
            self.assertIsNone(lin.edit_penalty(name))
        self.assertGreater(lin.edit_penalty("moss"), scan_editing.DEFAULT_EDIT_PENALTY)


class LineagePlastid(unittest.TestCase):
    """A lineage's own plastid references, used for its species on top of
    the land-plant plastid set: RNA-edited genes of DToL assemblies (fern,
    horsetail, lycophyte) or RefSeq plastomes (moss, liverwort)."""

    def test_where_they_live(self):
        from poa.pipeline.trans_splicing import reference_genes as rg

        self.assertEqual(set(rg.LINEAGE_PLASTID), {"fern", "horsetail", "lycophyte", "moss", "liverwort"})
        self.assertEqual(set(rg.LINEAGE_PLASTOMES), {"moss", "liverwort"})
        shared = PIPELINE / "trans_splicing" / "reference" / "plastid"
        self.assertEqual(rg.reference_dir(PIPELINE, "pltd"), shared)
        self.assertEqual(rg.reference_dir(PIPELINE, "pltd", "angiosperm"), shared)
        for name in rg.LINEAGE_PLASTID:
            own = rg.reference_dir(PIPELINE, "pltd", name)
            self.assertEqual(own, lin.trans_splicing_reference(PIPELINE, name) / "plastid")
            self.assertTrue((own / "profiles" / "rpl2" / "whole_gene.pssm").exists(), name)
            # RefSeq loci in raw/, RNA-edited ones in raw_rna/
            source = "raw" if name in rg.LINEAGE_PLASTOMES else "raw_rna"
            self.assertTrue((own / source / "junctions.tsv").exists(), name)
            self.assertFalse((own / "genbank").exists() and name not in rg.LINEAGE_PLASTOMES)

    def test_profiles_and_proteins(self):
        from poa.pipeline.editing import lineage_profiles as lp

        have = {p.stem for p in (lin.editing_profiles(PIPELINE, "fern") / "pltd").glob("*.pssm")}
        self.assertTrue({"rbcL", "atpA", "ndhD"} <= have)
        proteins = lp.rna_proteins("fern", organelle="pltd")
        self.assertIn("RNA_Onoclea_sensibilis", proteins["rbcL"])
        self.assertNotIn("rbcL", lp.rna_proteins("fern"))  # the mito set is a separate file
        self.assertEqual(lp.rna_proteins("moss", organelle="pltd"), {})
        with self.assertRaises(ValueError):
            lp.build_plastid_profiles("auto", Path("."))

    def test_refseq_plastid_proteins(self):
        from poa.pipeline.editing import lineage_profiles as lp
        from poa.pipeline.trans_splicing import reference_genes as rg

        self.assertEqual(lp.refseq_plastid_proteins("fern"), {})
        for name, plastomes in rg.LINEAGE_PLASTOMES.items():
            proteins = lp.refseq_plastid_proteins(name)
            self.assertTrue(set(proteins) <= set(lp.plastid_profile_genes()))
            self.assertTrue(set(proteins["rbcL"]) <= set(plastomes))
            self.assertGreater(len(proteins["rbcL"]), len(plastomes) // 2)
            acc = next(iter(proteins["rbcL"]))
            self.assertNotIn(acc, lp.refseq_plastid_proteins(name, (acc,))["rbcL"])
            have = {p.stem for p in (lin.editing_profiles(PIPELINE, name) / "pltd").glob("*.pssm")}
            self.assertTrue({"rbcL", "atpA", "psbB"} <= have, name)

    def test_reconstruct_uses_them_for_that_lineage_only(self):
        from poa.pipeline.trans_splicing import reconstruct
        from poa.pipeline.trans_splicing import reference_genes as rg

        # reconstruct's own pipeline dir: with poa installed, another test
        # module may have imported it from site-packages rather than src/
        own = rg.reference_dir(reconstruct.PIPELINE_DIR, "pltd", "fern") / "profiles" / "ndhB"
        shared = rg.reference_dir(reconstruct.PIPELINE_DIR, "pltd") / "profiles" / "ndhB"
        self.assertTrue(reconstruct.on_lineage_plastid_templates({"template_dir": str(own)}, "fern"))
        self.assertFalse(reconstruct.on_lineage_plastid_templates({"template_dir": str(shared)}, "fern"))
        self.assertFalse(reconstruct.on_lineage_plastid_templates({"template_dir": str(own)}, "angiosperm"))

    def test_no_recurrence_count_on_own_plastid_profiles(self):
        shifts = pd.DataFrame(
            [("Fern", "pltd", "rpoC2", 100), ("Fern", "pltd", "rpoB", 50)],
            columns=["species", "organelle", "gene", "profile_column"],
        )
        out = scan_editing.add_recurrence(shifts, {"Fern": "fern"}, own_plastid={("Fern", "rpoC2")})
        self.assertEqual(out["dtol_assemblies"].iloc[0], "")
        self.assertNotEqual(out["dtol_assemblies"].iloc[1], "")


class References(unittest.TestCase):
    """The bundled moss references exist where the stages look for them."""

    def test_bundled(self):
        ts = lin.trans_splicing_reference(PIPELINE, "moss") / "profiles"
        self.assertTrue((ts / "nad1" / "exon3.pssm").exists())
        self.assertFalse((ts / "nad1" / "exon4.pssm").exists())  # moss nad1 has 3 exons
        self.assertFalse((ts / "rps10").exists())  # lost from moss mitogenomes
        ed = lin.editing_profiles(PIPELINE, "moss") / "mito"
        self.assertEqual(sorted(p.stem for p in ed.glob("*.pssm")), sorted(GENES))
        self.assertEqual(
            lin.trans_splicing_reference(PIPELINE, "angiosperm"), PIPELINE / "trans_splicing" / "reference"
        )


if __name__ == "__main__":
    unittest.main()
