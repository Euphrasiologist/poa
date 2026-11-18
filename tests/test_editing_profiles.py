"""RefSeq profile fallback in editing.scan and the pinned lineage gene lists."""

import tempfile
import unittest
from pathlib import Path

from poa.pipeline.common import lineage as lin
from poa.pipeline.editing import lineage_profiles as lp
from poa.pipeline.editing.scan import profile_paths


class ProfilePaths(unittest.TestCase):
    def test_first_directory_wins_and_others_fill_in(self):
        with tempfile.TemporaryDirectory() as tmp:
            dataset, refseq = Path(tmp) / "dataset", Path(tmp) / "refseq"
            for d, genes in ((dataset, ["atp4", "nad9"]), (refseq, ["atp4", "rps12"])):
                (d / "mito").mkdir(parents=True)
                for g in genes:
                    (d / "mito" / f"{g}.pssm").write_text("")
            found = profile_paths([dataset, refseq], "mito")
            self.assertEqual(sorted(found), ["atp4", "nad9", "rps12"])
            self.assertEqual(found["atp4"], dataset / "mito" / "atp4.pssm")
            self.assertEqual(found["rps12"], refseq / "mito" / "rps12.pssm")
            self.assertEqual(profile_paths([dataset, refseq], "pltd"), {})


class BundledProfiles(unittest.TestCase):
    def test_moss_genes_pinned(self):
        # the lineage call was calibrated on these; adding genes changes it
        self.assertEqual(lp.MOSS_GENES, ["atp4", "ccmB", "ccmC", "mttB", "nad5", "nad9", "rps14"])
        moss = sorted(p.stem for p in (lin.editing_profiles(lp.PIPELINE_DIR, "moss") / "mito").glob("*.pssm"))
        self.assertEqual(moss, ["atp4", "ccmB", "ccmC", "mttB", "rps14"])

    def test_refseq_profiles_bundled_for_single_exon_genes(self):
        bundled = {p.stem for p in (lp.refseq_profiles_dir() / "mito").glob("*.pssm")}
        self.assertEqual(len(bundled), 26)
        self.assertTrue({"atp6", "atp9", "cox1", "rps12"} <= bundled)
        self.assertTrue(bundled <= set(lp.REFSEQ_GENES))


class CopiesToCorrect(unittest.TestCase):
    def calls(self, rows):
        import pandas as pd

        cols = ["species", "organelle", "contig_id", "start", "end", "gene", "score", "hmm_from", "hmm_to", "model_len"]
        return pd.DataFrame(rows, columns=cols)

    def test_best_plus_substantial_non_overlapping_copies(self):
        from poa.pipeline.editing.scan import hits_to_correct

        # Arabidopsis atp6: two full copies and a 21% fragment
        calls = self.calls(
            [
                ["Ath", "mito", "c1", 115747, 116672, "atp6", 1039.5, 157, 1026, 1532],
                ["Ath", "mito", "c1", 240429, 241852, "atp6", 1072.8, 60, 1026, 1532],
                ["Ath", "mito", "c1", 292280, 292630, "atp6", 88.4, 1108, 1424, 1532],
                ["Ath", "mito", "c1", 240500, 241800, "atp6", 900.0, 70, 1000, 1532],  # overlaps the best
                ["Ath", "pltd", "c2", 1, 900, "atp6", 999.0, 1, 1500, 1532],  # other organelle
            ]
        )
        got = hits_to_correct(calls, "mito", ["atp6"])
        self.assertEqual(sorted(got.start), [115747, 240429])

    def test_paralogue_locus_not_a_copy(self):
        from poa.pipeline.editing.scan import hits_to_correct

        # plastid psaA's HMM also hits psaB, near full length but far weaker
        calls = self.calls(
            [
                ["Ath", "mito", "c1", 37467, 39672, "psaB", 2724.4, 1, 2205, 2205],
                ["Ath", "mito", "c1", 37604, 39658, "psaA", 434.0, 107, 2111, 2253],
                ["Ath", "mito", "c1", 39697, 41950, "psaA", 2802.8, 1, 2253, 2253],
            ]
        )
        self.assertEqual(list(hits_to_correct(calls, "mito", ["psaA"]).start), [39697])

    def test_plastid_best_copy_only(self):
        from poa.pipeline.editing.scan import hits_to_correct

        calls = self.calls(
            [
                ["Ath", "pltd", "c1", 97548, 98016, "rps7", 599.5, 1, 468, 468],
                ["Ath", "pltd", "c1", 140775, 141243, "rps7", 599.4, 1, 468, 468],
            ]
        )
        self.assertEqual(list(hits_to_correct(calls, "pltd", ["rps7"]).start), [97548])

    def test_weak_extra_copy_dropped_with_its_edits(self):
        import pandas as pd

        from poa.pipeline.editing.scan import drop_weak_copies

        base = {"species": "Ath", "organelle": "mito", "gene": "atp6", "contig": "c1", "n_columns_total": "100"}
        calls = pd.DataFrame(
            [
                {**base, "hit_start": "1", "hit_end": "900", "score": "700", "n_columns_aligned": "60"},
                {**base, "hit_start": "2000", "hit_end": "2900", "score": "650", "n_columns_aligned": "95"},
                {**base, "hit_start": "5000", "hit_end": "5900", "score": "300", "n_columns_aligned": "70"},
            ]
        )
        edits = calls[["species", "organelle", "gene", "contig", "hit_start", "hit_end"]].assign(genomic_pos="5")
        kept, kept_edits = drop_weak_copies(calls, edits)
        self.assertEqual(list(kept.hit_start), ["1", "2000"])  # best kept despite low coverage
        self.assertEqual(sorted(kept_edits.hit_start), ["1", "2000"])

    def test_short_best_hit_still_kept(self):
        from poa.pipeline.editing.scan import hits_to_correct

        calls = self.calls([["Ath", "mito", "c1", 1, 300, "rps12", 50.0, 1, 100, 1000]])
        self.assertEqual(list(hits_to_correct(calls, "mito", ["rps12"]).start), [1])

    def test_copy_labels_round_trip(self):
        import pandas as pd

        from poa.pipeline.editing.scan import copy_label, split_copy_labels

        df = pd.DataFrame(
            {
                "species": ["Ath"],
                "organelle": ["mito"],
                "gene": [copy_label("atp6", "ctg000001c", 240429, 241852)],
                "score": ["688.4"],
            }
        )
        out = split_copy_labels(df)
        self.assertEqual(list(out.columns), ["species", "organelle", "gene", "contig", "hit_start", "hit_end", "score"])
        self.assertEqual(
            list(out.iloc[0][["gene", "contig", "hit_start", "hit_end"]]), ["atp6", "ctg000001c", "240429", "241852"]
        )
