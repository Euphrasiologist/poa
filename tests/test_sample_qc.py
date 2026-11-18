"""Single-sample QC: the dataset gate's rules on one assembly."""

import tempfile
import unittest
from pathlib import Path

from poa.pipeline.qc_basic_stats import sample_qc as q
from poa.pipeline.qc_basic_stats.qc_reference import poa_name

GFATK_TWO_SUBGRAPHS = """Subgraph 1:
	Number of nodes/segments: 3
	Number of edges/links: 4
	Circular: false
	Dead-end nodes: 1

	Total sequence length:	129887
	Average coverage of total segments:	5543
Subgraph 2:
	Number of nodes/segments: 1
	Number of edges/links: 1
	Circular: true
	Dead-end nodes: 0

	Total sequence length:	20000
	Average coverage of total segments:	4000
Total number of subgraphs: 2
"""


class Status(unittest.TestCase):
    def test_clean_assembly_passes(self):
        row = {"organelle": "mito", "parse_status": "ok", "dead_end_nodes": 0, "n_fragment_subgraphs": 1,
               "any_circular": True, "core_gene_pct": 1.0, "size_mad_z": 0.5, "has_ctg_fasta": True}  # fmt: skip
        self.assertEqual(q.build_status(row), ("pass", ""))

    def test_thresholds_and_reason_strings(self):
        row = {"organelle": "mito", "parse_status": "ok", "dead_end_nodes": 3, "n_fragment_subgraphs": 5,
               "any_circular": False, "core_gene_pct": 0.7, "size_mad_z": -4.5, "has_ctg_fasta": True,
               "cross_organelle_aligned_frac": 0.2}  # fmt: skip
        self.assertEqual(
            q.build_status(row),
            (
                "fail",
                "dead_end_nodes>2;fragmented(n_fragment_subgraphs>=5);non_circular;size_outlier(|z|>4);"
                "core_gene_pct<0.80;elevated_cross_organelle_alignment(>15%)",
            ),
        )

    def test_complete_plastid_exempt_from_cross_organelle(self):
        row = {"organelle": "pltd", "any_circular": True, "core_gene_pct": 0.95, "cross_organelle_aligned_frac": 0.5}
        self.assertEqual(q.build_status(row)[0], "pass")
        row["core_gene_pct"] = 0.85
        self.assertEqual(q.build_status(row), ("fail", "high_cross_organelle_alignment(>35%)"))

    def test_unjoined_and_empty(self):
        self.assertEqual(q.build_status({"unjoined": True}), ("flag", "no_resolved_ctg_fasta(unjoined)"))
        self.assertEqual(q.build_status({"parse_status": "empty"}), ("fail", "empty_assembly"))

    def test_nothing_measured_is_not_a_fail(self):
        # poa run without a GFA: the graph checks are simply absent
        self.assertEqual(q.build_status({"organelle": "mito", "has_ctg_fasta": True})[0], "pass")


class Measurements(unittest.TestCase):
    def test_closed_circle_at_genome_depth_is_not_a_fragment(self):
        m = q.graph_metrics(*q.parse_gfatk_stats(GFATK_TWO_SUBGRAPHS))
        self.assertEqual((m["n_subgraphs"], m["n_fragment_subgraphs"], m["dead_end_nodes"]), (2, 1, 1))
        self.assertEqual(m["total_length"], 149887)

    def test_contig_headers(self):
        with tempfile.TemporaryDirectory() as tmp:
            fa = Path(tmp) / "x.fa"
            fa.write_text("".join(f">ctg{i} length=4 wlength=1 nv=1 circular=false\nACGT\n" for i in range(5)))
            m = q.contig_metrics(fa)
            self.assertEqual((m["n_contigs"], m["assembly_length"], m["unjoined"]), (5, 20, True))
            fa.write_text(">a\nACGT\n>b\nAC\n")  # not oatk: circularity unknown
            m = q.contig_metrics(fa)
            self.assertEqual((m["assembly_length"], m["any_circular"], m["unjoined"]), (6, None, False))

    def test_core_genes_count_split_ccmfn(self):
        pct, missing = q.core_gene_completeness({"ccmFn1", "ccmFn2", "cox1"}, {"ccmFn", "cox1", "atp1"})
        self.assertAlmostEqual(pct, 2 / 3)
        self.assertEqual(missing, ["atp1"])

    def test_size_z_genus_then_dataset(self):
        ref = {
            ("mito", "*"): {"n": 1000, "median": 400000.0, "mad": 100000.0},
            ("mito", "Silene"): {"n": 7, "median": 1200000.0, "mad": 0.0},
            ("mito", "Carlina"): {"n": 1, "median": 327000.0, "mad": 0.0},
        }
        z, basis = q.size_z(1200000 + 1.4826 * 100000, "mito", "Silene", ref)  # genus MAD 0 -> dataset MAD
        self.assertAlmostEqual(z, 1.0)
        self.assertEqual(basis, "genus Silene (n=7)")
        z, basis = q.size_z(400000, "mito", "Carlina", ref)  # 1 congener: whole dataset
        self.assertEqual((z, basis), (0.0, "dataset"))

    def test_oatk_names(self):
        self.assertEqual([poa_name(g) for g in ("trnA_UGC", "rrn4", "atp1")], ["trnA-UGC", "rrn4.5", "atp1"])


class Bundled(unittest.TestCase):
    def test_reference_loads(self):
        core = q.load_core_genes()
        self.assertEqual((len(core["mito"]), len(core["pltd"])), (38, 108))
        self.assertIn("trnA-UGC", core["pltd"])
        ref = q.load_size_reference()
        self.assertGreater(ref[("mito", "*")]["n"], 1000)
