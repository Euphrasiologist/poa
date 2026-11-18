"""Write sample_qc's reference data from a dataset's QC gate results.

    python -m poa.pipeline.qc_basic_stats.qc_reference RESULTS_DIR

RESULTS_DIR is the dataset repo's analysis/qc_basic_stats/results/
(core_genes_{mito,pltd}.tsv, contig_stats.tsv, qc_summary.tsv). Writes,
into this package's qc_basic_stats/reference/:

  core_genes.tsv      genes present in >= 90% of assemblies (the dataset
                      gate's QC core set), renamed to poa's gene names:
                      oatk's plastid tRNAs trnA_UGC -> trnA-UGC, rrn4 ->
                      rrn4.5.
  size_reference.tsv  per organelle, per genus and over the whole
                      dataset ("*"): the number of assemblies, median and
                      MAD of genome size, size being the resolved contigs'
                      total length where the assembly has resolved contigs,
                      else the graph's (05_qc_summary.size_mad_z).
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

from poa.pipeline.qc_basic_stats.sample_qc import QC_CORE_GENE_PRESENCE_THRESHOLD, REFERENCE_DIR

RENAME = {"rrn4": "rrn4.5"}


def poa_name(gene: str) -> str:
    return RENAME.get(gene, gene.replace("_", "-") if gene.startswith("trn") else gene)


def build_qc_reference(results_dir: Path, out_dir: Path = REFERENCE_DIR) -> None:
    results_dir, out_dir = Path(results_dir), Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    core = []
    for organelle in ("mito", "pltd"):
        cg = pd.read_csv(results_dir / f"core_genes_{organelle}.tsv", sep="\t")
        genes = cg.loc[cg["pct_present"] >= QC_CORE_GENE_PRESENCE_THRESHOLD, "gene"]
        core += [{"organelle": organelle, "gene": poa_name(g)} for g in sorted(genes)]
    pd.DataFrame(core).to_csv(out_dir / "core_genes.tsv", sep="\t", index=False)

    qc = pd.read_csv(results_dir / "qc_summary.tsv", sep="\t")
    cs = pd.read_csv(results_dir / "contig_stats.tsv", sep="\t")
    resolved = (
        cs[cs["contig_source"] == "ctg_fasta"]
        .groupby(["species", "organelle"])["length"]
        .sum()
        .rename("resolved_length")
        .reset_index()
    )
    qc = qc.merge(resolved, on=["species", "organelle"], how="left")
    has_ctg = qc["has_ctg_fasta"].astype(str).str.lower() == "true"
    qc["size"] = qc["total_length"].where(~(has_ctg & qc["resolved_length"].notna()), qc["resolved_length"])
    qc = qc[qc["size"].notna()]
    qc["genus"] = qc["species"].str.split("_").str[0]

    def stats(sizes: pd.Series) -> dict:
        med = sizes.median()
        return {"n": len(sizes), "median": med, "mad": (sizes - med).abs().median()}

    rows = []
    for organelle, g in qc.groupby("organelle"):
        rows.append({"organelle": organelle, "genus": "*", **stats(g["size"])})
        for genus, gg in g.groupby("genus"):
            rows.append({"organelle": organelle, "genus": genus, **stats(gg["size"])})
    pd.DataFrame(rows).to_csv(out_dir / "size_reference.tsv", sep="\t", index=False, float_format="%.1f")
    n_core = pd.DataFrame(core).groupby("organelle").size().to_dict()
    print(f"[info] qc_reference: core genes {n_core}, {len(rows)} size rows -> {out_dir}", file=sys.stderr)


if __name__ == "__main__":
    build_qc_reference(Path(sys.argv[1]))
