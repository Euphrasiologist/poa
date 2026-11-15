"""Run `transsplice` over every species' calls for the genes in
reference_genes.GENES - trans-spliced (nad1/nad2/nad5/rps3) and, since
2026-09-29, cis-spliced ones too (transsplice classifies each junction
cis/trans per species) - reconstructing full-length genes from their per-exon
`.ctg.bed` hits.

Unlike editing.scan (best-scoring hit only per gene), every
raw `.ctg.bed` row for these genes is a candidate exon fragment here - all
of them go into the manifest, and `transsplice` itself does the
merging/assignment (see its README).

Each species uses the exon templates of its lineage (common/lineage.py), as
recorded by editing.scan in editing/results/lineage.tsv:
mosses get reference/lineages/moss/profiles/, whose exon counts are the
moss ones (nad1 3 exons, not 5). A gene with no template in the species'
lineage (moss rps10, lost from moss mitogenomes) isn't reconstructed.
Without lineage.tsv, or with --lineage, one lineage is used for all.

Output: analysis/trans_splicing/results/reconstructed_genes.tsv
        analysis/trans_splicing/results/reconstructed_exons.tsv
        analysis/trans_splicing/results/reconstructed_junctions.tsv
        analysis/trans_splicing/results/reconstructed_edits.tsv

Library use: reconstruct(root, ...) - everything passed in, StageError on
failure. src/01_reconstruct.py is the command-line wrapper (main() here).
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from poa import PIPELINE_DIR
from poa.pipeline.common import io_utils
from poa.pipeline.common import lineage as lin
from poa.pipeline.common import species_discovery as sd
from poa.pipeline.common.errors import StageError
from poa.pipeline.common.tools import cargo_tool
from poa.pipeline.trans_splicing.reference_genes import GENES

DEFAULT_THREADS = 4
DEFAULT_MERGE_DISTANCE = 5000
DEFAULT_JUNCTION_DISTANCE_THRESHOLD = 10000
DEFAULT_EXON_FLANK = 60
OUT_NAMES = ["reconstructed_genes", "reconstructed_exons", "reconstructed_junctions", "reconstructed_edits"]


@dataclass(frozen=True)
class ReconstructResult:
    tables: dict[str, Path]  # OUT_NAMES -> results/<name>.tsv


def build_manifest(
    root: Path, gene_calls: pd.DataFrame, species_filter: set | None, lineage_of: dict[str, str], default_lineage: str
) -> list[dict]:
    rows = []
    templates = {}
    for name in lin.LINEAGES:
        profiles_dir = lin.trans_splicing_reference(PIPELINE_DIR, name) / "profiles"
        templates[name] = {g: profiles_dir / g for g in GENES if (profiles_dir / g / "whole_gene.pssm").exists()}
    missing = set(GENES) - set(templates[lin.DEFAULT_LINEAGE])
    if missing:
        print(
            f"[warn] no template for {sorted(missing)} - run trans_splicing/src/00_build_exon_profiles.py first",
            file=sys.stderr,
        )

    sub = gene_calls[(gene_calls["organelle"] == "mito") & (gene_calls["gene"].isin(GENES))]
    data_root = sd.repo_data_root(root, "mito")
    fasta_cache: dict[str, str | None] = {}
    for row in sub.itertuples():
        if species_filter is not None and row.species not in species_filter:
            continue
        template = templates[lineage_of.get(row.species, default_lineage)].get(row.gene)
        if template is None:
            continue
        if row.species not in fasta_cache:
            r = sd.resolve_species(data_root / row.species, "mito")
            fasta_cache[row.species] = r.ctg_fasta if r.status == "ok" else None
        fasta_path = fasta_cache[row.species]
        if fasta_path is None:
            continue
        rows.append(
            {
                "species": row.species,
                "organelle": "mito",
                "gene": row.gene,
                "fasta": fasta_path,
                "template_dir": str(template),
                "contig": row.contig_id,
                "start": int(row.start),
                "end": int(row.end),
                "strand": row.strand,
            }
        )
    return rows


def write_manifest(rows: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    cols = ["species", "organelle", "gene", "fasta", "template_dir", "contig", "start", "end", "strand"]
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")


def reconstruct(
    root: Path,
    *,
    gene_calls: Path | None = None,
    species: set[str] | None = None,
    threads: int = DEFAULT_THREADS,
    merge_distance: int = DEFAULT_MERGE_DISTANCE,
    junction_distance_threshold: int = DEFAULT_JUNCTION_DISTANCE_THRESHOLD,
    exon_flank: int = DEFAULT_EXON_FLANK,
    lineage: str = "auto",
) -> ReconstructResult:
    """Reconstruct the spliced mito genes under data root `root`.

    gene_calls: a gene_calls.tsv-shaped table (default: the data root's
    denovo_annotation/results/gene_calls.tsv). species: only these species,
    merged into the existing result tables. lineage: "auto" takes each
    species' templates from editing's lineage.tsv; any of
    common/lineage.LINEAGES forces it."""
    root = Path(root)
    analysis_dir = root / "analysis"
    if lineage not in lin.CHOICES:
        raise ValueError(f"unknown lineage {lineage!r}; expected one of {lin.CHOICES}")
    transsplice = cargo_tool("transsplice")

    gene_calls_path = (
        Path(gene_calls) if gene_calls else analysis_dir / "denovo_annotation" / "results" / "gene_calls.tsv"
    )
    if not gene_calls_path.exists():
        raise StageError(
            f"missing {gene_calls_path} - run denovo_annotation/src/04_build_gene_calls.py first (or pass --gene-calls)"
        )
    calls = pd.read_csv(gene_calls_path, sep="\t")

    lineage_of = lin.load(analysis_dir) if lineage == "auto" else {}
    if lineage_of is None:
        print(
            f"[warn] no {lin.table_path(analysis_dir)} (run editing/src/01_scan_editing.py first) - "
            f"using {lin.DEFAULT_LINEAGE} templates for every species",
            file=sys.stderr,
        )
        lineage_of = {}
    default_lineage = lin.DEFAULT_LINEAGE if lineage == "auto" else lineage

    rows = build_manifest(root, calls, species, lineage_of, default_lineage)
    if not rows:
        raise StageError("empty manifest - nothing to reconstruct")

    work_dir = analysis_dir / "trans_splicing" / "work"
    manifest_path = work_dir / "manifest.tsv"
    write_manifest(rows, manifest_path)
    print(f"[info] reconstruct: manifest has {len(rows)} fragment rows -> {manifest_path}", file=sys.stderr)

    results_dir = analysis_dir / "trans_splicing" / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    # transsplice writes these paths itself (truncating), so a species-scoped
    # run must not point it straight at results_dir - it'd truncate the
    # shared tables down to just the scoped species. Write to work_dir
    # instead, then merge_species_scoped_tsv decides how that relates to
    # what's already in results_dir (see its docstring).
    work_outs = {name: work_dir / f"{name}.new.tsv" for name in OUT_NAMES}
    cmd = [
        transsplice,
        "scan-batch",
        "--manifest",
        str(manifest_path),
        "--genes-out",
        str(work_outs["reconstructed_genes"]),
        "--exons-out",
        str(work_outs["reconstructed_exons"]),
        "--junctions-out",
        str(work_outs["reconstructed_junctions"]),
        "--edits-out",
        str(work_outs["reconstructed_edits"]),
        "--threads",
        str(threads),
        "--merge-distance",
        str(merge_distance),
        "--junction-distance-threshold",
        str(junction_distance_threshold),
        "--exon-flank",
        str(exon_flank),
    ]
    proc = subprocess.run(cmd)
    if proc.returncode != 0:
        raise StageError(f"transsplice scan-batch failed (exit {proc.returncode})")

    tables = {name: results_dir / f"{name}.tsv" for name in OUT_NAMES}
    for name in OUT_NAMES:
        io_utils.merge_species_scoped_tsv(tables[name], work_outs[name], species)
    print(
        f"[info] reconstruct: merged into {results_dir} ({'species-scoped' if species else 'full rewrite'})",
        file=sys.stderr,
    )
    return ReconstructResult(tables=tables)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species-list", default=None)
    ap.add_argument("--threads", type=int, default=DEFAULT_THREADS)
    ap.add_argument("--merge-distance", type=int, default=DEFAULT_MERGE_DISTANCE)
    ap.add_argument("--junction-distance-threshold", type=int, default=DEFAULT_JUNCTION_DISTANCE_THRESHOLD)
    ap.add_argument("--exon-flank", type=int, default=DEFAULT_EXON_FLANK)
    ap.add_argument(
        "--gene-calls",
        default=None,
        help="path to a gene_calls.tsv-shaped table to anchor "
        "on (default: denovo_annotation/results/gene_calls.tsv under the data root)",
    )
    ap.add_argument(
        "--lineage",
        default="auto",
        choices=lin.CHOICES,
        help="exon-template lineage: auto (default) takes each species' from editing's lineage.tsv; "
        "or force one for all",
    )
    ap.add_argument(
        "--force",
        action="store_true",
        help="accepted for CLI consistency; transsplice "
        "always does a full recompute over the manifest (no incremental skip yet)",
    )
    args = ap.parse_args(argv)
    root = Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT") or os.getcwd())
    try:
        reconstruct(
            root,
            gene_calls=Path(args.gene_calls) if args.gene_calls else None,
            species=sd.load_species_list(Path(args.species_list)) if args.species_list else None,
            threads=args.threads,
            merge_distance=args.merge_distance,
            junction_distance_threshold=args.junction_distance_threshold,
            exon_flank=args.exon_flank,
            lineage=args.lineage,
        )
    except StageError as e:
        print(f"[err] {e}", file=sys.stderr)
        return 1
    return 0
