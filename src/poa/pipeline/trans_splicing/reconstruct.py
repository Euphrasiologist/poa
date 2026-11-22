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

Plastid genes with introns (reference_genes.PLASTID_SPLICED: atpF, clpP,
ndhA, ndhB, rpl2, rpoC1, rps12, rps16, ycf3) are reconstructed the same
way from one template set built from 13 RefSeq plastomes across land
plants (reference/plastid/).

With transsplice >= 0.3.0 the exons and the joined CDS may change reading
frame (--frameshift-penalty), and each frameshift is a row of
reconstructed_frameshifts.tsv: exon slot, the skipped bases' genomic
start/end, n_bases (1 or 2) and the whole-gene profile column.

Output: analysis/trans_splicing/results/reconstructed_genes.tsv
        analysis/trans_splicing/results/reconstructed_exons.tsv
        analysis/trans_splicing/results/reconstructed_junctions.tsv
        analysis/trans_splicing/results/reconstructed_edits.tsv
        analysis/trans_splicing/results/reconstructed_frameshifts.tsv

Library use: reconstruct(root, ...) - everything passed in, StageError on
failure. src/01_reconstruct.py is the command-line wrapper (main() here).
"""

from __future__ import annotations

import argparse
import logging
import math
import os
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from poa import PIPELINE_DIR
from poa.pipeline.common import io_utils
from poa.pipeline.common import lineage as lin
from poa.pipeline.common import species_discovery as sd
from poa.pipeline.common.errors import StageError
from poa.pipeline.common.tools import cargo_tool
from poa.pipeline.editing.scan import add_recurrence, supports_option
from poa.pipeline.trans_splicing.reference_genes import GENES, LINEAGE_MITOGENOMES, PLASTID_SPLICED, reference_dir
from poa.pipeline.common.log import setup_logging

log = logging.getLogger(__name__)

DEFAULT_THREADS = 4
DEFAULT_MERGE_DISTANCE = 5000
DEFAULT_JUNCTION_DISTANCE_THRESHOLD = 10000
DEFAULT_EXON_FLANK = 60
OUT_NAMES = [
    "reconstructed_genes",
    "reconstructed_exons",
    "reconstructed_junctions",
    "reconstructed_edits",
    "reconstructed_frameshifts",
]
ORGANELLES = ("mito", "pltd")
# transsplice's cost of a frame change inside a gene (as editing's orfedit)
DEFAULT_FRAMESHIFT_PENALTY = 20.0
FRAMESHIFT_COLUMNS = [
    "species", "organelle", "gene", "slot", "contig", "start", "end", "n_bases", "profile_column", "window_start", "window_end"
]  # fmt: skip


@dataclass(frozen=True)
class ReconstructResult:
    tables: dict[str, Path]  # OUT_NAMES -> results/<name>.tsv


def template_lineage(name: str) -> str:
    """The lineage whose exon templates `name` uses: its own if it has any
    (fern has, though it borrows the angiosperm editing profiles), else
    those of the lineage it borrows from."""
    own = lin.trans_splicing_reference(PIPELINE_DIR, name) / "profiles"
    return name if name in LINEAGE_MITOGENOMES and any(own.glob("*/whole_gene.pssm")) else lin.reference(name)


def build_manifest(
    root: Path,
    gene_calls: pd.DataFrame,
    species_filter: set | None,
    lineage_of: dict[str, str],
    default_lineage: str,
    organelles: list[str] | tuple[str, ...] = ("mito",),
) -> list[dict]:
    """organelles: mito genes (GENES) take their species' lineage's
    templates, plastid ones (PLASTID_SPLICED) the one plastid set."""
    rows = []
    templates = {}
    # every lineage with reference genomes; one that borrows another's
    # references (fern, for editing) still has its own exon templates
    for name in LINEAGE_MITOGENOMES:
        profiles_dir = lin.trans_splicing_reference(PIPELINE_DIR, name) / "profiles"
        templates[name] = {g: profiles_dir / g for g in GENES if (profiles_dir / g / "whole_gene.pssm").exists()}
    plastid_dir = reference_dir(PIPELINE_DIR, "pltd") / "profiles"
    plastid = {g: plastid_dir / g for g in PLASTID_SPLICED if (plastid_dir / g / "whole_gene.pssm").exists()}
    missing = set(GENES) - set(templates[lin.DEFAULT_LINEAGE])
    if "pltd" in organelles:
        missing |= {f"pltd {g}" for g in set(PLASTID_SPLICED) - set(plastid)}
    if missing:
        log.warning(f"no template for {sorted(missing)} - run trans_splicing/src/00_build_exon_profiles.py first")

    for organelle in organelles:
        genes = PLASTID_SPLICED if organelle == "pltd" else GENES
        sub = gene_calls[(gene_calls["organelle"] == organelle) & (gene_calls["gene"].isin(genes))]
        data_root = sd.repo_data_root(root, organelle)
        if not data_root.is_dir():
            continue
        rows += manifest_rows(sub, organelle, data_root, species_filter, lambda species, gene: (
            plastid.get(gene)
            if organelle == "pltd"
            else templates[template_lineage(lineage_of.get(species, default_lineage))].get(gene)
        ))
    return rows


def manifest_rows(sub: pd.DataFrame, organelle: str, data_root: Path, species_filter: set | None, template_of) -> list:
    rows = []
    fasta_cache: dict[str, str | None] = {}
    for row in sub.itertuples():
        if species_filter is not None and row.species not in species_filter:
            continue
        template = template_of(row.species, row.gene)
        if template is None:
            continue
        if row.species not in fasta_cache:
            r = sd.resolve_species(data_root / row.species, organelle)
            fasta_cache[row.species] = r.ctg_fasta if r.status == "ok" else None
        fasta_path = fasta_cache[row.species]
        if fasta_path is None:
            continue
        rows.append(
            {
                "species": row.species,
                "organelle": organelle,
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


def supports_frameshifts(transsplice: str) -> bool:
    """transsplice >= 0.3.0 (frame changes, --frameshifts-out)."""
    out = subprocess.run([transsplice, "--version"], capture_output=True, text=True).stdout
    m = re.search(r"(\d+)\.(\d+)\.(\d+)", out)
    return m is not None and tuple(map(int, m.groups())) >= (0, 3, 0)


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
    organelles: list[str] | tuple[str, ...] = ORGANELLES,
    frameshift_penalty: float = DEFAULT_FRAMESHIFT_PENALTY,
) -> ReconstructResult:
    """Reconstruct the spliced genes under data root `root` (mito GENES,
    plastid PLASTID_SPLICED, for `organelles`).

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
        log.warning(
            f"no {lin.table_path(analysis_dir)} (run editing/src/01_scan_editing.py first) - "
            f"using {lin.DEFAULT_LINEAGE} templates for every species"
        )
        lineage_of = {}
    default_lineage = lin.DEFAULT_LINEAGE if lineage == "auto" else lineage

    rows = build_manifest(root, calls, species, lineage_of, default_lineage, organelles)
    if not rows:
        raise StageError("empty manifest - nothing to reconstruct")

    work_dir = analysis_dir / "trans_splicing" / "work"
    results_dir = analysis_dir / "trans_splicing" / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    # transsplice writes these paths itself (truncating), so a species-scoped
    # run must not point it straight at results_dir - it'd truncate the
    # shared tables down to just the scoped species. Write to work_dir
    # instead, then merge_species_scoped_tsv decides how that relates to
    # what's already in results_dir (see its docstring).
    work_outs = {name: work_dir / f"{name}.new.tsv" for name in OUT_NAMES}
    if math.isfinite(frameshift_penalty) and not supports_frameshifts(transsplice):
        log.warning(f"reconstruct: {transsplice} is older than 0.3.0 - no frameshift calls")
        frameshift_penalty = math.inf

    # One transsplice run per editing model: species of a lineage with a
    # stop read-through cost (fern) apart from the rest, as the cost is an
    # option of the whole batch.
    # Other U-to-C edits (lin.reverse_edit_penalty) only for mito genes on
    # the lineage's own templates, so those go in a batch apart from the
    # species' plastid genes.
    groups: dict[tuple[float | None, float | None], list[dict]] = {}
    for row in rows:
        name = lineage_of.get(row["species"], default_lineage)
        own = row["organelle"] == "mito" and template_lineage(name) == name
        groups.setdefault((lin.stop_edit_penalty(name), lin.reverse_edit_penalty(name) if own else None), []).append(row)
    for i, option in enumerate(("--stop-edit-penalty", "--reverse-edit-penalty")):
        if any(key[i] is not None for key in groups) and not supports_option(transsplice, option):
            raise StageError(f"the fern lineage needs a transsplice with {option}; {transsplice} has none")
    parts: dict[str, list[pd.DataFrame]] = {name: [] for name in OUT_NAMES}
    order = sorted(groups.items(), key=lambda g: tuple((v is not None, v or 0) for v in g[0]))
    for n, ((stop, reverse), group_rows) in enumerate(order):
        tag = "" if len(groups) == 1 else f".{n}"
        manifest_path = work_dir / f"manifest{tag}.tsv"
        write_manifest(group_rows, manifest_path)
        log.info(
            f"reconstruct: manifest has {len(group_rows)} fragment rows -> {manifest_path}"
            + (f" (stop read-through at {stop})" if stop is not None else "")
            + (f" (other U-to-C edits at {reverse})" if reverse is not None else "")
        )
        outs = {name: work_dir / f"{name}.new{tag}.tsv" for name in OUT_NAMES}
        cmd = [
            transsplice,
            "scan-batch",
            "--manifest",
            str(manifest_path),
            "--genes-out",
            str(outs["reconstructed_genes"]),
            "--exons-out",
            str(outs["reconstructed_exons"]),
            "--junctions-out",
            str(outs["reconstructed_junctions"]),
            "--edits-out",
            str(outs["reconstructed_edits"]),
            "--threads",
            str(threads),
            "--merge-distance",
            str(merge_distance),
            "--junction-distance-threshold",
            str(junction_distance_threshold),
            "--exon-flank",
            str(exon_flank),
        ]
        if math.isfinite(frameshift_penalty):
            cmd += [
                "--frameshift-penalty",
                str(frameshift_penalty),
                "--frameshifts-out",
                str(outs["reconstructed_frameshifts"]),
            ]
        else:
            pd.DataFrame(columns=FRAMESHIFT_COLUMNS).to_csv(outs["reconstructed_frameshifts"], sep="\t", index=False)
        if stop is not None:
            cmd += ["--stop-edit-penalty", str(stop)]
        if reverse is not None:
            cmd += ["--reverse-edit-penalty", str(reverse)]
        proc = subprocess.run(cmd)
        if proc.returncode != 0:
            raise StageError(f"transsplice scan-batch failed (exit {proc.returncode})")
        for name in OUT_NAMES:
            parts[name].append(pd.read_csv(outs[name], sep="\t", dtype=str))
    if len(groups) > 1:
        for name in OUT_NAMES:
            pd.concat(parts[name], ignore_index=True).to_csv(work_outs[name], sep="\t", index=False)

    shifts = pd.read_csv(work_outs["reconstructed_frameshifts"], sep="\t")
    lineages = {sp: lineage_of.get(sp, default_lineage) for sp in shifts["species"].unique()}
    add_recurrence(shifts, lineages, source="reconstruction", reference_of=template_lineage).to_csv(
        work_outs["reconstructed_frameshifts"], sep="\t", index=False
    )
    tables = {name: results_dir / f"{name}.tsv" for name in OUT_NAMES}
    for name in OUT_NAMES:
        io_utils.merge_species_scoped_tsv(tables[name], work_outs[name], species)
    log.info(f"reconstruct: merged into {results_dir} ({'species-scoped' if species else 'full rewrite'})")
    return ReconstructResult(tables=tables)


def main(argv: list[str] | None = None) -> int:
    setup_logging()
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
        log.error(e)
        return 1
    return 0
