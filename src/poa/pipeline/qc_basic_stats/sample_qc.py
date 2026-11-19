"""Pass/flag/fail QC for one sample's organelle assemblies, with the rules
the dataset gate (qc_basic_stats 01-05 in the dataset repo) applies, and
the dataset-derived parts of those rules shipped as fixed reference data.

The dataset gate judges an assembly partly against the rest of the
dataset: its core-gene set is the genes present in >= 90% of assemblies,
and genome size is compared with the same genus. Neither needs the
dataset at run time once it is computed, so qc_reference.py writes both
from a dataset's QC results into reference/:

  core_genes.tsv      organelle, gene (poa's names; 38 mito, 108 plastid)
  size_reference.tsv  organelle, genus ("*" = whole dataset), n, median, mad

Checks, thresholds and reason strings are the dataset gate's
(05_qc_summary.py), so a status here means the same as in qc_summary.tsv:

  graph (needs the GFA)
    dead-end nodes                          flag >0, fail >2
    fragment subgraphs (not counting closed
      circles near the median depth)        flag >=5, fail >=10
    empty graph                             fail
  contigs
    no circular contig                      flag (header circular=)
    every contig an unjoined graph segment
      (>=5 contigs, all linear nv=1)        flag no_resolved_ctg_fasta(unjoined)
  genes
    core genes present                      flag <0.80, fail <0.50
  size
    |MAD z| vs the genus (>=3 in the
      reference) or the whole dataset       flag >4
  both organelles given
    share of the assembly aligning to the
      other organelle (minimap2 -x asm5)    flag >0.15, fail >0.35
      (not applied to a complete plastid - circular, core genes >= 90%)

Differences from the dataset gate, all from having one sample:
  - No GFA is not a fail (`poa run` takes contigs alone); the graph checks
    are reported as not assessed.
  - Core genes are counted from poa's own gene calls (any source), not
    oatk's .ctg.bed; ccmFn1/ccmFn2 count as ccmFn.
  - The genus is the first word of the sample name unless given. A genus
    with fewer than 3 assemblies in the reference uses the whole-dataset
    median and MAD, as the dataset gate does.
  - A check whose tool (gfatk, minimap2) isn't installed is reported as
    not assessed instead of failing the run.

Reported but not part of the status: frameshifts found in edited genes
(editing/results/frameshifts.tsv) and reconstructed spliced genes
(trans_splicing/results/reconstructed_frameshifts.tsv) - n_frameshifts, and the genes
carrying them as gene:N, N being how many DToL assemblies have a
frameshift at the same place in that gene. A frameshift is either an
assembly error worth checking against the reads or real biology (a
pseudogene, a gene read in two frames), which no rule here can tell apart
for certain; one seen in at most NEW_FRAMESHIFT_MAX DToL assemblies counts
towards n_new_frameshifts, the ones most worth checking.

Output: analysis/qc_basic_stats/results/sample_qc.tsv, one row per
species x organelle (COLUMNS).

Library use: sample_qc(root, organelles, gene_calls=None, genus=None).
"""

from __future__ import annotations

import csv
import logging
import math
import re
import shutil
import subprocess
from pathlib import Path

import pandas as pd

from poa.pipeline.common import species_discovery as sd
from poa.pipeline.common.errors import StageError

log = logging.getLogger(__name__)

REFERENCE_DIR = Path(__file__).resolve().parent / "reference"
ORGANELLES = ("mito", "pltd")

# 05_qc_summary.py's thresholds
CROSS_ORG_FLAG_THRESHOLD = 0.15
CROSS_ORG_FAIL_THRESHOLD = 0.35
SIZE_Z_FLAG_THRESHOLD = 4.0
DEAD_END_FLAG_THRESHOLD = 0
DEAD_END_FAIL_THRESHOLD = 2
FRAGMENTATION_FLAG_THRESHOLD = 5
FRAGMENTATION_FAIL_THRESHOLD = 10
CHROMOSOME_DEPTH_BAND = 3.0
CORE_GENE_FLAG_THRESHOLD = 0.80
CORE_GENE_FAIL_THRESHOLD = 0.50
QC_CORE_GENE_PRESENCE_THRESHOLD = 0.90
UNJOINED_MIN_CONTIGS = 5
MIN_GENUS_N = 3
GENE_ALIASES = {"ccmFn1": "ccmFn", "ccmFn2": "ccmFn"}
NEW_FRAMESHIFT_MAX = 2

COLUMNS = [
    "species",
    "organelle",
    "status",
    "status_reasons",
    "not_assessed",
    "assembly_length",
    "n_contigs",
    "any_circular",
    "pct_unjoined_contigs",
    "n_subgraphs",
    "n_fragment_subgraphs",
    "dead_end_nodes",
    "core_gene_pct",
    "core_genes_missing",
    "size_mad_z",
    "size_reference",
    "cross_organelle_aligned_frac",
    "n_frameshifts",
    "n_new_frameshifts",
    "frameshifted_genes",
]

HEADER_RE = re.compile(r"^>(\S+).*?length=(\d+).*?nv=(\d+).*?circular=(true|false)")


# ---------------------------------------------------------------- reference


def load_core_genes(path: Path = REFERENCE_DIR / "core_genes.tsv") -> dict[str, set[str]]:
    core: dict[str, set[str]] = {o: set() for o in ORGANELLES}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            core[r["organelle"]].add(r["gene"])
    return core


def load_size_reference(path: Path = REFERENCE_DIR / "size_reference.tsv") -> dict[tuple[str, str], dict]:
    with open(path) as fh:
        return {
            (r["organelle"], r["genus"]): {"n": int(r["n"]), "median": float(r["median"]), "mad": float(r["mad"])}
            for r in csv.DictReader(fh, delimiter="\t")
        }


# ---------------------------------------------------------------- measurements


def parse_gfatk_stats(text: str) -> tuple[int | None, list[dict]]:
    """`gfatk stats` output -> (number of subgraphs, per-subgraph metrics)."""
    total = re.search(r"Total number of subgraphs:\s*(\d+)", text)
    subgraphs = []
    for block in re.split(r"Subgraph \d+:", text)[1:]:

        def grab(pattern, cast):
            m = re.search(pattern, block)
            return cast(m.group(1)) if m else None

        subgraphs.append(
            {
                "circular": grab(r"Circular:\s*(\w+)", lambda x: x.lower() == "true"),
                "dead_end_nodes": grab(r"Dead-end nodes:\s*(\d+)", int),
                "total_length": grab(r"Total sequence length:\s*(\d+)", int),
                "avg_coverage": grab(r"Average coverage of total segments:\s*([\d.eE+-]+)", float),
            }
        )
    return (int(total.group(1)) if total else None), subgraphs


def graph_metrics(n_subgraphs: int | None, subgraphs: list[dict]) -> dict:
    """The dataset gate's per-assembly graph aggregate: total length, worst
    dead-end count, and subgraphs that count as fragments (anything but a
    closed circle within CHROMOSOME_DEPTH_BAND of the length-weighted
    median depth - a chromosome of a multichromosomal genome)."""
    if n_subgraphs is None:
        return {"parse_status": "error"}
    if n_subgraphs == 0:
        return {"parse_status": "empty", "n_subgraphs": 0}
    out = {
        "parse_status": "ok",
        "n_subgraphs": n_subgraphs,
        "total_length": sum(s["total_length"] or 0 for s in subgraphs),
        "dead_end_nodes": max(
            (s["dead_end_nodes"] for s in subgraphs if s["dead_end_nodes"] is not None), default=None
        ),
        "n_fragment_subgraphs": None,
    }
    depth = sorted(
        (s for s in subgraphs if s["avg_coverage"] is not None and s["total_length"] is not None),
        key=lambda s: s["avg_coverage"],
    )
    if depth:
        half, cum, ref = sum(s["total_length"] for s in depth) / 2, 0, depth[-1]["avg_coverage"]
        for s in depth:
            cum += s["total_length"]
            if cum >= half:
                ref = s["avg_coverage"]
                break
        chromosomes = sum(
            1
            for s in subgraphs
            if s["circular"]
            and s["avg_coverage"] is not None
            and ref / CHROMOSOME_DEPTH_BAND <= s["avg_coverage"] <= ref * CHROMOSOME_DEPTH_BAND
        )
        out["n_fragment_subgraphs"] = len(subgraphs) - chromosomes
    return out


def gfa_stats(gfa: Path, gfatk: str) -> dict:
    proc = subprocess.run([gfatk, "stats", str(gfa)], capture_output=True, text=True, timeout=300)
    if proc.returncode != 0:
        return {"parse_status": "error"}
    return graph_metrics(*parse_gfatk_stats(proc.stdout))


def contig_metrics(fasta: Path) -> dict:
    """Contig count, circularity and unjoined segments from oatk-style
    headers (length= nv= circular=); sequence length counted directly.
    Headers without those fields leave circularity/unjoined unknown."""
    contigs: list[dict] = []
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                m = HEADER_RE.match(line)
                contigs.append(
                    {"nv": int(m.group(3)), "circular": m.group(4) == "true", "length": 0}
                    if m
                    else {"nv": None, "circular": None, "length": 0}
                )
            elif contigs:
                contigs[-1]["length"] += len(line.strip())
    n = len(contigs)
    known = [c for c in contigs if c["nv"] is not None]
    out = {
        "n_contigs": n,
        "assembly_length": sum(c["length"] for c in contigs),
        "any_circular": any(c["circular"] for c in known) if known and len(known) == n else None,
        "pct_unjoined_contigs": None,
        "unjoined": False,
    }
    if known and len(known) == n:
        unjoined = sum(1 for c in contigs if c["nv"] == 1 and not c["circular"])
        out["pct_unjoined_contigs"] = unjoined / n
        out["unjoined"] = n >= UNJOINED_MIN_CONTIGS and unjoined == n
    return out


def core_gene_completeness(present: set[str], core: set[str]) -> tuple[float | None, list[str]]:
    if not core:
        return None, []
    present = {GENE_ALIASES.get(g, g) for g in present}
    missing = sorted(core - present)
    return (len(core) - len(missing)) / len(core), missing


def size_z(length: float | None, organelle: str, genus: str | None, ref: dict) -> tuple[float | None, str]:
    if not length:
        return None, ""
    overall = ref.get((organelle, "*"))
    entry, basis = overall, "dataset"
    g = ref.get((organelle, genus)) if genus else None
    if g and g["n"] >= MIN_GENUS_N:
        entry, basis = g, f"genus {genus} (n={g['n']})"
        if entry["mad"] == 0:
            entry = {**entry, "mad": overall["mad"]}
    if not entry or not entry["mad"]:
        return 0.0, basis
    return (length - entry["median"]) / (1.4826 * entry["mad"]), basis


def merged_aligned_length(intervals: dict[str, list[tuple[int, int]]]) -> int:
    total = 0
    for spans in intervals.values():
        spans.sort()
        cur_start, cur_end = spans[0]
        for start, end in spans[1:]:
            if start <= cur_end:
                cur_end = max(cur_end, end)
            else:
                total += cur_end - cur_start
                cur_start, cur_end = start, end
        total += cur_end - cur_start
    return total


def cross_organelle(mito_fasta: Path, pltd_fasta: Path, minimap2: str) -> tuple[float, float] | None:
    """(share of the mito assembly aligned to the plastid, share of the
    plastid aligned to the mito), minimap2 -x asm5 with the mito as
    reference, overlapping blocks merged per sequence."""
    proc = subprocess.run(
        [minimap2, "-x", "asm5", str(mito_fasta), str(pltd_fasta)], capture_output=True, text=True, timeout=600
    )
    if proc.returncode != 0:
        return None
    query: dict[str, list] = {}
    target: dict[str, list] = {}
    for line in proc.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 12:
            continue
        query.setdefault(f[0], []).append((int(f[2]), int(f[3])))
        target.setdefault(f[5], []).append((int(f[7]), int(f[8])))
    mito_len = contig_metrics(mito_fasta)["assembly_length"]
    pltd_len = contig_metrics(pltd_fasta)["assembly_length"]
    mito_aligned = merged_aligned_length(target) if target else 0
    pltd_aligned = merged_aligned_length(query) if query else 0
    return (mito_aligned / mito_len if mito_len else 0.0), (pltd_aligned / pltd_len if pltd_len else 0.0)


# ---------------------------------------------------------------- status


def build_status(row: dict) -> tuple[str, str]:
    """05_qc_summary.py's build_status, except that a missing GFA leaves the
    graph checks unassessed rather than failing."""
    fail, flag = [], []

    def known(x):
        return x is not None and not (isinstance(x, float) and math.isnan(x))

    if row.get("parse_status") == "empty":
        fail.append("empty_assembly")
    else:
        den = row.get("dead_end_nodes")
        if known(den):
            if den > DEAD_END_FAIL_THRESHOLD:
                fail.append("dead_end_nodes>2")
            elif den > DEAD_END_FLAG_THRESHOLD:
                flag.append("dead_end_nodes>0")
        n_sub = row.get("n_fragment_subgraphs")
        if not known(n_sub):
            n_sub = row.get("n_subgraphs")
        if known(n_sub):
            if n_sub >= FRAGMENTATION_FAIL_THRESHOLD:
                fail.append(f"fragmented(n_fragment_subgraphs>={FRAGMENTATION_FAIL_THRESHOLD})")
            elif n_sub >= FRAGMENTATION_FLAG_THRESHOLD:
                flag.append(f"fragmented(n_fragment_subgraphs>={FRAGMENTATION_FLAG_THRESHOLD})")
        circ = row.get("any_circular")
        if known(circ) and not circ:
            flag.append("non_circular")
        z = row.get("size_mad_z")
        if known(z) and abs(z) > SIZE_Z_FLAG_THRESHOLD:
            flag.append("size_outlier(|z|>4)")
        cgp = row.get("core_gene_pct")
        if known(cgp):
            if cgp < CORE_GENE_FAIL_THRESHOLD:
                fail.append("core_gene_pct<0.50")
            elif cgp < CORE_GENE_FLAG_THRESHOLD:
                flag.append("core_gene_pct<0.80")
        if row.get("unjoined"):
            flag.append("no_resolved_ctg_fasta(unjoined)")
        elif row.get("has_ctg_fasta") is False:
            flag.append("no_resolved_ctg_fasta")
        cof = row.get("cross_organelle_aligned_frac")
        complete_plastid = (
            row.get("organelle") == "pltd"
            and bool(row.get("any_circular"))
            and known(cgp)
            and cgp >= QC_CORE_GENE_PRESENCE_THRESHOLD
        )
        if known(cof) and not complete_plastid:
            if cof > CROSS_ORG_FAIL_THRESHOLD:
                fail.append("high_cross_organelle_alignment(>35%)")
            elif cof > CROSS_ORG_FLAG_THRESHOLD:
                flag.append("elevated_cross_organelle_alignment(>15%)")

    if fail:
        return "fail", ";".join(fail + flag)
    if flag:
        return "flag", ";".join(flag)
    return "pass", ""


def assess(
    species: str,
    organelle: str,
    *,
    ctg_fasta: Path | None,
    gfa: Path | None,
    genes_present: set[str] | None,
    genus: str | None,
    cross_frac: float | None,
    core: dict[str, set[str]],
    size_ref: dict,
    gfatk: str | None,
) -> dict:
    """One organelle assembly's QC row (COLUMNS, plus the internal fields
    build_status reads)."""
    row: dict = {"species": species, "organelle": organelle, "cross_organelle_aligned_frac": cross_frac}
    not_assessed = []
    if gfa is None:
        not_assessed.append("graph(no GFA)")
    elif gfatk is None:
        not_assessed.append("graph(gfatk not found)")
    else:
        row.update(gfa_stats(gfa, gfatk))
        if row.get("parse_status") == "error":
            not_assessed.append("graph(gfatk stats failed)")
    if ctg_fasta is not None:
        row.update(contig_metrics(ctg_fasta))
        row["has_ctg_fasta"] = not row["unjoined"]
        length = row["assembly_length"] if row["has_ctg_fasta"] else row.get("total_length")
    else:
        row["has_ctg_fasta"] = False
        length = row.get("total_length")
    if genes_present is None:
        not_assessed.append("core_genes(no gene calls)")
    else:
        row["core_gene_pct"], missing = core_gene_completeness(genes_present, core.get(organelle, set()))
        row["core_genes_missing"] = ",".join(missing)
    row["size_mad_z"], row["size_reference"] = size_z(length, organelle, genus, size_ref)
    if cross_frac is None:
        not_assessed.append("cross_organelle(needs both organelles)")
    row["not_assessed"] = ";".join(not_assessed)
    row["status"], row["status_reasons"] = build_status(row)
    return row


# ---------------------------------------------------------------- stage


def present_genes(calls: pd.DataFrame, species: str, organelle: str) -> set[str]:
    sub = calls[(calls["species"] == species) & (calls["organelle"] == organelle)]
    return set(sub["gene"])


def fmt(v) -> str:
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return ""
    if isinstance(v, bool):
        return str(v).lower()
    if isinstance(v, float):
        return f"{v:.4f}"
    return str(v)


def sample_qc(
    root: Path,
    organelles: list[str] | tuple[str, ...] = ORGANELLES,
    *,
    gene_calls: Path | None = None,
    genus: str | None = None,
) -> Path:
    """QC every assembly under data root `root` (module docstring). genus:
    for the size check (default: the first word of each species name).
    Returns the output table's path."""
    root = Path(root)
    analysis_dir = root / "analysis"
    calls_path = Path(gene_calls) if gene_calls else analysis_dir / "denovo_annotation" / "results" / "gene_calls.tsv"
    calls = pd.read_csv(calls_path, sep="\t", usecols=["species", "organelle", "gene"]) if calls_path.exists() else None
    if calls is None:
        log.warning(f"sample_qc: no {calls_path} - core genes not assessed")
    core, size_ref = load_core_genes(), load_size_reference()
    shift_tables = [
        pd.read_csv(p, sep="\t")
        for p in (
            analysis_dir / "editing" / "results" / "frameshifts.tsv",
            analysis_dir / "trans_splicing" / "results" / "reconstructed_frameshifts.tsv",
        )
        if p.exists()
    ]
    shifts = pd.concat(shift_tables, ignore_index=True) if shift_tables else None
    gfatk, minimap2 = shutil.which("gfatk"), shutil.which("minimap2")

    found = {
        o: {r.species: r for r in sd.discover_all(sd.repo_data_root(root, o), o) if r.ctg_fasta or r.gfa}
        for o in organelles
        if sd.repo_data_root(root, o).is_dir()
    }
    cross: dict[tuple[str, str], float] = {}
    if {"mito", "pltd"} <= set(found):
        for species in sorted(set(found["mito"]) & set(found["pltd"])):
            m, p = found["mito"][species], found["pltd"][species]
            if not (m.ctg_fasta and p.ctg_fasta) or minimap2 is None:
                continue
            if contig_metrics(Path(m.ctg_fasta))["unjoined"] or contig_metrics(Path(p.ctg_fasta))["unjoined"]:
                continue  # the dataset gate only aligns resolved assemblies
            fracs = cross_organelle(Path(m.ctg_fasta), Path(p.ctg_fasta), minimap2)
            if fracs is not None:
                cross[(species, "mito")], cross[(species, "pltd")] = fracs

    rows = []
    for organelle, by_species in found.items():
        for species, r in sorted(by_species.items()):
            row = assess(
                species,
                organelle,
                ctg_fasta=Path(r.ctg_fasta) if r.ctg_fasta else None,
                gfa=Path(r.gfa) if r.gfa else None,
                genes_present=present_genes(calls, species, organelle) if calls is not None else None,
                genus=genus or species.split("_")[0],
                cross_frac=cross.get((species, organelle)),
                core=core,
                size_ref=size_ref,
                gfatk=gfatk,
            )
            if shifts is not None:
                fs = shifts[(shifts.species == species) & (shifts.organelle == organelle)]
                seen = (
                    pd.to_numeric(fs["dtol_assemblies"], errors="coerce")
                    if "dtol_assemblies" in fs.columns
                    else pd.Series(float("nan"), index=fs.index)
                )
                row["n_frameshifts"] = len(fs)
                row["n_new_frameshifts"] = int((seen.fillna(0) <= NEW_FRAMESHIFT_MAX).sum())
                row["frameshifted_genes"] = ",".join(
                    f"{g}:{'' if pd.isna(n) else int(n)}" for g, n in sorted(zip(fs.gene, seen), key=lambda x: x[0])
                )
            rows.append(row)
            log.info(
                f"sample_qc: {species} {organelle}: {row['status']}"
                + (f" ({row['status_reasons']})" if row["status_reasons"] else "")
                + (f"; not assessed: {row['not_assessed']}" if row["not_assessed"] else "")
            )
    if not rows:
        raise StageError(f"sample_qc: no assemblies under {root / 'data'}")

    out = analysis_dir / "qc_basic_stats" / "results" / "sample_qc.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        for row in rows:
            fh.write("\t".join(fmt(row.get(c)) for c in COLUMNS) + "\n")
    log.info(f"sample_qc: {len(rows)} rows -> {out}")
    return out
