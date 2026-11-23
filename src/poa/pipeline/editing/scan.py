"""Run `orfedit` over every (species, gene) anchor that has a reference
profile, correcting gene boundaries for C->U RNA editing.

Builds one manifest TSV from `annotation/results/gene_calls.tsv` (already
QC-gated - see `annotation/src/01_gene_table.py`) joined against the
profiles `00_build_reference_profiles.py` built, then makes a single call
to `orfedit scan-batch`, which rayon-parallelizes across every row itself -
deliberately not the `--jobs`/`--cpu` multi-process pattern `orf_scan` uses
for `hmmscan`. That pattern exists because `hmmscan` isn't internally
parallel; `orfedit` is (via rayon), so one process using `--threads` is the
right shape here, and it sidesteps the LSF-array-of-many-small-jobs
scaling pain `orf_scan/README.md` already documents for its own heavy step.

Mito genes are also scanned against each other lineage's bundled profiles
(00b_build_lineage_profiles.py - moss), and each species keeps the run of
the lineage its genes fit best (common/lineage.py; --lineage to force
one). A species' genes with no profile in its lineage (moss nad5/nad9,
multi-exon there) are left uncorrected.

Every substantial mito copy of a gene is corrected, not just the
best-scoring hit: the best hit, plus any other hit covering at least
MIN_COPY_SPAN of the HMM that overlaps neither a hit already taken nor a
stronger hit of another gene (Arabidopsis has two atp6 genes with
different N-termini; 41-129 species per gene in the DToL dataset have two
or more such copies). An extra copy keeps its correction only if orfedit
aligns at least MIN_COPY_COVERAGE of the profile. Plastids keep the best
hit only (see COPY_ORGANELLES). Each copy goes to orfedit as gene
`<gene>|<contig>:<start>-<end>`, which orfedit passes through untouched,
and both output tables get the copy back as `contig`, `hit_start` and
`hit_end`. The lineage call uses each gene's best-scoring copy.

Angiosperm mito genes with no dataset-built (or bundled dataset) profile -
everything but atp4 ccmB ccmC mttB nad5 nad9 rps14, which are the genes
single-copy across the dataset - use the bundled RefSeq profiles
(00b_build_lineage_profiles.py --lineage angiosperm) instead. The lineage
call only compares genes both lineages have profiles for, so this doesn't
change it. The call is recorded in
lineage.tsv, which trans_splicing reads to pick its exon templates.

Frameshifts: orfedit may change reading frame inside a gene, skipping 1
or 2 bases at a cost (DEFAULT_FRAMESHIFT_PENALTY), the way a
frameshifting indel would - an assembly error (HiFi homopolymer length) or
a real pseudogene. Each one is a row of frameshifts.tsv: the skipped
bases' genomic [start, end), n_bases (1: the frame wants one base fewer,
e.g. a 1-bp insertion; 2: one base more, e.g. a 1-bp deletion) and the
profile column it falls before; edited_gene_calls.tsv counts them per copy
(n_frameshifts). dtol_assemblies says how many DToL assemblies have one in
the same gene within 5 profile columns (add_recurrence), which separates
a gene's shared frameshifts (rpl16, infA, ndhK) from this assembly's own.

Output: analysis/editing/results/edited_gene_calls.tsv
        analysis/editing/results/predicted_edits.tsv
        analysis/editing/results/frameshifts.tsv
        analysis/editing/results/lineage.tsv

Library use - everything is passed in, nothing is read from the
environment, and failures raise StageError:

    from poa.pipeline.editing.scan import scan
    result = scan(root, ["mito"], gene_calls=..., threads=8)

src/01_scan_editing.py is the command-line wrapper (main() here), which
takes the data root from PLANT_ORGANELLE_DATA_ROOT or the current directory.
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
from poa.pipeline.editing.lineage_profiles import refseq_profiles_dir
from poa.pipeline.trans_splicing.reference_genes import LINEAGE_EXTRA_SPLICED
from poa.pipeline.common import species_discovery as sd
from poa.pipeline.common.errors import StageError
from poa.pipeline.common.tools import cargo_tool
from poa.pipeline.common.log import setup_logging

log = logging.getLogger(__name__)

ORGANELLES = ("mito", "pltd")
DEFAULT_THREADS = 4
DEFAULT_FLANK = 300
DEFAULT_EDIT_PENALTY = 2.5
# orfedit's cost of a frame change inside a gene (see FRAMESHIFTS below)
DEFAULT_FRAMESHIFT_PENALTY = 20.0
# a hit other than a gene's best is corrected as a copy of its own when it
# covers at least this much of the HMM (a short fragment - Arabidopsis's
# third atp6 hit spans 21% - stays raw)
MIN_COPY_SPAN = 0.5
# ...and keeps its correction only when orfedit aligns at least this much of
# the profile (a gene's best copy is always kept, as before). The 5th
# percentile for single-copy mito genes on the DToL dataset is 0.88-0.99.
MIN_COPY_COVERAGE = 0.8
# Extra copies are mito only. In plastids the candidates are mostly the
# truncated rps19 at the inverted-repeat boundary (404 of 2328 across the
# dataset) and duplicates from redundant contigs; the genes the inverted
# repeat really duplicates (rps7, ycf2, ndhB, rpl2, rpl23) have no profile.
COPY_ORGANELLES = ("mito",)
COPY_SEP = "|"


@dataclass(frozen=True)
class ScanResult:
    edited_gene_calls: Path
    predicted_edits: Path
    frameshifts: Path
    lineage_table: Path
    lineages: dict[str, str]  # species -> lineage, for the species scanned


def hits_to_correct(gene_calls: pd.DataFrame, organelle: str, genes: list[str]) -> pd.DataFrame:
    """Per species and gene: the best-scoring hit, plus every other hit
    covering >= MIN_COPY_SPAN of the HMM that overlaps neither a hit already
    taken on the same contig nor a stronger hit of a different gene (a
    paralogue's locus: plastid psaA's HMM also hits psaB at 89% of its
    length, scoring 434 against psaB's 2724). Only the best hit outside
    COPY_ORGANELLES, or without hmm_from/hmm_to/model_len columns."""
    org = gene_calls[gene_calls["organelle"] == organelle]
    sub = org[org["gene"].isin(genes)]
    if sub.empty:
        return sub
    has_span = {"hmm_from", "hmm_to", "model_len"} <= set(sub.columns) and organelle in COPY_ORGANELLES
    by_contig = {k: g for k, g in org.groupby(["species", "contig_id"])}
    keep = []
    for (species, _), g in sub.sort_values("score", ascending=False, kind="stable").groupby(
        ["species", "gene"], sort=False
    ):
        taken: list = []
        for i, h in g.iterrows():
            if taken:
                if not has_span or not (h.hmm_to - h.hmm_from + 1) >= MIN_COPY_SPAN * h.model_len:
                    continue
                if any(t.contig_id == h.contig_id and t.start <= h.end and h.start <= t.end for t in taken):
                    continue
                here = by_contig[(species, h.contig_id)]
                rivals = here[(here.gene != h.gene) & (here.start <= h.end) & (h.start <= here.end)]
                if (rivals.score > h.score).any():
                    continue
            taken.append(h)
            keep.append(i)
    return sub.loc[keep]


def copy_label(gene: str, contig: str, start: int, end: int) -> str:
    return f"{gene}{COPY_SEP}{contig}:{start}-{end}"


def split_copy_labels(df: pd.DataFrame) -> pd.DataFrame:
    """orfedit output with copy-labelled genes -> gene, contig, hit_start, hit_end."""
    if df.empty:
        return df.assign(contig=pd.Series(dtype=str), hit_start=pd.Series(dtype=str), hit_end=pd.Series(dtype=str))
    parts = df["gene"].str.extract(r"^(?P<gene>[^|]+)\|(?P<contig>.+):(?P<hit_start>\d+)-(?P<hit_end>\d+)$")
    if parts["gene"].isna().any():
        raise StageError(
            f"orfedit returned gene labels without a copy: {df.loc[parts.gene.isna(), 'gene'].unique()[:5]}"
        )
    df = df.drop(columns="gene")
    at = list(df.columns).index("organelle") + 1
    for j, col in enumerate(["gene", "contig", "hit_start", "hit_end"]):
        df.insert(at + j, col, parts[col])
    return df


def profile_paths(profiles_dirs: list[Path], organelle: str) -> dict[str, Path]:
    """{gene: profile} for one organelle; each gene from the first of
    profiles_dirs (each holding <organelle>/<gene>.pssm) that has it."""
    found: dict[str, Path] = {}
    for d in profiles_dirs:
        for p in sorted((Path(d) / organelle).glob("*.pssm")):
            found.setdefault(p.stem, p)
    return found


def build_manifest(
    root: Path,
    gene_calls: pd.DataFrame,
    organelles: list[str],
    species_filter: set | None,
    profiles_dir: Path | list[Path],
) -> list[dict]:
    """profiles_dir: one directory of <organelle>/<gene>.pssm, or several in
    priority order."""
    dirs = [profiles_dir] if isinstance(profiles_dir, (str, Path)) else list(profiles_dir)
    rows = []
    for organelle in organelles:
        profiles = profile_paths(dirs, organelle)
        genes = sorted(profiles)
        if not genes:
            log.warning(
                f"{organelle}: no profiles in {', '.join(str(Path(d) / organelle) for d in dirs)} - run "
                f"editing/src/00_build_reference_profiles.py first"
            )
            continue
        hits = hits_to_correct(gene_calls, organelle, genes)
        if hits.empty:
            continue
        data_root = sd.repo_data_root(root, organelle)
        # Resolve each species' ctg.fasta path once, reused across all its genes.
        fasta_cache: dict[str, str | None] = {}
        for row in hits.itertuples():
            if species_filter is not None and row.species not in species_filter:
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
                    "gene": copy_label(row.gene, row.contig_id, int(row.start), int(row.end)),
                    "fasta": fasta_path,
                    "profile": str(profiles[row.gene]),
                    "contig": row.contig_id,
                    "start": int(row.start),
                    "end": int(row.end),
                    "strand": row.strand,
                }
            )
    return rows


def lineage_spliced(lineage: str) -> set[str]:
    """Mito genes single-exon in angiosperms that `lineage` splices and has
    exon templates for (reference_genes.LINEAGE_EXTRA_SPLICED)."""
    templates = lin.trans_splicing_reference(PIPELINE_DIR, lineage) / "profiles"
    return {g for g in LINEAGE_EXTRA_SPLICED.get(lineage, ()) if (templates / g / "whole_gene.pssm").exists()}


def without_genes(rows: list[dict], genes: set[str]) -> list[dict]:
    """Manifest rows minus the mito ones for `genes`."""
    if not genes:
        return rows
    return [r for r in rows if not (r["organelle"] == "mito" and r["gene"].split(COPY_SEP)[0] in genes)]


def write_manifest(rows: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    cols = ["species", "organelle", "gene", "fasta", "profile", "contig", "start", "end", "strand"]
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")


def supports_option(tool: str, option: str) -> bool:
    """Whether `tool scan-batch` takes `option`."""
    out = subprocess.run([tool, "scan-batch", "--help"], capture_output=True, text=True).stdout
    return option in out


def supports_frameshifts(orfedit: str) -> bool:
    """orfedit >= 0.2.0 (frame changes, --frameshifts-out)."""
    out = subprocess.run([orfedit, "--version"], capture_output=True, text=True).stdout
    m = re.search(r"(\d+)\.(\d+)\.(\d+)", out)
    return m is not None and tuple(map(int, m.groups())) >= (0, 2, 0)


FRAMESHIFT_COLUMNS = ["species", "organelle", "gene", "start", "end", "n_bases", "profile_column", "window_start", "window_end"]
# frameshifts orfedit found across the DToL dataset (frameshift_reference.py)
FRAMESHIFT_REFERENCE = PIPELINE_DIR / "editing" / "reference" / "frameshifts_dtol.tsv"
RECURRENCE_WINDOW = 5  # profile columns


def add_recurrence(
    shifts: pd.DataFrame,
    lineages: dict[str, str],
    source: str = "editing",
    reference_of=lin.reference,
    own_plastid: set[tuple[str, str]] | frozenset = frozenset(),
) -> pd.DataFrame:
    """Add dtol_assemblies: how many DToL assemblies have a frameshift in
    the same gene within RECURRENCE_WINDOW profile columns - many means a
    shared feature of the gene (a pseudogene, an alternative frame), none
    a frameshift of this assembly's own. Empty for a mito gene aligned to
    another lineage's profile, whose columns the reference doesn't use.
    source: "editing" (orfedit's profiles) or "reconstruction"
    (transsplice's whole-gene profiles), which number columns differently.
    reference_of: lineage -> the lineage whose profiles its genes were
    aligned to (for reconstruction, its exon templates'). own_plastid:
    (species, gene) pairs of plastid genes aligned to a lineage's own
    plastid references, which get no count either."""
    shifts = shifts.copy()
    if not FRAMESHIFT_REFERENCE.exists() or shifts.empty:
        shifts["dtol_assemblies"] = pd.Series(dtype=object)
        return shifts
    ref = pd.read_csv(FRAMESHIFT_REFERENCE, sep="\t")
    ref = ref[ref["source"] == source]
    by_gene = {k: g for k, g in ref.groupby(["organelle", "gene"])}
    out = []
    for r in shifts.itertuples():
        if r.organelle == "mito" and reference_of(lineages.get(r.species, lin.DEFAULT_LINEAGE)) != lin.DEFAULT_LINEAGE:
            out.append("")
            continue
        if r.organelle == "pltd" and (r.species, r.gene) in own_plastid:
            out.append("")
            continue
        g = by_gene.get((r.organelle, r.gene))
        col = int(r.profile_column)
        near = 0 if g is None else int(g.loc[(g.profile_column - col).abs() <= RECURRENCE_WINDOW, "n_assemblies"].sum())
        out.append(near)
    shifts["dtol_assemblies"] = out
    return shifts


def run_orfedit(
    orfedit: str,
    rows: list[dict],
    tag: str,
    work_dir: Path,
    threads: int,
    flank: int,
    edit_penalty: float,
    frameshift_penalty: float,
    stop_edit_penalty: float | None = None,
    reverse_edit_penalty: float | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """One orfedit scan-batch over `rows`; (gene calls, edits, frameshifts)
    as all-string frames. stop_edit_penalty: read in-frame stop codons
    through as U-to-C edits at this cost (common/lineage.py, fern);
    reverse_edit_penalty: any other U-to-C edit at this cost."""
    manifest_path = work_dir / f"manifest.{tag}.tsv"
    write_manifest(rows, manifest_path)
    log.info(f"scan_editing: {tag} manifest has {len(rows)} rows -> {manifest_path}")
    calls_out, edits_out = work_dir / f"edited_gene_calls.{tag}.tsv", work_dir / f"predicted_edits.{tag}.tsv"
    shifts_out = work_dir / f"frameshifts.{tag}.tsv"
    cmd = [
        orfedit,
        "scan-batch",
        "--manifest",
        str(manifest_path),
        "--gene-calls-out",
        str(calls_out),
        "--edits-out",
        str(edits_out),
        "--threads",
        str(threads),
        "--flank",
        str(flank),
        "--edit-penalty",
        str(edit_penalty),
    ]
    if math.isfinite(frameshift_penalty):
        cmd += ["--frameshifts-out", str(shifts_out), "--frameshift-penalty", str(frameshift_penalty)]
    else:
        pd.DataFrame(columns=FRAMESHIFT_COLUMNS).to_csv(shifts_out, sep="\t", index=False)
    if stop_edit_penalty is not None:
        cmd += ["--stop-edit-penalty", str(stop_edit_penalty)]
    if reverse_edit_penalty is not None:
        cmd += ["--reverse-edit-penalty", str(reverse_edit_penalty)]
    proc = subprocess.run(cmd)
    if proc.returncode != 0:
        raise StageError(f"orfedit scan-batch failed (exit {proc.returncode})")
    calls, edits, shifts = (split_copy_labels(pd.read_csv(p, sep="\t", dtype=str)) for p in (calls_out, edits_out, shifts_out))
    return drop_weak_copies(calls, edits, shifts)


def drop_weak_copies(calls: pd.DataFrame, *per_gene: pd.DataFrame) -> tuple[pd.DataFrame, ...]:
    """Keep each species x gene's best-scoring corrected copy, and any other
    copy orfedit aligned to at least MIN_COPY_COVERAGE of the profile;
    `per_gene` tables (edits, frameshifts) keep the kept copies' rows."""
    if calls.empty:
        return (calls, *per_gene)
    score = calls["score"].astype(float)
    frac = calls["n_columns_aligned"].astype(float) / calls["n_columns_total"].astype(float)
    best = score.groupby([calls.species, calls.organelle, calls.gene]).transform("max") == score
    kept = calls[best | (frac >= MIN_COPY_COVERAGE)]
    key = ["species", "organelle", "gene", "contig", "hit_start", "hit_end"]
    kept_keys = kept[key].drop_duplicates()
    return (kept, *(t.merge(kept_keys, on=key, how="inner") if not t.empty else t for t in per_gene))


def assign_lineages(runs: dict[str, pd.DataFrame], forced: str | None) -> pd.DataFrame:
    """One row per species with mito calls: the lineage it takes and why."""
    scores: dict[str, dict[str, dict[str, float]]] = {}
    for name, calls in runs.items():
        mito = calls[calls["organelle"] == "mito"]
        for row in mito.itertuples():
            # a gene's best copy speaks for it (one copy each when calibrated)
            genes = scores.setdefault(row.species, {}).setdefault(name, {})
            genes[row.gene] = max(genes.get(row.gene, float("-inf")), float(row.score))
    out = []
    for species in sorted(scores):
        chosen, n, totals = lin.choose(scores[species])
        out.append(
            {
                "species": species,
                "lineage": forced or chosen,
                "decided_by": "forced" if forced else "auto",
                "n_genes": n,
                **{f"{x}_score": f"{totals[x]:.3f}" if x in totals else "" for x in lin.LINEAGES},
            }
        )
    return pd.DataFrame(out, columns=lin.COLUMNS)


def keep_lineage(df: pd.DataFrame, lineages: dict[str, str], name: str) -> pd.DataFrame:
    """Rows of lineage `name`'s run that belong in the output: its own mito
    species, plus (for the default run) everything that isn't mito."""
    if df.empty:
        return df
    takes = df["species"].map(lineages).fillna(lin.DEFAULT_LINEAGE).map(lin.reference) == name
    if name == lin.DEFAULT_LINEAGE:
        return df[(df["organelle"] != "mito") | takes]
    return df[(df["organelle"] == "mito") & takes]


def scan(
    root: Path,
    organelles: list[str] | tuple[str, ...] = ORGANELLES,
    *,
    gene_calls: Path | None = None,
    species: set[str] | None = None,
    threads: int = DEFAULT_THREADS,
    flank: int = DEFAULT_FLANK,
    edit_penalty: float = DEFAULT_EDIT_PENALTY,
    frameshift_penalty: float = DEFAULT_FRAMESHIFT_PENALTY,
    lineage: str = "auto",
) -> ScanResult:
    """Correct the gene calls under data root `root` for RNA editing.

    gene_calls: a gene_calls.tsv-shaped table to anchor on (default: the
    data root's denovo_annotation/results/gene_calls.tsv). species: only
    these species, merged into the existing result tables rather than
    replacing them. lineage: "auto" picks each species' mito reference
    lineage, any of common/lineage.LINEAGES forces it."""
    root = Path(root)
    analysis_dir = root / "analysis"
    organelles = list(organelles)
    if unknown := set(organelles) - set(ORGANELLES):
        raise ValueError(f"unknown organelle(s) {sorted(unknown)}; expected {ORGANELLES}")
    if lineage not in lin.CHOICES:
        raise ValueError(f"unknown lineage {lineage!r}; expected one of {lin.CHOICES}")
    orfedit = cargo_tool("orfedit")
    if math.isfinite(frameshift_penalty) and not supports_frameshifts(orfedit):
        log.warning(f"scan_editing: {orfedit} is older than 0.2.0 - no frameshift calls")
        frameshift_penalty = math.inf

    gene_calls_path = (
        Path(gene_calls) if gene_calls else analysis_dir / "denovo_annotation" / "results" / "gene_calls.tsv"
    )
    if not gene_calls_path.exists():
        raise StageError(
            f"missing {gene_calls_path} - run denovo_annotation/src/04_build_gene_calls.py first (or pass --gene-calls)"
        )
    calls = pd.read_csv(gene_calls_path, sep="\t")

    profiles_dir = analysis_dir / "editing" / "results" / "profiles"
    if not any(profiles_dir.glob("*/*.pssm")):
        # a dataset too small for 00_build_reference_profiles.py (>=5 species
        # per gene) uses the bundled ones, built on the 1250-species dataset
        profiles_dir = PIPELINE_DIR / "editing" / "reference" / "profiles"
        log.info(f"scan_editing: no dataset-built profiles, using bundled {profiles_dir}")

    refseq_dir = refseq_profiles_dir(PIPELINE_DIR)
    rows = build_manifest(root, calls, organelles, species, [profiles_dir, refseq_dir])
    # genes with introns in the forced lineage (fern atp1, cox1, rps14):
    # trans_splicing reconstructs them from the lineage's exon templates
    spliced_here = lineage_spliced(lineage)
    rows = without_genes(rows, spliced_here)
    if not rows:
        raise StageError("empty manifest - nothing to scan")

    work_dir = analysis_dir / "editing" / "work"
    work_dir.mkdir(parents=True, exist_ok=True)
    results_dir = analysis_dir / "editing" / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    # orfedit writes its output paths itself (truncating), so a
    # species-scoped run must not point it straight at results_dir - it'd
    # truncate the shared table down to just the scoped species. It writes
    # to work_dir instead, then merge_species_scoped_tsv decides how that
    # relates to what's already in results_dir (see its docstring).
    stop = lin.stop_edit_penalty(lineage)
    if stop is not None and not supports_option(orfedit, "--stop-edit-penalty"):
        raise StageError(f"--lineage {lineage} needs an orfedit with --stop-edit-penalty; {orfedit} has none")
    opts = dict(threads=threads, flank=flank, edit_penalty=edit_penalty, frameshift_penalty=frameshift_penalty)
    runs = {
        lin.DEFAULT_LINEAGE: run_orfedit(orfedit, rows, lin.DEFAULT_LINEAGE, work_dir, stop_edit_penalty=stop, **opts)
    }
    # a forced lineage that borrows the default references (fern) has no
    # other lineage's profiles to try
    if "mito" in organelles and lin.reference(lineage) != lin.DEFAULT_LINEAGE:
        for name in lin.LINEAGES[1:]:
            lrows = build_manifest(root, calls, ["mito"], species, lin.editing_profiles(PIPELINE_DIR, name))
            if lrows:
                runs[name] = run_orfedit(orfedit, lrows, name, work_dir, **opts)

    # A forced lineage that borrows the default references but has some
    # editing profiles of its own (fern): its genes with such a profile are
    # corrected against it, the rest against the default ones as above.
    own_run = None
    own_dir = lin.editing_profiles(PIPELINE_DIR, lineage)
    # mito, and plastid where the lineage has plastid profiles too (fern)
    own_orgs = [o for o in organelles if lineage in lin.BORROWS and any((own_dir / o).glob("*.pssm"))]
    own_parts = []
    for o in own_orgs:
        lrows = without_genes(build_manifest(root, calls, [o], species, own_dir), spliced_here)
        if not lrows:
            continue
        # other U-to-C edits too, on the lineage's own mito profiles only:
        # in plastid genes most such calls are not in the RNA (see
        # trans_splicing/reconstruct.py)
        rev = lin.reverse_edit_penalty(lineage) if o == "mito" else None
        if rev is not None and not supports_option(orfedit, "--reverse-edit-penalty"):
            raise StageError(f"--lineage {lineage} needs an orfedit with --reverse-edit-penalty; {orfedit} has none")
        tag = lineage if o == "mito" else f"{lineage}.{o}"
        part = run_orfedit(orfedit, lrows, tag, work_dir, stop_edit_penalty=stop, reverse_edit_penalty=rev, **opts)
        own_parts.append(part)
        log.info(f"scan_editing: {part[0]['gene'].nunique()} {o} genes use {lineage} profiles")
    if own_parts:
        own_run = tuple(pd.concat([part[i] for part in own_parts], ignore_index=True) for i in range(3))

    forced = None if lineage == "auto" else lineage
    lineage_table = assign_lineages({k: run[0] for k, run in runs.items()}, forced)
    lineages = dict(zip(lineage_table["species"], lineage_table["lineage"]))
    for name, n in lineage_table["lineage"].value_counts().items():
        log.info(f"scan_editing: {n} species use {name} mito profiles")
    work_gene_calls_out = work_dir / "edited_gene_calls.new.tsv"
    work_edits_out = work_dir / "predicted_edits.new.tsv"
    work_shifts_out = work_dir / "frameshifts.new.tsv"
    for i, out in enumerate((work_gene_calls_out, work_edits_out, work_shifts_out)):
        parts = [keep_lineage(run[i], lineages, name) for name, run in runs.items()]
        if own_run is not None:
            covered = set(zip(own_run[0]["organelle"], own_run[0]["gene"]))
            parts = [p[[pair not in covered for pair in zip(p["organelle"], p["gene"])]] for p in parts]
        table = pd.concat(parts, ignore_index=True)
        if out == work_shifts_out:
            table = add_recurrence(table, lineages)
        if own_run is not None:
            # (no DToL recurrence for these: the reference counts frameshifts
            # by the default profiles' columns)
            own = own_run[i].assign(dtol_assemblies="") if out == work_shifts_out else own_run[i]
            table = pd.concat([table, own], ignore_index=True)
        table.to_csv(out, sep="\t", index=False)
    work_lineage_out = work_dir / "lineage.new.tsv"
    lineage_table.to_csv(work_lineage_out, sep="\t", index=False)

    result = ScanResult(
        edited_gene_calls=results_dir / "edited_gene_calls.tsv",
        predicted_edits=results_dir / "predicted_edits.tsv",
        frameshifts=results_dir / "frameshifts.tsv",
        lineage_table=lin.table_path(analysis_dir),
        lineages=lineages,
    )
    io_utils.merge_species_scoped_tsv(result.edited_gene_calls, work_gene_calls_out, species)
    io_utils.merge_species_scoped_tsv(result.predicted_edits, work_edits_out, species)
    io_utils.merge_species_scoped_tsv(result.frameshifts, work_shifts_out, species)
    io_utils.merge_species_scoped_tsv(result.lineage_table, work_lineage_out, species)
    log.info(f"scan_editing: merged into {results_dir} ({'species-scoped' if species else 'full rewrite'})")
    return result


def main(argv: list[str] | None = None) -> int:
    setup_logging()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--organelle", default="both", choices=["mito", "pltd", "both"])
    ap.add_argument("--species-list", default=None)
    ap.add_argument("--threads", type=int, default=DEFAULT_THREADS)
    ap.add_argument("--flank", type=int, default=DEFAULT_FLANK)
    ap.add_argument("--edit-penalty", type=float, default=DEFAULT_EDIT_PENALTY)
    ap.add_argument(
        "--frameshift-penalty",
        type=float,
        default=DEFAULT_FRAMESHIFT_PENALTY,
        help="cost of a frame change inside a gene (inf: single-frame alignment, no frameshift calls)",
    )
    ap.add_argument(
        "--lineage",
        default="auto",
        choices=lin.CHOICES,
        help="mito reference lineage: auto picks per species (default), or force one for all",
    )
    ap.add_argument(
        "--gene-calls",
        default=None,
        help="path to a gene_calls.tsv-shaped table to anchor "
        "on (default: denovo_annotation/results/gene_calls.tsv under the data root)",
    )
    ap.add_argument(
        "--force",
        action="store_true",
        help="accepted for CLI consistency; orfedit "
        "always does a full recompute over the manifest (no incremental skip yet)",
    )
    args = ap.parse_args(argv)
    root = Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT") or os.getcwd())
    try:
        scan(
            root,
            list(ORGANELLES) if args.organelle == "both" else [args.organelle],
            gene_calls=Path(args.gene_calls) if args.gene_calls else None,
            species=sd.load_species_list(Path(args.species_list)) if args.species_list else None,
            threads=args.threads,
            flank=args.flank,
            edit_penalty=args.edit_penalty,
            frameshift_penalty=args.frameshift_penalty,
            lineage=args.lineage,
        )
    except StageError as e:
        log.error(e)
        return 1
    return 0
