#!/usr/bin/env python3
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
multi-exon there) are left uncorrected. The call is recorded in
lineage.tsv, which trans_splicing reads to pick its exon templates.

Output: analysis/editing/results/edited_gene_calls.tsv
        analysis/editing/results/predicted_edits.tsv
        analysis/editing/results/lineage.tsv
"""
from __future__ import annotations

import os
import argparse
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

# CODE_DIR is the installed poa/pipeline/ (shared modules, bundled reference
# data); data/ and every module's results/ and work/ live under the data
# root - PLANT_ORGANELLE_DATA_ROOT, or the current directory if unset.
CODE_DIR = Path(__file__).resolve().parents[2]
ROOT_DIR = Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT") or os.getcwd())
ANALYSIS_DIR = ROOT_DIR / "analysis"
sys.path.insert(0, str(CODE_DIR / "common"))

import species_discovery as sd  # noqa: E402
import io_utils  # noqa: E402
import lineage as lin  # noqa: E402

ORFEDIT = shutil.which("orfedit") or str(Path.home() / ".cargo" / "bin" / "orfedit")


def best_hits_per_species(gene_calls: pd.DataFrame, organelle: str, genes: list[str]) -> pd.DataFrame:
    sub = gene_calls[(gene_calls["organelle"] == organelle) & (gene_calls["gene"].isin(genes))]
    if sub.empty:
        return sub
    idx = sub.groupby(["species", "gene"])["score"].idxmax()
    return sub.loc[idx]


def build_manifest(gene_calls: pd.DataFrame, organelles: list[str], species_filter: set | None,
                    profiles_dir: Path) -> list[dict]:
    rows = []
    for organelle in organelles:
        organelle_profiles_dir = profiles_dir / organelle
        genes = sorted(p.stem for p in organelle_profiles_dir.glob("*.pssm"))
        if not genes:
            print(f"[warn] {organelle}: no profiles in {organelle_profiles_dir} - run "
                  f"00_build_reference_profiles.py first", file=sys.stderr)
            continue
        hits = best_hits_per_species(gene_calls, organelle, genes)
        if hits.empty:
            continue
        data_root = sd.repo_data_root(ROOT_DIR, organelle)
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
            rows.append({
                "species": row.species, "organelle": organelle, "gene": row.gene,
                "fasta": fasta_path, "profile": str(organelle_profiles_dir / f"{row.gene}.pssm"),
                "contig": row.contig_id, "start": int(row.start), "end": int(row.end),
                "strand": row.strand,
            })
    return rows


def write_manifest(rows: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    cols = ["species", "organelle", "gene", "fasta", "profile", "contig", "start", "end", "strand"]
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")


def run_orfedit(rows: list[dict], tag: str, work_dir: Path, args) -> tuple[pd.DataFrame, pd.DataFrame]:
    """One orfedit scan-batch over `rows`; (gene calls, edits) as all-string frames."""
    manifest_path = work_dir / f"manifest.{tag}.tsv"
    write_manifest(rows, manifest_path)
    print(f"[info] scan_editing: {tag} manifest has {len(rows)} rows -> {manifest_path}", file=sys.stderr)
    calls_out, edits_out = work_dir / f"edited_gene_calls.{tag}.tsv", work_dir / f"predicted_edits.{tag}.tsv"
    cmd = [
        ORFEDIT, "scan-batch",
        "--manifest", str(manifest_path),
        "--gene-calls-out", str(calls_out),
        "--edits-out", str(edits_out),
        "--threads", str(args.threads),
        "--flank", str(args.flank),
        "--edit-penalty", str(args.edit_penalty),
    ]
    proc = subprocess.run(cmd)
    if proc.returncode != 0:
        print(f"[err] orfedit scan-batch failed (exit {proc.returncode})", file=sys.stderr)
        sys.exit(1)
    return pd.read_csv(calls_out, sep="\t", dtype=str), pd.read_csv(edits_out, sep="\t", dtype=str)


def assign_lineages(runs: dict[str, pd.DataFrame], forced: str | None) -> pd.DataFrame:
    """One row per species with mito calls: the lineage it takes and why."""
    scores: dict[str, dict[str, dict[str, float]]] = {}
    for name, calls in runs.items():
        mito = calls[calls["organelle"] == "mito"]
        for row in mito.itertuples():
            scores.setdefault(row.species, {}).setdefault(name, {})[row.gene] = float(row.score)
    out = []
    for species in sorted(scores):
        chosen, n, totals = lin.choose(scores[species])
        out.append({"species": species, "lineage": forced or chosen, "decided_by": "forced" if forced else "auto",
                    "n_genes": n, **{f"{x}_score": f"{totals[x]:.3f}" if x in totals else "" for x in lin.LINEAGES}})
    return pd.DataFrame(out, columns=lin.COLUMNS)


def keep_lineage(df: pd.DataFrame, lineages: dict[str, str], name: str) -> pd.DataFrame:
    """Rows of lineage `name`'s run that belong in the output: its own mito
    species, plus (for the default run) everything that isn't mito."""
    if df.empty:
        return df
    takes = df["species"].map(lineages).fillna(lin.DEFAULT_LINEAGE) == name
    if name == lin.DEFAULT_LINEAGE:
        return df[(df["organelle"] != "mito") | takes]
    return df[(df["organelle"] == "mito") & takes]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organelle", default="both", choices=["mito", "pltd", "both"])
    ap.add_argument("--species-list", default=None)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--flank", type=int, default=300)
    ap.add_argument("--edit-penalty", type=float, default=2.5)
    ap.add_argument("--lineage", default="auto", choices=lin.CHOICES,
                    help="mito reference lineage: auto picks per species (default), or force one for all")
    ap.add_argument("--gene-calls", default=None, help="path to a gene_calls.tsv-shaped table to anchor "
                     "on (default: annotation/results/gene_calls.tsv, oatk's own calls) - e.g. point this "
                     "at denovo_annotation/results/gene_calls.tsv to use the oatk-independent "
                     "nhmmscan-based calls instead")
    ap.add_argument("--force", action="store_true", help="accepted for CLI consistency; orfedit "
                     "always does a full recompute over the manifest (no incremental skip yet)")
    args = ap.parse_args()

    if not Path(ORFEDIT).exists():
        print(f"[err] orfedit binary not found at {ORFEDIT} - see orfedit's own README "
              f"(`cargo install --path .`)", file=sys.stderr)
        sys.exit(1)

    gene_calls_path = Path(args.gene_calls) if args.gene_calls else ANALYSIS_DIR / "denovo_annotation" / "results" / "gene_calls.tsv"
    if not gene_calls_path.exists():
        print(f"[err] missing {gene_calls_path} - run denovo_annotation/src/04_build_gene_calls.py first "
              f"(or pass --gene-calls)", file=sys.stderr)
        sys.exit(1)
    gene_calls = pd.read_csv(gene_calls_path, sep="\t")

    organelles = ["mito", "pltd"] if args.organelle == "both" else [args.organelle]
    species_filter = sd.load_species_list(Path(args.species_list)) if args.species_list else None
    profiles_dir = ANALYSIS_DIR / "editing" / "results" / "profiles"
    if not any(profiles_dir.glob("*/*.pssm")):
        # a dataset too small for 00_build_reference_profiles.py (>=5 species
        # per gene) uses the bundled ones, built on the 1250-species dataset
        profiles_dir = CODE_DIR / "editing" / "reference" / "profiles"
        print(f"[info] scan_editing: no dataset-built profiles, using bundled {profiles_dir}", file=sys.stderr)

    rows = build_manifest(gene_calls, organelles, species_filter, profiles_dir)
    if not rows:
        print("[err] empty manifest - nothing to scan", file=sys.stderr)
        sys.exit(1)

    work_dir = ANALYSIS_DIR / "editing" / "work"
    work_dir.mkdir(parents=True, exist_ok=True)
    results_dir = ANALYSIS_DIR / "editing" / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    # orfedit writes its output paths itself (truncating), so a
    # --species-list-scoped run must not point it straight at results_dir -
    # it'd truncate the shared table down to just the scoped species. It
    # writes to work_dir instead, then merge_species_scoped_tsv decides how
    # that relates to what's already in results_dir (see its docstring).
    runs = {lin.DEFAULT_LINEAGE: run_orfedit(rows, lin.DEFAULT_LINEAGE, work_dir, args)}
    if "mito" in organelles and args.lineage != lin.DEFAULT_LINEAGE:
        for name in lin.LINEAGES[1:]:
            lrows = build_manifest(gene_calls, ["mito"], species_filter, lin.editing_profiles(CODE_DIR, name))
            if lrows:
                runs[name] = run_orfedit(lrows, name, work_dir, args)

    forced = None if args.lineage == "auto" else args.lineage
    lineage_table = assign_lineages({k: calls for k, (calls, _) in runs.items()}, forced)
    lineages = dict(zip(lineage_table["species"], lineage_table["lineage"]))
    for name, n in lineage_table["lineage"].value_counts().items():
        print(f"[info] scan_editing: {n} species use {name} mito profiles", file=sys.stderr)
    work_gene_calls_out = work_dir / "edited_gene_calls.new.tsv"
    work_edits_out = work_dir / "predicted_edits.new.tsv"
    for i, out in enumerate((work_gene_calls_out, work_edits_out)):
        parts = [keep_lineage(run[i], lineages, name) for name, run in runs.items()]
        pd.concat(parts, ignore_index=True).to_csv(out, sep="\t", index=False)
    work_lineage_out = work_dir / "lineage.new.tsv"
    lineage_table.to_csv(work_lineage_out, sep="\t", index=False)

    io_utils.merge_species_scoped_tsv(results_dir / "edited_gene_calls.tsv", work_gene_calls_out, species_filter)
    io_utils.merge_species_scoped_tsv(results_dir / "predicted_edits.tsv", work_edits_out, species_filter)
    io_utils.merge_species_scoped_tsv(lin.table_path(ANALYSIS_DIR), work_lineage_out, species_filter)
    print(f"[info] scan_editing: merged into {results_dir} "
          f"({'species-scoped' if species_filter else 'full rewrite'})", file=sys.stderr)


if __name__ == "__main__":
    main()
