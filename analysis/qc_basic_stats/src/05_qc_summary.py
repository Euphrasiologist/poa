#!/usr/bin/env python3
"""Combine gfa_stats/contig_stats/gene_matrix/core_genes into one QC gate table.

One row per (species, organelle). Every downstream module (annotation,
repeats, synteny, phylogeny) filters its species list against this table's
`status` column (pass/flag/fail) - default is to use `pass` only, with an
explicit override so "what does bad QC actually look like" can itself be a
teaching exercise.

Checks (all computed from data already on disk, or one cheap minimap2 call
per species - no heavy new compute):
  - assembly presence (has_gfa)                          -> fail if absent
  - non-empty graph (parse_status == 'empty')             -> fail
  - dead-end nodes                                        -> flag >0, fail >2
  - fragmentation (n_subgraphs, i.e. disconnected
    components in the assembly graph)                      -> flag >=5, fail >=10
  - circularity (any_circular, read from the RESOLVED
    .ctg.fasta header, not gfatk's raw-graph-level flag -
    see note below)                                        -> flag only if false
                                                              (many real plant
                                                              mitos are legitimately
                                                              non-circular/multipartite)
  - genome size vs genus-level median/MAD (falls back to
    dataset-wide median/MAD when <3 congeners exist)       -> flag if |z|>4
  - core-gene completeness (data-driven core-gene set,
    see 04_core_gene_list.py)                              -> flag <0.80, fail <0.50
  - resolved contig fasta available (vs gfa_fallback only) -> flag if not
  - cross-organelle alignment fraction (see note below)    -> flag >0.15, fail >0.35

Note on cross-contamination: an earlier design used a hardcoded list of
organelle-exclusive gene names and looked for "wrong-organelle" hits in a
species' own annotation. That was verified empirically to never fire:
oatk's mito/plastid annotation each only ever searches its own family HMM
database (acrogymnospermae_mito.fam / angiosperm_pltd.fam), so a plastid
gene name can never appear in a .mito.ctg.bed file regardless of whether
the underlying sequence is contaminated. Instead this script runs one
`minimap2 -x asm5` alignment between a species' own mito and plastid
resolved contigs (when both exist) and flags when an implausibly large
fraction of one organelle's assembly aligns to the other - genuine MTPT
transfer is a small, isolated insert (a few percent); bulk alignment
covering >15-35% of a genome indicates likely mis-binned/contaminated
contigs rather than real biology. This is a coarse QC pre-screen; the
synteny module computes the full per-species MTPT picture (PAFs +
dotplots) for actual study.

Note on circularity: gfatk's "Circular" field describes the raw unitig
GRAPH before path resolution, not the final assembled genome - verified
empirically on real data (e.g. Achillea_maritima plastid: gfatk reports
Circular: false for the 3-unitig graph, yet the resolved contig's own
path visits one repeat unitig twice (u2-,u1+,u3+,u1-) to close a genuine
circle, and its .ctg.fasta header says circular=true). Using gfatk's flag
would flag nearly every plastid genome as "non_circular" - meaningless
noise. any_circular here is read from contig_stats' ctg.fasta-header
value instead, which reflects the structure actually delivered to
students.

Note on unjoined contigs (a has_ctg_fasta blind spot, found via manual
inspection of Azolla_filiculoides prompted by a user question about ferns):
a non-empty .ctg.fasta is not the same as Pathfinder actually resolving
anything - oatk can emit one contig per raw graph segment (nv=1, i.e. zero
nodes joined into any path) when it can't make sense of a graph, and
that file is just as "non-empty" as a real resolved genome. Verified by
comparing contig count against nv=1 count across the whole dataset: 78
species (mostly fern/pteridophyte mitogenomes - Vandenboschia_speciosa,
Asplenium spp., Polystichum spp., Dryopteris spp., Pteridium_aquilinum,
etc., known for repeat-rich/complex mitogenomes that trip up graph
heuristics - plus a scatter of others e.g. Phacelia_tanacetifolia with
353 singleton contigs) have contig_source=="ctg_fasta" but EVERY contig
is nv=1, with 5+ such contigs (the >=5 floor matches the fragmentation
threshold's precedent: real multipartite plant mitogenomes are 2-4
subgenomic circles, so 5+ never-joined pieces isn't plausible biology).
has_ctg_fasta is corrected to treat this case as NOT resolved - it's
functionally identical to Pathfinder producing nothing, just disguised
by a non-empty file - which automatically folds these species into the
existing no_resolved_ctg_fasta reason (tagged "(unjoined)" to distinguish
from a genuinely-empty ctg_fasta) and therefore into every downstream
process that already keys off that reason (revision_plan_v2, and
analysis/linearize's gfatk-resolve fallback tier) with no further changes
needed there.

Note on fragmentation: n_subgraphs distribution (mito, full dataset) shows a
clean natural break - 1154/1250 species have 1-4 subgraphs (consistent with
real single-chromosome or modest multipartite biology, matches the literature),
then a long tail of 46 species with 5-20. Real multipartite plant
mitogenomes are essentially always described as 2-4 sub-genomic circles;
5+ disconnected pieces is fragmentation, not biology (e.g. Clematis_viticella
at 10 subgraphs visually confirmed as "lots of small linear segments" via
the topology_plots module - clearly a failed assembly, not previously
caught since nothing checked n_subgraphs directly before this).

Output: analysis/qc_basic_stats/results/qc_summary.tsv
Columns: species organelle has_gfa has_ctg_fasta n_subgraphs any_circular
         dead_end_nodes total_length size_mad_z core_gene_pct
         cross_organelle_aligned_frac status status_reasons
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Defaults to the old same-tree layout (code and data siblings under one
# repo root); set PLANT_ORGANELLE_DATA_ROOT to point at a separate data
# checkout instead (e.g. when this code is installed from its own repo).
ROOT_DIR = Path(os.environ.get("PLANT_ORGANELLE_DATA_ROOT", str(Path(__file__).resolve().parents[3])))
ANALYSIS_DIR = ROOT_DIR / "analysis"
sys.path.insert(0, str(ANALYSIS_DIR / "common"))

import minimap_utils  # noqa: E402
import species_discovery as sd  # noqa: E402

RESULTS_DIR = ANALYSIS_DIR / "qc_basic_stats" / "results"
COLUMNS = ["species", "organelle", "has_gfa", "has_ctg_fasta", "n_subgraphs", "any_circular",
           "dead_end_nodes", "total_length", "size_mad_z", "core_gene_pct", "n_contigs",
           "pct_unjoined_contigs", "cross_organelle_aligned_frac", "status", "status_reasons"]

CROSS_ORG_FLAG_THRESHOLD = 0.15
CROSS_ORG_FAIL_THRESHOLD = 0.35
SIZE_Z_FLAG_THRESHOLD = 4.0
DEAD_END_FLAG_THRESHOLD = 0
DEAD_END_FAIL_THRESHOLD = 2
FRAGMENTATION_FLAG_THRESHOLD = 5
FRAGMENTATION_FAIL_THRESHOLD = 10
CORE_GENE_FLAG_THRESHOLD = 0.80
CORE_GENE_FAIL_THRESHOLD = 0.50
QC_CORE_GENE_PRESENCE_THRESHOLD = 0.90
# A ctg.fasta where every single contig is nv=1 (no graph segments ever
# joined into a path) with this many or more contigs is Pathfinder having
# produced nothing usable, not a real multipartite genome - see the
# "Note on unjoined contigs" module docstring above.
UNJOINED_MIN_CONTIGS = 5
UNJOINED_PCT_THRESHOLD = 1.0


def aggregate_gfa_stats(gfa_stats: pd.DataFrame, organelle: str) -> pd.DataFrame:
    sub = gfa_stats[gfa_stats["organelle"] == organelle].copy()
    for col in ("n_subgraphs", "dead_end_nodes", "total_length"):
        sub[col] = pd.to_numeric(sub[col], errors="coerce")

    def agg(group):
        parse_status = group["parse_status"].iloc[0]
        has_gfa = parse_status != "missing_gfa"
        if parse_status == "ok":
            total_length = group["total_length"].sum()
            dead_end_nodes = group["dead_end_nodes"].max()
        else:
            total_length = np.nan
            dead_end_nodes = np.nan
        return pd.Series({
            "has_gfa": has_gfa,
            "parse_status": parse_status,
            "n_subgraphs": group["n_subgraphs"].iloc[0],
            "dead_end_nodes": dead_end_nodes,
            "total_length": total_length,
        })

    return sub.groupby("species").apply(agg, include_groups=False).reset_index()


def aggregate_contig_stats(contig_stats: pd.DataFrame, organelle: str) -> pd.DataFrame:
    sub = contig_stats[contig_stats["organelle"] == organelle].copy()
    sub["nv"] = pd.to_numeric(sub["nv"], errors="coerce")

    def agg(group):
        if (group["contig_source"] == "ctg_fasta").any():
            contig_source = "ctg_fasta"
        elif (group["contig_source"] == "gfa_fallback").any():
            contig_source = "gfa_fallback"
        else:
            contig_source = "no_assembly"
        # circularity/joined-ness only meaningful for real resolved contigs, not gfa_fallback unitig soup
        resolved = group[group["contig_source"] == "ctg_fasta"]
        any_circular = (resolved["circular"] == "true").any() if not resolved.empty else np.nan
        n_contigs = len(resolved) if not resolved.empty else np.nan
        pct_unjoined = (resolved["nv"] == 1).mean() if not resolved.empty else np.nan
        return pd.Series({"contig_source": contig_source, "any_circular": any_circular,
                           "n_contigs": n_contigs, "pct_unjoined_contigs": pct_unjoined})

    return sub.groupby("species").apply(agg, include_groups=False).reset_index()


def size_mad_z(df: pd.DataFrame) -> pd.Series:
    """MAD-based z-score of total_length, genus-relative (>=3 congeners) else dataset-wide."""
    df = df.copy()
    df["genus"] = df["species"].str.split("_").str[0]
    valid = df["total_length"].notna()

    def dataset_stats():
        vals = df.loc[valid, "total_length"]
        med = vals.median()
        mad = (vals - med).abs().median()
        return med, mad

    global_med, global_mad = dataset_stats()
    z = pd.Series(np.nan, index=df.index)

    for genus, idx in df.groupby("genus").groups.items():
        idx = [i for i in idx if valid.loc[i]]
        if len(idx) >= 3:
            vals = df.loc[idx, "total_length"]
            med = vals.median()
            mad = (vals - med).abs().median()
            if mad == 0:
                mad = global_mad
        else:
            med, mad = global_med, global_mad
        if not mad or mad == 0:
            z.loc[idx] = 0.0
            continue
        z.loc[idx] = (df.loc[idx, "total_length"] - med) / (1.4826 * mad)
    return z


def core_gene_completeness(gene_matrix: pd.DataFrame, core_genes: pd.DataFrame) -> pd.Series:
    qc_core = set(core_genes.loc[core_genes["pct_present"] >= QC_CORE_GENE_PRESENCE_THRESHOLD, "gene"])
    if not qc_core:
        return pd.Series(dtype=float)
    present = gene_matrix[gene_matrix["gene"].isin(qc_core) & (gene_matrix["n_hits"] >= 1)]
    n_present = present.groupby("species")["gene"].nunique()
    return (n_present / len(qc_core)).rename("core_gene_pct")


def fasta_total_length(fasta_path: str) -> int:
    total = 0
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                total += 0
            else:
                total += len(line.strip())
    return total


def build_status(row: pd.Series) -> tuple[str, str]:
    reasons_fail, reasons_flag = [], []

    if not row["has_gfa"]:
        reasons_fail.append("missing_gfa")
    elif row.get("n_subgraphs") == 0:
        reasons_fail.append("empty_assembly")
    else:
        den = row.get("dead_end_nodes")
        if pd.notna(den):
            if den > DEAD_END_FAIL_THRESHOLD:
                reasons_fail.append("dead_end_nodes>2")
            elif den > DEAD_END_FLAG_THRESHOLD:
                reasons_flag.append("dead_end_nodes>0")
        n_sub = row.get("n_subgraphs")
        if pd.notna(n_sub):
            if n_sub >= FRAGMENTATION_FAIL_THRESHOLD:
                reasons_fail.append(f"fragmented(n_subgraphs>={FRAGMENTATION_FAIL_THRESHOLD})")
            elif n_sub >= FRAGMENTATION_FLAG_THRESHOLD:
                reasons_flag.append(f"fragmented(n_subgraphs>={FRAGMENTATION_FLAG_THRESHOLD})")
        circ = row.get("any_circular")
        if pd.notna(circ) and not bool(circ):
            reasons_flag.append("non_circular")
        z = row.get("size_mad_z")
        if pd.notna(z) and abs(z) > SIZE_Z_FLAG_THRESHOLD:
            reasons_flag.append("size_outlier(|z|>4)")
        cgp = row.get("core_gene_pct")
        if pd.notna(cgp):
            if cgp < CORE_GENE_FAIL_THRESHOLD:
                reasons_fail.append("core_gene_pct<0.50")
            elif cgp < CORE_GENE_FLAG_THRESHOLD:
                reasons_flag.append("core_gene_pct<0.80")
        has_ctg = row.get("has_ctg_fasta")
        if pd.notna(has_ctg) and not bool(has_ctg):
            # distinguish "ctg.fasta exists but every contig is an unjoined
            # singleton" (Pathfinder ran, produced nothing usable) from a
            # genuinely absent/empty ctg.fasta - same underlying failure
            # (nothing resolved), different visible symptom. Substring match
            # on "no_resolved_ctg_fasta" still catches both for every
            # downstream consumer of status_reasons.
            if row.get("contig_source") == "ctg_fasta":
                reasons_flag.append("no_resolved_ctg_fasta(unjoined)")
            else:
                reasons_flag.append("no_resolved_ctg_fasta")
        cof = row.get("cross_organelle_aligned_frac")
        if pd.notna(cof):
            if cof > CROSS_ORG_FAIL_THRESHOLD:
                reasons_fail.append("high_cross_organelle_alignment(>35%)")
            elif cof > CROSS_ORG_FLAG_THRESHOLD:
                reasons_flag.append("elevated_cross_organelle_alignment(>15%)")

    if reasons_fail:
        return "fail", ";".join(reasons_fail + reasons_flag)
    if reasons_flag:
        return "flag", ";".join(reasons_flag)
    return "pass", ""


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organelle", default="both", choices=["mito", "pltd", "both"])
    args = ap.parse_args()

    gfa_stats = pd.read_csv(RESULTS_DIR / "gfa_stats.tsv", sep="\t", dtype=str)
    contig_stats = pd.read_csv(RESULTS_DIR / "contig_stats.tsv", sep="\t", dtype=str)

    organelles = ["mito", "pltd"] if args.organelle == "both" else [args.organelle]
    per_organelle = {}
    for organelle in organelles:
        gm_path = RESULTS_DIR / f"gene_matrix_{organelle}.tsv"
        cg_path = RESULTS_DIR / f"core_genes_{organelle}.tsv"
        gene_matrix = pd.read_csv(gm_path, sep="\t")
        core_genes = pd.read_csv(cg_path, sep="\t")

        gfa_agg = aggregate_gfa_stats(gfa_stats, organelle)
        contig_agg = aggregate_contig_stats(contig_stats, organelle)
        cgp = core_gene_completeness(gene_matrix, core_genes)

        df = gfa_agg.merge(contig_agg, on="species", how="left")
        unjoined_failure = (df["n_contigs"] >= UNJOINED_MIN_CONTIGS) & \
                            (df["pct_unjoined_contigs"] >= UNJOINED_PCT_THRESHOLD)
        df["has_ctg_fasta"] = (df["contig_source"] == "ctg_fasta") & ~unjoined_failure.fillna(False)
        df = df.merge(cgp, on="species", how="left")
        df["size_mad_z"] = size_mad_z(df)
        df["organelle"] = organelle
        per_organelle[organelle] = df

    # cross-organelle alignment: only meaningful where BOTH organelles present for a species
    if set(organelles) == {"mito", "pltd"}:
        common_species = sorted(set(per_organelle["mito"]["species"]) & set(per_organelle["pltd"]["species"]))
        mito_by_sp = per_organelle["mito"].set_index("species")
        pltd_by_sp = per_organelle["pltd"].set_index("species")
        eligible = [s for s in common_species
                    if mito_by_sp.loc[s, "has_ctg_fasta"] and pltd_by_sp.loc[s, "has_ctg_fasta"]]
        mito_root = sd.repo_data_root(ROOT_DIR, "mito")
        pltd_root = sd.repo_data_root(ROOT_DIR, "pltd")
        mito_frac, pltd_frac = {}, {}
        for species in eligible:
            mito_r = sd.resolve_species(mito_root / species, "mito")
            pltd_r = sd.resolve_species(pltd_root / species, "pltd")
            if mito_r.status != "ok" or pltd_r.status != "ok":
                continue
            try:
                records = minimap_utils.run_asm5(mito_r.ctg_fasta, pltd_r.ctg_fasta)
            except Exception:  # noqa: BLE001
                continue
            mito_len = fasta_total_length(mito_r.ctg_fasta)
            pltd_len = fasta_total_length(pltd_r.ctg_fasta)
            mito_aligned = minimap_utils.total_aligned_length(records, side="target")
            pltd_aligned = minimap_utils.total_aligned_length(records, side="query")
            if mito_len:
                mito_frac[species] = mito_aligned / mito_len
            if pltd_len:
                pltd_frac[species] = pltd_aligned / pltd_len
        per_organelle["mito"]["cross_organelle_aligned_frac"] = per_organelle["mito"]["species"].map(mito_frac)
        per_organelle["pltd"]["cross_organelle_aligned_frac"] = per_organelle["pltd"]["species"].map(pltd_frac)
    else:
        for organelle in organelles:
            per_organelle[organelle]["cross_organelle_aligned_frac"] = np.nan

    all_rows = []
    for organelle in organelles:
        df = per_organelle[organelle]
        statuses = df.apply(build_status, axis=1)
        df["status"] = statuses.apply(lambda t: t[0])
        df["status_reasons"] = statuses.apply(lambda t: t[1])
        all_rows.append(df[COLUMNS])

    out_df = pd.concat(all_rows, ignore_index=True)
    out_path = RESULTS_DIR / "qc_summary.tsv"
    out_df.to_csv(out_path, sep="\t", index=False)

    counts = out_df.groupby(["organelle", "status"]).size()
    print(f"[info] qc_summary: {len(out_df)} rows -> {out_path}\n{counts}", file=sys.stderr)


if __name__ == "__main__":
    main()
