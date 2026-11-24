"""detect_cal.py OUT THREADS [TAG...] : scores for calibrating automatic lineage detection, over the whole dataset.
Every dataset species' mito gene calls (the dataset's gene_calls.tsv) are aligned with orfedit to the default
(angiosperm) profiles and to each other lineage's own profiles, plain (C-to-U only) and with the lineage's own
editing options. Writes OUT: species, run, gene, score (a gene's best copy), n_columns_aligned, n_columns_total.
TAG...: pilot runs (pilot/TAG/<species>/tables/gene_calls.tsv) whose fresh gene calls replace the dataset's stored ones
for their species (nine fern assemblies were replaced after the stored calls were made). Plain C-to-U runs only.
Reads the dataset repo only; work files go to pilot/detect_cal.work (removed)."""
import sys, shutil, math
from pathlib import Path
import pandas as pd
R = Path("/lustre/scratch122/tol/teams/blaxter/users/mb39/plant_organelle_annotator"); H = R / "build/rnaedit"
DB = Path("/lustre/scratch122/tol/teams/blaxter/users/mb39/plant_organellar_database")
from poa.pipeline.common import lineage as lin
from poa.pipeline.common.tools import cargo_tool
from poa.pipeline.editing import scan
from poa.pipeline.editing.lineage_profiles import refseq_profiles_dir
out, threads, tags = sys.argv[1], int(sys.argv[2]), sys.argv[3:]
P = R / "src/poa/pipeline"; orfedit = cargo_tool("orfedit")
calls = pd.read_csv(DB / "analysis/denovo_annotation/results/gene_calls.tsv", sep="\t"); calls = calls[calls.organelle == "mito"]
fresh = [pd.read_csv(f, sep="\t") for tag in tags for f in sorted((H / "pilot" / tag).glob("*/tables/gene_calls.tsv"))]
if fresh:
    fresh = pd.concat(fresh); fresh = fresh[fresh.organelle == "mito"].drop_duplicates()
    calls = pd.concat([calls[~calls.species.isin(set(fresh.species))], fresh], ignore_index=True)
W = H / "pilot" / "detect_cal.work"; shutil.rmtree(W, ignore_errors=True); W.mkdir(parents=True)
others = ("moss", "liverwort", "fern", "horsetail", "lycophyte")
own = {L: lin.editing_profiles(P, L) for L in others}
genes = {L: {p.stem for p in (own[L] / "mito").glob("*.pssm")} - scan.lineage_spliced(L) for L in others}
union = set().union(*genes.values())
drows = [r for r in scan.build_manifest(DB, calls, ["mito"], None, [refseq_profiles_dir(P), P / "editing/reference/profiles"]) if r["gene"].split("|")[0] in union]
opts = dict(threads=threads, flank=scan.DEFAULT_FLANK, edit_penalty=scan.DEFAULT_EDIT_PENALTY, frameshift_penalty=scan.DEFAULT_FRAMESHIFT_PENALTY)
res = []
def go(tag, rows, **kw):
    df = scan.run_orfedit(orfedit, rows, tag, W, **{**opts, **kw})[0]
    df = df.assign(score=df.score.astype(float)).sort_values("score", ascending=False).drop_duplicates(["species", "gene"])
    res.append(df[["species", "gene", "score", "n_columns_aligned", "n_columns_total"]].assign(run=tag)); print(tag, len(rows), "rows", flush=True)
go("default.plain", drows)
for L in others:
    rows = [r for r in scan.build_manifest(DB, calls, ["mito"], None, own[L]) if r["gene"].split("|")[0] in genes[L]]
    go(f"{L}.plain", rows, **({"edit_penalty": lin.edit_penalty(L)} if lin.edit_penalty(L) else {}))  # at the lineage's own edit cost, as scan.py runs it
pd.concat(res).to_csv(out, sep="\t", index=False); shutil.rmtree(W, ignore_errors=True)
