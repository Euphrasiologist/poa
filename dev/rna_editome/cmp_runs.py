"""cmp_runs.py TAG_A TAG_B : poa runs of the same species (pilot/TAG/<species>/tables) before and after a change: lineage
calls that differ, and per taxonomic group and organelle the single-exon genes aligned >= 99%, spliced genes complete,
edits called and frameshift calls. Writes pilot/TAG_A.vs.TAG_B.tsv (per species)."""
import sys
from pathlib import Path
import pandas as pd
A, B = sys.argv[1], sys.argv[2]
T = pd.read_csv("pilot/dataset_taxa.tsv", sep="\t", header=None, names=["species", "taxid", "cls", "order", "family"]).set_index("species")
def group(s):
    if s not in T.index: return "?"
    c, o = T.cls[s], T.order[s]
    if c == "Magnoliopsida": return "angiosperm"
    if c in ("Bryopsida", "Sphagnopsida", "Polytrichopsida", "Andreaeopsida", "Oedipodiopsida"): return "moss"
    if c in ("Jungermanniopsida", "Marchantiopsida"): return "liverwort"
    if o == "Equisetales": return "horsetail"
    if c == "Polypodiopsida": return "fern"
    if o == "Lycopodiales": return "clubmoss"
    return "other"
def read(tag, sp):
    t = Path("pilot") / tag / sp / "tables"
    if not (t / "edited_gene_calls.tsv").exists(): return None
    r = {}
    try: L = pd.read_csv(t / "lineage.tsv", sep="\t", dtype=str).fillna(""); r["lineage"] = L.lineage.iloc[0]; r["sites"] = L.edit_sites.iloc[0] if "edit_sites" in L else ""
    except Exception: r["lineage"] = r["sites"] = ""
    C = pd.read_csv(t / "edited_gene_calls.tsv", sep="\t"); E = pd.read_csv(t / "predicted_edits.tsv", sep="\t"); F = pd.read_csv(t / "frameshifts.tsv", sep="\t")
    G = pd.read_csv(t / "reconstructed_genes.tsv", sep="\t"); X = pd.read_csv(t / "reconstructed_edits.tsv", sep="\t"); RF = pd.read_csv(t / "reconstructed_frameshifts.tsv", sep="\t")
    for o in ("mito", "pltd"):
        c = C[C.organelle == o]; g = G[G.organelle == o]
        r[f"{o}.single"] = len(c); r[f"{o}.aligned"] = int((c.n_columns_aligned / c.n_columns_total >= 0.99).sum()); r[f"{o}.edits"] = int((E.organelle == o).sum()); r[f"{o}.fs"] = int((F.organelle == o).sum()) if len(F) else 0
        r[f"{o}.spliced"] = len(g); r[f"{o}.complete"] = int((g.complete.astype(str).str.lower() == "true").sum()); r[f"{o}.sp_edits"] = int((X.organelle == o).sum()) if len(X) else 0; r[f"{o}.sp_fs"] = int((RF.organelle == o).sum()) if len(RF) else 0
    return r
rows = []
for d in sorted(p.name for p in (Path("pilot") / B).iterdir() if p.is_dir()):
    a, b = read(A, d), read(B, d)
    if a is None or b is None: continue
    rows.append({"species": d, "group": group(d), **{f"A.{k}": v for k, v in a.items()}, **{f"B.{k}": v for k, v in b.items()}})
D = pd.DataFrame(rows); D.to_csv(f"pilot/{A}.vs.{B}.tsv", sep="\t", index=False)
pd.set_option("display.width", 320); pd.set_option("display.max_rows", 200)
ch = D[D["A.lineage"] != D["B.lineage"]]
print(f"{len(D)} species in both; lineage differs in {len(ch)}:"); print(ch[["species", "group", "A.lineage", "B.lineage"]].to_string(index=False) if len(ch) else "")
print("lineage by group (B):"); print(D.groupby(["group", "B.lineage"]).size().unstack(fill_value=0).to_string())
print("edit sites by group (B):"); print(D.groupby(["group", "B.sites"]).size().unstack(fill_value=0).to_string())
for o in ("mito", "pltd"):
    cols = [f"{o}.{k}" for k in ("single", "aligned", "edits", "fs", "spliced", "complete", "sp_edits", "sp_fs")]
    g = D.groupby("group")[[f"{x}.{c}" for c in cols for x in "AB"]].sum(); g.columns = [c.replace(f"{o}.", "") for c in g.columns]
    print(o); print(g.to_string())
