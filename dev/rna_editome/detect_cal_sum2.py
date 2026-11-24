"""detect_cal_sum2.py SCORES : per truth group (taxonomy) and lineage, the mean per-gene score advantage of the lineage's
profiles over the default ones, over the genes scored both ways (n >= 4). Writes pilot/detect_cal2.margins.tsv."""
import pandas as pd, sys
S = pd.read_csv(sys.argv[1], sep="\t")
T = pd.read_csv("pilot/dataset_taxa.tsv", sep="\t", header=None, names=["species", "taxid", "cls", "order", "family"])
def truth(r):
    if r.cls == "Magnoliopsida": return "angiosperm"
    if r.cls in ("Bryopsida", "Sphagnopsida", "Polytrichopsida", "Andreaeopsida", "Oedipodiopsida"): return "moss"
    if r.order == "Equisetales": return "horsetail"
    if r.cls == "Polypodiopsida": return "fern"
    if r.order == "Lycopodiales": return "lycophyte"
    if r.cls in ("Jungermanniopsida", "Marchantiopsida"): return "liverwort"
    return "other:" + (r.order if r.cls == "Lycopodiopsida" else r.cls)
T["truth"] = T.apply(truth, axis=1); tr = dict(zip(T.species, T.truth))
P = S.pivot_table(index=["species", "gene"], columns="run", values="score")
L = [c[:-6] for c in P.columns if c != "default.plain"]; rows = []
for l in L:
    x = P[[f"{l}.plain", "default.plain"]].dropna(); m = (x[f"{l}.plain"] - x["default.plain"]).groupby("species").agg(["mean", "size"])
    for sp, r in m.iterrows(): rows.append((sp, tr.get(sp, "?"), l, r["mean"], int(r["size"])))
M = pd.DataFrame(rows, columns=["species", "truth", "lineage", "margin", "n"]); M.to_csv("pilot/detect_cal2.margins.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 300)
g = M[M.n >= 4].groupby(["truth", "lineage"]).margin.agg(["size", "min", "median", "max"]).round(0).unstack("lineage")
print(g.to_string())
W = M[M.n >= 4].pivot(index="species", columns="lineage", values="margin").round(0); W["truth"] = W.index.map(lambda s: tr.get(s, "?"))
print(W[~W.truth.isin(["angiosperm", "moss", "liverwort"])].sort_values(["truth", "liverwort"]).to_string())
for l in L:
    x = W[W.truth != l].sort_values(l, ascending=False).head(5); print(l, "top non-members:", [(i, r[l], r.truth) for i, r in x.iterrows()])
