"""scope_sum.py TAG : poa runs in pilot/TAG/<species>/tables summarised by taxonomic group (pilot/dataset_taxa.tsv) and organelle,
no truth needed: single-exon genes aligned >= 99% of their profile, spliced genes complete, edits called and frameshift calls
per species. Writes pilot/TAG.scope.tsv (per species x organelle) and prints the group table."""
import sys
from pathlib import Path
import pandas as pd
tag = sys.argv[1]; P = Path("pilot") / tag
T = pd.read_csv("pilot/dataset_taxa.tsv", sep="\t", header=None, names=["species", "taxid", "cls", "order", "family"])
def group(r):
    if r.cls == "Magnoliopsida": return "angiosperm"
    if r.cls in ("Bryopsida", "Sphagnopsida", "Polytrichopsida", "Andreaeopsida", "Oedipodiopsida"): return "moss"
    if r.cls == "Marchantiopsida": return "liverwort, complex thalloid"
    if r.cls == "Jungermanniopsida": return "liverwort, leafy and simple thalloid"
    if r.order == "Equisetales": return "horsetail"
    if r.cls == "Polypodiopsida": return "fern"
    if r.cls == "Lycopodiopsida": return {"Lycopodiales": "clubmoss"}.get(r.order, r.order)
    return {"Anthocerotopsida": "hornwort"}.get(r.cls, "gymnosperm")
grp = dict(zip(T.species, T.apply(group, axis=1))); rows = []
for d in sorted(P.glob("*/tables")):
    sp = d.parent.name
    def rd(n):
        f = d / f"{n}.tsv"; return pd.read_csv(f, sep="\t") if f.exists() and f.stat().st_size > 1 else pd.DataFrame(columns=["organelle"])
    E, FS, G, RF, L = rd("edited_gene_calls"), rd("frameshifts"), rd("reconstructed_genes"), rd("reconstructed_frameshifts"), rd("lineage")
    for o in ("mito", "pltd"):
        e, g = E[E.organelle == o], G[G.organelle == o]
        rows.append(dict(species=sp, group=grp.get(sp, "?"), organelle=o, lineage=L.lineage.iloc[0] if len(L) else "",
                         single=len(e), single_full=int((e.n_columns_aligned / e.n_columns_total >= 0.99).sum()) if len(e) else 0, single_edits=int(e.n_edits.sum()) if len(e) else 0,
                         single_fs=int((FS.organelle == o).sum()), spliced=len(g), spliced_complete=int((g.complete.astype(str).str.lower() == "true").sum()) if len(g) else 0,
                         spliced_edits=int(g.whole_gene_n_edits.sum()) if len(g) else 0, spliced_fs=int((RF.organelle == o).sum()) if len(RF) else 0))
D = pd.DataFrame(rows); D.to_csv(f"pilot/{tag}.scope.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
for o in ("mito", "pltd"):
    d = D[D.organelle == o]; g = d.groupby("group")
    out = pd.DataFrame({"species": g.size(), "single-exon aligned": (100 * g.single_full.sum() / g.single.sum()).round(0), "spliced complete": (100 * g.spliced_complete.sum() / g.spliced.sum()).round(0),
                        "edits/species single": (g.single_edits.sum() / g.size()).round(0), "edits/species spliced": (g.spliced_edits.sum() / g.size()).round(0),
                        "fs/species single": (g.single_fs.sum() / g.size()).round(1), "fs/species spliced": (g.spliced_fs.sum() / g.size()).round(1)})
    print("==", o); print(out.to_string())
print(D[D.organelle == "mito"].groupby(["group", "lineage"]).size().to_string())
