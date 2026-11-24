"""rna_sum.py TAG... : RNA-scored gene tables (pilot/TAG.genes.tsv mito, TAG.pltd.genes.tsv) summed by group (moss; liverwort,
complex thalloid; liverwort, leafy and simple thalloid) and kind: genes with >= 90% RNA coverage, aligned >= 99%, frameshift
genes, RNA edits, edits called, calls right."""
import sys, pandas as pd
T = pd.read_csv("pilot/dataset_taxa.tsv", sep="\t", header=None, names=["species", "taxid", "cls", "order", "family"]).set_index("species")
grp = lambda s: "moss" if T.cls.get(s) in ("Bryopsida", "Sphagnopsida", "Polytrichopsida", "Andreaeopsida") else {"Marchantiopsida": "liverwort, complex thalloid", "Jungermanniopsida": "liverwort, leafy"}.get(T.cls.get(s), T.cls.get(s))
rows = []
for tag in sys.argv[1:]:
    for org, f in (("mito", f"pilot/{tag}.genes.tsv"), ("pltd", f"pilot/{tag}.pltd.genes.tsv")):
        try: d = pd.read_csv(f, sep="\t")
        except FileNotFoundError: continue
        d["group"] = d.species.map(grp); c = d[d.covered >= 0.9]
        for (g, k), x in c.groupby(["group", "kind"]):
            rna = x.rna_CU.sum() + x.rna_UC.sum()
            rows.append(dict(organelle=org, group=g, kind=k, tag=tag, species=x.species.nunique(), genes=len(x), aligned=int((x.frac >= 0.99).sum()), fs_genes=int((x.n_fs > 0).sum()),
                             rna_edits=int(rna), called=int(x.called.sum()), right=int(x.right.sum()), found=f"{100 * x.right.sum() / max(rna, 1):.0f}%", precision=f"{100 * x.right.sum() / max(x.called.sum(), 1):.0f}%",
                             protein_ok=int((x.rna_protein == "ok").sum()) if "rna_protein" in x else ""))
pd.set_option("display.width", 250)
print(pd.DataFrame(rows).sort_values(["organelle", "group", "kind", "tag"]).to_string(index=False))
