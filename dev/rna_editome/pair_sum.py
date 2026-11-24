"""pair_sum.py TAG TAG... : RNA-scored gene tables (pilot/TAG.genes.tsv mito, TAG.pltd.genes.tsv) of the same species under
different references, summed over the genes (species, kind, gene) every TAG has with >= 90% RNA coverage: genes, aligned
>= 99%, frameshift genes, RNA edits, edits called, calls right, RNA-edited call reads cleanly."""
import sys, pandas as pd
rows = []
for org, ext in (("mito", "genes.tsv"), ("pltd", "pltd.genes.tsv")):
    try: D = {t: pd.read_csv(f"pilot/{t}.{ext}", sep="\t").drop_duplicates(["species", "kind", "gene"]).set_index(["species", "kind", "gene"]) for t in sys.argv[1:]}
    except FileNotFoundError: continue
    keep = None
    for t, d in D.items():
        k = set(d.index[d.covered >= 0.9]); keep = k if keep is None else keep & k
    for t, d in D.items():
        d = d.loc[sorted(keep)].reset_index()
        for k, x in d.groupby("kind"):
            rna = x.rna_CU.sum() + x.rna_UC.sum()
            rows.append(dict(organelle=org, kind=k, tag=t, species=x.species.nunique(), genes=len(x), aligned=int((x.frac >= 0.99).sum()), fs_genes=int((x.n_fs > 0).sum()), rna_edits=int(rna),
                             called=int(x.called.sum()), right=int(x.right.sum()), found=f"{100 * x.right.sum() / max(rna, 1):.0f}%", precision=f"{100 * x.right.sum() / max(x.called.sum(), 1):.0f}%",
                             protein_ok=int((x.rna_protein == "ok").sum()) if "rna_protein" in x else ""))
pd.set_option("display.width", 250)
print(pd.DataFrame(rows).sort_values(["organelle", "kind", "tag"]).to_string(index=False))
