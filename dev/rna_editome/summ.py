"""summ.py TAG SPECIES... : gene-level summary of poa runs (no RNA needed): per gene, copies, aligned >= 99% / complete, frameshift calls"""
import sys
from pathlib import Path
import pandas as pd
tag, species = sys.argv[1], sys.argv[2:]; rows = []
for sp in species:
    T = Path(tag) / sp / "tables"
    if not T.exists(): continue
    E = pd.read_csv(T / "edited_gene_calls.tsv", sep="\t"); E = E[E.organelle == "mito"]
    FS = pd.read_csv(T / "frameshifts.tsv", sep="\t"); FS = FS[FS.organelle == "mito"]
    for r in E.itertuples(): rows.append((sp, "single", r.gene, r.n_columns_aligned / r.n_columns_total >= 0.99, int(((FS.gene == r.gene) & (FS.hit_start == r.hit_start)).sum()), int(r.n_edits)))
    G = pd.read_csv(T / "reconstructed_genes.tsv", sep="\t"); G = G[G.organelle == "mito"]
    RF = pd.read_csv(T / "reconstructed_frameshifts.tsv", sep="\t")
    for g in G.itertuples(): rows.append((sp, "spliced", g.gene, str(g.complete).lower() == "true", int(((RF.gene == g.gene) & (RF.organelle == "mito")).sum()) if len(RF) else 0, int(g.whole_gene_n_edits)))
D = pd.DataFrame(rows, columns=["species", "kind", "gene", "full", "n_fs", "n_edits"])
for kind, d in D.groupby("kind"):
    print(f"{tag} {kind}: {d.species.nunique()} species, {len(d)} genes, full {d.full.sum()}, with frameshift calls {(d.n_fs > 0).sum()}, edits called {d.n_edits.sum()}")
if "-v" in sys.argv or True:
    x = D.groupby(["kind", "gene"]).agg(n=("gene", "size"), full=("full", "sum"), fs=("n_fs", lambda v: (v > 0).sum()))
    bad = x[(x.full < x.n) | (x.fs > 0)]
    print("  not full / frameshift:", "; ".join(f"{g} {r.full}/{r.n} fs{r.fs}" for (k, g), r in bad.iterrows()))
