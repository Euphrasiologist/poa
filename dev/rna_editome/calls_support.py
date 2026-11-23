"""calls_support.py TAG ORG SPECIES... : every edit poa called (pilot/TAG), by how much of the RNA shows it"""
import sys, glob
from collections import Counter
from pathlib import Path
import pandas as pd
H = Path(__file__).resolve().parent.parent; sys.path.insert(0, str(H)); import rna, boot
tag, org, species = sys.argv[1], sys.argv[2], sys.argv[3:]
tot = Counter()
for sp in species:
    T = H / "pilot" / tag / sp / "tables"
    pile = boot.load_piles(sp) if org == "mito" else rna.load_pile(H / f"tol.{sp}.pltd.pile.tsv.gz")
    gc = pd.read_csv(T / "gene_calls.tsv", sep="\t"); gc = gc[gc.organelle == org]
    E = pd.read_csv(T / "edited_gene_calls.tsv", sep="\t"); E = E[E.organelle == org]
    PE = pd.read_csv(T / "predicted_edits.tsv", sep="\t"); PE = PE[PE.organelle == org]
    X = pd.read_csv(T / "reconstructed_exons.tsv", sep="\t"); X = X[X.organelle == org]
    RE = pd.read_csv(T / "reconstructed_edits.tsv", sep="\t")
    RE = RE[RE.organelle == org] if "organelle" in RE else RE
    sites = []
    for r in E.itertuples():
        m = gc[(gc.gene == r.gene) & (gc.contig_id == r.contig) & (gc.start < r.hit_end + 300) & (gc.end > r.hit_start - 300)]
        if m.empty: continue
        for e in PE[(PE.gene == r.gene) & (PE.hit_start == r.hit_start)].itertuples(): sites.append(("single", r.contig, int(e.genomic_pos) + 1, m.iloc[0].strand, e.edit))
    for g, xs in X.groupby("gene"):
        for e in RE[RE.gene == g].itertuples(): sites.append(("spliced", xs.contig.iloc[0], int(e.genomic_pos) + 1, xs.strand.iloc[0], e.edit))
    for kind, c, p, strand, edit in sites:
        if (c, p) not in pile: continue
        ref, n, cnt = pile[(c, p)]
        if n < 10: tot[(kind, edit, "depth < 10")] += 1; continue
        alt = {"C>U": "T", "U>C": "C"}[edit]; alt = alt if strand == "+" else rna.comp[alt]; na = cnt[alt]
        if na >= 3 and na / n >= 0.2: k = "in the RNA (>= 20%)"
        elif na >= 2 and na / n >= 0.05: k = "partly (5-20%)"
        elif na >= 1: k = "a read or more (< 5%)"
        else: k = "none"
        tot[(kind, edit, k)] += 1
D = pd.Series(tot).unstack(fill_value=0); print(tag, org); print(D.to_string())
