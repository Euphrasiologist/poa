"""score_gb.py TAG NAME... : poa tables (pilot/TAG/NAME) for a GenBank genome (extra_genomes/accessions.tsv) against
the genome's own annotation. Per annotated protein gene: exons placed (both ends within 3 bp of an annotated exon;
the first exon's start and the last exon's end may differ by the start codon chosen or the stop codon), genes with
every exon placed, frameshift calls (the annotated genes are intact, so any is false)."""
import sys, warnings
from pathlib import Path
import pandas as pd
warnings.simplefilter("ignore")
from Bio import SeqIO
H = Path(__file__).resolve().parent.parent; R = H.parent.parent
GB = R / "src/poa/pipeline/trans_splicing/reference/genbank"
ACC = dict(l.split() for l in open(H / "extra_genomes" / "accessions.tsv"))
SYN = {"tatc": "mttB", "mat-r": "matR"}
tag, names = sys.argv[1], sys.argv[2:]; rows = []
for name in names:
    T = H / "pilot" / tag / name / "tables"
    if not T.exists(): print("skip", name); continue
    rec = SeqIO.read(GB / f"{ACC[name]}.gb", "genbank"); ann = {}
    for f in rec.features:
        if f.type != "CDS" or "gene" not in f.qualifiers or "translation" not in f.qualifiers: continue
        g = f.qualifiers["gene"][0]; g = SYN.get(g.lower(), g)
        parts = list(f.location.parts)
        if "".join(str(p.extract(rec.seq)) for p in parts) != str(f.extract(rec.seq)): parts = parts[::-1]
        ann.setdefault(g, []).append([(int(p.start), int(p.end)) for p in parts])
    E = pd.read_csv(T / "edited_gene_calls.tsv", sep="\t"); E = E[E.organelle == "mito"]
    FS = pd.read_csv(T / "frameshifts.tsv", sep="\t"); FS = FS[FS.organelle == "mito"]
    X = pd.read_csv(T / "reconstructed_exons.tsv", sep="\t"); X = X[X.organelle == "mito"]
    G = pd.read_csv(T / "reconstructed_genes.tsv", sep="\t"); G = G[G.organelle == "mito"]
    RF = pd.read_csv(T / "reconstructed_frameshifts.tsv", sep="\t"); RF = RF[RF.organelle == "mito"] if len(RF) else RF
    pred = {}
    for g in G.gene: pred[g] = ("spliced", [tuple(sorted((int(x.start), int(x.end)))) for x in X[X.gene == g].sort_values("slot").itertuples()], int((RF.gene == g).sum()) if len(RF) else 0)
    for r in E.sort_values("score", ascending=False).itertuples():
        if r.gene in pred: continue
        lo, hi = sorted((int(r.corrected_start), int(r.corrected_end)))
        pred[r.gene] = ("single", [(lo, hi)], int(((FS.gene == r.gene) & (FS.hit_start == r.hit_start)).sum()))
    for g, loci in ann.items():
        if g.startswith(("ORF", "orf", "mat-")): continue
        if g not in pred: rows.append(dict(genome=name, gene=g, kind="none", n_exons=len(loci[0]), placed=0, all=False, inner=False, n_fs=0, pred_exons=0)); continue
        kind, px, nfs = pred[g]
        best = None
        for ex in loci:  # the annotated copy the prediction fits best
            k = len(ex); placed = inner = 0
            for i, (a, b) in enumerate(ex):
                # ends at an intron must match; the gene's own ends may differ by a start codon choice (<= 60 bp) or the stop codon
                tol_a, tol_b = (60 if i == 0 else 3), (6 if i == k - 1 else 3)
                hit = any(abs(p[0] - a) <= 3 and abs(p[1] - b) <= 3 for p in px)
                fwd = any(abs(p[0] - a) <= tol_a and abs(p[1] - b) <= tol_b for p in px); rev = any(abs(p[0] - a) <= tol_b and abs(p[1] - b) <= tol_a for p in px)
                placed += hit; inner += (fwd or rev)
            if best is None or inner > best[1]: best = (placed, inner, k)
        placed, inner, k = best
        rows.append(dict(genome=name, gene=g, kind=kind, n_exons=k, placed=placed, all=placed == k, inner=inner == k, n_fs=nfs, pred_exons=len(px)))
D = pd.DataFrame(rows); D.to_csv(H / "pilot" / f"{tag}.gb.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
for name, d in D.groupby("genome"):
    print(f"{tag} {name}: {len(d)} annotated genes, {d.n_exons.sum()} exons; exons placed within 3 bp {d.placed.sum()}; genes with every exon placed {d['all'].sum()} "
          f"(allowing for start/stop codon choice {d.inner.sum()}); not called {(d.kind == 'none').sum()}; genes with a frameshift call {(d.n_fs > 0).sum()}")
bad = D[~D.inner | (D.n_fs > 0)]
print(bad.groupby("gene").agg(n=("gene", "size"), exons=("n_exons", "first"), pred=("pred_exons", lambda x: ",".join(map(str, x))), placed=("placed", lambda x: ",".join(map(str, x))), fs=("n_fs", "sum")).to_string())
