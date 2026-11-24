"""score.py BASE TAG... : gene ends of orfedit's calls (ends/TAG.calls.tsv, from run.sh) against the GenBank annotation of the test
genomes (ends/genomes.tsv): per single-exon annotated CDS the call of that gene overlapping it most; start and stop within 3 bp
(the stop codon in or out). Each TAG is compared with BASE: ends right in one and not in the other. Writes ends/TAG.ends.tsv."""
import sys, warnings
from pathlib import Path
import pandas as pd
warnings.simplefilter("ignore")
from Bio import SeqIO
H = Path(__file__).resolve().parent; R = H.parents[2]
GB = [R / "src/poa/pipeline/trans_splicing/reference/genbank", H.parent / "extra_genomes" / "gb", R / "src/poa/pipeline/trans_splicing/reference/plastid/genbank"]
SYN = {"tatc": "mttb", "mat-r": "matr", "orfx": "mttb", "orf25": "atp4", "orfb": "atp8", "ccb206": "ccmb", "ccb256": "ccmc", "ccb452": "ccmfc", "ccb382": "ccmfn1", "ccb203": "ccmfn2", "cox-2": "cox2", "coxii": "cox2", "yejv": "ccmfc", "yejr": "ccmfn"}
MITO = set("nad1 nad2 nad3 nad4 nad4l nad5 nad6 nad7 nad9 cob cox1 cox2 cox3 atp1 atp4 atp6 atp8 atp9 ccmb ccmc ccmfc ccmfn ccmfn1 ccmfn2 matr mttb rpl2 rpl5 rpl10 rpl16 rps1 rps2 rps3 rps4 rps7 rps10 rps11 rps12 rps13 rps14 rps19 sdh3 sdh4".split())
norm = lambda g: SYN.get(g.lower(), g.lower())
G = pd.read_csv(H / "genomes.tsv", sep="\t", header=None, names=["name", "org", "fasta", "acc"])
truth = []
for g in G.itertuples():
    rec = SeqIO.read(next(d / f"{g.acc}.gb" for d in GB if (d / f"{g.acc}.gb").exists()), "genbank")
    for f in rec.features:
        if f.type != "CDS" or "gene" not in f.qualifiers or "translation" not in f.qualifiers or "pseudo" in f.qualifiers or len(f.location.parts) != 1: continue
        gene = norm(f.qualifiers["gene"][0])
        if g.org == "mito" and gene not in MITO: continue
        truth.append(dict(genome=g.name, org=g.org, gene=gene, a=int(f.location.start), b=int(f.location.end), strand=f.location.strand))
T = pd.DataFrame(truth)
def score(tag):
    C = pd.read_csv(H / f"{tag}.calls.tsv", sep="\t"); C["g"] = C.gene.str.split("|").str[0].map(norm)
    by = {k: v for k, v in C.groupby(["species", "g"])}; genes = set(C.g); rows = []
    for t in T.itertuples():
        if t.gene not in genes and not (t.gene in ("ccmfn1", "ccmfn2")): continue  # no profile for it at all
        cand = [by[k] for k in ((t.genome, t.gene), *(((t.genome, x) for x in ("ccmfn", "ccmfn1", "ccmfn2")) if t.gene.startswith("ccmfn") else ())) if k in by]
        best = None
        for r in (r for c in cand for r in c.itertuples()):
            pa, pb = sorted((int(r.corrected_start), int(r.corrected_end))); ov = min(t.b, pb) - max(t.a, pa)
            if ov > 0 and (best is None or ov > best[0]): best = (ov, pa, pb, r)
        if best is None: rows.append(dict(genome=t.genome, org=t.org, gene=t.gene, a=t.a, called=False, s=False, e=False)); continue
        _, pa, pb, r = best
        ds, de = (pa - t.a, pb - t.b) if t.strand == 1 else (t.b - pb, t.a - pa); de = min((de, de + 3, de - 3), key=abs)
        rows.append(dict(genome=t.genome, org=t.org, gene=t.gene, a=t.a, b=t.b, ann_strand=t.strand, call=r.gene, called=True, s=abs(ds) <= 3, e=abs(de) <= 3, ds=ds, de=de, aligned=f"{r.n_columns_aligned}/{r.n_columns_total}", stop=r.includes_stop, fs=r.n_frameshifts, istop=getattr(r, "n_internal_stops", 0)))
    D = pd.DataFrame(rows); D.to_csv(H / f"{tag}.ends.tsv", sep="\t", index=False); return D.set_index(["genome", "gene", "a"])
base = score(sys.argv[1])
def line(tag, D):
    out = []
    for org in ("mito", "pltd"):
        x = D[D.org == org]; out.append(f"{org}: {len(x)} CDS, called {int(x.called.sum())}, stop {int(x.e.sum())}, start {int(x.s.sum())}, both {int((x.s & x.e).sum())}, with a stop codon {int(x.get('stop', pd.Series(dtype=bool)).fillna(False).astype(bool).sum())}, frame changes {int(x.get('fs', pd.Series(dtype=float)).fillna(0).sum())}, internal stops {int(x.get('istop', pd.Series(dtype=float)).fillna(0).sum())}")
    print(f"{tag:14s} " + " | ".join(out))
line(sys.argv[1], base)
for tag in sys.argv[2:]:
    D = score(tag); line(tag, D); j = base.join(D, lsuffix="_0", how="outer")
    for end, lab in (("e", "stop"), ("s", "start")):
        lost = j[j[f"{end}_0"].fillna(False).astype(bool) & ~j[end].fillna(False).astype(bool)]; won = j[~j[f"{end}_0"].fillna(False).astype(bool) & j[end].fillna(False).astype(bool)]
        print(f"   {lab}: gained {len(won)}, lost {len(lost)}" + (": " + ", ".join(f"{g[0].split('_')[0][:5]}.{g[1]}({int(r['d' + end]) if pd.notna(r['d' + end]) else 'nc'})" for g, r in lost.iterrows()) if len(lost) else ""))
