"""score_models.py LABEL SOURCE NAME... : protein-gene models of a GenBank genome (extra_genomes/accessions.tsv) against
the genome's own annotation, the same way for poa and for PMGA.
SOURCE is poa:TAG (tables in pilot/TAG/NAME) or pmga:DIR (build/pmga/DIR/<name>/<name>.gff, <name> being NAME without
_refseq). One row per annotated CDS (every copy): the prediction of that gene (any copy) that fits it best; exons placed
(both ends within 3 bp; abutting predicted exons also count as the one exon they make), every exon placed, the same
allowing the gene's own ends to differ by a start codon choice (<= 60 bp) or the stop codon (<= 6 bp), and exact (every
exon end to the base, stop codon in or out). Writes pilot/LABEL.models.tsv and prints totals."""
import sys, warnings, re
from pathlib import Path
import pandas as pd
warnings.simplefilter("ignore")
from Bio import SeqIO
H = Path(__file__).resolve().parent.parent; R = H.parent.parent
GB0 = R / "src/poa/pipeline/trans_splicing/reference/genbank"
ACC = dict(l.split() for l in open(H / "extra_genomes" / "accessions.tsv"))
SYN = {"tatc": "mttb", "mat-r": "matr", "orfx": "mttb", "orf25": "atp4", "orfb": "atp8", "ccb206": "ccmb", "ccb256": "ccmc", "ccb452": "ccmfc", "ccb382": "ccmfn1", "ccb203": "ccmfn2",
       "ccmfn1": "ccmfn1", "ccmfn2": "ccmfn2", "cox-2": "cox2", "coxii": "cox2", "yejv": "ccmfc", "yejr": "ccmfn"}
# the mitochondrial protein genes (plastid-derived CDS and ORFs annotated in some records are left out)
MITO = set("nad1 nad2 nad3 nad4 nad4l nad5 nad6 nad7 nad9 cob cox1 cox2 cox3 atp1 atp4 atp6 atp8 atp9 ccmb ccmc ccmfc ccmfn ccmfn1 ccmfn2 matr mttb "
           "rpl2 rpl5 rpl10 rpl16 rps1 rps2 rps3 rps4 rps7 rps10 rps11 rps12 rps13 rps14 rps19 sdh3 sdh4".split())
norm = lambda g: SYN.get(g.lower(), g.lower())
label, source, names = sys.argv[1], sys.argv[2], sys.argv[3:]; kind, where = source.split(":"); rows = []; extra = []
def poa_models(name):
    T = H / "pilot" / where / name / "tables"; out = []
    E = pd.read_csv(T / "edited_gene_calls.tsv", sep="\t"); E = E[E.organelle == "mito"]
    X = pd.read_csv(T / "reconstructed_exons.tsv", sep="\t"); X = X[X.organelle == "mito"]
    G = pd.read_csv(T / "reconstructed_genes.tsv", sep="\t"); G = G[G.organelle == "mito"]
    for g in G.gene: out.append((norm(g), [tuple(sorted((int(x.start), int(x.end)))) for x in X[X.gene == g].sort_values("slot").itertuples()]))
    done = {g for g, _ in out}
    for r in E.itertuples():
        if norm(r.gene) not in done: out.append((norm(r.gene), [tuple(sorted((int(r.corrected_start), int(r.corrected_end))))]))
    return out
def pmga_models(name):
    n = name.replace("_refseq", ""); gff = R / "build/pmga" / where / n / f"{n}.gff"; mrna, cds = {}, {}
    for line in open(gff):
        if line.startswith("##FASTA"): break
        f = line.rstrip("\n").split("\t")
        if len(f) != 9: continue
        a = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x)
        if f[2] == "mRNA": mrna[a["ID"]] = a.get("Name", "")
        elif f[2] == "CDS": cds.setdefault(a["Parent"], []).append((int(f[3]) - 1, int(f[4])))
    return [(norm(mrna[m]), ex) for m, ex in cds.items() if m in mrna]
for name in names:
    try: models = poa_models(name) if kind == "poa" else pmga_models(name)
    except FileNotFoundError: print("skip", name); continue
    GB = GB0 if (GB0 / f"{ACC[name]}.gb").exists() else H / "extra_genomes" / "gb"; rec = SeqIO.read(GB / f"{ACC[name]}.gb", "genbank"); used = set()
    for f in rec.features:
        if f.type != "CDS" or "gene" not in f.qualifiers or "translation" not in f.qualifiers: continue
        g = norm(f.qualifiers["gene"][0])
        if g not in MITO or "pseudo" in f.qualifiers: continue
        parts = list(f.location.parts)
        if "".join(str(p.extract(rec.seq)) for p in parts) != str(f.extract(rec.seq)): parts = parts[::-1]
        ex = [(int(p.start), int(p.end)) for p in parts]; k = len(ex); best = None
        # a split ccmFn's pieces may be named ccmFn1/ccmFn2 or both plain ccmFn
        cand = [(i, m) for i, m in enumerate(models) if m[0] == g or (g in ("ccmfn1", "ccmfn2", "ccmfn") and m[0] in ("ccmfn", "ccmfn1", "ccmfn2"))]
        for i, (_, px) in cand:
            run = []
            for q in sorted(px):
                if run and q[0] - run[-1][1] <= 1: run[-1] = (run[-1][0], q[1], run[-1][2] + 1)
                else: run.append((q[0], q[1], 1))
            pxx = px + [(a, b) for a, b, n in run if n > 1]; placed = inner = exact = 0
            for j, (a, b) in enumerate(ex):
                ta, tb = (60 if j == 0 else 3), (6 if j == k - 1 else 3)
                placed += any(abs(p[0] - a) <= 3 and abs(p[1] - b) <= 3 for p in pxx)
                inner += any(abs(p[0] - a) <= ta and abs(p[1] - b) <= tb for p in pxx) or any(abs(p[0] - a) <= tb and abs(p[1] - b) <= ta for p in pxx)
                # to the base: the stop codon may be in or out of the last exon (either end, by strand)
                exact += any((p[0] == a or (j in (0, k - 1) and abs(p[0] - a) == 3)) and (p[1] == b or (j in (0, k - 1) and abs(p[1] - b) == 3)) for p in pxx)
            if best is None or (inner, placed, exact) > best[:3]: best = (inner, placed, exact, i, len(px))
        if best is None: rows.append(dict(genome=name, gene=g, n_exons=k, called=False, placed=0, all=False, inner=False, exact=False, pred_exons=0)); continue
        used.add(best[3])
        rows.append(dict(genome=name, gene=g, n_exons=k, called=True, placed=best[1], all=best[1] == k, inner=best[0] == k, exact=best[2] == k, pred_exons=best[4]))
    extra += [(name, m[0]) for i, m in enumerate(models) if i not in used and m[0] in MITO]
D = pd.DataFrame(rows); D.to_csv(H / "pilot" / f"{label}.models.tsv", sep="\t", index=False)
s = D[D.n_exons > 1]; o = D[D.n_exons == 1]
print(f"{label}: {D.genome.nunique()} genomes, {len(D)} annotated CDS, {D.n_exons.sum()} exons | exons placed within 3 bp {D.placed.sum()} | CDS with every exon placed {int(D['all'].sum())} "
      f"(spliced {int(s['all'].sum())}/{len(s)}, single-exon {int(o['all'].sum())}/{len(o)}) | allowing start/stop choice {int(D.inner.sum())} | exact to the base {int(D.exact.sum())} | not called {int((~D.called).sum())} | "
      f"extra models {len(extra)}")
bad = D[~D["all"]]; print("  not fully placed:", bad.gene.value_counts().head(16).to_dict())
print("  extra models:", pd.Series([g for _, g in extra]).value_counts().head(14).to_dict() if extra else {})
