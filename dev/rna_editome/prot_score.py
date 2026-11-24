"""prot_score.py LABEL SOURCE NAME... : each tool's protein for the annotated mitochondrial CDS of a GenBank genome whose
annotated proteins are the RNA-edited ones, against those proteins. SOURCE is poa:TAG or pmga:DIR (as score_models.py).
Per annotated CDS (canonical mito genes): the tool's protein of that gene nearest in length; identical or not; and, where
the lengths agree, the residues RNA editing changes (annotated protein != genomic translation): how many the tool has
(edits found), and how many other residues differ (wrong residues). Writes pilot/LABEL.prot.tsv."""
import sys, warnings
from pathlib import Path
import pandas as pd
warnings.simplefilter("ignore")
from Bio import SeqIO
H = Path(__file__).resolve().parent.parent; R = H.parent.parent; sys.path.insert(0, str(H / "pilot"))
GB0 = R / "src/poa/pipeline/trans_splicing/reference/genbank"
ACC = dict(l.split() for l in open(H / "extra_genomes" / "accessions.tsv"))
SYN = {"tatc": "mttb", "mat-r": "matr", "orfx": "mttb", "orf25": "atp4", "orfb": "atp8", "ccb206": "ccmb", "ccb256": "ccmc", "ccb452": "ccmfc", "ccb382": "ccmfn1", "ccb203": "ccmfn2", "cox-2": "cox2", "coxii": "cox2", "yejv": "ccmfc", "yejr": "ccmfn"}
MITO = set("nad1 nad2 nad3 nad4 nad4l nad5 nad6 nad7 nad9 cob cox1 cox2 cox3 atp1 atp4 atp6 atp8 atp9 ccmb ccmc ccmfc ccmfn ccmfn1 ccmfn2 matr mttb "
           "rpl2 rpl5 rpl10 rpl16 rps1 rps2 rps3 rps4 rps7 rps10 rps11 rps12 rps13 rps14 rps19 sdh3 sdh4".split())
norm = lambda g: SYN.get(g.lower(), g.lower())
label, source, names = sys.argv[1], sys.argv[2], sys.argv[3:]; kind, where = source.split(":"); rows = []
def poa_proteins(name):
    T = H / "pilot" / where / name / "tables"; out = []
    for f in ("reconstructed_genes", "edited_gene_calls"):
        d = pd.read_csv(T / f"{f}.tsv", sep="\t"); d = d[d.organelle == "mito"]
        out += [(norm(r.gene), str(r.protein).rstrip("*")) for r in d.itertuples() if isinstance(r.protein, str)]
    return out
def pmga_proteins(name):
    n = name.replace("_refseq", ""); out = []
    for rec in SeqIO.parse(R / "build/pmga" / where / n / f"{n}.Protein.fasta", "fasta"): out.append((norm(rec.id.split("_[")[0]), str(rec.seq).rstrip("*")))
    return out
for name in names:
    try: prots = poa_proteins(name) if kind == "poa" else pmga_proteins(name)
    except FileNotFoundError: print("skip", name); continue
    GB = GB0 if (GB0 / f"{ACC[name]}.gb").exists() else H / "extra_genomes" / "gb"; rec = SeqIO.read(GB / f"{ACC[name]}.gb", "genbank")
    for f in rec.features:
        if f.type != "CDS" or "gene" not in f.qualifiers or "translation" not in f.qualifiers or "pseudo" in f.qualifiers: continue
        g = norm(f.qualifiers["gene"][0])
        if g not in MITO: continue
        p = f.qualifiers["translation"][0]; geno = str(f.extract(rec.seq).translate(table=1)).rstrip("*")
        edited = [i for i in range(min(len(p), len(geno))) if p[i] != geno[i]] if len(p) == len(geno) else None
        cand = [q for n_, q in prots if n_ == g or (g.startswith("ccmfn") and n_.startswith("ccmfn"))]
        if not cand: rows.append(dict(genome=name, gene=g, aa=len(p), called=False, same_length=False, identical=False, edits=len(edited or []), found=0, wrong=0)); continue
        # the tool's copy nearest the annotated protein: identical first, then same length with fewest differences, then nearest in length
        q = min(cand, key=lambda q: (q != p, len(q) != len(p), sum(a != b for a, b in zip(q, p)) if len(q) == len(p) else abs(len(q) - len(p))))
        # the same length and nearly the same protein (a model shifted by a residue can have the length by chance)
        same = len(q) == len(p) and sum(x == y for x, y in zip(q, p)) >= 0.8 * len(p)
        found = sum(q[i] == p[i] for i in edited) if same and edited is not None else 0
        wrong = sum(q[i] != p[i] for i in range(len(p)) if i not in set(edited or [])) if same else 0
        rows.append(dict(genome=name, gene=g, aa=len(p), called=True, same_length=same, identical=q == p, edits=len(edited or []), found=found, wrong=wrong))
D = pd.DataFrame(rows); D.to_csv(H / "pilot" / f"{label}.prot.tsv", sep="\t", index=False)
for name, d in D.groupby("genome"):
    s = d[d.same_length]
    print(f"{label} {name}: {len(d)} annotated CDS; protein identical to the annotated (edited) one {int(d.identical.sum())}; same length {len(s)}: edited residues {int(s.edits.sum())}, found {int(s.found.sum())}, other residues wrong {int(s.wrong.sum())}; not called {int((~d.called).sum())}")
s = D[D.same_length]; print(f"{label} TOTAL: {len(D)} CDS; identical {int(D.identical.sum())}; same length {len(s)}: edited residues {int(s.edits.sum())}, found {int(s.found.sum())}, wrong {int(s.wrong.sum())}; not called {int((~D.called).sum())}")
