"""starts.py TAG : starts of orfedit's single-exon calls (TAG.calls.tsv, run.sh) against the GenBank annotation of the test genomes,
and what start rules would do: per annotated CDS (ends/score.py's match), the codon each call starts at, and where in-frame start
codons lie around it. Writes TAG.starts.tsv; prints how many starts each rule gets within 3 bp, by organelle and gene."""
import sys, warnings
from pathlib import Path
import pandas as pd
warnings.simplefilter("ignore")
from Bio import SeqIO
from Bio.Seq import Seq
H = Path(__file__).resolve().parent; R = H.parents[2]
GB = [R / "src/poa/pipeline/trans_splicing/reference/genbank", H.parent / "extra_genomes" / "gb", R / "src/poa/pipeline/trans_splicing/reference/plastid/genbank"]
STOP = {"TAA", "TAG", "TGA"}
tag = sys.argv[1]
E = pd.read_csv(H / f"{tag}.ends.tsv", sep="\t"); E = E[E.called == True]
C = pd.read_csv(H / f"{tag}.calls.tsv", sep="\t")
M = pd.concat([pd.read_csv(H / f, sep="\t") for f in ("sites.manifest.tsv", "flat.manifest.tsv")])
strand = {(r.species, r.gene): r.strand for r in M.itertuples()}
G = pd.read_csv(H / "genomes.tsv", sep="\t", header=None, names=["name", "org", "fasta", "acc"]).set_index("name")
seqs = {}
def genome(name):
    if name not in seqs: seqs[name] = {r.id: str(r.seq).upper() for r in SeqIO.parse(G.loc[name, "fasta"], "fasta")}
    return seqs[name]
SYN = {"tatc": "mttb", "mat-r": "matr", "orfx": "mttb", "orf25": "atp4", "orfb": "atp8", "ccb206": "ccmb", "ccb256": "ccmc", "ccb452": "ccmfc", "ccb382": "ccmfn1", "ccb203": "ccmfn2", "cox-2": "cox2", "coxii": "cox2", "yejv": "ccmfc", "yejr": "ccmfn"}
norm = lambda g: SYN.get(g.lower(), g.lower())
C["g"] = C.gene.str.split("|").str[0].map(norm)
rows = []
for e in E.itertuples():
    r = C[(C.species == e.genome) & (C.gene == e.call)].iloc[0]
    st = strand[(r.species, r.gene)]
    cs = int(r.corrected_start) if st == "+" else int(r.corrected_end)
    ans = int(e.a) if e.ann_strand == 1 else int(e.b)
    contig = r.gene.split("|")[1].split(":")[0]; s = genome(r.species)[contig]; n = len(s)
    def cod(p):  # codon starting at sense position p (genomic start for +, end for -)
        if st == "+": return s[p:p + 3] if 0 <= p <= n - 3 else "NNN"
        return str(Seq(s[p - 3:p]).reverse_complement()) if 3 <= p <= n else "NNN"
    step = 3 if st == "+" else -3
    up = []; p = cs - step
    while cod(p) not in STOP and cod(p) != "NNN" and len(up) < 300: up.append(p); p -= step
    down = [cs + step * i for i in range(1, 41)]
    first_stop_down = next((i for i, q in enumerate(down) if cod(q) in STOP), None)
    if first_stop_down is not None: down = down[:first_stop_down]
    aa = r.protein.lstrip("-")[:1] if isinstance(r.protein, str) else ""
    lead = len(r.protein) - len(r.protein.lstrip("-")) if isinstance(r.protein, str) else 0
    rows.append(dict(genome=e.genome, org=e.org, gene=e.gene, strand=st, call_start=cs, ann_start=ans, ds=int((cs - ans) if st == '+' else (ans - cs)), ok=abs(e.ds) <= 3,
                     cur=cod(cs), ann=cod(ans), aa=aa, lead_unaligned=lead,
                     up_atg=[abs(cs - q) // 3 for q in up if cod(q) == "ATG"], up_acg=[abs(cs - q) // 3 for q in up if cod(q) == "ACG"],
                     down_atg=[abs(q - cs) // 3 for q in down if cod(q) == "ATG"], down_acg=[abs(q - cs) // 3 for q in down if cod(q) == "ACG"],
                     open_up=len(up)))
D = pd.DataFrame(rows); D.to_csv(H / f"{tag}.starts.tsv", sep="\t", index=False)
def rule(r, kind):
    """codons to move: negative = upstream"""
    starts = ("ATG", "ACG") if kind != "atg_only" else ("ATG",)
    if kind == "most_upstream_atg":
        return -max(r.up_atg) if r.up_atg else (0 if r.cur in ("ATG", "ACG") else (r.down_atg[0] if r.down_atg else 0))
    if r.cur in starts: return 0
    ups = r.up_atg + (r.up_acg if "ACG" in starts else []); dns = r.down_atg + (r.down_acg if "ACG" in starts else [])
    c = [(u, -u) for u in ups if u <= 20] + [(d, d) for d in dns if d <= 20]
    if kind == "nearest_atg_first":
        ca = [(u, -u) for u in r.up_atg if u <= 20] + [(d, d) for d in r.down_atg if d <= 20]
        if ca: return min(ca)[1]
    return min(c)[1] if c else 0
for kind in ("none", "nearest", "nearest_atg_first", "atg_only", "most_upstream_atg"):
    ok = D.apply(lambda r: r.ok if kind == "none" else abs(r.ds + 3 * rule(r, kind)) <= 3, axis=1)
    D["ok_" + kind] = ok
print(D.groupby("org")[[c for c in D.columns if c.startswith("ok")]].sum().to_string())
print("current start codon:", D.groupby("org").cur.apply(lambda x: x.map(lambda c: c if c in ("ATG", "ACG", "GTG", "TTG") else "other").value_counts().to_dict()).to_dict())
print("annotated start codon:", D.groupby("org").ann.apply(lambda x: x.map(lambda c: c if c in ("ATG", "ACG", "GTG", "TTG") else "other").value_counts().to_dict()).to_dict())
bad = D[~D.ok]
print(bad.groupby(["org", "gene"]).agg(n=("ok", "size"), ds=("ds", lambda x: sorted(set(int(v // 3) for v in x))[:8]), cur=("cur", lambda x: sorted(set(x))), ann=("ann", lambda x: sorted(set(x)))).to_string())
