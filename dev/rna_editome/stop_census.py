"""stop_census.py TAG... : in poa's single-exon calls (pilot/TAG/<species>/tables), literal in-frame stop codons inside the corrected span
(an insertion that hops a stop), and where frame changes sit relative to the profile's ends. Writes pilot/stop_census.tsv."""
import sys, glob
from pathlib import Path
import pandas as pd
H = Path(__file__).resolve().parent.parent
D = Path("/lustre/scratch122/tol/teams/blaxter/users/mb39/plant_organellar_database/data")
COMP = str.maketrans("ACGTN", "TGCAN"); STOP = {"TAA", "TAG", "TGA"}
def fasta(path):
    out, name = {}, None
    for line in open(path):
        if line.startswith(">"): name = line[1:].split()[0]; out[name] = []
        else: out[name].append(line.strip().upper())
    return {k: "".join(v) for k, v in out.items()}
rows = []
for tag in sys.argv[1:]:
    for t in sorted(glob.glob(str(H / "pilot" / tag / "*" / "tables"))):
        sp = t.split("/")[-2]
        try:
            E = pd.read_csv(f"{t}/edited_gene_calls.tsv", sep="\t"); G = pd.read_csv(f"{t}/gene_calls.tsv", sep="\t"); L = pd.read_csv(f"{t}/lineage.tsv", sep="\t").lineage.iloc[0]
            F = pd.read_csv(f"{t}/frameshifts.tsv", sep="\t")
        except Exception as e: print("skip", sp, e); continue
        strand = {(r.organelle, r.gene, r.contig_id, r.start, r.end): r.strand for r in G.itertuples()}
        seqs = {}
        for org, sub, ext in (("mito", "mito", "mito"), ("pltd", "plastid", "pltd")):
            fs = sorted(glob.glob(str(D / sub / sp / f"*.{ext}.ctg.fasta")))
            if fs: seqs[org] = fasta(fs[0])
        for r in E.itertuples():
            st = strand.get((r.organelle, r.gene, r.contig, r.hit_start, r.hit_end)); s = seqs.get(r.organelle, {}).get(r.contig)
            if st is None or s is None: continue
            a, b = int(r.corrected_start), int(r.corrected_end); cds = s[a:b] if st == "+" else s[a:b][::-1].translate(COMP)
            cod = [cds[i:i + 3] for i in range(0, len(cds) - 2, 3)]
            inner = [i for i, c in enumerate(cod[:-1]) if c in STOP]
            f = F[(F.gene == r.gene) & (F.organelle == r.organelle) & (F.hit_start == r.hit_start)] if len(F) else F
            cols = sorted(int(x) for x in f.profile_column) if len(f) else []
            rows.append(dict(tag=tag, species=sp, lineage=L, organelle=r.organelle, gene=r.gene, aligned=r.n_columns_aligned, total=r.n_columns_total, n_fs=int(r.n_frameshifts), bp=b - a,
                             inner_stops=len(inner) if int(r.n_frameshifts) == 0 else -1, first_stop_codon=inner[0] if inner and int(r.n_frameshifts) == 0 else -1, n_codons=len(cod),
                             fs_cols=",".join(map(str, cols)), fs_from_end=min((r.n_columns_total - c for c in cols), default=-1), fs_from_start=min(cols, default=-1), stop_in=r.includes_stop))
T = pd.DataFrame(rows); T.to_csv(H / "pilot" / "stop_census.tsv", sep="\t", index=False)
N = T[T.n_fs == 0]; S = N[N.inner_stops > 0]
print(f"{len(T)} calls; without frame changes {len(N)}; with a literal stop inside the span {len(S)}")
print("by lineage:", S.groupby("lineage").size().to_dict(), "of", N.groupby("lineage").size().to_dict())
print(S.groupby(["organelle", "gene"]).agg(n=("gene", "size"), stops=("inner_stops", "median"), codons_after_first=("first_stop_codon", lambda x: 0)).sort_values("n", ascending=False).head(40).to_string())
S = S.assign(tail=S.n_codons - S.first_stop_codon)
print("codons from the first inner stop to the span's end:", pd.cut(S.tail, [0, 5, 10, 20, 50, 100, 1e5]).value_counts().sort_index().to_dict())
print("span length over 3 x aligned columns (an intron bridged?):", pd.cut((S.bp / 3 - S.aligned), [-1e5, 5, 20, 50, 100, 1e5]).value_counts().sort_index().to_dict())
X = T[T.n_fs > 0]
print(f"calls with frame changes {len(X)}; nearest to the profile's end:", pd.cut(X.fs_from_end, [-1, 5, 10, 15, 20, 30, 1e5]).value_counts().sort_index().to_dict(), "| to its start:", pd.cut(X.fs_from_start, [-1, 5, 10, 15, 20, 30, 1e5]).value_counts().sort_index().to_dict())
