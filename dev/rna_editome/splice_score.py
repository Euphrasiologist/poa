"""splice_score.py PREFIX [LABEL] : RNA-seq truth for spliced mito genes. PREFIX.{exons,edits,genes}.tsv are transsplice
tables; piles are tol.<species>.pile.tsv.gz. Per exon: coverage and RNA edits (as rebuild.py), and how transsplice's
edit calls compare. Writes PREFIX.rna_exons.tsv and prints a per-gene and overall summary."""
import gzip, sys, glob, re
from pathlib import Path
import pandas as pd
H = Path(__file__).resolve().parent; MIN_DEPTH, MIN_FRAC = 10, 0.2
pre = sys.argv[1]; label = sys.argv[2] if len(sys.argv) > 2 else pre
comp = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
species = sorted(re.sub(r"^tol\.|\.pile\.tsv\.gz$", "", Path(p).name) for p in glob.glob(str(H / "tol.*.pile.tsv.gz")))
X = pd.read_csv(f"{pre}.exons.tsv", sep="\t"); X = X[X.species.isin(species) & (X.organelle == "mito")]
E = pd.read_csv(f"{pre}.edits.tsv", sep="\t"); G = pd.read_csv(f"{pre}.genes.tsv", sep="\t")
if "edit" not in E.columns: E["edit"] = "C>U"
rows = []
for sp, xs in X.groupby("species"):
    pile = {}
    with gzip.open(H / f"tol.{sp}.pile.tsv.gz", "rt") as fh:
        for line in fh:
            c, p, ref, n, a, cc, g, t = line.split("\t")
            pile[(c, int(p))] = (ref, int(n), {"A": int(a), "C": int(cc), "G": int(g), "T": int(t)})
    es = E[E.species == sp]
    for x in xs.itertuples():
        lo, hi = sorted((int(x.start), int(x.end))); truth = {}; ncov = 0; other = 0
        for p in range(lo + 1, hi + 1):
            ref, n, cnt = pile.get((x.contig, p), ("N", 0, {}))
            if n < MIN_DEPTH: continue
            ncov += 1
            sref = ref if x.strand == "+" else comp[ref]
            for alt in "ACGT":
                galt = alt if x.strand == "+" else comp[alt]
                if alt == sref or cnt[galt] < MIN_FRAC * n or cnt[galt] < 3: continue
                if (sref, alt) == ("C", "T"): truth[p] = "C>U"
                elif (sref, alt) == ("T", "C"): truth[p] = "U>C"
                else: other += 1
        r = {"species": sp, "gene": x.gene, "slot": x.slot, "len": hi - lo, "covered": round(ncov / max(hi - lo, 1), 2), "rna_other": other}
        e = es[(es.gene == x.gene) & (es.slot == x.slot)]
        for kind, lab in (("C>U", "CU"), ("U>C", "UC")):
            tk = {p for p, k in truth.items() if k == kind}; pk = {int(v) + 1 for v in e[e.edit == kind].genomic_pos}
            pk = {p for p in pk if pile.get((x.contig, p), ("N", 0))[1] >= MIN_DEPTH}
            r.update({f"{lab}_true": len(tk), f"{lab}_called": len(pk), f"{lab}_right": len(tk & pk)})
        rows.append(r)
T = pd.DataFrame(rows); T.to_csv(f"{pre}.rna_exons.tsv", sep="\t", index=False)
ok = T[T.covered >= 0.9]
def summ(d, by):
    s = d.groupby(by).agg(exons=("slot", "size"), bp=("len", "sum"), CU_true=("CU_true", "sum"), CU_called=("CU_called", "sum"), CU_right=("CU_right", "sum"),
                          UC_true=("UC_true", "sum"), UC_called=("UC_called", "sum"), UC_right=("UC_right", "sum"), other=("rna_other", "sum")).reset_index()
    for lab in ("CU", "UC"):
        s[f"{lab}_found%"] = (100 * s[f"{lab}_right"] / s[f"{lab}_true"].clip(lower=1)).round(0); s[f"{lab}_right%"] = (100 * s[f"{lab}_right"] / s[f"{lab}_called"].clip(lower=1)).round(0)
    return s
pd.set_option("display.width", 250)
Gs = G[G.species.isin(species) & (G.organelle == "mito")]
print(f"{label}: {len(species)} species, {len(Gs)} spliced mito genes, complete {int((Gs.complete.astype(str).str.lower() == 'true').sum())}; exons placed {len(T)}, covered >=90% at depth {MIN_DEPTH}: {len(ok)}")
print(summ(ok, "gene").to_string(index=False))
print(summ(ok.assign(all="all"), "all").to_string(index=False))
