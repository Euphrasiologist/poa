"""assess_pltd.py TAG SPECIES... : plastid calls of pilot/TAG/SPECIES against the species' RNA-seq plastid pileup.
Per gene (single-exon corrections and reconstructed spliced genes): profile columns aligned, frameshift calls, RNA
coverage, RNA edits found / calls right, and whether the RNA-edited call translates cleanly. Writes TAG.pltd.genes.tsv."""
import sys, glob
from pathlib import Path
import pandas as pd
H = Path(__file__).resolve().parent.parent; sys.path.insert(0, str(H)); import rna
D = Path("/lustre/scratch122/tol/teams/blaxter/users/mb39/plant_organellar_database/data")
tag, species = sys.argv[1], sys.argv[2:]
def read_fasta(path):
    out, name = {}, None
    for line in open(path):
        if line.startswith(">"): name = line[1:].split()[0]; out[name] = []
        else: out[name].append(line.strip().upper())
    return {k: "".join(v) for k, v in out.items()}
rows = []
for sp in species:
    T = H / "pilot" / tag / sp / "tables"
    # every run mapped to the species: the ToL library and public ones (map4.sh), base counts summed
    pfs = sorted(glob.glob(str(H / f"tol.{sp}.pltd.pile.tsv.gz")) + glob.glob(str(H / f"pub.{sp}.*[0-9].pltd.pile.tsv.gz")))
    if not pfs or not T.exists(): print("skip", sp); continue
    pile = rna.load_pile(pfs[0])
    for pf in pfs[1:]:
        for k, (ref, n, cnt) in rna.load_pile(pf).items():
            if k in pile:
                for b in cnt: pile[k][2][b] += cnt[b]
                pile[k] = (ref, pile[k][1] + n, pile[k][2])
            else: pile[k] = (ref, n, cnt)
    genome = read_fasta(sorted(glob.glob(str(D / "plastid" / sp / "*.pltd.ctg.fasta")))[0])
    gc = pd.read_csv(T / "gene_calls.tsv", sep="\t"); gc = gc[gc.organelle == "pltd"]
    E = pd.read_csv(T / "edited_gene_calls.tsv", sep="\t"); E = E[E.organelle == "pltd"]
    PE = pd.read_csv(T / "predicted_edits.tsv", sep="\t"); PE = PE[PE.organelle == "pltd"]
    FS = pd.read_csv(T / "frameshifts.tsv", sep="\t"); FS = FS[FS.organelle == "pltd"]
    def score(tr, called):
        tru = rna.truth(tr); pos = {k: n for k, b, n, na in tr}
        ok = {k for k in called if pos.get(k, 0) >= rna.MIN_DEPTH}
        return sum(v == "C>U" for v in tru.values()), sum(v == "U>C" for v in tru.values()), len(ok), sum(k in tru for k in ok), sum(tru.get(k) == "U>C" for k in ok)
    for r in E.itertuples():
        m = gc[(gc.gene == r.gene) & (gc.contig_id == r.contig) & (gc.start < r.hit_end + 300) & (gc.end > r.hit_start - 300)]
        if m.empty: continue
        lo, hi = sorted((int(r.corrected_start), int(r.corrected_end))); lo, hi = max(lo, 0), min(hi, len(genome[r.contig]))
        tr = rna.walk(pile, genome, r.contig, lo, hi, m.iloc[0].strand); d = rna.derive(tr)
        e = PE[(PE.gene == r.gene) & (PE.hit_start == r.hit_start)]
        cu, uc, nc, nr, nru = score(tr, {(r.contig, int(p) + 1) for p in e.genomic_pos})
        rows.append(dict(species=sp, kind="single", gene=r.gene, frac=round(r.n_columns_aligned / r.n_columns_total, 2), n_fs=int(((FS.gene == r.gene) & (FS.hit_start == r.hit_start)).sum()),
                         covered=d["covered"], rna_CU=cu, rna_UC=uc, called=nc, right=nr, right_UC=nru, rna_protein=d["status"], over_ext=d["over_extended"], start_moved=d["start"], aa=len(d["protein"])))
    X = pd.read_csv(T / "reconstructed_exons.tsv", sep="\t"); X = X[X.organelle == "pltd"]
    G = pd.read_csv(T / "reconstructed_genes.tsv", sep="\t"); G = G[G.organelle == "pltd"]
    RE = pd.read_csv(T / "reconstructed_edits.tsv", sep="\t"); RE = RE[RE.organelle == "pltd"]
    RF = pd.read_csv(T / "reconstructed_frameshifts.tsv", sep="\t"); RF = RF[RF.organelle == "pltd"] if len(RF) else RF
    for g in G.itertuples():
        xs = X[X.gene == g.gene].sort_values("slot"); tr = []
        for x in xs.itertuples():
            lo, hi = sorted((int(x.start), int(x.end))); tr += rna.walk(pile, genome, x.contig, max(lo, 0), min(hi, len(genome[x.contig])), x.strand)
        d = rna.derive(tr) if tr else dict(status="none", covered=0, over_extended=0, start=0, protein="")
        c0 = xs.contig.iloc[0] if len(xs) else ""
        cu, uc, nc, nr, nru = score(tr, {(c0, int(p) + 1) for p in RE[RE.gene == g.gene].genomic_pos})
        rows.append(dict(species=sp, kind="spliced", gene=g.gene, frac=round(g.n_filled / g.n_slots, 2), n_fs=int((RF.gene == g.gene).sum()) if len(RF) else 0, covered=d["covered"],
                         rna_CU=cu, rna_UC=uc, called=nc, right=nr, right_UC=nru, rna_protein=d["status"], over_ext=d["over_extended"], start_moved=d["start"], aa=len(d["protein"])))
Gd = pd.DataFrame(rows); Gd.to_csv(H / "pilot" / f"{tag}.pltd.genes.tsv", sep="\t", index=False)
ok = Gd[Gd.covered >= 0.9]
for kind, g in ok.groupby("kind"):
    print(f"{tag} pltd {kind}: {len(g)} of {(Gd.kind == kind).sum()} genes covered >= 90% in {g.species.nunique()} species; aligned >= 99%: {(g.frac >= 0.99).sum()}; frameshift genes: {(g.n_fs > 0).sum()}; "
          f"RNA C>U {g.rna_CU.sum()} U>C {g.rna_UC.sum()}; called {g.called.sum()} right {g.right.sum()} (U>C right {g.right_UC.sum()}) -> found {100 * g.right.sum() / max(g.rna_CU.sum() + g.rna_UC.sum(), 1):.0f}%, right {100 * g.right.sum() / max(g.called.sum(), 1):.0f}%; "
          f"RNA protein ok {(g.rna_protein == 'ok').sum()}; over-extended {(g.over_ext > 0).sum()}")
print(Gd[Gd.covered >= 0.9].groupby("rna_protein").size().to_string())
