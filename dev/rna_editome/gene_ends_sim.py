"""gene_ends_sim.py TAG NAME... : projected gene ends if (1) a gene ended at the first literal in-frame stop after its start unless it now ends
earlier at an edit-made stop, and (2) its start were moved to a start codon by one rule per gene. Also how consistent RefSeq's starts are."""
import sys, warnings
from pathlib import Path
import pandas as pd
warnings.simplefilter("ignore")
from Bio import SeqIO
R = Path("/lustre/scratch122/tol/teams/blaxter/users/mb39/plant_organelle_annotator"); H = R / "build/rnaedit"
GB0 = R / "src/poa/pipeline/trans_splicing/reference/genbank"
ACC = dict(l.split() for l in open(H / "extra_genomes" / "accessions.tsv"))
SYN = {"tatc": "mttb", "mat-r": "matr", "orfx": "mttb", "orf25": "atp4", "orfb": "atp8", "ccb206": "ccmb", "ccb256": "ccmc", "ccb452": "ccmfc", "ccb382": "ccmfn1", "ccb203": "ccmfn2", "cox-2": "cox2", "coxii": "cox2", "yejv": "ccmfc", "yejr": "ccmfn"}
MITO = set("nad1 nad2 nad3 nad4 nad4l nad5 nad6 nad7 nad9 cob cox1 cox2 cox3 atp1 atp4 atp6 atp8 atp9 ccmb ccmc ccmfc ccmfn ccmfn1 ccmfn2 matr mttb rpl2 rpl5 rpl10 rpl16 rps1 rps2 rps3 rps4 rps7 rps10 rps11 rps12 rps13 rps14 rps19 sdh3 sdh4".split())
norm = lambda g: SYN.get(g.lower(), g.lower()); STOP = ("TAA", "TAG", "TGA")
tag, names = sys.argv[1], sys.argv[2:]; rows = []; P10 = set(open(H / "pilot/ang_gb.list").read().split())
for name in names:
    T = H / "pilot" / tag / name / "tables"
    E = pd.read_csv(T / "edited_gene_calls.tsv", sep="\t"); E = E[E.organelle == "mito"]
    GB = GB0 if (GB0 / f"{ACC[name]}.gb").exists() else H / "extra_genomes" / "gb"; rec = SeqIO.read(GB / f"{ACC[name]}.gb", "genbank"); seq = rec.seq; n = len(seq)
    for f in rec.features:
        if f.type != "CDS" or "gene" not in f.qualifiers or "translation" not in f.qualifiers: continue
        g = norm(f.qualifiers["gene"][0])
        if g not in MITO or "pseudo" in f.qualifiers or len(f.location.parts) != 1: continue
        a, b, st = int(f.location.start), int(f.location.end), f.location.strand
        cand = E[(E.gene.str.lower().map(norm) == g) | ((g in ("ccmfn1", "ccmfn2", "ccmfn")) & E.gene.str.lower().isin(["ccmfn", "ccmfn1", "ccmfn2"]))]
        best = None
        for r in cand.itertuples():
            pa, pb = sorted((int(r.corrected_start), int(r.corrected_end))); ov = min(b, pb) - max(a, pa)
            if ov > 0 and (best is None or ov > best[0]): best = (ov, pa, pb, r)
        if best is None: rows.append(dict(genome=name, gene=g, called=False, p10=name in P10)); continue
        _, pa, pb, r = best
        def cod(p):
            if st == 1: return str(seq[p:p + 3]).upper() if 0 <= p <= n - 3 else "NNN"
            return str(seq[p - 3:p].reverse_complement()).upper() if 3 <= p <= n else "NNN"
        step = 3 * st; ann_s, ann_e = (a, b) if st == 1 else (b, a); poa_s, poa_e = (pa, pb) if st == 1 else (pb, pa); dist = lambda p, q: (p - q) * st
        up = []; p = poa_s - step
        while cod(p) not in STOP and cod(p) != "NNN" and len(up) < 400: up.append(p); p -= step
        down = [poa_s + step * i for i in range(1, 16)]; cur = cod(poa_s)
        def nearest(codons):
            if cur in codons: return poa_s
            c = [(abs(dist(q, poa_s)), 0 if dist(q, poa_s) < 0 else 1, q) for q in up + down if cod(q) in codons]
            return min(c)[2] if c else None
        A = nearest(("ATG",)); A = A if A is not None else (nearest(("ACG",)) or poa_s)
        c = [q for q in up if cod(q) == "ATG"]; B = c[-1] if c else (poa_s if cur == "ATG" else A)
        okS = lambda q: abs(dist(q, ann_s)) <= 3
        # stop: first literal stop in frame from poa's start; kept at poa's present end if that is earlier and ends on a stop made by an edit
        q = poa_s; k = 0
        while cod(q) not in STOP and cod(q) != "NNN" and k < 1500: q += step; k += 1
        orf_end = q + step
        last = cod(poa_e - step)
        edited_stop = bool(r.includes_stop) and last in ("CAA", "CAG", "CGA") and dist(poa_e, orf_end) < 0 and dist(poa_e, poa_s) % 3 == 0
        new_e = poa_e if edited_stop else orf_end
        nd = lambda e: min((dist(e, ann_e), dist(e, ann_e) + 3, dist(e, ann_e) - 3), key=abs)
        rows.append(dict(genome=name, gene=g, called=True, p10=name in P10, s0=okS(poa_s), sA=okS(A), sB=okS(B), e0=abs(nd(poa_e)) <= 3, e1=abs(nd(new_e)) <= 3,
                         ds=round(dist(ann_s, poa_s) / 3), ann_start=cod(ann_s), cur=cur))
D = pd.DataFrame(rows); D.to_csv(H / "pilot" / f"{tag}.sim2.tsv", sep="\t", index=False)
for lab, X in (("19 genomes", D), ("the ten of the PMGA comparison", D[D.p10])):
    C = X[X.called == True].astype({c: bool for c in ("s0", "sA", "sB", "e0", "e1")})
    use_b = {g: x.sB.sum() > x.sA.sum() for g, x in C.groupby("gene")}; sR = C.apply(lambda r: r.sB if use_b[r.gene] else r.sA, axis=1)
    print(f"{lab}: {len(X)} annotated single-exon CDS, {len(C)} called")
    print(f"  stop within 3 bp: now {C.e0.sum()} -> ends at its first stop {C.e1.sum()} (made wrong: {(C.e0 & ~C.e1).sum()})")
    print(f"  start within 3 bp: now {C.s0.sum()} -> nearest ATG {C.sA.sum()} | one rule per gene {sR.sum()} (upstream-most ATG for: {sorted(g for g, v in use_b.items() if v)})")
    print(f"  both ends: now {(C.s0 & C.e0).sum()} -> stops only {(C.s0 & C.e1).sum()} | stops + nearest ATG {(C.sA & C.e1).sum()} | stops + one start rule per gene {(sR & C.e1).sum()}")
C = D[D.called == True]
# how far RefSeq's own starts agree: per gene, the share of annotated starts at the gene's commonest position relative to poa's start
m = C.groupby("gene").ds.agg(lambda s: (s == s.mode().iloc[0]).sum()); t = C.groupby("gene").size()
print(f"annotated starts at their gene's commonest position (relative to the profile's first column): {m.sum()} of {t.sum()}")
print("  genes below 80%:", {g: f"{m[g]}/{t[g]}" for g in t.index if m[g] / t[g] < 0.8})
print("  start codon of poa's call where the annotation has ATG/ACG and poa's is neither:", int(((~C.cur.isin(["ATG", "ACG"])) & C.ann_start.isin(["ATG", "ACG"])).sum()))
