"""assess.py LINEAGE SPECIES... : a pilot poa run (pilot/LINEAGE/SPECIES) against the species' RNA-seq.
Per mito gene: profile columns aligned, frameshift calls, RNA coverage, RNA edits found / calls right, and whether
the RNA-edited call translates cleanly (rna.derive). Then every intron the RNA shows (split reads, >= MINJ reads)
against the predicted exon junctions. Writes pilot/LINEAGE.genes.tsv and pilot/LINEAGE.introns.tsv."""
import sys, re, subprocess, glob
from collections import Counter
from pathlib import Path
import pandas as pd
H = Path(__file__).resolve().parent.parent; sys.path.insert(0, str(H)); import rna, boot
D = Path("/lustre/scratch122/tol/teams/blaxter/users/mb39/plant_organellar_database/data")
lin, species = sys.argv[1], sys.argv[2:]; MINJ = 5
CIG = re.compile(r"(\d+)([MIDNSHP=X])")
def span(pos, cigar):
    ops = [(int(n), o) for n, o in CIG.findall(cigar)]
    ref = sum(n for n, o in ops if o in "MDN=X"); q = sum(n for n, o in ops if o in "MI=X")
    left = ops[0][0] if ops[0][1] in "SH" else 0
    return pos, pos + ref - 1, left, left + q
def read_fasta(path):
    out, name = {}, None
    for line in open(path):
        if line.startswith(">"): name = line[1:].split()[0]; out[name] = []
        else: out[name].append(line.strip().upper())
    return {k: "".join(v) for k, v in out.items()}
def introns(bams):
    obs = Counter()
    for line in (l for bam in bams for l in subprocess.Popen(["samtools", "view", "-F", "0x900", bam], stdout=subprocess.PIPE, text=True).stdout):
        if "\tSA:Z:" not in line: continue
        f = line.split("\t", 12); flag = int(f[1]); contig = f[2]; rev = bool(flag & 16)
        s1, e1, qs1, qe1 = span(int(f[3]), f[5])
        sa = line.split("\tSA:Z:", 1)[1].split("\t", 1)[0].rstrip("\n")
        for part in sa.rstrip(";").split(";"):
            c2, pos2, st2, cig2 = part.split(",")[:4]
            if c2 != contig or (st2 == "-") != rev: continue
            s2, e2, qs2, qe2 = span(int(pos2), cig2)
            (ls, le, lqs, lqe), (rs, re_, rqs, rqe) = sorted([(s1, e1, qs1, qe1), (s2, e2, qs2, qe2)])
            if lqs < rqs and 20 < rs - le < 20000 and abs(rqs - lqe) <= 15: obs[(contig.replace("mito_", ""), le, rs)] += 1
    # merge junctions within 6 bp at both ends (edited bases clipped, microhomology)
    out = []
    for (c, le, rs), n in sorted(obs.items(), key=lambda kv: -kv[1]):
        for o in out:
            if o[0] == c and abs(o[1] - le) <= 6 and abs(o[2] - rs) <= 6: o[3] += n; break
        else: out.append([c, le, rs, n])
    return [o for o in out if o[3] >= MINJ]
grows, irows = [], []
for sp in species:
    T = H / "pilot" / lin / sp / "tables"
    pile = boot.load_piles(sp)
    if not pile or not T.exists(): print("skip", sp); continue
    genome = read_fasta(boot.genome_fasta(sp))
    gc = pd.read_csv(T / "gene_calls.tsv", sep="\t"); gc = gc[gc.organelle == "mito"]
    E = pd.read_csv(T / "edited_gene_calls.tsv", sep="\t"); E = E[E.organelle == "mito"]
    PE = pd.read_csv(T / "predicted_edits.tsv", sep="\t"); PE = PE[PE.organelle == "mito"]
    FS = pd.read_csv(T / "frameshifts.tsv", sep="\t"); FS = FS[FS.organelle == "mito"]
    def score(tr, called):
        tru = rna.truth(tr); pos = {k: (b, n) for k, b, n, na in tr}
        ok = {k for k in called if k in pos and pos[k][1] >= rna.MIN_DEPTH}
        return sum(v == "C>U" for v in tru.values()), sum(v == "U>C" for v in tru.values()), len(ok), sum(k in tru for k in ok)
    for r in E.itertuples():
        m = gc[(gc.gene == r.gene) & (gc.contig_id == r.contig) & (gc.start < r.hit_end + 300) & (gc.end > r.hit_start - 300)]
        if m.empty: continue
        strand = m.iloc[0].strand; lo, hi = sorted((int(r.corrected_start), int(r.corrected_end))); lo, hi = max(lo, 0), min(hi, len(genome[r.contig]))
        tr = rna.walk(pile, genome, r.contig, lo, hi, strand); d = rna.derive(tr)
        e = PE[(PE.gene == r.gene) & (PE.hit_start == r.hit_start)]
        cu, uc, nc, nr = score(tr, {(r.contig, int(p) + 1) for p in e.genomic_pos})
        grows.append(dict(species=sp, kind="single", gene=r.gene, aligned=f"{r.n_columns_aligned}/{r.n_columns_total}", frac=round(r.n_columns_aligned / r.n_columns_total, 2),
                          n_fs=int(((FS.gene == r.gene) & (FS.hit_start == r.hit_start)).sum()), covered=d["covered"], rna_CU=cu, rna_UC=uc, called=nc, right=nr, rna_protein=d["status"], over_ext=d["over_extended"], start_moved=d["start"]))
    X = pd.read_csv(T / "reconstructed_exons.tsv", sep="\t"); X = X[X.organelle == "mito"]
    G = pd.read_csv(T / "reconstructed_genes.tsv", sep="\t"); G = G[G.organelle == "mito"]
    RE = pd.read_csv(T / "reconstructed_edits.tsv", sep="\t"); RF = pd.read_csv(T / "reconstructed_frameshifts.tsv", sep="\t")
    pred = []
    for g in G.itertuples():
        xs = X[X.gene == g.gene].sort_values("slot"); tr = []
        for x in xs.itertuples():
            lo, hi = sorted((int(x.start), int(x.end))); lo, hi = max(lo, 0), min(hi, len(genome[x.contig])); tr += rna.walk(pile, genome, x.contig, lo, hi, x.strand)
        xl = xs.to_dict("records")
        for a, b in zip(xl, xl[1:]):
            if a["contig"] == b["contig"] and a["strand"] == b["strand"]:
                lo_a, hi_a = sorted((a["start"], a["end"])); lo_b, hi_b = sorted((b["start"], b["end"]))
                pred.append((g.gene, f"{a['slot']}-{b['slot']}", a["contig"], *((hi_a, lo_b + 1) if a["strand"] == "+" else (hi_b, lo_a + 1))))
        d = rna.derive(tr) if tr else dict(status="none", covered=0, over_extended=0, start=0)
        e = RE[(RE.gene == g.gene) & (RE.organelle == "mito")] if "organelle" in RE else RE[RE.gene == g.gene]
        c0 = xs.contig.iloc[0] if len(xs) else ""
        cu, uc, nc, nr = score(tr, {(c0, int(p) + 1) for p in e.genomic_pos})
        nfs = int(((RF.gene == g.gene) & (RF.organelle == "mito")).sum()) if len(RF) else 0
        grows.append(dict(species=sp, kind="spliced", gene=g.gene, aligned=f"{g.n_filled}/{g.n_slots} exons", frac=round(g.n_filled / g.n_slots, 2), n_fs=nfs,
                          covered=d["covered"], rna_CU=cu, rna_UC=uc, called=nc, right=nr, rna_protein=d["status"], over_ext=d["over_extended"], start_moved=d["start"]))
    # introns the RNA shows, against predicted junctions and gene calls
    spans = [(r.gene, r.contig, *sorted((int(r.corrected_start), int(r.corrected_end)))) for r in E.itertuples()] + [(x.gene + f".e{x.slot}", x.contig, *sorted((int(x.start), int(x.end)))) for x in X.itertuples()]
    used = set()
    for c, le, rs, n in introns(boot.bams(sp)):
        hit = [p for p in pred if p[2] == c and abs(p[3] - le) <= 60 and abs(p[4] - rs) <= 60]
        exact = [p for p in hit if abs(p[3] - le) <= 6 and abs(p[4] - rs) <= 6]
        left = [s[0] for s in spans if s[1] == c and s[2] - 50 <= le <= s[3] + 50]; right = [s[0] for s in spans if s[1] == c and s[2] - 50 <= rs <= s[3] + 50]
        for p in hit: used.add(p[:2])
        irows.append(dict(species=sp, contig=c, last_before=le, first_after=rs, intron_bp=rs - le - 1, reads=n,
                          predicted=(f"{exact[0][0]} {exact[0][1]}" if exact else (f"near {hit[0][0]} {hit[0][1]} ({le - hit[0][3]},{rs - hit[0][4]})" if hit else "")),
                          left_in=",".join(sorted(set(left))), right_in=",".join(sorted(set(right)))))
    for p in pred:
        if p[:2] not in used: irows.append(dict(species=sp, contig=p[2], last_before=p[3], first_after=p[4], intron_bp=p[4] - p[3] - 1, reads=0, predicted=f"{p[0]} {p[1]} (no reads)", left_in="", right_in=""))
Gd = pd.DataFrame(grows, columns=["species", "kind", "gene", "aligned", "frac", "n_fs", "covered", "rna_CU", "rna_UC", "called", "right", "rna_protein", "over_ext", "start_moved"]); Gd.to_csv(H / "pilot" / f"{lin}.genes.tsv", sep="\t", index=False)
I = pd.DataFrame(irows); I.to_csv(H / "pilot" / f"{lin}.introns.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 500)
ok = Gd[Gd.covered >= 0.9]
for kind, g in ok.groupby("kind"):
    print(f"{lin} {kind}: {len(g)} genes covered >= 90% in {g.species.nunique()} species; aligned >= 99%: {(g.frac >= 0.99).sum()}; frameshift genes: {(g.n_fs > 0).sum()}; "
          f"RNA C>U {g.rna_CU.sum()} U>C {g.rna_UC.sum()}; called {g.called.sum()} right {g.right.sum()} -> found {100 * g.right.sum() / max(g.rna_CU.sum() + g.rna_UC.sum(), 1):.0f}%, right {100 * g.right.sum() / max(g.called.sum(), 1):.0f}%; RNA protein ok {(g.rna_protein == 'ok').sum()}")
