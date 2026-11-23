"""pltd_refs.py TAG OUT SPECIES... : RNA-edited plastid references from poa runs (pilot/TAG/SPECIES) and the species'
plastid RNA pileups, by the rules of rna.derive (as for the fern mitochondrial references, 2026-11).
Single-exon genes -> OUT/rna_proteins.pltd.fasta (species|gene), one per species and gene (the best-scoring copy;
the inverted repeat's second copy is the same protein). Spliced genes transsplice reconstructed completely, without a
frameshift call -> OUT/raw_rna/RNA_<species>.fasta + junctions.tsv, exons cut to the protein's span."""
import sys, glob, csv
from pathlib import Path
import pandas as pd
H = Path(__file__).resolve().parent; sys.path.insert(0, str(H)); import rna
D = Path("/lustre/scratch122/tol/teams/blaxter/users/mb39/plant_organellar_database/data")
tag, out, species = sys.argv[1], Path(sys.argv[2]), sys.argv[3:]
comp = rna.comp
def read_fasta(path):
    o, name = {}, None
    for line in open(path):
        if line.startswith(">"): name = line[1:].split()[0]; o[name] = []
        else: o[name].append(line.strip().upper())
    return {k: "".join(v) for k, v in o.items()}
(out / "raw_rna").mkdir(parents=True, exist_ok=True)
for f in (out / "raw_rna").glob("RNA_*.fasta"): f.unlink()
jcols = ["accession", "gene", "locus_tag", "exon_num", "n_exons", "start", "end", "strand", "junction_to_next", "gap_bp", "notes"]
jrows, prots, stat = [], [], []
for sp in species:
    T = H / "pilot" / tag / sp / "tables"; pf = H / f"tol.{sp}.pltd.pile.tsv.gz"
    if not pf.exists() or not T.exists(): print("skip", sp); continue
    pile = rna.load_pile(pf); genome = read_fasta(sorted(glob.glob(str(D / "plastid" / sp / "*.pltd.ctg.fasta")))[0])
    gc = pd.read_csv(T / "gene_calls.tsv", sep="\t"); gc = gc[gc.organelle == "pltd"]
    E = pd.read_csv(T / "edited_gene_calls.tsv", sep="\t"); E = E[E.organelle == "pltd"].sort_values("score", ascending=False)
    done = set()
    for r in E.itertuples():
        if r.gene in done: continue
        m = gc[(gc.gene == r.gene) & (gc.contig_id == r.contig) & (gc.start < r.hit_end + 300) & (gc.end > r.hit_start - 300)]
        if m.empty: continue
        lo, hi = sorted((int(r.corrected_start), int(r.corrected_end))); lo, hi = max(lo, 0), min(hi, len(genome[r.contig]))
        d = rna.derive(rna.walk(pile, genome, r.contig, lo, hi, m.iloc[0].strand))
        why = d["status"] if d["status"] != "ok" else ("length not a multiple of 3" if d["mod3"] else "ok")
        stat.append((sp, r.gene, "single", why, d["start"], d["over_extended"]))
        if why == "ok": prots.append((sp, r.gene, d["protein"])); done.add(r.gene)
    X = pd.read_csv(T / "reconstructed_exons.tsv", sep="\t"); X = X[X.organelle == "pltd"]
    G = pd.read_csv(T / "reconstructed_genes.tsv", sep="\t"); G = G[G.organelle == "pltd"]
    RF = pd.read_csv(T / "reconstructed_frameshifts.tsv", sep="\t"); fs = set(RF[RF.organelle == "pltd"].gene) if len(RF) else set()
    with open(out / "raw_rna" / f"RNA_{sp}.fasta", "w") as fh:
        for g in G.itertuples():
            xs = X[X.gene == g.gene].sort_values("slot"); why, d = "", None
            if str(g.complete).lower() != "true": why = "incomplete"
            elif g.gene in fs: why = "frameshift call"
            if not why:
                tr, owner = [], []
                for x in xs.itertuples():
                    lo, hi = sorted((int(x.start), int(x.end))); w = rna.walk(pile, genome, x.contig, max(lo, 0), min(hi, len(genome[x.contig])), x.strand)
                    tr += w; owner += [x.slot] * len(w)
                d = rna.derive(tr); why = "" if d["status"] == "ok" else d["status"]
                if not why and d["mod3"]: why = "length not a multiple of 3"
            if not why:
                first, last = 3 * d["start"], 3 * d["end"] + (3 if d["has_stop"] else 0); exons = []
                for x in xs.itertuples():
                    idx = [i for i in range(first, min(last, len(tr))) if owner[i] == x.slot]
                    if not idx: continue
                    ps = [tr[i][0][1] for i in idx]; lo, hi = min(ps) - 1, max(ps); raw = genome[x.contig][lo:hi]
                    if x.strand == "-": raw = "".join(comp[b] for b in reversed(raw))
                    exons.append((x, lo, hi, raw))
                if len(exons) != len(xs): why = "exon lost in trimming"
            stat.append((sp, g.gene, "spliced", why or "ok", d["start"] if d else 0, d["over_extended"] if d else 0))
            if why: continue
            acc, k = f"RNA_{sp}", len(exons)
            for i, (x, lo, hi, raw) in enumerate(exons, 1):
                fh.write(f">{acc}|{g.gene}|{g.gene}_1|exon{i}of{k}|{lo + 1}-{hi}|{x.strand}\n{raw}\n")
                nx = exons[i][0] if i < k else None
                jt = "" if nx is None else ("cis" if nx.contig == x.contig and nx.strand == x.strand else "trans")
                jrows.append(dict(accession=acc, gene=g.gene, locus_tag=f"{g.gene}_1", exon_num=i, n_exons=k, start=lo + 1, end=hi, strand=x.strand, junction_to_next=jt, gap_bp="", notes=f"{sp}: RNA-seq edited"))
            fh.write(f">{acc}|{g.gene}|{g.gene}_1|translation_REFSEQ_ANNOTATED|len={len(d['protein'])}\n{d['protein']}\n")
with open(out / "raw_rna" / "junctions.tsv", "w") as fh:
    w = csv.DictWriter(fh, jcols, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(jrows)
with open(out / "rna_proteins.pltd.fasta", "w") as fh:
    for sp, gene, p in prots: fh.write(f">{sp}|{gene}\n{p}\n")
S = pd.DataFrame(stat, columns=["species", "gene", "kind", "status", "start_trim", "over_extended"]); S.to_csv(out / "loci.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
print(S.groupby(["kind", "status"]).size().to_string()); print("single-exon genes:", S[(S.kind == "single") & (S.status == "ok")].gene.nunique(), "| proteins:", len(prots), "| spliced loci:", (S[S.kind == "spliced"].status == "ok").sum())
print(S[(S.status == "ok")].groupby("kind").agg(start_trimmed=("start_trim", lambda x: (x > 0).sum()), over_extended=("over_extended", lambda x: (x > 0).sum())).to_string())
print(S[S.kind == "spliced"].pivot_table(index="gene", columns="status", values="species", aggfunc="count", fill_value=0).to_string())
