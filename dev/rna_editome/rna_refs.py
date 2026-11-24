"""rna_refs.py ORG TAG OUT SPECIES... : RNA-edited references for one organelle (mito or pltd) from poa runs
(pilot/TAG/SPECIES) and the species' RNA pileups of that organelle (tol.*, pub.*: every run summed), by the rules of
rna.derive. pltd_refs.py made general: the angiosperm references, 2026-11.
Single-exon genes -> OUT/rna_proteins.fasta (mito) or rna_proteins.pltd.fasta (species|gene), one per species and gene
(the best-scoring copy). Spliced genes transsplice reconstructed completely, without a frameshift call ->
OUT/raw_rna/RNA_<species>.fasta + junctions.tsv, exons cut to the protein's span. A gene is used when at least 90% of
it has 10 or more reads and, with the edits the reads show applied, it translates cleanly. OUT/loci.tsv says what was.

RNA_REVERSE=0 in the environment, for a lineage without U-to-C editing (angiosperms): no T is read as C, and the sites
where the reads do show a C for the genome's T measure how far the reads are from the assembly's own plant (public
RNA-seq is often another individual's, now and then another species'). OUT/species.tsv has the counts. A species with
more than max(3, 5% of its C-to-U sites) such sites is "noisy": its proteins are used, but not as evidence of where
editing happens (single-exon proteins without lower case, spliced loci noted "sites not used"). One with more than
max(50, 50%) is left out (the reads are another species'). RNA_MIN_ALIGNED=0.95: a single-exon gene is used only if the
call it is read from aligned to at least that share of its profile, or runs from a Met to a stop codon. The proteins of the others have the residues a C-to-U edit made in lower case."""
import sys, glob, csv, os
from pathlib import Path
import pandas as pd
H = Path(__file__).resolve().parent; sys.path.insert(0, str(H)); import rna
D = Path("/lustre/scratch122/tol/teams/blaxter/users/mb39/plant_organellar_database/data")
org, tag, out, species = sys.argv[1], sys.argv[2], Path(sys.argv[3]), sys.argv[4:]
PILE = "pile" if org == "mito" else "pltd.pile"; DIR = "mito" if org == "mito" else "plastid"; EXT = "mito" if org == "mito" else "pltd"
comp = rna.comp
REV = os.environ.get("RNA_REVERSE", "1") != "0"
MIN_ALIGNED = float(os.environ.get("RNA_MIN_ALIGNED", 0))  # single-exon genes: least share of the profile the seed call aligned
def read_fasta(path):
    o, name = {}, None
    for line in open(path):
        if line.startswith(">"): name = line[1:].split()[0]; o[name] = []
        else: o[name].append(line.strip().upper())
    return {k: "".join(v) for k, v in o.items()}
(out / "raw_rna").mkdir(parents=True, exist_ok=True)
for f in (out / "raw_rna").glob("RNA_*.fasta"): f.unlink()
jcols = ["accession", "gene", "locus_tag", "exon_num", "n_exons", "start", "end", "strand", "junction_to_next", "gap_bp", "notes"]
jrows, prots, stat, noise = [], [], [], []
for sp in species:
    sp_prots, sp_jrows, sp_fasta, sp_stat, sites = [], [], [], [], {}
    T = H / "pilot" / tag / sp / "tables"
    # every run mapped to the species: the ToL library and public ones (map4.sh), base counts summed
    pfs = sorted(glob.glob(str(H / f"tol.{sp}.{PILE}.tsv.gz")) + glob.glob(str(H / f"pub.{sp}.*[0-9].{PILE}.tsv.gz")))
    if not pfs or not T.exists(): print("skip", sp); continue
    pile = rna.load_pile(pfs[0])
    for pf in pfs[1:]:
        for k, (ref, n, cnt) in rna.load_pile(pf).items():
            if k in pile:
                for b in cnt: pile[k][2][b] += cnt[b]
                pile[k] = (ref, pile[k][1] + n, pile[k][2])
            else: pile[k] = (ref, n, cnt)
    genome = read_fasta(sorted(glob.glob(str(D / DIR / sp / f"*.{EXT}.ctg.fasta")))[0])
    gc = pd.read_csv(T / "gene_calls.tsv", sep="\t"); gc = gc[gc.organelle == org]
    E = pd.read_csv(T / "edited_gene_calls.tsv", sep="\t"); E = E[E.organelle == org].sort_values("score", ascending=False)
    done = set()
    for r in E.itertuples():
        if r.gene in done: continue
        m = gc[(gc.gene == r.gene) & (gc.contig_id == r.contig) & (gc.start < r.hit_end + 300) & (gc.end > r.hit_start - 300)]
        if m.empty: continue
        lo, hi = sorted((int(r.corrected_start), int(r.corrected_end))); lo, hi = max(lo, 0), min(hi, len(genome[r.contig]))
        tr = rna.walk(pile, genome, r.contig, lo, hi, m.iloc[0].strand); sites.update(rna.truth(tr))
        d = rna.derive(tr, REV)
        why = d["status"] if d["status"] != "ok" else ("length not a multiple of 3" if d["mod3"] else "ok")
        # a call that did not reach both ends of its profile can be part of the protein (rbcL: the HMM hit is half the
        # gene and the window round it does not hold the rest), and a profile built from such parts would hide that.
        # It is the whole protein all the same if it runs from a Met to a stop codon (the profile was the longer one).
        if why == "ok" and r.n_columns_aligned < MIN_ALIGNED * r.n_columns_total and not (d["has_stop"] and d["protein"][:1] == "M"): why = "partial call"
        sp_stat.append((sp, r.gene, "single", why, d["start"], d["over_extended"]))
        if why == "ok": sp_prots.append((sp, r.gene, d["marked"])); done.add(r.gene)
    X = pd.read_csv(T / "reconstructed_exons.tsv", sep="\t"); X = X[X.organelle == org]
    G = pd.read_csv(T / "reconstructed_genes.tsv", sep="\t"); G = G[G.organelle == org]
    RF = pd.read_csv(T / "reconstructed_frameshifts.tsv", sep="\t"); fs = set(RF[RF.organelle == org].gene) if len(RF) else set()
    if True:
        for g in G.itertuples():
            xs = X[X.gene == g.gene].sort_values("slot"); why, d = "", None
            if str(g.complete).lower() != "true": why = "incomplete"
            elif g.gene in fs: why = "frameshift call"
            if not why:
                tr, owner = [], []
                for x in xs.itertuples():
                    lo, hi = sorted((int(x.start), int(x.end))); w = rna.walk(pile, genome, x.contig, max(lo, 0), min(hi, len(genome[x.contig])), x.strand)
                    tr += w; owner += [x.slot] * len(w)
                sites.update(rna.truth(tr)); d = rna.derive(tr, REV); why = "" if d["status"] == "ok" else d["status"]
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
            sp_stat.append((sp, g.gene, "spliced", why or "ok", d["start"] if d else 0, d["over_extended"] if d else 0))
            if why: continue
            acc, k = f"RNA_{sp}", len(exons)
            for i, (x, lo, hi, raw) in enumerate(exons, 1):
                sp_fasta.append(f">{acc}|{g.gene}|{g.gene}_1|exon{i}of{k}|{lo + 1}-{hi}|{x.strand}\n{raw}\n")
                nx = exons[i][0] if i < k else None
                jt = "" if nx is None else ("cis" if nx.contig == x.contig and nx.strand == x.strand else "trans")
                sp_jrows.append(dict(accession=acc, gene=g.gene, locus_tag=f"{g.gene}_1", exon_num=i, n_exons=k, start=lo + 1, end=hi, strand=x.strand, junction_to_next=jt, gap_bp="", notes=f"{sp}: RNA-seq edited"))
            sp_fasta.append(f">{acc}|{g.gene}|{g.gene}_1|translation_REFSEQ_ANNOTATED|len={len(d['protein'])}\n{d['protein']}\n")
    n_cu, n_uc = sum(v == "C>U" for v in sites.values()), sum(v == "U>C" for v in sites.values())
    cls = "ok" if REV or n_uc <= max(3, 0.05 * n_cu) else "left out" if n_uc > max(50, 0.5 * n_cu) else "noisy"
    noise.append((sp, n_cu, n_uc, cls))
    if cls == "left out": sp_stat = [(a, b, c, "reads not the assembly's" if w == "ok" else w, e, f) for a, b, c, w, e, f in sp_stat]; sp_prots, sp_jrows, sp_fasta = [], [], []
    if cls == "noisy":
        sp_prots = [(a, b, p.upper()) for a, b, p in sp_prots]
        for j in sp_jrows: j["notes"] += "; sites not used"
    stat += sp_stat; prots += sp_prots; jrows += sp_jrows
    if sp_fasta: open(out / "raw_rna" / f"RNA_{sp}.fasta", "w").write("".join(sp_fasta))
with open(out / "raw_rna" / "junctions.tsv", "w") as fh:
    w = csv.DictWriter(fh, jcols, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(jrows)
with open(out / ("rna_proteins.fasta" if org == "mito" else "rna_proteins.pltd.fasta"), "w") as fh:
    for sp, gene, p in prots: fh.write(f">{sp}|{gene}\n{p}\n")
S = pd.DataFrame(stat, columns=["species", "gene", "kind", "status", "start_trim", "over_extended"]); S.to_csv(out / "loci.tsv", sep="\t", index=False)
N = pd.DataFrame(noise, columns=["species", "c_to_u_sites", "u_to_c_sites", "reads"]); N.to_csv(out / "species.tsv", sep="\t", index=False)
print("species:", N.reads.value_counts().to_dict(), "| proteins with sites:", sum(p != p.upper() for _, _, p in prots))
pd.set_option("display.width", 250)
print(S.groupby(["kind", "status"]).size().to_string()); print("single-exon genes:", S[(S.kind == "single") & (S.status == "ok")].gene.nunique(), "| proteins:", len(prots), "| spliced loci:", (S[S.kind == "spliced"].status == "ok").sum())
print(S[(S.status == "ok")].groupby("kind").agg(start_trimmed=("start_trim", lambda x: (x > 0).sum()), over_extended=("over_extended", lambda x: (x > 0).sum())).to_string())
print(S[S.kind == "spliced"].pivot_table(index="gene", columns="status", values="species", aggfunc="count", fill_value=0).to_string())
