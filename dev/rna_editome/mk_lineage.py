"""mk_lineage.py LINEAGE BOOTDIR OUT [--only SP,SP] [--genes G,G] : a lineage's reference files from boot.py gene models.

Takes each species' models with status ok, one locus per species and gene, whose exon count is the gene's most
common one and whose protein is within 10% of the gene's median length. Genes with more than one exon (or in the
package's spliced list, GENES) go to OUT/raw_rna (RNA_<species>.fasta, junctions.tsv: the raw/ format of
trans_splicing), the rest to OUT/rna_proteins.fasta (species|gene). OUT/used.tsv lists what went where."""
import sys, glob, csv, re
from pathlib import Path
import pandas as pd
H = Path(__file__).resolve().parent
lin, boot, out = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
opt = dict(zip(sys.argv[4::2], sys.argv[5::2]))
only = set(opt["--only"].split(",")) if "--only" in opt else None
genes = set(opt["--genes"].split(",")) if "--genes" in opt else None
ALWAYS = ["nad1", "nad2", "nad5", "rps3", "nad4", "nad7", "ccmFc", "cox2", "rpl2", "rps10"]
G = pd.concat(pd.read_csv(f, sep="\t") for f in sorted(boot.glob("*.genes.tsv")))
if only: G = G[G.species.isin(only)]
if genes: G = G[G.gene.isin(genes)]
G = G[G.status == "ok"].drop_duplicates(["species", "gene"])
keep = []
for gene, g in G.groupby("gene"):
    k = g.n_exons.mode().iloc[0]; g = g[g.n_exons == k]; med = g.aa.median()
    if len(g) == 2 and abs(g.aa.iloc[0] - g.aa.iloc[1]) > 0.1 * med:
        # two models that disagree: the one a homologue (Marchantia) covers better
        miss = g.march_missing.map(lambda x: sum(map(int, x.split("/"))) if isinstance(x, str) else 10**6)
        g = g[miss == miss.min()].head(1); med = g.aa.median()
    keep.append(g[(g.aa - med).abs() <= 0.1 * med])
G = pd.concat(keep)
seqs = {}
for sp in G.species.unique():
    name = None
    for line in open(boot / f"{sp}.fasta"):
        if line.startswith(">"): name = line[1:].strip(); seqs[name] = ""
        else: seqs[name] += line.strip()
(out / "raw_rna").mkdir(parents=True, exist_ok=True)
for f in (out / "raw_rna").glob("RNA_*.fasta"): f.unlink()
jcols = ["accession", "gene", "locus_tag", "exon_num", "n_exons", "start", "end", "strand", "junction_to_next", "gap_bp", "notes"]
jrows, prot, used = [], [], []
for r in G.sort_values(["species", "gene"]).itertuples():
    acc = f"RNA_{r.species}"; spliced = r.n_exons > 1 or r.gene in ALWAYS
    used.append((r.species, r.gene, "template" if spliced else "protein", r.n_exons, r.aa, r.note if isinstance(r.note, str) else "", r.junction_reads if isinstance(r.junction_reads, str) else ""))
    if not spliced: prot.append((r.species, r.gene, r.protein)); continue
    recs = sorted(((int(re.search(r"exon(\d+)of", h).group(1)), h, s) for h, s in seqs.items() if h.startswith(f"{acc}|{r.gene}|{r.locus}|exon")))
    with open(out / "raw_rna" / f"{acc}.fasta", "a") as fh:
        for i, h, s in recs:
            f_ = h.split("|"); lo, hi = f_[4].split("-"); contig = f_[6]
            fh.write(f">{acc}|{r.gene}|{r.gene}_1|exon{i}of{len(recs)}|{lo}-{hi}|{f_[5]}\n{s}\n")
            nxt = recs[i][1].split("|") if i < len(recs) else None
            jt = "" if nxt is None else ("cis" if (nxt[6] == contig and nxt[5] == f_[5]) else "trans")
            jrows.append(dict(accession=acc, gene=r.gene, locus_tag=f"{r.gene}_1", exon_num=i, n_exons=len(recs), start=lo, end=hi, strand=f_[5], junction_to_next=jt, gap_bp="", notes=f"{r.species}: RNA-seq gene model"))
        fh.write(f">{acc}|{r.gene}|{r.gene}_1|translation_REFSEQ_ANNOTATED|len={len(r.protein)}\n{r.protein}\n")
with open(out / "raw_rna" / "junctions.tsv", "w") as fh:
    w = csv.DictWriter(fh, jcols, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(jrows)
with open(out / "rna_proteins.fasta", "w") as fh:
    for sp, gene, p in prot: fh.write(f">{sp}|{gene}\n{p}\n")
U = pd.DataFrame(used, columns=["species", "gene", "use", "n_exons", "aa", "start_note", "junction_reads"]); U.to_csv(out / "used.tsv", sep="\t", index=False)
print(U.pivot_table(index="gene", columns="use", values="species", aggfunc="count", fill_value=0).T.to_string())
