"""Gene ends against RNA: how often a call runs past a stop codon the RNA shows is never edited, on the committed and
the new fern profiles/templates. Single-exon: wire/single2{.old,}.calls.tsv; spliced: ts.cur and wire/spliced."""
import sys, glob, re
from pathlib import Path
import pandas as pd
from Bio import SeqIO
R = Path("/lustre/scratch122/tol/teams/blaxter/users/mb39/plant_organelle_annotator"); H = R / "build/rnaedit"; sys.path.insert(0, str(H)); import rna
species = sorted(re.sub(r"^tol\.|\.pile\.tsv\.gz$", "", Path(p).name) for p in glob.glob(str(H / "tol.*.pile.tsv.gz")))
man = pd.read_csv(H / "wire/single2.tsv", sep="\t"); strand_of = {(r.species, r.gene): r.strand for r in man.itertuples()}
tm = pd.read_csv(R / "build/fern_tpl/m.fern.tsv", sep="\t"); fasta_of = dict(zip(tm.species, tm.fasta))
S = {k: pd.read_csv(H / f"wire/{v}.calls.tsv", sep="\t") for k, v in (("committed", "single2.old"), ("new", "single2"))}
X = {"committed": pd.read_csv(H / "ts.cur.exons.tsv", sep="\t"), "new": pd.read_csv(H / "wire/spliced.exons.tsv", sep="\t")}
rows = []
for sp in species:
    pile = rna.load_pile(H / f"tol.{sp}.pile.tsv.gz"); genome = {r.id: str(r.seq).upper() for r in SeqIO.parse(fasta_of[sp], "fasta")}
    for k, c in S.items():
        for r in c[c.species == sp].itertuples():
            contig = r.gene.split("|")[1].split(":")[0]; lo, hi = sorted((int(r.corrected_start), int(r.corrected_end)))
            if hi > len(genome[contig]): continue
            d = rna.derive(rna.walk(pile, genome, contig, lo, hi, strand_of[(sp, r.gene)]))
            rows.append(("single-exon", k, sp, r.gene.split("|")[0], d["status"], d["over_extended"], d["start"], d["has_stop"]))
    for k, x in X.items():
        for gene, xs in x[(x.species == sp) & (x.organelle == "mito")].groupby("gene"):
            tr = []
            for e in xs.sort_values("slot").itertuples():
                lo, hi = sorted((int(e.start), int(e.end))); tr += rna.walk(pile, genome, e.contig, lo, hi, e.strand)
            d = rna.derive(tr); rows.append(("spliced", k, sp, gene, d["status"], d["over_extended"], d["start"], d["has_stop"]))
T = pd.DataFrame(rows, columns=["kind", "set", "species", "gene", "status", "extra_codons", "start_trim", "real_stop_seen"]); T.to_csv(H / "ends.tsv", sep="\t", index=False)
pd.set_option("display.width", 220)
print(T.groupby(["kind", "set"]).agg(genes=("gene", "size"), real_stop_seen=("real_stop_seen", "sum"), run_past_stop=("extra_codons", lambda x: (x > 0).sum()),
      median_extra_codons=("extra_codons", lambda x: int(x[x > 0].median()) if (x > 0).any() else 0), starts_upstream=("start_trim", lambda x: (x > 0).sum()), rna_protein_ok=("status", lambda x: (x == "ok").sum())).to_string())
print(T[T.extra_codons > 0].groupby(["kind", "set", "gene"]).extra_codons.agg(["size", "median"]).unstack("set").fillna(0).astype(int).to_string())
