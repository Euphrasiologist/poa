"""angio_readme.py MITO_STAGE PLTD_STAGE : the README of the angiosperm RNA-derived references, from the stage directories
rna_refs.py wrote (species.tsv, loci.tsv) and the run tables (pub3.tsv, pub4.tsv). Prints it."""
import sys, pandas as pd
m, p = sys.argv[1], sys.argv[2]
runs = {}
for f in ("pub3.tsv", "pub4.tsv"):
    for line in open(f):
        sp, run = line.split("\t")[:2]; runs.setdefault(sp, []).append(run)
S = {o: pd.read_csv(f"stage/{d}/species.tsv", sep="\t").set_index("species") for o, d in (("mito", m), ("pltd", p))}
L = {o: pd.read_csv(f"stage/{d}/loci.tsv", sep="\t") for o, d in (("mito", m), ("pltd", p))}
ok = {o: L[o][L[o].status == "ok"] for o in L}
n = {o: ok[o].groupby(["species", "kind"]).size().unstack(fill_value=0) for o in ok}
def count(o, kind): return int(n[o][kind].sum()) if kind in n[o] else 0
def cls(o): return ", ".join(f"{v} {k}" for k, v in S[o].reads.value_counts().items())
species = sorted(set(S["mito"].index) | set(S["pltd"].index))
used = {o: sorted(set(ok[o].species)) for o in ok}
print(f"""# Angiosperm editing references from RNA-seq

`rna_proteins.fasta` and `rna_proteins.pltd.fasta` (headers `species|gene`) are single-exon mitochondrial and plastid
protein genes of Darwin Tree of Life angiosperm assemblies, each the gene poa called with the C-to-U editing the
species' own RNA-seq shows applied: {count("mito", "single")} mitochondrial proteins of {ok["mito"][ok["mito"].kind == "single"].gene.nunique()} genes from {len(used["mito"])} species, and {count("pltd", "single")} plastid proteins of
{ok["pltd"][ok["pltd"].kind == "single"].gene.nunique()} genes from {len(used["pltd"])}. The spliced genes of the same species ({count("mito", "spliced")} mitochondrial and {count("pltd", "spliced")} plastid loci) are in
`trans_splicing/reference/raw_rna/` and `trans_splicing/reference/plastid/raw_rna/`. With the RefSeq mitogenomes
annotated with edited proteins (`reference_genes.ANGIOSPERM_EDITED_MITOGENOMES`) they are what
`00b_build_lineage_profiles.py --lineage angiosperm` builds `profiles_refseq/mito/` from, and alone what
`--organelle pltd` builds `pltd/` from.

Residues in lower case are those a C-to-U edit made. From them each profile takes its editing sites (the `edit`
column; `editing/README.md`).

## Method

`dev/rna_editome/rna_refs.py`, run with `RNA_REVERSE=0`. The RNA-seq is public (ENA; the runs are below), one paired
run a species, mapped with bwa-mem2 to the species' own mitochondrial and plastid contigs
(`dev/rna_editome/map4.sh`). A site is edited where at least two reads and a fifth of the reads show U for the
genome's C. A gene is used when at least 90% of it has 10 or more reads and, with those edits applied, it
translates without an internal stop: for a single-exon gene the reading frame poa corrected, for a spliced gene the
exons transsplice placed, complete and without a frameshift call.

The reads are often not from the plant that was assembled, and where the two differ by a C/T the reads look edited.
No angiosperm edits U to C, so the sites where the reads show C for the genome's T count those differences. A species
with more than max(3, 5% of its C-to-U sites) of them is "noisy": its proteins are used, in upper case throughout,
and its loci are noted `sites not used` in `junctions.tsv`, so that they shape the profile but are no evidence of
where editing happens. One with more than max(50, half) is left out: the reads are another species'.
Mitochondrion: {cls("mito")}. Plastid: {cls("pltd")}.

## Species and runs

| species | run | mito C-to-U / U-to-C sites | mito | plastid C-to-U / U-to-C sites | plastid |
|---|---|---|---|---|---|""")
for sp in species:
    def cell(o):
        if sp not in S[o].index: return "", ""
        r = S[o].loc[sp]; return f"{r.c_to_u_sites} / {r.u_to_c_sites}", r.reads
    a, b = cell("mito"), cell("pltd")
    print(f"| {sp.replace('_', ' ')} | {', '.join(runs.get(sp, ['?']))} | {a[0]} | {a[1]} | {b[0]} | {b[1]} |")
