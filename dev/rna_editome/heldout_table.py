"""heldout_table.py BEFORE FLAT SITES : the held-out table of reports/pmga_comparison (rows from the RNA-scored tables
pilot/TAG.genes.tsv and TAG.pltd.genes.tsv, over the genes all three have with >= 90% RNA coverage). The RefSeq rows are
given as arguments after the tags: three "placed found% right%" triples."""
import sys, pandas as pd
tags = sys.argv[1:4]; gb = sys.argv[4:]
names = ["edited RefSeq references", "RNA-derived references", "with editing sites"]
import os
ST = os.environ.get("STAGE", "f")  # species.tsv of stage/angio_{mito,pltd}_<ST>: whose reads are the assembly's own ("ok")
clean = {o: set(pd.read_csv(f"stage/angio_{o}_{ST}/species.tsv", sep="\t").query("reads == 'ok'").species) for o in ("mito", "pltd")}
out = {n: {} for n in names}; order = []
for org, ext in (("mito", "genes.tsv"), ("pltd", "pltd.genes.tsv")):
    D = {t: pd.read_csv(f"pilot/{t}.{ext}", sep="\t").drop_duplicates(["species", "kind", "gene"]).set_index(["species", "kind", "gene"]) for t in tags}
    keep = None
    for d in D.values():
        k = set(d.index[d.covered >= 0.9]); keep = k if keep is None else keep & k
    for n, t in zip(names, tags):
        d = D[t].loc[sorted(keep)].reset_index()
        for kind, x in d.groupby("kind"):
            rna = int(x.rna_CU.sum()); lab = f"{'mitochondrial' if org == 'mito' else 'plastid'} {'single-exon' if kind == 'single' else 'spliced'} genes"
            c = x[x.species.isin(clean[org])]; crna = int(c.rna_CU.sum())
            rows = [(f"{lab}: aligned end to end (of {len(x):,})", f"{int((x.frac >= 0.99).sum()):,}"),
                    ("... C-to-U sites in the RNA found", f"{100 * x.right.sum() / max(rna, 1):.0f}% (of {rna:,})"),
                    (f"... the same, {c.species.nunique()} species whose reads are the assembly's own", f"{100 * c.right.sum() / max(crna, 1):.0f}% (of {crna:,})"),
                    ("... edit calls in the RNA", f"{100 * x.right.sum() / max(x.called.sum(), 1):.0f}% ({int(x.right.sum()):,} of {int(x.called.sum()):,})")]
            for i, (k, v) in enumerate(rows):
                key = (lab, i, k)
                if key not in order: order.append(key)
                out[n][key] = v
print("\t".join(["", *names]))
for key in order: print("\t".join([key[2] if key[1] == 0 else key[2], *[out[n].get(key, "") for n in names]]))
if gb:
    labs = ["RefSeq coding sequences with every exon placed (of 655)", "RefSeq annotated edits found", "... edit calls that are annotated"]
    for i, lab in enumerate(labs): print("\t".join([lab, *[gb[3 * j + i] for j in range(3)]]))
