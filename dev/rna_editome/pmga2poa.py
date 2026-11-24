"""pmga2poa.py DIR TAG SPECIES... : PMGA's protein-gene models (build/pmga/DIR/<sp>/<sp>.gff) as poa-shaped tables in
pilot/TAG/<sp>/tables, so assess.py scores them against RNA-seq exactly as it scores poa's. PMGA works on the contigs
joined end to end; coordinates go back to each contig (<sp>_merged_stat.txt). A gene's exons are its CDS rows; a
single-exon gene becomes an edited_gene_calls row (no edits, PMGA's container predicts none), a multi-exon one
reconstructed exons. A second copy of a gene name is called <gene>#2."""
import sys
from pathlib import Path
import pandas as pd
R = Path("/lustre/scratch122/tol/teams/blaxter/users/mb39/plant_organelle_annotator"); H = R / "build/rnaedit"
d, tag, species = sys.argv[1], sys.argv[2], sys.argv[3:]
# PMGA's predicted editing sites in genome coordinates (build/pmga/scripts/map_editing_sites.py), if the editing step was run
SITES = pd.read_csv(R / "build/pmga/editing_sites.tsv", sep="\t") if (R / "build/pmga/editing_sites.tsv").exists() and d == "out" else None
for sp in species:
    P = R / "build/pmga" / d / sp
    if not (P / f"{sp}.gff").exists(): print("skip", sp); continue
    ctg = []
    for line in open(P / f"{sp}_merged_stat.txt"):
        f = line.split()
        if len(f) >= 4 and f[0].startswith("contig"): ctg.append((int(f[1]), int(f[2]), f[3]))
    def local(s, e):  # 1-based inclusive on the joined sequence -> contig, 0-based half-open
        for a, b, name in ctg:
            if a <= s <= b: return name, s - a, e - a + 1
        raise ValueError((s, e))
    mrna, cds = {}, {}
    for line in open(P / f"{sp}.gff"):
        if line.startswith("##FASTA"): break
        f = line.rstrip("\n").split("\t")
        if len(f) != 9: continue
        a = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x)
        if f[2] == "mRNA": mrna[a["ID"]] = a.get("Name", "")
        elif f[2] == "CDS": cds.setdefault(a["Parent"], []).append((*local(int(f[3]), int(f[4])), f[6]))
    gc, E, X, G, seen, PE, RE = [], [], [], [], {}, [], []
    sites = {} if SITES is None else {c: sorted(g.position - 1) for c, g in SITES[SITES.genome == sp].groupby("contig")}
    inside = lambda c, s, e: [q for q in sites.get(c, []) if s <= q < e]
    for m, ex in cds.items():
        g = mrna.get(m, ""); seen[g] = seen.get(g, 0) + 1; name = g if seen[g] == 1 else f"{g}#{seen[g]}"
        if len(ex) == 1:
            c, s, e, st = ex[0]
            gc.append(dict(species=sp, organelle="mito", contig_id=c, start=s, end=e, gene=name, score=1, strand=st))
            E.append(dict(species=sp, organelle="mito", gene=name, contig=c, hit_start=s, hit_end=e, score=1, corrected_start=s, corrected_end=e, n_edits=len(inside(c, s, e)), n_columns_aligned=1, n_columns_total=1))
            PE += [dict(species=sp, organelle="mito", gene=name, contig=c, hit_start=s, hit_end=e, genomic_pos=q) for q in inside(c, s, e)]
        else:
            for i, (c, s, e, st) in enumerate(ex, 1):
                X.append(dict(species=sp, organelle="mito", gene=name, slot=i, contig=c, start=s, end=e, strand=st))
                RE += [dict(species=sp, organelle="mito", gene=name, genomic_pos=q) for q in inside(c, s, e)]
            G.append(dict(species=sp, organelle="mito", gene=name, n_slots=len(ex), n_filled=len(ex), complete=True, whole_gene_n_edits=0))
    T = H / "pilot" / tag / sp / "tables"; T.mkdir(parents=True, exist_ok=True)
    w = lambda rows, cols, n: pd.DataFrame(rows, columns=cols).to_csv(T / f"{n}.tsv", sep="\t", index=False)
    w(gc, ["species", "organelle", "contig_id", "start", "end", "gene", "score", "strand"], "gene_calls")
    w(E, ["species", "organelle", "gene", "contig", "hit_start", "hit_end", "score", "corrected_start", "corrected_end", "n_edits", "n_columns_aligned", "n_columns_total"], "edited_gene_calls")
    w(PE, ["species", "organelle", "gene", "contig", "hit_start", "hit_end", "genomic_pos"], "predicted_edits")
    w([], ["species", "organelle", "gene", "contig", "hit_start", "hit_end"], "frameshifts")
    w(X, ["species", "organelle", "gene", "slot", "contig", "start", "end", "strand"], "reconstructed_exons")
    w(G, ["species", "organelle", "gene", "n_slots", "n_filled", "complete", "whole_gene_n_edits"], "reconstructed_genes")
    w(RE, ["species", "organelle", "gene", "genomic_pos"], "reconstructed_edits"); w([], ["species", "organelle", "gene"], "reconstructed_frameshifts")
    print(sp, len(E), "single-exon,", len(G), "multi-exon,", len(ctg), "contig(s),", len(PE) + len(RE), "predicted edits")
