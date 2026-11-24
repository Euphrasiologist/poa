"""gene_ends.py TAG NAME... : single-exon annotated CDS of RefSeq test genomes against poa's call of the gene: offsets of start and stop
(bp along the gene; + = poa's end lies downstream of the annotated one), start codons, and what lies between."""
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
norm = lambda g: SYN.get(g.lower(), g.lower())
tag, names = sys.argv[1], sys.argv[2:]; rows = []
for name in names:
    T = H / "pilot" / tag / name / "tables"
    E = pd.read_csv(T / "edited_gene_calls.tsv", sep="\t"); E = E[E.organelle == "mito"]
    GB = GB0 if (GB0 / f"{ACC[name]}.gb").exists() else H / "extra_genomes" / "gb"; rec = SeqIO.read(GB / f"{ACC[name]}.gb", "genbank"); seq = rec.seq
    for f in rec.features:
        if f.type != "CDS" or "gene" not in f.qualifiers or "translation" not in f.qualifiers: continue
        g = norm(f.qualifiers["gene"][0])
        if g not in MITO or "pseudo" in f.qualifiers or len(f.location.parts) != 1: continue
        a, b, st = int(f.location.start), int(f.location.end), f.location.strand
        cds = str(f.extract(seq)).upper(); has_stop = cds[-3:] in ("TAA", "TAG", "TGA")
        cand = E[(E.gene.str.lower().map(norm) == g) | ((g in ("ccmfn1", "ccmfn2", "ccmfn")) & E.gene.str.lower().isin(["ccmfn", "ccmfn1", "ccmfn2"]))]
        best = None
        for r in cand.itertuples():
            pa, pb = sorted((int(r.corrected_start), int(r.corrected_end)))
            ov = min(b, pb) - max(a, pa)
            if ov > 0 and (best is None or ov > best[0]): best = (ov, pa, pb, r)
        if best is None: rows.append(dict(genome=name, gene=g, called=False)); continue
        _, pa, pb, r = best
        if st == 1: ds, de = pa - a, pb - b
        else: ds, de = b - pb, a - pa
        # normalise the stop: annotation with the stop codon, poa without (or with): take the offset nearest zero of de, de+3, de-3
        de_n = min((de, de + 3, de - 3), key=abs)
        def codon_at(off):  # codon starting off bp downstream of the annotated start
            if st == 1: p = a + off; c = seq[p:p + 3]
            else: p = b - off; c = seq[p - 3:p].reverse_complement()
            return str(c).upper()
        between = ""
        if ds % 3 == 0 and ds != 0:
            lo, hi = sorted((0, ds)); cods = [codon_at(o) for o in range(lo, hi, 3)]
            between = f"{sum(c in ('TAA','TAG','TGA') for c in cods)} stops, {sum(c == 'ATG' for c in cods)} ATG in {len(cods)} codons"
        rows.append(dict(genome=name.replace("_refseq", ""), gene=g, called=True, len_aa=len(cds) // 3, ds=ds, de=de_n, ann_start=cds[:3], poa_start=codon_at(ds), between=between,
                         aligned=f"{r.n_columns_aligned}/{r.n_columns_total}", stop_in=r.includes_stop))
D = pd.DataFrame(rows); D.to_csv(H / "pilot" / f"{tag}.ends.tsv", sep="\t", index=False)
C = D[D.called]
ok_s, ok_e = C.ds.abs() <= 3, C.de.abs() <= 3
print(f"{tag}: {len(D)} single-exon CDS, called {len(C)} | both ends within 3 bp {int((ok_s & ok_e).sum())} | start off only {int((~ok_s & ok_e).sum())} | stop off only {int((ok_s & ~ok_e).sum())} | both off {int((~ok_s & ~ok_e).sum())}")
S = C[~ok_s]
print("start off: poa starts downstream (shorter)", int((S.ds > 0).sum()), "upstream (longer)", int((S.ds < 0).sum()), "| in frame", int((S.ds % 3 == 0).sum()))
print("  |offset| codons:", pd.cut(S.ds.abs() / 3, [0, 2, 5, 10, 20, 50, 1e4]).value_counts().sort_index().to_dict())
print("  annotated start codon:", S.ann_start.value_counts().head(8).to_dict()); print("  poa start codon:", S.poa_start.value_counts().head(8).to_dict())
print("  all annotated start codons:", C.ann_start.value_counts().head(8).to_dict())
print("  by gene:", S.gene.value_counts().to_dict())
X = C[~ok_e]
print("stop off: poa ends downstream (longer)", int((X.de > 0).sum()), "upstream (shorter)", int((X.de < 0).sum())); print("  |offset| codons:", pd.cut(X.de.abs() / 3, [0, 2, 5, 10, 20, 50, 1e4]).value_counts().sort_index().to_dict()); print("  by gene:", X.gene.value_counts().to_dict())
