"""gb_edit.py OUT NAME... : a GenBank genome's annotated genes with the RNA editing its own RNA-seq shows applied.
NAME is in extra_genomes/accessions.tsv, with pub.NAME.*.pile.tsv.gz from map5.sh. Each annotated CDS (its exons as
annotated) is walked through the pileup and edited by the rules of rna.derive; a gene is used if its protein comes
out clean (no stop left, >= 90% of it with 10 or more reads, a multiple of three long). Appends to OUT/raw_rna
(spliced genes: RNA_<NAME>.fasta, junctions.tsv) and OUT/rna_proteins.fasta (single-exon genes), the formats of
mk_lineage.py, and prints what was used."""
import sys, csv, warnings
from pathlib import Path
import pandas as pd
warnings.simplefilter("ignore")
from Bio import SeqIO
H = Path(__file__).resolve().parent; sys.path.insert(0, str(H)); import rna, boot
R = H.parent.parent; GB = R / "src/poa/pipeline/trans_splicing/reference/genbank"
ACC = dict(l.split() for l in open(H / "extra_genomes" / "accessions.tsv"))
SYN = {"tatc": "mttB", "mat-r": "matR"}
out, names = Path(sys.argv[1]), sys.argv[2:]
(out / "raw_rna").mkdir(parents=True, exist_ok=True)
jcols = ["accession", "gene", "locus_tag", "exon_num", "n_exons", "start", "end", "strand", "junction_to_next", "gap_bp", "notes"]
jp = out / "raw_rna" / "junctions.tsv"
jrows = list(csv.DictReader(open(jp), delimiter="\t")) if jp.exists() else []
stat = []
for name in names:
    acc = f"RNA_{name}"; jrows = [r for r in jrows if r["accession"] != acc]
    rec = SeqIO.read(GB / f"{ACC[name]}.gb", "genbank"); genome = {rec.id: str(rec.seq).upper()}; pile = boot.load_piles(name)
    pile = {(rec.id if c != rec.id else c, p): v for (c, p), v in pile.items()}
    seen, prots = set(), []
    with open(out / "raw_rna" / f"{acc}.fasta", "w") as fh:
        for f in rec.features:
            if f.type != "CDS" or "gene" not in f.qualifiers or "translation" not in f.qualifiers: continue
            g = f.qualifiers["gene"][0]; g = SYN.get(g.lower(), g)
            if g in seen or g.startswith(("ORF", "orf", "mat-")): continue
            parts = list(f.location.parts)
            if "".join(str(p.extract(rec.seq)) for p in parts) != str(f.extract(rec.seq)): parts = parts[::-1]
            tr, owner = [], []
            for i, p in enumerate(parts):
                w = rna.walk(pile, genome, rec.id, int(p.start), int(p.end), "+" if p.strand == 1 else "-"); tr += w; owner += [i] * len(w)
            d = rna.derive(tr); why = d["status"] if d["status"] != "ok" else ("length not a multiple of 3" if d["mod3"] else "ok")
            tru = rna.truth(tr)
            stat.append((name, g, len(parts), why, len(d["protein"]), len(f.qualifiers["translation"][0]), sum(v == "C>U" for v in tru.values()), sum(v == "U>C" for v in tru.values()), d["covered"]))
            if why != "ok": continue
            seen.add(g)
            if len(parts) == 1: prots.append((name, g, d["protein"])); continue
            first, last = 3 * d["start"], 3 * d["end"] + (3 if d["has_stop"] else 0); ex = []
            for i, p in enumerate(parts):
                ps = [tr[k][0][1] for k in range(first, min(last, len(tr))) if owner[k] == i]
                if not ps: continue
                lo, hi = min(ps) - 1, max(ps); raw = genome[rec.id][lo:hi]; s = "+" if p.strand == 1 else "-"
                if s == "-": raw = "".join(rna.comp[b] for b in reversed(raw))
                ex.append((lo, hi, s, raw))
            for i, (lo, hi, s, raw) in enumerate(ex, 1):
                fh.write(f">{acc}|{g}|{g}_1|exon{i}of{len(ex)}|{lo + 1}-{hi}|{s}\n{raw}\n")
                jt = "" if i == len(ex) else ("cis" if ex[i][2] == s and 0 < ((ex[i][0] - hi) if s == "+" else (lo - ex[i][1])) <= 10000 else "trans")
                jrows.append(dict(accession=acc, gene=g, locus_tag=f"{g}_1", exon_num=i, n_exons=len(ex), start=lo + 1, end=hi, strand=s, junction_to_next=jt, gap_bp="", notes=f"{name} {ACC[name]}: annotated exons, RNA-seq edited"))
            fh.write(f">{acc}|{g}|{g}_1|translation_REFSEQ_ANNOTATED|len={len(d['protein'])}\n{d['protein']}\n")
    pf = out / "rna_proteins.fasta"; keep = []
    if pf.exists():
        lines = pf.read_text().split(">")[1:]; keep = [">" + x for x in lines if not x.startswith(name + "|")]
    pf.write_text("".join(keep) + "".join(f">{sp}|{g}\n{p}\n" for sp, g, p in prots))
with open(jp, "w") as fh:
    w = csv.DictWriter(fh, jcols, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(jrows)
S = pd.DataFrame(stat, columns=["genome", "gene", "exons", "status", "aa", "aa_annotated", "rna_CU", "rna_UC", "covered"]); S.to_csv(out / "gb_edit.tsv", sep="\t", index=False)
pd.set_option("display.width", 220); print(S.groupby(["genome", "status"]).size().to_string()); print(S[S.status == "ok"].groupby("genome").agg(genes=("gene", "size"), CU=("rna_CU", "sum"), UC=("rna_UC", "sum")).to_string()); print(S[S.status != "ok"].to_string(index=False))
