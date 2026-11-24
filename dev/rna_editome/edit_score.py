"""edit_score.py LABEL SOURCE NAME... : predicted RNA edits against the edits a GenBank record's annotation applies,
codon by codon, the same way for poa and PMGA. SOURCE is poa:TAG (predicted_edits.tsv and reconstructed_edits.tsv of
pilot/TAG/NAME) or pmga:FILE (build/pmga/FILE, map_editing_sites.py's table; NAME without _refseq there).
Truth: in each annotated CDS of a canonical mitochondrial gene, the codons whose annotated residue differs from the
genome's own translation (the record must annotate edited proteins). A prediction that changes the amino acid and
falls in an annotated CDS is right when it is in such a codon and gives the annotated residue. Predictions outside
the annotated CDS, and silent ones, cannot be judged and are counted apart. Writes pilot/LABEL.edits.tsv."""
import sys, warnings
from pathlib import Path
import pandas as pd
warnings.simplefilter("ignore")
from Bio import SeqIO
from Bio.Seq import Seq
H = Path(__file__).resolve().parent.parent; R = H.parent.parent
GB0 = R / "src/poa/pipeline/trans_splicing/reference/genbank"
ACC = dict(l.split() for l in open(H / "extra_genomes" / "accessions.tsv"))
SYN = {"tatc": "mttb", "mat-r": "matr", "orfx": "mttb", "orf25": "atp4", "orfb": "atp8", "ccb206": "ccmb", "ccb256": "ccmc", "ccb452": "ccmfc", "ccb382": "ccmfn1", "ccb203": "ccmfn2", "cox-2": "cox2", "coxii": "cox2", "yejv": "ccmfc", "yejr": "ccmfn"}
MITO = set("nad1 nad2 nad3 nad4 nad4l nad5 nad6 nad7 nad9 cob cox1 cox2 cox3 atp1 atp4 atp6 atp8 atp9 ccmb ccmc ccmfc ccmfn ccmfn1 ccmfn2 matr mttb "
           "rpl2 rpl5 rpl10 rpl16 rps1 rps2 rps3 rps4 rps7 rps10 rps11 rps12 rps13 rps14 rps19 sdh3 sdh4".split())
AA3 = dict(Ala="A", Arg="R", Asn="N", Asp="D", Cys="C", Gln="Q", Glu="E", Gly="G", His="H", Ile="I", Leu="L", Lys="K", Met="M", Phe="F", Pro="P", Ser="S", Thr="T", Trp="W", Tyr="Y", Val="V", Stop="*", Ter="*", End="*")
label, source, names = sys.argv[1], sys.argv[2], sys.argv[3:]; kind, where = source.split(":"); rows = []
PM = pd.read_csv(R / "build/pmga" / where, sep="\t") if kind == "pmga" else None
def predictions(name):
    """[(0-based genome position, resulting amino acid, changes the amino acid)]"""
    if kind == "pmga":
        d = PM[PM.genome == name.replace("_refseq", "")]; out = []
        for r in d.itertuples():
            a, b = r.aa_change.split("->"); out.append((int(r.position) - 1, AA3.get(b, "?"), a != b))
        return out
    T = H / "pilot" / where / name / "tables"; out = []
    for f in ("predicted_edits", "reconstructed_edits"):
        d = pd.read_csv(T / f"{f}.tsv", sep="\t"); d = d[d.organelle == "mito"] if "organelle" in d else d
        out += [(int(r.genomic_pos), str(r.resulting_aa), True) for r in d.itertuples()]
    return list(set(out))
for name in names:
    try: pred = predictions(name)
    except FileNotFoundError: print("skip", name); continue
    GB = GB0 if (GB0 / f"{ACC[name]}.gb").exists() else H / "extra_genomes" / "gb"; rec = SeqIO.read(GB / f"{ACC[name]}.gb", "genbank")
    truth, cds_pos = {}, {}   # genome position -> (codon id, annotated residue) for edited codons; every CDS position -> codon id
    for fi, f in enumerate(rec.features):
        if f.type != "CDS" or "gene" not in f.qualifiers or "translation" not in f.qualifiers or "pseudo" in f.qualifiers: continue
        g = SYN.get(f.qualifiers["gene"][0].lower(), f.qualifiers["gene"][0].lower())
        if g not in MITO: continue
        parts = list(f.location.parts)
        if "".join(str(p.extract(rec.seq)) for p in parts) != str(f.extract(rec.seq)): parts = parts[::-1]
        pos = []
        for p in parts: pos += list(range(int(p.start), int(p.end))) if p.strand != -1 else list(range(int(p.end) - 1, int(p.start) - 1, -1))
        cds = str(f.extract(rec.seq)).upper(); prot = f.qualifiers["translation"][0]; geno = str(Seq(cds[: len(cds) - len(cds) % 3]).translate())
        # a CDS whose annotated protein is not the genome's translation plus edits (a different frame or exon set
        # from the one the location gives) is left out
        if abs(len(geno.rstrip("*")) - len(prot)) > 1 or sum(x != y for x, y in zip(geno, prot)) > 0.25 * len(prot): continue
        for i in range(min(len(prot), len(geno))):
            cid = (fi, i)
            for q in pos[3 * i: 3 * i + 3]:
                cds_pos.setdefault(q, cid)
                if prot[i] != geno[i]: truth.setdefault(q, (cid, prot[i]))
    codons = {c for c, _ in truth.values()}; hit = set(); n_in = n_right = n_out = n_silent = 0
    for q, aa, changes in pred:
        if not changes: n_silent += 1; continue
        if q not in cds_pos: n_out += 1; continue
        n_in += 1
        if q in truth and truth[q][1] == aa: n_right += 1; hit.add(truth[q][0])
    rows.append(dict(genome=name, edited_codons=len(codons), found=len(hit), calls_in_cds=n_in, calls_right=n_right, calls_outside_annotated_cds=n_out, silent_calls=n_silent))
D = pd.DataFrame(rows); D.to_csv(H / "pilot" / f"{label}.edits.tsv", sep="\t", index=False)
print(D.to_string(index=False))
t = D.sum(numeric_only=True)
print(f"{label} TOTAL: edited codons {int(t.edited_codons)}, found {int(t.found)} ({100 * t.found / max(t.edited_codons, 1):.0f}%); calls in annotated CDS {int(t.calls_in_cds)}, right {int(t.calls_right)} ({100 * t.calls_right / max(t.calls_in_cds, 1):.0f}%); outside {int(t.calls_outside_annotated_cds)}; silent {int(t.silent_calls)}")
