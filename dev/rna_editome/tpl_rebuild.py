"""tpl_rebuild.py : fern exon templates with RNA-derived loci from the dataset's ferns, cross-validated by genus.

RNA-derived locus: a spliced mito gene transsplice reconstructed completely (committed fern templates, ts.cur.*)
with no frameshift call and every exon covered >= 90% by RNA-seq; its exons with RNA edits at >= 0.5 applied are
joined and translated, and the locus is kept if that protein has no internal stop. Written in the raw/ format of
trans_splicing/reference (accession RNA_<species>) next to the GenBank fern loci, in tpl/raw.
Then per genus: templates built without that genus's RNA loci (unmasked), transsplice on its species, with stop
read-through only and with other U-to-C edits at cost REV. Outputs ts.cv.* and ts.cvrev.* for splice_score.py."""
import gzip, glob, re, shutil, subprocess, sys, csv
from pathlib import Path
import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO
R = Path("/lustre/scratch122/tol/teams/blaxter/users/mb39/plant_organelle_annotator"); sys.path.insert(0, str(R / "src"))
from poa.pipeline.trans_splicing import exon_profiles as ep
H = R / "build/rnaedit"; TS = str(R / "build/rev_edit/bin2/transsplice"); REV = sys.argv[1] if len(sys.argv) > 1 else "5"
VARIANT = sys.argv[2] if len(sys.argv) > 2 else "all"  # "edited": leave out the GenBank ferns whose proteins only have stops edited out
from poa.pipeline.trans_splicing import reference_genes as rg
STOP_ONLY = tuple(a for a, n in rg.FERN_MITOGENOMES.items() if n != "Haplopteris ensiformis") if VARIANT == "edited" else ()
SRC = sys.argv[3] if len(sys.argv) > 3 else "ts.cur"  # transsplice tables the RNA loci are cut from
SUF = ("" if VARIANT == "all" else "E") + (sys.argv[4] if len(sys.argv) > 4 else "")
MIN_DEPTH, MIN_FRAC = 10, 0.2
comp = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
species = sorted(re.sub(r"^tol\.|\.pile\.tsv\.gz$", "", Path(p).name) for p in glob.glob(str(H / "tol.*.pile.tsv.gz")))
man = pd.read_csv(R / "build/fern_tpl/m.fern.tsv", sep="\t"); man = man[man.species.isin(species) & (man.organelle == "mito")]
X = pd.read_csv(H / f"{SRC}.exons.tsv", sep="\t"); G = pd.read_csv(H / f"{SRC}.genes.tsv", sep="\t"); F = pd.read_csv(H / f"{SRC}.fs.tsv", sep="\t")
fs = set(zip(F.species, F.gene)); fasta_of = dict(zip(man.species, man.fasta))
pkg_raw = R / "src/poa/pipeline/trans_splicing/reference/lineages/fern/raw"
T = H / ("tpl" + SUF); shutil.rmtree(T, ignore_errors=True); (T / "raw").mkdir(parents=True)
for f in pkg_raw.glob("*.fasta"): shutil.copy(f, T / "raw" / f.name)
jrows = list(csv.DictReader(open(pkg_raw / "junctions.tsv"), delimiter="\t")); jcols = list(jrows[0].keys())
stat = []
sys.path.insert(0, str(H)); import rna
# exon boundaries the RNA disagrees with (junctions.py on SRC): moved to the junction most intron-spanning reads show
JC = {}
if (H / f"{SRC}.junction_check.tsv").exists() and len(sys.argv) > 5 and sys.argv[5] == "fixjunctions":
    jc = pd.read_csv(H / f"{SRC}.junction_check.tsv", sep="\t")
    jc = jc[(jc.status == "other junction") & (jc.reads_other >= 5) & (jc.reads_other >= 3 * jc.reads_agree)]
    for r in jc.itertuples(): JC[(r.species, r.gene, r.junction)] = (int(r.other_le), int(r.other_rs), int(r.other_gap))
n_fixed = 0
def fix_boundaries(sp, gene, xs, genome):
    """exon rows (dicts, slot order) with RNA-supported junctions applied. Read bases the aligner left between the
    two pieces (or gave to both) go to the side that leaves the intron starting most like a group II intron (GTGCG)."""
    global n_fixed
    xs = [dict(x._asdict()) for x in xs.itertuples()]
    for a, b in zip(xs, xs[1:]):
        key = (sp, gene, f"{a['slot']}-{b['slot']}")
        if key not in JC or a["contig"] != b["contig"] or a["strand"] != b["strand"]: continue
        le, rs, gap = JC[key]; best = None
        for t in range(0, abs(gap) + 1):
            l, r = (le + t, rs - (gap - t)) if gap >= 0 else (le - t, rs + (-gap - t))
            intron = genome[a["contig"]][l:r - 1]
            if a["strand"] == "-": intron = str(Seq(intron).reverse_complement())
            score = sum(x == y for x, y in zip(intron[:5], "GTGCG"))
            if best is None or score > best[0]: best = (score, l, r)
        _, l, r = best  # 1-based last exonic base before the intron, first after (genomic order)
        up, down = (a, b) if a["strand"] == "+" else (b, a)  # genomic-left exon, genomic-right exon
        lo_u, hi_u = sorted((up["start"], up["end"])); lo_d, hi_d = sorted((down["start"], down["end"]))
        if not (lo_u < l and r - 1 < hi_d): continue
        up["start"], up["end"] = lo_u, l; down["start"], down["end"] = r - 1, hi_d
        n_fixed += 1
    return pd.DataFrame(xs)

for sp in species:
    pile = rna.load_pile(H / f"tol.{sp}.pile.tsv.gz")
    genome = {r.id: str(r.seq).upper() for r in SeqIO.parse(fasta_of[sp], "fasta")}
    for g in G[(G.species == sp) & (G.organelle == "mito")].itertuples():
        xs = X[(X.species == sp) & (X.gene == g.gene)].sort_values("slot")
        if JC and len(xs): xs = fix_boundaries(sp, g.gene, xs, genome)
        why, d = "", None
        if str(g.complete).lower() != "true": why = "incomplete"
        elif (sp, g.gene) in fs: why = "frameshift call"
        if not why:
            tr, owner = [], []
            for x in xs.itertuples():
                lo, hi = sorted((int(x.start), int(x.end))); w = rna.walk(pile, genome, x.contig, lo, hi, x.strand)
                tr += w; owner += [x.slot] * len(w)
            d = rna.derive(tr)
            why = "" if d["status"] == "ok" else d["status"]
            if not why and d["mod3"]: why = "length not a multiple of 3"
        stat.append((sp, g.gene, len(xs), why or "ok", d["start"] if d else 0, d["over_extended"] if d else 0))
        if why: continue
        # the exons, cut to the protein's span (stop codon kept on the last one when there is one)
        first, last = 3 * d["start"], 3 * d["end"] + (3 if d["has_stop"] else 0)
        exons = []
        for x in xs.itertuples():
            idx = [i for i in range(first, min(last, len(tr))) if owner[i] == x.slot]
            if not idx: continue
            ps = [tr[i][0][1] for i in idx]; lo, hi = min(ps) - 1, max(ps)
            raw = genome[x.contig][lo:hi]
            if x.strand == "-": raw = str(Seq(raw).reverse_complement())
            exons.append((x, lo, hi, raw))
        if len(exons) != len(xs): stat[-1] = (sp, g.gene, len(xs), "exon lost in trimming", d["start"], d["over_extended"]); continue
        body = d["protein"]
        acc, locus, k = f"RNA_{sp}", f"{g.gene}_1", len(exons)
        with open(T / "raw" / f"{acc}_{g.gene}.fasta", "w") as fh:
            for i, (x, lo, hi, raw) in enumerate(exons, 1):
                fh.write(f">{acc}|{g.gene}|{locus}|exon{i}of{k}|{lo + 1}-{hi}|{x.strand}\n{raw}\n")
                jrows.append({**{c: "" for c in jcols}, "accession": acc, "gene": g.gene, "locus_tag": locus, "exon_num": i, "n_exons": k, "start": lo + 1, "end": hi,
                              "strand": x.strand, "junction_to_next": "cis" if i < k else "", "notes": f"{sp}: RNA-seq edited"})
            fh.write(f">{acc}|{g.gene}|{locus}|translation_REFSEQ_ANNOTATED|len={len(body)}\n{body}\n")
with open(T / "raw" / "junctions.tsv", "w") as fh:
    w = csv.DictWriter(fh, jcols, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(jrows)
S = pd.DataFrame(stat, columns=["species", "gene", "exons", "status", "start_trim", "over_extended"]); S.to_csv(H / f"tpl_loci{SUF}.tsv", sep="\t", index=False)
print(S.groupby(["gene", "status"]).size().unstack(fill_value=0).to_string()); print("exon junctions moved to the RNA-supported position:", n_fixed); print("RNA loci:", (S.status == "ok").sum(), "of", len(S)); print(S[S.status == "ok"].groupby("gene").agg(loci=("exons", "size"), start_trimmed=("start_trim", lambda x: (x > 0).sum()), over_extended=("over_extended", lambda x: (x > 0).sum()), median_extra=("over_extended", lambda x: int(x[x > 0].median()) if (x > 0).any() else 0)).T.to_string())

ep.reference_dir = lambda *a, **k: T  # raw/ with the RNA loci
parts = {"cv": [], "cvrev": []}
for genus in sorted({s.split("_")[0] for s in species}):
    test = [s for s in species if s.startswith(genus + "_")]
    out = T / genus; ep.build_exon_profiles(lineage="fern", exclude=tuple(f"RNA_{s}" for s in test) + STOP_ONLY, out_dir=out / "profiles", work_root=out / "work", use_mask=False)
    shutil.rmtree(out / "work", ignore_errors=True)
    m = man[man.species.isin(test)].copy(); m["template_dir"] = [str(out / "profiles" / Path(t).name) for t in m.template_dir]
    m = m[[Path(t).exists() for t in m.template_dir]]; m.to_csv(out / "manifest.tsv", sep="\t", index=False)
    for tag, extra in (("cv", []), ("cvrev", ["--reverse-edit-penalty", REV])):
        o = str(out / tag)
        subprocess.run([TS, "scan-batch", "--manifest", str(out / "manifest.tsv"), "--genes-out", o + ".genes.tsv", "--exons-out", o + ".exons.tsv", "--junctions-out", o + ".junctions.tsv",
                        "--edits-out", o + ".edits.tsv", "--frameshifts-out", o + ".fs.tsv", "--threads", "4", "--stop-edit-penalty", "2.5"] + extra, check=True, capture_output=True)
        parts[tag].append(o)
# templates from every locus, nothing held out
ep.build_exon_profiles(lineage="fern", exclude=STOP_ONLY, out_dir=T / "all" / "profiles", work_root=T / "all" / "work", use_mask=False); shutil.rmtree(T / "all" / "work", ignore_errors=True)
for tag, ps in parts.items():
    for kind in ("genes", "exons", "edits", "fs"):
        pd.concat([pd.read_csv(f"{p}.{kind}.tsv", sep="\t") for p in ps]).to_csv(H / f"ts.{tag}{SUF}.{kind}.tsv", sep="\t", index=False)
    gg = pd.read_csv(H / f"ts.{tag}{SUF}.genes.tsv", sep="\t"); print(tag, "genes", len(gg), "complete", int((gg.complete.astype(str).str.lower() == "true").sum()))
    print(tag, "frameshift calls", len(pd.read_csv(H / f"ts.{tag}{SUF}.fs.tsv", sep="\t")), "(committed templates:", len(F), ")")
    subprocess.run([sys.executable, str(H / "splice_score.py"), str(H / f"ts.{tag}{SUF}"), {"cv": "RNA + reference templates, genus held out", "cvrev": f"same, other U>C at {REV}"}[tag] + (" (edited references only)" if SUF else "")], check=True)
