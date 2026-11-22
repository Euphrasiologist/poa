"""rebuild.py : fern editing profiles from RNA-seq editomes of the dataset's own ferns, cross-validated.

For every species with a pileup (tol.<species>.pile.tsv.gz from map.sh) and every single-exon mito gene copy poa
called there (spans from build/fern_edit/s.masked.calls.tsv; nad9 and matR from d.angio.calls.tsv):
  truth   = RNA edits in the span (C>U / U>C on the gene strand, >= MIN_FRAC of >= MIN_DEPTH reads)
  protein = the span with edits at >= 0.5 applied, translated; kept if covered >= 90%, a multiple of three,
            and free of internal stops.
Writes proteins.fasta (species|gene) and truth.tsv. Then, holding out one genus at a time, builds profiles from
the five GenBank species with editomes (rg.FERN_EDITED_MITOGENOMES) plus the other genera's RNA proteins, runs
orfedit on the held-out species and scores its calls against their RNA truth, next to the committed fern profiles
and the angiosperm ones. Writes cv.tsv (per fold and profile set), cv_genes.tsv, and final profiles from everything
in profiles_all/."""
import gzip, glob, re, shutil, subprocess, sys
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

R = Path("/lustre/scratch122/tol/teams/blaxter/users/mb39/plant_organelle_annotator"); sys.path.insert(0, str(R / "src"))
from poa.pipeline.common import pssm
from poa.pipeline.editing import lineage_profiles as lp
from poa.pipeline.trans_splicing import reference_genes as rg

H = R / "build/rnaedit"; sys.path.insert(0, str(H)); FE = R / "build/fern_edit"; ORF = str(R / "build/rev_edit/bin2/orfedit")
MIN_DEPTH, MIN_FRAC = 10, 0.2
comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
species = sorted(re.sub(r"^tol\.|\.pile\.tsv\.gz$", "", Path(p).name) for p in glob.glob(str(H / "tol.*.pile.tsv.gz")))
print("species with RNA-seq:", species)

man = pd.concat([pd.read_csv(FE / "s.masked.tsv", sep="\t"), pd.read_csv(R / "build/rev_edit/fern.fs_val_manifest.tsv", sep="\t")]).drop_duplicates(["species", "gene"])
calls = pd.read_csv(FE / "s.masked.calls.tsv", sep="\t")
extra = pd.read_csv(FE / "d.angio.calls.tsv", sep="\t")
extra = extra[extra.gene.str.split("|").str[0].isin(["nad9", "matR"])]
calls = pd.concat([calls, extra])
# assemblies replaced (gfatk re-resolved, 2026-11-19) after the gene calls these manifests come from: their coordinates
# no longer point at the genes
STALE = {l.split(".")[0] for l in open(R / "build/gfatk063/pass.txt") if l.strip()}
species = [s for s in species if s not in STALE]; print("single-exon species:", species)
calls = calls[calls.species.isin(species)].copy(); calls["g"] = calls.gene.str.split("|").str[0]
strand_of = {(r.species, r.gene): r.strand for r in man.itertuples()}
fasta_of = {r.species: r.fasta for r in man[man.organelle == "mito"].itertuples()}

truth_rows, prot = [], {}  # prot[(species, gene)] = (score, protein)
info = []
import rna
for sp in species:
    pile = rna.load_pile(H / f"tol.{sp}.pile.tsv.gz")
    genome = {r.id: str(r.seq).upper() for r in SeqIO.parse(fasta_of[sp], "fasta")}
    for r in calls[calls.species == sp].itertuples():
        contig = r.gene.split("|")[1].split(":")[0]; strand = strand_of[(sp, r.gene)]
        lo, hi = sorted((int(r.corrected_start), int(r.corrected_end)))
        if hi > len(genome[contig]) or lo < 0: continue  # runs off the contig end
        tr = rna.walk(pile, genome, contig, lo, hi, strand)
        tt = rna.truth(tr)
        for (c, p), kind in tt.items(): truth_rows.append((sp, r.gene, r.g, p, kind, 0, 0))
        d = rna.derive(tr)
        covered = sum(t[2] >= MIN_DEPTH for t in tr) / max(len(tr), 1)
        ok = d["status"] == "ok" and d["mod3"] == 0 and int(r.n_frameshifts) == 0
        info.append((sp, r.g, r.gene, hi - lo, round(covered, 2), sum(v == "C>U" for v in tt.values()), sum(v == "U>C" for v in tt.values()), d["mod3"], d["status"], int(r.n_frameshifts), ok, d["start"], d["over_extended"], d["has_stop"]))
        if ok and ((sp, r.g) not in prot or r.score > prot[(sp, r.g)][0]): prot[(sp, r.g)] = (r.score, d["protein"])
T = pd.DataFrame(truth_rows, columns=["species", "label", "gene", "pos", "edit", "frac", "depth"]); T.to_csv(H / "truth.tsv", sep="\t", index=False)
I = pd.DataFrame(info, columns=["species", "gene", "label", "len", "covered", "rna_CU", "rna_UC", "len_mod3", "status", "frameshifts", "protein_ok", "start_trim", "over_extended", "real_stop"])
I.to_csv(H / "copies.tsv", sep="\t", index=False)
with open(H / "proteins.fasta", "w") as fh:
    for (sp, g), (_, s) in sorted(prot.items()): fh.write(f">{sp}|{g}\n{s}\n")
print(f"gene copies {len(I)}, covered>=90% {(I.covered >= 0.9).sum()}, usable proteins (one per species and gene) {len(prot)}")
print(I.groupby("gene").agg(copies=("label", "size"), covered=("covered", lambda x: (x >= 0.9).sum()), ok=("protein_ok", "sum"), start_trimmed=("start_trim", lambda x: (x > 0).sum()),
      over_extended=("over_extended", lambda x: (x > 0).sum()), median_extra_codons=("over_extended", lambda x: int(x[x > 0].median()) if (x > 0).any() else 0)).T.to_string())
print(I.status.value_counts().to_string())

# ---- profiles
GENES = sorted(set(calls.g))
REF = lp.one_per_species(lp.collect_proteins(rg.FERN_EDITED_MITOGENOMES, GENES, lp.LINEAGE_SYNONYMS["fern"], (), reverse=True)[0], rg.FERN_EDITED_MITOGENOMES)

def build(out, leave_out=(), with_ref=True, min_n=4, sharp=False):
    shutil.rmtree(out, ignore_errors=True); out.mkdir(parents=True); work = out / "w"; work.mkdir()
    for g in GENES:
        seqs = {f"ref_{a}": s for a, s in REF.get(g, {}).items()} if with_ref else {}
        seqs.update({f"rna_{sp}": s for (sp, gg), (_, s) in prot.items() if gg == g and sp.split("_")[0] not in leave_out})
        if len(seqs) < min_n: continue
        # sharp: as the exon-slot profiles - 2 pseudo-observations per column, shared out by the proteins' composition
        kw = dict(background=pssm.composition(seqs.values()), pseudocount_mass=2.0) if sharp else {}
        cols = pssm.align_and_build(seqs, work, g, 0.5, **kw)
        if cols: pssm.write_profile(cols, out / f"{g}.pssm")
    shutil.rmtree(work)
    return out

E = R / "src/poa/pipeline/editing/reference"
def run(profile_dirs, sp_list, tag, rev=None):
    m = man[man.species.isin(sp_list) & (man.organelle == "mito")].copy(); m["g"] = m.gene.str.split("|").str[0]
    m = m[m.g.isin(GENES)]
    def pick(g):
        return next((str(d / f"{g}.pssm") for d in profile_dirs if (d / f"{g}.pssm").exists()), None)
    m["profile"] = m.g.map(pick); m = m[m.profile.notna()].drop(columns="g")
    mp = H / f"cv/{tag}.manifest.tsv"; mp.parent.mkdir(exist_ok=True); m.to_csv(mp, sep="\t", index=False)
    o = H / f"cv/{tag}"
    subprocess.run([ORF, "scan-batch", "--manifest", str(mp), "--gene-calls-out", f"{o}.calls.tsv", "--edits-out", f"{o}.edits.tsv", "--frameshifts-out", f"{o}.fs.tsv",
                    "--frameshift-penalty", "20", "--frameshift-flank", "0", "--stop-edit-penalty", "2.5", "--threads", "4"] + (["--reverse-edit-penalty", str(rev)] if rev else []), check=True, capture_output=True)
    return pd.read_csv(f"{o}.calls.tsv", sep="\t"), pd.read_csv(f"{o}.edits.tsv", sep="\t")

okcopies = set(zip(I[I.covered >= 0.9].species, I[I.covered >= 0.9].label))
def score(c, e, fold, name):
    rows = []
    c = c[pd.Series([(s, l) in okcopies for s, l in zip(c.species, c.gene)], index=c.index, dtype=bool)]
    for g, cg in c.groupby(c.gene.str.split("|").str[0]):
        keys = set(zip(cg.species, cg.gene))
        t = T[pd.Series([(s, l) in keys for s, l in zip(T.species, T.label)], index=T.index, dtype=bool)]; ee = e[pd.Series([(s, l) in keys for s, l in zip(e.species, e.gene)], index=e.index, dtype=bool)]
        r = {"fold": fold, "profiles": name, "gene": g, "copies": len(cg), "aligned99": int((cg.n_columns_aligned >= 0.99 * cg.n_columns_total).sum()), "frameshifts": int(cg.n_frameshifts.sum())}
        for kind, lab in (("C>U", "CU"), ("U>C", "UC")):
            tk = set(zip(t[t.edit == kind].species, t[t.edit == kind].label, t[t.edit == kind].pos))
            pk = set(zip(ee[ee.edit == kind].species, ee[ee.edit == kind].gene, ee[ee.edit == kind].genomic_pos + 1))
            r.update({f"{lab}_true": len(tk), f"{lab}_called": len(pk), f"{lab}_right": len(tk & pk)})
        rows.append(r)
    return rows

rows = []
genera = sorted({s.split("_")[0] for s in species})
for genus in genera:
    test = [s for s in species if s.startswith(genus + "_")]
    d = H / f"cv/{genus}"
    sets = {
        "angiosperm": [E / "profiles/mito", E / "profiles_refseq/mito"],
        "fern (committed)": [E / "lineages/fern/mito"],
        "RNA only": [build(d / "rna", leave_out=(genus,), with_ref=False)],
        "RNA + references": [build(d / "rnaref", leave_out=(genus,))],
        "RNA + references, sharp": [build(d / "rnarefsharp", leave_out=(genus,), sharp=True)],
        "RNA only, sharp": [build(d / "rnasharp", leave_out=(genus,), with_ref=False, sharp=True)],
    }
    # the same profiles with other U-to-C edits (not only stop removal) allowed, at a cost
    for cost in (3, 4, 5, 6, 8):
        sets[f"RNA + references, U>C at {cost}"] = (sets["RNA + references"], cost)
    sets["fern (committed), U>C at 5"] = (sets["fern (committed)"], 5)
    for cost in (3, 4, 5, 6):
        sets[f"RNA + references, sharp, U>C at {cost}"] = (sets["RNA + references, sharp"], cost)
    for name, dirs in sets.items():
        dirs, rev = dirs if isinstance(dirs, tuple) else (dirs, None)
        c, e = run(dirs, test, f"{genus}." + re.sub(r"[^A-Za-z0-9]+", "_", name), rev)
        rows += score(c, e, genus, name)
G = pd.DataFrame(rows); G.to_csv(H / "cv_genes.tsv", sep="\t", index=False)
# only genes every profile set covers in a fold
n_sets = G.groupby(["fold", "gene"]).profiles.nunique(); common = n_sets[n_sets == G.profiles.nunique()].index
Gc = G[pd.Series([x in common for x in zip(G.fold, G.gene)], index=G.index, dtype=bool)]
def summ(d):
    s = d.drop(columns=["gene"]).groupby(["fold", "profiles"], sort=False).sum(numeric_only=True).reset_index()
    for lab in ("CU", "UC"):
        s[f"{lab}_found%"] = (100 * s[f"{lab}_right"] / s[f"{lab}_true"].clip(lower=1)).round(0); s[f"{lab}_right%"] = (100 * s[f"{lab}_right"] / s[f"{lab}_called"].clip(lower=1)).round(0)
    return s
S = summ(Gc); S.to_csv(H / "cv.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
print(S[["fold", "profiles", "copies", "aligned99", "frameshifts", "CU_true", "CU_found%", "CU_right%", "UC_true", "UC_found%", "UC_right%"]].to_string(index=False))
tot = summ(Gc.assign(fold="all")); print(tot[["profiles", "copies", "aligned99", "frameshifts", "CU_true", "CU_called", "CU_found%", "CU_right%", "UC_true", "UC_called", "UC_found%", "UC_right%"]].to_string(index=False))
build(H / "profiles_all"); build(H / "profiles_all_sharp", sharp=True)
print("final:", sorted(p.stem for p in (H / "profiles_all").glob("*.pssm")))
