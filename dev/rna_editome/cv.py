"""cv.py LINEAGE SEED TAG BY SPECIES... : a lineage's RNA-derived references, held out.

For each group of SPECIES (BY = species, or genus) the lineage's exon templates and editing profiles are rebuilt
from the package's reference files without that group's RNA_<species> loci, and its species are annotated with
them: transsplice on the spliced genes, orfedit on the single-exon ones, with the options poa gives that lineage,
from the gene calls of the poa run in pilot/SEED/<species>. Results go to pilot/TAG/<species>/tables in poa's
layout, for assess.py. Genes the lineage has no profile or template for are not in them.
Run with the package's Python (PYTHONPATH=src) and orfedit/transsplice on PATH."""
import sys, shutil, subprocess, math
from pathlib import Path
import pandas as pd
R = Path("/lustre/scratch122/tol/teams/blaxter/users/mb39/plant_organelle_annotator"); H = R / "build/rnaedit"
D = Path("/lustre/scratch122/tol/teams/blaxter/users/mb39/plant_organellar_database/data")
from poa.pipeline.common import lineage as lin
from poa.pipeline.common.tools import cargo_tool
from poa.pipeline.editing import scan, lineage_profiles as lp
from poa.pipeline.trans_splicing import reconstruct as rc, exon_profiles as ep, reference_genes as rg

import os
NOREV = bool(os.environ.get("CV_PLASTID_NOREV"))  # plastid genes without other U-to-C edits (stop read-through only)
EDITP = float(os.environ.get("CV_EDIT_PENALTY", scan.DEFAULT_EDIT_PENALTY))  # cost of a C-to-U edit, both tools, both organelles
if os.environ.get("CV_PSEUDOCOUNT_MASS"):  # pseudo-observations per profile column (pssm.PSEUDOCOUNT_MASS)
    from poa.pipeline.common import pssm as _pssm; _pssm.PSEUDOCOUNT_MASS = float(os.environ["CV_PSEUDOCOUNT_MASS"])
PEDITP = float(os.environ.get("CV_PLTD_EDIT_PENALTY", 0)) or None  # cost of a C-to-U edit in plastid genes, when not EDITP
# editing sites: cost of a C-to-U edit at a profile column not edited in the references (mito, plastid); genera needed for a site
NOVEL = float(os.environ.get("CV_NOVEL", 0)) or None; PNOVEL = float(os.environ.get("CV_PLTD_NOVEL", 0)) or NOVEL
KNOWN = float(os.environ.get("CV_KNOWN", 0)) or None; PKNOWN = float(os.environ.get("CV_PLTD_KNOWN", 0)) or KNOWN  # and at one that is
if os.environ.get("CV_SITE_GROUPS"):
    from poa.pipeline.common import pssm as _pssm2; _pssm2.EDIT_SITE_MIN_GROUPS = int(os.environ["CV_SITE_GROUPS"])
NOPLASTID = bool(os.environ.get("CV_NO_PLASTID"))  # mito only (GenBank test genomes have no plastid here)
def same_group(organism, grp):  # a RefSeq reference genome of the group's genus (or species)
    w = organism.split(); return (w[0] if by != "species" else "_".join(w[:2])) == grp
lineage, seed, tag, by, species = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5:]
orfedit, transsplice = cargo_tool("orfedit"), cargo_tool("transsplice")
ACC = dict(l.split() for l in open(H / "extra_genomes" / "accessions.tsv")) if (H / "extra_genomes" / "accessions.tsv").exists() else {}
groups = {}
for sp in species: groups.setdefault(sp if by == "species" else sp.split("_")[0], []).append(sp)
TS = lin.trans_splicing_reference(Path(os.environ.get("CV_SRC", R / "src")) / "poa/pipeline", lineage)  # CV_SRC: another copy of the package (with PYTHONPATH)
mito_names = set(n for by_name in lp.rna_proteins(lineage).values() for n in by_name)
have = {p.stem for d in (TS / "raw_rna", TS / "plastid" / "raw_rna") for p in d.glob("RNA_*.fasta")} | mito_names | set(
    n for by_name in lp.rna_proteins(lineage, organelle="pltd").values() for n in by_name)
for grp, sps in groups.items():
    W = H / "pilot" / tag / f"_work.{grp}.{os.environ.get('LSB_JOBID', 'x')}"; shutil.rmtree(W, ignore_errors=True); W.mkdir(parents=True)  # (two jobs can hold species of one genus)
    excl = tuple(f"RNA_{s}" for s in species if (s if by == "species" else s.split("_")[0]) == grp and f"RNA_{s}" in have)
    # every dataset species of the group is left out, whether or not it was asked for
    excl = tuple(sorted(n for n in have if (n[4:] if by == "species" else n[4:].split("_")[0]) == grp))
    # a GenBank reference genome annotated as a test (extra_genomes/accessions.tsv): its own accession is left out
    excl += tuple(ACC[s] for s in sps if s in ACC)
    # and the lineage's RefSeq reference genomes of the group, mito and plastid
    for refs in (rg.LINEAGE_MITOGENOMES.get(lineage, {}), rg.LINEAGE_PLASTOMES.get(lineage, {})):
        excl += tuple(a for a, organism in refs.items() if same_group(organism, grp) and a not in excl)
    tpl = ep.build_exon_profiles(lineage, None, W, "mito", excl, W / "tpl")
    gb = rg.LINEAGE_EDITED_MITOGENOMES.get(lineage, rg.LINEAGE_MITOGENOMES[lineage])
    ed = lp.build_lineage_profiles(lineage, W, exclude=tuple(e for e in excl if e in mito_names or e in gb or (lineage == lin.DEFAULT_LINEAGE and e in lp.SPLIT_CCMFN_MITOGENOMES)), out_dir=W / "ed" / "mito")
    spliced = [g for g in rg.mito_spliced(lineage) if (tpl / g / "whole_gene.pssm").exists()]
    plastid = (lineage in rg.LINEAGE_PLASTID or lineage == lin.DEFAULT_LINEAGE) and not NOPLASTID
    if plastid:
        ptpl = ep.build_exon_profiles(lineage, None, W, "pltd", excl, W / "ptpl")
        ped = lp.build_plastid_profiles(lineage, W, exclude=excl, out_dir=W / "ed" / "pltd")
        pspliced = [g for g in rg.PLASTID_SPLICED if (ptpl / g / "whole_gene.pssm").exists()]
    for sp in sps:
        src = H / "pilot" / seed / sp / "tables"; out = H / "pilot" / tag / sp / "tables"; out.mkdir(parents=True, exist_ok=True)
        calls = pd.read_csv(src / "gene_calls.tsv", sep="\t"); shutil.copy(src / "gene_calls.tsv", out / "gene_calls.tsv")
        x = H / "extra_genomes" / f"{sp}.fasta"
        fasta = str(x) if x.exists() and not os.environ.get("USE_DATASET") else str(sorted((D / "mito" / sp).glob("*.mito.ctg.fasta"))[0])
        mito = calls[calls.organelle == "mito"]
        rows = [dict(species=sp, organelle="mito", gene=r.gene, fasta=fasta, template_dir=str(tpl / r.gene), contig=r.contig_id, start=int(r.start), end=int(r.end), strand=r.strand)
                for r in mito[mito.gene.isin(spliced)].itertuples()]
        rc.write_manifest(rows, W / f"{sp}.ts.tsv")
        outs = {n: out / f"{n}.tsv" for n in rc.OUT_NAMES}
        cmd = [transsplice, "scan-batch", "--manifest", str(W / f"{sp}.ts.tsv"), "--genes-out", str(outs["reconstructed_genes"]), "--exons-out", str(outs["reconstructed_exons"]),
               "--junctions-out", str(outs["reconstructed_junctions"]), "--edits-out", str(outs["reconstructed_edits"]), "--threads", "2",
               "--merge-distance", str(rc.DEFAULT_MERGE_DISTANCE), "--junction-distance-threshold", str(rc.DEFAULT_JUNCTION_DISTANCE_THRESHOLD), "--exon-flank", str(rc.DEFAULT_EXON_FLANK),
               "--frameshift-penalty", str(rc.DEFAULT_FRAMESHIFT_PENALTY), "--frameshifts-out", str(outs["reconstructed_frameshifts"]), "--edit-penalty", str(EDITP),
               "--internal-stop-penalty", str(scan.DEFAULT_INTERNAL_STOP_PENALTY)]
        for opt, v in (("--stop-edit-penalty", lin.stop_edit_penalty(lineage)), ("--reverse-edit-penalty", lin.reverse_edit_penalty(lineage)), ("--junction-slide", lin.junction_slide(lineage))):
            if v is not None: cmd += [opt, str(v)]
        if NOVEL: cmd += ["--novel-edit-penalty", str(NOVEL)]
        if KNOWN: cmd += ["--known-edit-penalty", str(KNOWN)]
        if rows: subprocess.run(cmd, check=True, stderr=subprocess.DEVNULL)
        if plastid:
            # plastid spliced genes on the held-out plastid templates: a run of their own (no junction slide), tables appended
            pf = str(sorted((D / "plastid" / sp).glob("*.pltd.ctg.fasta"))[0]); pl = calls[calls.organelle == "pltd"]
            prow = [dict(species=sp, organelle="pltd", gene=r.gene, fasta=pf, template_dir=str(ptpl / r.gene), contig=r.contig_id, start=int(r.start), end=int(r.end), strand=r.strand)
                    for r in pl[pl.gene.isin(pspliced)].itertuples()]
            rc.write_manifest(prow, W / f"{sp}.pts.tsv"); pouts = {n: W / f"{sp}.p.{n}.tsv" for n in rc.OUT_NAMES}
            pcmd = [transsplice, "scan-batch", "--manifest", str(W / f"{sp}.pts.tsv"), "--genes-out", str(pouts["reconstructed_genes"]), "--exons-out", str(pouts["reconstructed_exons"]),
                    "--junctions-out", str(pouts["reconstructed_junctions"]), "--edits-out", str(pouts["reconstructed_edits"]), "--threads", "2",
                    "--merge-distance", str(rc.DEFAULT_MERGE_DISTANCE), "--junction-distance-threshold", str(rc.DEFAULT_JUNCTION_DISTANCE_THRESHOLD), "--exon-flank", str(rc.DEFAULT_EXON_FLANK),
                    "--frameshift-penalty", str(rc.DEFAULT_FRAMESHIFT_PENALTY), "--frameshifts-out", str(pouts["reconstructed_frameshifts"]), "--edit-penalty", str(PEDITP or EDITP),
                    "--internal-stop-penalty", str(scan.DEFAULT_INTERNAL_STOP_PENALTY)]
            for opt, v in (("--stop-edit-penalty", lin.stop_edit_penalty(lineage)), ("--reverse-edit-penalty", None if NOREV else lin.reverse_edit_penalty(lineage))):
                if v is not None: pcmd += [opt, str(v)]
            if PNOVEL: pcmd += ["--novel-edit-penalty", str(PNOVEL)]
            if PKNOWN: pcmd += ["--known-edit-penalty", str(PKNOWN)]
            subprocess.run(pcmd, check=True, stderr=subprocess.DEVNULL)
            for n in rc.OUT_NAMES: pd.concat([pd.read_csv(outs[n], sep="\t", dtype=str), pd.read_csv(pouts[n], sep="\t", dtype=str)]).to_csv(outs[n], sep="\t", index=False)
        # single-exon genes on the held-out profiles (not the lineage's spliced ones)
        hits = scan.hits_to_correct(calls, "mito", sorted(p.stem for p in ed.glob("*.pssm") if p.stem not in spliced))
        erows = [dict(species=sp, organelle="mito", gene=scan.copy_label(r.gene, r.contig_id, int(r.start), int(r.end)), fasta=fasta, profile=str(ed / f"{r.gene}.pssm"),
                      contig=r.contig_id, start=int(r.start), end=int(r.end), strand=r.strand) for r in hits.itertuples()]
        res = scan.run_orfedit(orfedit, erows, sp, W, 2, scan.DEFAULT_FLANK, EDITP, scan.DEFAULT_FRAMESHIFT_PENALTY,
                               stop_edit_penalty=lin.stop_edit_penalty(lineage), reverse_edit_penalty=lin.reverse_edit_penalty(lineage), novel_edit_penalty=NOVEL, known_edit_penalty=KNOWN)
        if plastid:
            ph = scan.hits_to_correct(calls, "pltd", sorted(p.stem for p in ped.glob("*.pssm")))
            prows = [dict(species=sp, organelle="pltd", gene=scan.copy_label(r.gene, r.contig_id, int(r.start), int(r.end)), fasta=pf, profile=str(ped / f"{r.gene}.pssm"),
                          contig=r.contig_id, start=int(r.start), end=int(r.end), strand=r.strand) for r in ph.itertuples()]
            pres = scan.run_orfedit(orfedit, prows, sp + ".pltd", W, 2, scan.DEFAULT_FLANK, PEDITP or EDITP, scan.DEFAULT_FRAMESHIFT_PENALTY,
                                    stop_edit_penalty=lin.stop_edit_penalty(lineage), reverse_edit_penalty=None if NOREV else lin.reverse_edit_penalty(lineage), novel_edit_penalty=PNOVEL, known_edit_penalty=PKNOWN)
            res = tuple(pd.concat([a, b], ignore_index=True) for a, b in zip(res, pres))
        for df, name in zip(res, ("edited_gene_calls", "predicted_edits", "frameshifts")): df.to_csv(out / f"{name}.tsv", sep="\t", index=False)
        print(f"{grp} {sp}: held out {len(excl)} reference species; {len(rows)} spliced-gene hits, {len(erows)} single-exon copies", flush=True)
    shutil.rmtree(W, ignore_errors=True)
