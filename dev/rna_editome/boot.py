"""boot.py LINEAGE MODE SPECIES... : gene models of mito protein genes from RNA-seq, seeded by a poa run.

For a lineage without annotated references. poa (pilot/LINEAGE/SPECIES, any lineage's references) says roughly where
each gene is; the species' RNA-seq says what the transcript is:
- introns are the junctions intron-spanning reads show (bwa split alignments), kept with >= MINJ reads that are at
  least MIN_SPLICED of the reads at the donor;
- the transcript is walked from the middle of the first called exon in both directions, jumping every such intron
  (poa's own junction is used only where the reads show none within 60 bp, e.g. around a 22 bp exon no read piece
  aligns to), 1.5 kb past the called ends;
- RNA edits are applied (rna.py thresholds; MODE cu = C-to-U only, both = U-to-C as well, with stop read-through);
- the gene is the open reading frame holding most of the called bases, from its ATG nearest the called start (the
  first ATG of the frame when poa's call did not reach the start of its profile) to its stop.
Pieces on different contigs or strands (trans-splicing) are joined at poa's boundaries.
Writes boot/LINEAGE/SPECIES.genes.tsv and SPECIES.fasta (exons + edited protein, the raw/ format of trans_splicing)."""
import sys, re, subprocess, glob, gzip
from collections import Counter, defaultdict
from pathlib import Path
import pandas as pd
H = Path(__file__).resolve().parent; sys.path.insert(0, str(H)); import rna
D = Path("/lustre/scratch122/tol/teams/blaxter/users/mb39/plant_organellar_database/data")
MINJ, MIN_SPLICED, REACH, NEAR = 3, 0.1, 1500, 60
REF = {}
SPLICED = ["nad1", "nad2", "nad5", "rps3", "nad4", "nad7", "ccmFc", "cox2", "rpl2", "rps10"]
CIG = re.compile(r"(\d+)([MIDNSHP=X])")
comp = rna.comp

def read_fasta(path):
    out, name = {}, None
    for line in open(path):
        if line.startswith(">"): name = line[1:].split()[0]; out[name] = []
        else: out[name].append(line.strip().upper())
    return {k: "".join(v) for k, v in out.items()}

def read_fasta_full(path):
    out, name = {}, None
    for line in open(path):
        if line.startswith(">"): name = line[1:].strip(); out[name] = ""
        else: out[name] += line.strip()
    return out

def load_piles(sp):
    """summed base counts of every run mapped to the species (tol.SP.pile.tsv.gz, pub.SP.RUN.pile.tsv.gz)"""
    pile = {}
    for f in sorted(glob.glob(str(H / f"tol.{sp}.pile.tsv.gz")) + glob.glob(str(H / f"pub.{sp}.*[0-9].pile.tsv.gz"))):
        with gzip.open(f, "rt") as fh:
            for line in fh:
                c, p, ref, n, a, cc, g, t = line.split("\t"); k = (c, int(p))
                if k in pile:
                    o = pile[k][2]; o["A"] += int(a); o["C"] += int(cc); o["G"] += int(g); o["T"] += int(t); pile[k] = (ref, pile[k][1] + int(n), o)
                else: pile[k] = (ref, int(n), {"A": int(a), "C": int(cc), "G": int(g), "T": int(t)})
    return pile

def genome_fasta(sp):
    """the dataset's mito contigs, or extra_genomes/SP.fasta (a GenBank genome)"""
    x = H / "extra_genomes" / f"{sp}.fasta"
    return str(x) if x.exists() else sorted(glob.glob(str(D / "mito" / sp / "*.mito.ctg.fasta")))[0]

def bams(sp): return sorted(glob.glob(str(H / f"tol.{sp}.mito.bam")) + glob.glob(str(H / f"pub.{sp}.*.mito.bam")))

def span(pos, cigar):
    ops = [(int(n), o) for n, o in CIG.findall(cigar)]
    ref = sum(n for n, o in ops if o in "MDN=X"); q = sum(n for n, o in ops if o in "MI=X")
    left = ops[0][0] if ops[0][1] in "SH" else 0
    return pos, pos + ref - 1, left, left + q

def norm(g, a, b):
    """leftmost of the junctions (last base before, first base after; 1-based) that join the same transcript"""
    while a > 1 and g[a - 1] == g[b - 2]: a -= 1; b -= 1
    return a, b

def introns(sp, genome, pile):
    """[(contig, last_before, first_after, reads, exact)] - junctions left-normalised; reads whose pieces leave
    unaligned bases between them (an edited base next to the junction is clipped) count towards an exact junction
    they are compatible with"""
    exact, loose = Counter(), Counter()
    for bam in bams(sp):
        p = subprocess.Popen(["samtools", "view", "-F", "0x900", bam], stdout=subprocess.PIPE, text=True)
        for line in p.stdout:
            if "\tSA:Z:" not in line: continue
            f = line.split("\t", 12); flag = int(f[1]); contig = f[2]; rev = bool(flag & 16)
            s1, e1, qs1, qe1 = span(int(f[3]), f[5])
            sa = line.split("\tSA:Z:", 1)[1].split("\t", 1)[0].rstrip("\n")
            for part in sa.rstrip(";").split(";"):
                c2, pos2, st2, cig2 = part.split(",")[:4]
                if c2 != contig or (st2 == "-") != rev: continue
                s2, e2, qs2, qe2 = span(int(pos2), cig2)
                (ls, le, lqs, lqe), (rs, re_, rqs, rqe) = sorted([(s1, e1, qs1, qe1), (s2, e2, qs2, qe2)])
                gap = rqs - lqe
                if not (lqs < rqs and rs - le > 20 and abs(gap) <= 15): continue
                c = contig.replace("mito_", ""); g = genome[c]
                if gap <= 0: exact[(c, *norm(g, le + gap, rs))] += 1  # shared bases given to the right piece
                else: loose[(c, le, rs, gap)] += 1
        p.wait()
    out = dict(exact); inexact = Counter()
    for (c, le, rs, gap), n in loose.items():
        g = genome[c]; hits = [k for t in range(gap + 1) if (k := (c, *norm(g, le + t, rs - (gap - t)))) in exact]
        if hits: out[max(hits, key=lambda k: exact[k])] += n
        else: inexact[(c, le, rs, gap)] += n
    res = [(c, a, b, n, True) for (c, a, b), n in out.items()]
    # junctions only reads with clipped bases show: the unaligned bases go to the left exon
    seen = Counter()
    for (c, le, rs, gap), n in inexact.items(): seen[(c, *norm(genome[c], le + gap, rs))] += n
    res += [(c, a, b, n, False) for (c, a, b), n in seen.items()]
    keep = []
    for c, a, b, n, ex in res:
        if n < MINJ: continue
        # reads running on into the intron at either end
        un = max(pile.get((c, a + 3), ("", 0))[1], pile.get((c, b - 3), ("", 0))[1])
        if n / (n + un) >= MIN_SPLICED: keep.append((c, a, b, n, ex))
    return keep

def edited_codons(tr, frame, both):
    """codons of the edited transcript in `frame`: [(codon, is_stop)]; MODE both reads a stop through on U>C evidence"""
    seq = []
    for _, b, n, na in tr:
        ok = b in rna.FLIP and na >= 2 and n and na / n >= rna.MIN_FRAC and (both or b == "C")
        seq.append(rna.FLIP[b] if ok else b)
    out = []
    for i in range(frame, len(seq) - 2, 3):
        cod = seq[i:i + 3]
        if rna.CODE.get("".join(cod), "X") == "*" and both:
            for j in range(3):
                _, b, n, na = tr[i + j]
                if cod[j] == "T" and b == "T" and na >= 2 and n and na / n >= rna.STOP_FRAC:
                    trial = cod[:]; trial[j] = "C"
                    if rna.CODE.get("".join(trial), "X") != "*": cod = trial; break
        out.append("".join(cod))
    return out

def model(sp, gene, seeds, genome, pile, J, both, reached_start):
    """seeds: called exons in transcript order [(contig, lo, hi, strand)] (0-based half-open)"""
    # consecutive seeds on one contig and strand, in transcript order, form a piece
    pieces = []
    for s in seeds:
        if pieces and pieces[-1][-1][0] == s[0] and pieces[-1][-1][3] == s[3] and ((s[1] >= pieces[-1][-1][1]) == (s[3] == "+")): pieces[-1].append(s)
        else: pieces.append([s])
    tr, seedpos, junc = [], set(), []
    for pi, piece in enumerate(pieces):
        c, strand = piece[0][0], piece[0][3]; g = genome[c]; L = len(g); step = 1 if strand == "+" else -1
        def tx(a, b): return (a, b) if strand == "+" else (b, a)  # genomic (last before, first after) -> (donor, acceptor)
        def after(p, k=40):  # k transcript bases starting at p
            s_ = g[p - 1:p - 1 + k] if strand == "+" else "".join(comp[x] for x in reversed(g[max(p - k, 0):p]))
            return s_
        def before(p, k=40):  # k transcript bases ending at p
            s_ = g[max(p - k, 0):p] if strand == "+" else "".join(comp[x] for x in reversed(g[p - 1:p - 1 + k]))
            return s_
        predicted = []
        for x, y in zip(piece, piece[1:]):
            d, ac = (x[2], y[1] + 1) if strand == "+" else (x[1] + 1, y[2])
            if abs(ac - d) > 1: predicted.append((d, ac))
        def in_seed(p): return any(s_[1] + 1 - NEAR <= p <= s_[2] + NEAR for s_ in piece)
        cand = defaultdict(list)
        for cc, a, b, n, ex in (j for j in J if j[0] == c):
            d, ac = tx(a, b)
            # few reads only count next to a junction poa called
            if n < 5 and not any(abs(d - pd) <= NEAR and abs(ac - pa) <= NEAR for pd, pa in predicted): continue
            cand[d].append((ac, n, ex))
        don, acc, hops = {}, {}, []
        for d, lst in cand.items():
            top = max(n for _, n, _ in lst); lst = [t for t in lst if t[1] >= 0.3 * top]
            # a donor read into several acceptors (repeated exons): the one in a called exon, else the nearest
            ac, n, ex = min(lst, key=lambda t: (not in_seed(t[0]), abs(t[0] - d)))
            don[d] = (ac, n, ex)
        for d, (ac, n, ex) in sorted(don.items(), key=lambda kv: (in_seed(kv[0]), -abs(kv[1][0] - kv[0]), kv[1][1])): acc[ac] = (d, n, ex)
        don = {d: (ac, n, ex) for ac, (d, n, ex) in acc.items()}
        # a junction the reads place in another copy of the exon (repeat): moved to the copy poa called
        for pd, pa in predicted:
            for d in [d for d in don if abs(d - pd) <= NEAR]:
                ac, n, ex = don[d]
                if abs(ac - pa) <= NEAR: continue
                hit = [q for q in range(pa - NEAR, pa + NEAR + 1) if 1 <= q <= L and after(q) == after(ac)]
                if hit: acc.pop(ac, None); don[d] = (hit[0], n, ex); acc[hit[0]] = (d, n, ex); continue
                # the copy the reads join starts differently and then runs into a repeat shared with the called exon
                # (the genome recombines across it): follow it to the repeat and cross to the called copy there
                y = next(s_ for s_ in piece if (s_[1] + 1 if strand == "-" else s_[2]) and abs((s_[1] + 1 if strand == "+" else s_[2]) - pa) <= 1)
                A, B = after(ac, 1500), after(y[1] + 1 if strand == "+" else y[2], y[2] - y[1])
                from Bio import Align
                al = Align.PairwiseAligner(); al.mode = "local"; al.match_score = 1; al.mismatch_score = -3; al.open_gap_score = -6; al.extend_gap_score = -2
                x = al.align(A, B)[0]; (ia, ea), (ib, eb) = x.aligned[0][0], x.aligned[1][0]
                if ea - ia >= 100 and ia > 0:
                    pA, pB = ac + step * ia, (y[1] + 1 if strand == "+" else y[2]) + step * ib
                    hops.append((pA - step, pB))
            for ac in [a_ for a_ in acc if abs(a_ - pa) <= NEAR]:
                d, n, ex = acc[ac]
                if abs(d - pd) <= NEAR: continue
                hit = [q for q in range(pd - NEAR, pd + NEAR + 1) if 1 <= q <= L and before(q) == before(d)]
                if hit: don.pop(d, None); acc[ac] = (hit[0], n, ex); don[hit[0]] = (ac, n, ex)
        # poa's junction where the reads show none nearby
        for d, ac in predicted:
            if not any(abs(k - d) <= NEAR for k in don) and not any(abs(k - ac) <= NEAR for k in acc): don[d] = (ac, 0, False); acc[ac] = (d, 0, False)
        for d, ac in hops: don[d] = (ac, -1, False); acc[ac] = (d, -1, False)
        first, last = piece[0], piece[-1]
        anchor = (first[1] + first[2]) // 2 + 1
        end_seed = (last[2] if strand == "+" else last[1] + 1); start_seed = (first[1] + 1 if strand == "+" else first[2])
        fwd, p, past = [], anchor, None
        while 1 <= p <= L:
            fwd.append(p)
            if p == end_seed: past = 0
            if past is not None:
                past += 1
                if past > (REACH if pi == len(pieces) - 1 else 0): break
            if p in don: junc.append((gene, c, *sorted((p, don[p][0])), don[p][1], don[p][2])); p = don[p][0]
            else: p += step
            if len(fwd) > 30000: break
        back, p, past = [], anchor, None
        while True:
            if p == start_seed: past = 0
            if past is not None:
                past += 1
                if past > (REACH if pi == 0 else 0): break
            if p in acc: junc.append((gene, c, *sorted((p, acc[p][0])), acc[p][1], acc[p][2])); p = acc[p][0]
            else: p -= step
            if not 1 <= p <= L or len(back) > 30000: break
            back.append(p)
        path = back[::-1] + fwd
        for p in path:
            ref, n, cnt = pile.get((c, p), (g[p - 1], 0, None))
            sref = ref if strand == "+" else comp[ref]; alt = rna.FLIP.get(sref)
            na = cnt[alt if strand == "+" else comp[alt]] if (alt and cnt) else 0
            tr.append(((c, p), sref, n, na, strand))
        for s in piece: seedpos |= {(c, p) for p in range(s[1] + 1, s[2] + 1)}
    t4 = [t[:4] for t in tr]
    best = None
    for f in range(3):
        cods = edited_codons(t4, f, both); i = 0
        while i < len(cods):
            j = i
            while j < len(cods) and rna.CODE.get(cods[j], "X") != "*": j += 1
            n_seed = sum(tr[f + 3 * k + o][0] in seedpos for k in range(i, j) for o in range(3))
            if best is None or n_seed > best[0]: best = (n_seed, f, i, j, cods)
            i = j + 1
    n_seed, f, i, j, cods = best
    has_stop = j < len(cods)
    # start codon: the ATG nearest the called start when poa's call reached the start of its profile, else the most
    # upstream ATG of the frame before the called bases; with neither, the first called codon (non-ATG start)
    atgs = [k for k in range(i, j) if cods[k] == "ATG"]
    seed_first = next((k for k in range(i, j) if any(tr[f + 3 * k + o][0] in seedpos for o in range(3))), i)
    near = [k for k in atgs if -30 <= k - seed_first <= 5]
    n_all = sum(tr[f + 3 * k + o][0] in seedpos for k in range(i, j) for o in range(3))
    def kept(k0): return sum(tr[f + 3 * k + o][0] in seedpos for k in range(k0, j) for o in range(3)) / max(n_all, 1)
    # where a homologue's first residue falls in this frame (codon index), when the two align well
    est = None
    h = homology("".join(rna.CODE.get(c_, "X") for c_ in cods[i:j]), gene, full=True)
    if h and h[0] >= 0.4 and h[5] >= 0.3 and h[1] <= 10: est = i + h[3] - h[1]  # aligned from the homologue's start
    by_hom = [k for k in atgs if est is not None and abs(k - est) <= 15]
    if by_hom: start, note = min(by_hom, key=lambda k: abs(k - est)), "ATG by homology"
    elif est is not None and est >= i: start, note = est, "non-ATG start (homology)"
    elif reached_start and near: start, note = min(near, key=lambda k: abs(k - seed_first)), ""
    elif not reached_start and atgs and kept(atgs[0]) >= 0.9: start, note = atgs[0], "first ATG"
    else: start, note = seed_first, "non-ATG start"
    lo, hi = f + 3 * start, f + 3 * j + (3 if has_stop else 0)
    cds = tr[lo:hi]
    prot = "".join(rna.CODE.get(c, "X") for c in cods[start:j])
    exons = []
    for t in cds:
        (c, p), strand = t[0], t[4]
        if exons and exons[-1][0] == c and exons[-1][3] == strand and p == exons[-1][2] + (1 if strand == "+" else -1): exons[-1][2] = p
        else: exons.append([c, p, p, strand])
    body = cds[:3 * (j - start)]
    cov = sum(t[2] >= rna.MIN_DEPTH for t in body) / max(len(body), 1)
    tru = rna.truth([t[:4] for t in body])
    used = {(c, a, b): (n, ex) for (_, c, a, b, n, ex) in junc}
    jn = []
    for x, y in zip(exons, exons[1:]):
        if x[0] != y[0] or x[3] != y[3]: jn.append("trans"); continue
        a, b = sorted((x[2], y[1])); n, ex = used.get((x[0], a, b), (None, None))
        jn.append("?" if n is None else ("repeat" if n == -1 else "poa" if n == 0 else f"{n}{'' if ex else '~'}"))
    status = "ok"; hfin = homology(prot, gene)  # a call across an intron holds little of the frame: homology decides
    if not has_stop: status = "no stop"
    elif n_seed < 0.7 * len(seedpos) and not (hfin and hfin[0] >= 0.6 and hfin[1] + hfin[2] <= 0.15 * len(REF[gene])): status = "frame holds < 70% of the call"
    elif cov < 0.9: status = "low RNA coverage"
    return dict(species=sp, gene=gene, status=status, note=note, n_exons=len(exons), aa=len(prot), covered=round(cov, 2), rna_CU=sum(v == "C>U" for v in tru.values()),
                rna_UC=sum(v == "U>C" for v in tru.values()), seed_in_frame=round(n_seed / max(len(seedpos), 1), 2), junction_reads=",".join(jn), jn=jn,
                exons=";".join(f"{c}:{min(a, b)}-{max(a, b)}{s}" for c, a, b, s in exons), protein=prot), exons

def load_ref():
    """Marchantia polymorpha mito proteins (NC_001660; unedited, complete): {gene: protein}"""
    out, name = {}, None
    for line in open(H / "boot/ref/NC_001660.cds_aa.fasta"):
        if line.startswith(">"):
            m = re.search(r"gene=([^\]]+)", line); name = m.group(1) if m else None
            name = {"tatC": "mttB"}.get(name, name)
            if name in out: name = None  # the three ccmF pieces: keep the first
            if name: out[name] = ""
        elif name: out[name] += line.strip()
    return out

def homology(prot, gene, full=False):
    """(identity over aligned pairs, reference residues before / after the alignment, own residues before / after,
    fraction of the reference aligned)"""
    ref = REF.get(gene)
    if not ref or not prot: return None
    from Bio import Align
    from Bio.Align import substitution_matrices
    al = Align.PairwiseAligner(); al.substitution_matrix = substitution_matrices.load("BLOSUM62"); al.open_gap_score = -11; al.extend_gap_score = -1; al.mode = "local"
    a = al.align(prot.replace("X", "A").replace("*", "A"), ref)[0]
    q, r = a.aligned; ident = sum(prot[i] == ref[k] for (qs, qe), (rs, re_) in zip(q, r) for i, k in zip(range(qs, qe), range(rs, re_)))
    n = sum(qe - qs for qs, qe in q)
    return (round(ident / max(n, 1), 2), int(r[0][0]), len(ref) - int(r[-1][1]), int(q[0][0]), len(prot) - int(q[-1][1]), n / len(ref))

def main():
    lin, mode, species = sys.argv[1], sys.argv[2], sys.argv[3:]
    REF.update(load_ref())
    extra = H / f"boot/ref/{lin}_extra.fasta"  # the lineage's own proteins where it has some
    if extra.exists(): REF.update({k.split()[0]: v for k, v in read_fasta_full(extra).items()})
    out = H / "boot" / lin; out.mkdir(parents=True, exist_ok=True)
    for sp in species:
        T = H / "pilot" / lin / sp / "tables"
        genome = read_fasta(genome_fasta(sp)); pile = load_piles(sp)
        if not pile or not T.exists(): print("skip", sp); continue
        J = introns(sp, genome, pile)
        gc = pd.read_csv(T / "gene_calls.tsv", sep="\t"); gc = gc[gc.organelle == "mito"]
        E = pd.read_csv(T / "edited_gene_calls.tsv", sep="\t"); E = E[(E.organelle == "mito") & ~E.gene.isin(SPLICED)]
        X = pd.read_csv(T / "reconstructed_exons.tsv", sep="\t"); X = X[X.organelle == "mito"]
        loci = []
        for r in E.sort_values("score", ascending=False).itertuples():
            m = gc[(gc.gene == r.gene) & (gc.contig_id == r.contig) & (gc.start < r.hit_end + 300) & (gc.end > r.hit_start - 300)]
            if m.empty: continue
            lo, hi = sorted((int(r.corrected_start), int(r.corrected_end))); lo, hi = max(lo, 0), min(hi, len(genome[r.contig]))
            loci.append((r.gene, [(r.contig, lo, hi, m.iloc[0].strand)], r.n_columns_aligned >= 0.95 * r.n_columns_total, 1))
        for gene, xs in X.groupby("gene"):
            xs = xs.sort_values("slot")
            seeds = [(x.contig, max(min(int(x.start), int(x.end)), 0), min(max(int(x.start), int(x.end)), len(genome[x.contig])), x.strand) for x in xs.itertuples()]
            loci.append((gene, seeds, int(xs.slot.iloc[0]) == 1, len(seeds)))
        rows, seen = [], set()
        with open(out / f"{sp}.fasta", "w") as fh:
            for gene, seeds, reached, n_poa in loci:
                r, exons = model(sp, gene, seeds, genome, pile, J, mode == "both", reached); jn = r.pop("jn")
                if (gene, r["exons"]) in seen: continue
                seen.add((gene, r["exons"])); k = sum(1 for g, _ in seen if g == gene); r["locus"] = f"{gene}_{k}"; r["poa_exons"] = n_poa
                h = homology(r["protein"], gene); r["march_id"], r["march_missing"], r["own_extra"] = (h[0], f"{h[1]}/{h[2]}", f"{h[3]}/{h[4]}") if h else ("", "", "")
                rows.append(r)
                # exon records: pieces either side of a repeat crossing are one exon
                recs = []
                for i, (c, a, b, s) in enumerate(exons):
                    lo, hi = min(a, b) - 1, max(a, b); raw = genome[c][lo:hi]
                    if s == "-": raw = "".join(comp[x] for x in reversed(raw))
                    if i and jn[i - 1] == "repeat": recs[-1][0] += raw; recs[-1][1:] = [lo, hi, s, c] if hi - lo > recs[-1][2] - recs[-1][1] else recs[-1][1:]
                    else: recs.append([raw, lo, hi, s, c])
                r["n_exons"] = len(recs)
                for i, (raw, lo, hi, s, c) in enumerate(recs, 1): fh.write(f">RNA_{sp}|{gene}|{gene}_{k}|exon{i}of{len(recs)}|{lo + 1}-{hi}|{s}|{c}\n{raw}\n")
                fh.write(f">RNA_{sp}|{gene}|{gene}_{k}|translation_REFSEQ_ANNOTATED|len={len(r['protein'])}\n{r['protein']}\n")
        G = pd.DataFrame(rows); G.to_csv(out / f"{sp}.genes.tsv", sep="\t", index=False)
        pd.DataFrame(J, columns=["contig", "last_before", "first_after", "reads", "exact"]).sort_values(["contig", "last_before"]).to_csv(out / f"{sp}.introns.tsv", sep="\t", index=False)
        print(f"{sp}: {len(G)} loci, {(G.status == 'ok').sum()} ok; {len(J)} introns")

if __name__ == "__main__": main()
