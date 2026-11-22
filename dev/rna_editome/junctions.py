"""junctions.py PREFIX [LABEL] : exon boundaries of transsplice's spliced mito genes against intron-spanning RNA reads.
Reads tol.<species>.mito.bam through `samtools view` (bwa-mem2 alignments; a read across an intron is a primary
alignment plus a supplementary one, SA tag). An observed junction is (last base before the intron, first base after
it) in genomic order, same contig and strand. For each pair of consecutive exons on one contig and strand the
predicted junction is compared with the observed ones nearby; bwa places the break anywhere within a repeat shared
by the two sides, so a junction that differs by the same shift at both ends (<= 3 bp) counts as the same.
Writes PREFIX.junction_check.tsv."""
import re, subprocess, sys, glob
from collections import Counter
from pathlib import Path
import pandas as pd
H = Path(__file__).resolve().parent; pre = sys.argv[1]; label = sys.argv[2] if len(sys.argv) > 2 else pre
CIG = re.compile(r"(\d+)([MIDNSHP=X])")
def span(pos, cigar):
    """(ref start 1-based, ref end, query start, query end) - query coordinates in the orientation the SAM line is
    written in, hard and soft clips counted"""
    ops = [(int(n), o) for n, o in CIG.findall(cigar)]
    ref = sum(n for n, o in ops if o in "MDN=X"); q = sum(n for n, o in ops if o in "MI=X")
    left = ops[0][0] if ops[0][1] in "SH" else 0
    return pos, pos + ref - 1, left, left + q
def agrees(pred, le, rs, gap):
    """does the predicted (last base before, first base after) fit a split read whose left piece ends at le and
    right piece starts at rs, with `gap` read bases between the pieces (negative: the pieces share -gap bases)?
    Unaligned bases between the pieces (an edited base at the exon end is clipped, not mismatched) belong to
    either exon; shared bases to either as well."""
    a, b = pred[0] - le, rs - pred[1]
    if gap >= 0: return a >= 0 and b >= 0 and a + b == gap
    return a <= 0 and b <= 0 and a + b == gap
X = pd.read_csv(f"{pre}.exons.tsv", sep="\t"); X = X[X.organelle == "mito"]
_man = pd.read_csv(H.parent / "fern_tpl/m.fern.tsv", sep="\t"); FASTA = dict(zip(_man.species, _man.fasta))
def read_fasta(path):
    out, name = {}, None
    for line in open(path):
        if line.startswith(">"): name = line[1:].split()[0]; out[name] = []
        else: out[name].append(line.strip().upper())
    return {k: "".join(v) for k, v in out.items()}
def same_transcript(genome, pred, d):
    """a junction d bases further along at both ends joins the same bases if the d bases moved across are identical"""
    le, rs = pred  # 1-based last base before the intron, first after
    if d > 0: return genome[le:le + d] == genome[rs - 1:rs - 1 + d]
    return genome[le + d:le] == genome[rs - 1 + d:rs - 1]
rows = []
for bam in sorted(glob.glob(str(H / "tol.*.mito.bam"))):
    sp = re.sub(r"^tol\.|\.mito\.bam$", "", Path(bam).name)
    obs = Counter()
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
            # the left piece must come first in the read, and an intron lie between them
            if lqs < rqs and rs - le > 20 and abs(rqs - lqe) <= 15: obs[(contig.replace("mito_", ""), le, rs, rqs - lqe)] += 1
    p.wait()
    xs = X[X.species == sp]
    genome = read_fasta(FASTA[sp])
    for gene, g in xs.groupby("gene"):
        g = g.sort_values("slot").to_dict("records")
        for a, b in zip(g, g[1:]):
            if a["contig"] != b["contig"] or a["strand"] != b["strand"]: rows.append({"species": sp, "gene": gene, "junction": f"{a['slot']}-{b['slot']}", "status": "not cis"}); continue
            lo_a, hi_a = sorted((a["start"], a["end"])); lo_b, hi_b = sorted((b["start"], b["end"]))
            # 0-based half-open exons -> 1-based (last exonic base, first exonic base) in genomic order
            pred = (hi_a, lo_b + 1) if a["strand"] == "+" else (hi_b, lo_a + 1)
            near = {k: v for k, v in obs.items() if k[0] == a["contig"] and abs(k[1] - pred[0]) <= 60 and abs(k[2] - pred[1]) <= 60}
            g1 = genome[a["contig"]]
            def ok(k, W=15):
                """Does the predicted junction give the transcript this split read shows? Compared as sequence: the
                read's two pieces with their unaligned bases between them unknown (an edited base within a few bases
                of an exon end is clipped by the aligner, not mismatched) against the predicted exons joined. This
                also covers a junction placed anywhere within bases the two sides share."""
                _, le, rs, gap = k
                if pred[0] < le - W or pred[1] - 1 > rs - 1 + W: return False
                left, right = g1[le - W:le], g1[rs - 1:rs - 1 + W]
                t_read = left + "?" * gap + right if gap >= 0 else left + right[-gap:]
                t_pred = g1[le - W:pred[0]] + g1[pred[1] - 1:rs - 1 + W]
                return len(t_read) == len(t_pred) and all(x == "?" or x == y for x, y in zip(t_read, t_pred))
            exact = sum(v for k, v in near.items() if ok(k))
            # the others, as the junction they imply with the unaligned or shared bases given to the upstream side
            other = Counter()
            for k, v in near.items():
                if not ok(k): other[(k[1] + max(k[3], 0) - pred[0], k[2] - min(k[3], 0) - pred[1])] += v
            top = other.most_common(1)[0] if other else ((0, 0), 0)
            raw = Counter({k: v for k, v in near.items() if not ok(k)}).most_common(1)
            rle, rrs, rgap = (raw[0][0][1], raw[0][0][2], raw[0][0][3]) if raw else ("", "", "")
            status = "no spanning reads" if exact + top[1] < 3 else ("confirmed" if exact >= top[1] else "other junction")
            rows.append({"species": sp, "gene": gene, "junction": f"{a['slot']}-{b['slot']}", "status": status, "reads_agree": exact, "reads_other": top[1],
                         "other_shift": f"{top[0][0]},{top[0][1]}" if top[1] else "", "intron_bp": pred[1] - pred[0] - 1,
                         "other_le": rle, "other_rs": rrs, "other_gap": rgap, "strand": a["strand"], "contig": a["contig"]})
T = pd.DataFrame(rows); T.to_csv(f"{pre}.junction_check.tsv", sep="\t", index=False)
pd.set_option("display.width", 220)
print(f"{label}: {len(T)} exon junctions in {T.species.nunique()} species"); print(T.status.value_counts().to_string())
print(T.groupby(["gene", "junction", "status"]).size().unstack(fill_value=0).to_string())
print(T[T.status == "other junction"].groupby(["gene", "junction", "other_shift"]).size().sort_values(ascending=False).head(25).to_string())
