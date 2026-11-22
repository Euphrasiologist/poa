"""Shared pieces: pileups, and the RNA-edited protein of a called gene.

derive(tr) takes the called coding sequence in transcript order, one tuple per base:
(key, strand_base, depth, n_alt) where n_alt counts reads with the editing alternative (T for a C, C for a T).
- an edit is applied where n_alt >= 2 and n_alt/depth >= MIN_FRAC (0.2); poly-A libraries carry unprocessed
  transcripts, so sites are often edited in well under half the reads.
- a stop codon left after that is read through if any of its U>C sites has n_alt >= 2 at >= 0.05.
- a stop codon with depth >= MIN_DEPTH at all three bases and no such evidence is a real stop: the first one in
  the second half ends the protein (what follows is over-extension by stop read-through); one in the first 15%
  moves the start to after it (the call began upstream of the real start).
- any other remaining stop (real one mid-gene, or one without coverage) rejects the gene."""
import gzip
MIN_DEPTH, MIN_FRAC, STOP_FRAC = 10, 0.2, 0.05
# a stop codon counts as real (never edited) only with this depth at all three bases: at 20 reads a site edited in a
# fifth of transcripts shows fewer than two edited reads about 7% of the time, at 10 reads 38%
REAL_STOP_DEPTH = 20
comp = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
B = "TCAG"; AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE = {a + b + c: AA[16 * i + 4 * j + k] for i, a in enumerate(B) for j, b in enumerate(B) for k, c in enumerate(B)}
FLIP = {"C": "T", "T": "C"}

def load_pile(path):
    pile = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            c, p, ref, n, a, cc, g, t = line.split("\t")
            pile[(c, int(p))] = (ref, int(n), {"A": int(a), "C": int(cc), "G": int(g), "T": int(t)})
    return pile

def walk(pile, genome, contig, lo, hi, strand):
    """(key=(contig, 1-based pos), strand base, depth, n_alt) over the 0-based half-open span, 5'->3'"""
    out = []
    for p in (range(lo + 1, hi + 1) if strand == "+" else range(hi, lo, -1)):
        ref, n, cnt = pile.get((contig, p), (genome[contig][p - 1], 0, None))
        sref = ref if strand == "+" else comp[ref]; alt = FLIP.get(sref)
        na = cnt[alt if strand == "+" else comp[alt]] if (alt and cnt) else 0
        out.append(((contig, p), sref, n, na))
    return out

def truth(tr):
    """{key: 'C>U' | 'U>C'} at MIN_FRAC of >= MIN_DEPTH reads (>= 3 reads)"""
    return {k: ("C>U" if b == "C" else "U>C") for k, b, n, na in tr if b in FLIP and n >= MIN_DEPTH and na >= 3 and na / n >= MIN_FRAC}

def derive(tr):
    n_cod = len(tr) // 3
    seq = [FLIP[b] if (b in FLIP and na >= 2 and n and na / n >= MIN_FRAC) else b for _, b, n, na in tr]
    real, bad = [], []
    for i in range(n_cod):
        cod = "".join(seq[3 * i:3 * i + 3])
        if CODE.get(cod, "X") != "*": continue
        fixed = False
        for j in range(3 * i, 3 * i + 3):
            _, b, n, na = tr[j]
            if seq[j] == "T" and b == "T" and na >= 2 and n and na / n >= STOP_FRAC:
                trial = seq[3 * i:3 * i + 3]; trial[j - 3 * i] = "C"
                if CODE.get("".join(trial), "X") != "*": seq[j] = "C"; fixed = True; break
        if fixed: continue
        (real if all(tr[j][2] >= REAL_STOP_DEPTH for j in range(3 * i, 3 * i + 3)) else bad).append(i)
    start = max([i + 1 for i in real + bad if i < 0.15 * n_cod], default=0)
    ends = [i for i in real if i >= 0.5 * n_cod]
    end = ends[0] if ends else n_cod
    has_stop = bool(ends)
    if not ends and n_cod and CODE.get("".join(seq[3 * (n_cod - 1):3 * n_cod]), "X") == "*": end, has_stop = n_cod - 1, True
    inner = [i for i in real + bad if start <= i < end]
    prot = "".join(CODE.get("".join(seq[3 * i:3 * i + 3]), "X") for i in range(start, end))
    span = tr[3 * start:3 * end]
    cov = sum(t[2] >= MIN_DEPTH for t in span) / max(len(span), 1)
    status = "ok"
    if inner: status = f"internal stop ({'real' if any(i in real for i in inner) else 'no coverage'})"
    elif cov < 0.9: status = "low RNA coverage"
    elif len(prot) < 30: status = "short"
    return {"status": status, "protein": prot, "start": start, "end": end, "n_codons": n_cod, "has_stop": has_stop, "covered": round(cov, 2),
            "over_extended": n_cod - end - (1 if has_stop and end == n_cod - 1 else 0), "mod3": len(tr) % 3}
