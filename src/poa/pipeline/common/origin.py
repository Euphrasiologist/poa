"""Genes across the origin of a circular contig.

A circular genome given as linear sequence (oatk's circular contigs, or a
finished, linearised organelle FASTA) starts at an arbitrary base, and a
gene that runs through that base is called as two partial genes, one at
each end. poa run finds such genes after gene calling: the same gene with
a call within EDGE bp of a contig's start and another within EDGE bp of
its end, on the same strand, on a contig not marked circular=false. It then
moves that contig's origin to the middle of its longest gene-free stretch
(rotate_fasta), annotates again, and at the end maps every coordinate in
the GFF and the output tables back onto the input sequence
(unrotate_gff, unrotate_table).

Back on the input, a feature that crosses the origin keeps start <= end
with end past the contig's length, GFF3's convention for circular
sequences (the contig's region feature gets Is_circular=true); its child
features (CDS, edits, frameshifts) are placed the same way, so they stay
inside their gene. In the tables an interval crossing the origin likewise
ends past the length, while a single position (an edit's genomic_pos) is
always on the contig, 0 to length-1. Offsets are 0-based:
rotated[i] = original[(i + offset) % length].
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

EDGE = 10
CIRCULAR_RE = re.compile(r"\bcircular=(true|false)\b")
ATTR_RE = re.compile(r"(?:^|;)(ID|Parent)=([^;]+)")

# per output table: (contig column, [(start, end) interval columns], [point columns])
TABLE_COORDS: dict[str, tuple[str, list[tuple[str, str]], list[str]]] = {
    "gene_calls.tsv": ("contig_id", [("start", "end")], []),
    "edited_gene_calls.tsv": (
        "contig",
        [("hit_start", "hit_end"), ("raw_start", "raw_end"), ("corrected_start", "corrected_end")],
        [],
    ),
    "predicted_edits.tsv": ("contig", [("hit_start", "hit_end")], ["genomic_pos"]),
    "frameshifts.tsv": ("contig", [("hit_start", "hit_end"), ("start", "end")], []),
    "reconstructed_exons.tsv": ("contig", [("start", "end")], []),
}


def read_fasta(path: Path) -> list[tuple[str, str]]:
    """[(header line without '>', sequence)]."""
    records: list[tuple[str, list[str]]] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                records.append((line[1:], []))
            elif records:
                records[-1][1].append(line.strip())
    return [(h, "".join(s)) for h, s in records]


def contigs(path: Path) -> dict[str, tuple[int, bool | None]]:
    """contig -> (length, circular from an oatk-style header, else None)."""
    out = {}
    for header, seq in read_fasta(path):
        m = CIRCULAR_RE.search(header)
        out[header.split()[0]] = (len(seq), None if m is None else m.group(1) == "true")
    return out


def split_genes(calls: pd.DataFrame, lengths: dict[str, tuple[int, bool | None]]) -> dict[str, list[str]]:
    """contig -> genes called at both ends of it (one species, one organelle)."""
    out: dict[str, list[str]] = {}
    for contig, g in calls.groupby("contig_id"):
        length, circular = lengths.get(contig, (None, None))
        if length is None or circular is False:
            continue
        head = g[g.start <= EDGE]
        tail = g[g.end >= length - EDGE]
        genes = sorted(set(zip(head.gene, head.strand)) & set(zip(tail.gene, tail.strand)))
        if genes:
            out[contig] = sorted({gene for gene, _ in genes})
    return out


def new_origin(calls: pd.DataFrame, length: int) -> int:
    """The middle of the longest stretch no call covers, going round the
    circle (calls: one contig's)."""
    ivs = sorted(zip(calls.start.astype(int), calls.end.astype(int)))
    if not ivs:
        return 0
    merged = [list(ivs[0])]
    for s, e in ivs[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    best_len, best_mid = -1, 0
    for k, (_, e) in enumerate(merged):
        nxt = merged[k + 1][0] if k + 1 < len(merged) else merged[0][0] + length
        if nxt - e > best_len:
            best_len, best_mid = nxt - e, (e + (nxt - e) // 2) % length
    return best_mid


def rotate_fasta(path: Path, offsets: dict[str, int]) -> None:
    """Rewrite `path` with each contig in `offsets` starting at its offset."""
    records = read_fasta(path)
    with open(path, "w") as fh:
        for header, seq in records:
            k = offsets.get(header.split()[0], 0)
            seq = seq[k:] + seq[:k]
            fh.write(f">{header}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i : i + 80] + "\n")


def point(p: int, offset: int, length: int) -> int:
    return (p + offset) % length


def interval(s: int, e: int, offset: int, length: int) -> tuple[int, int]:
    """0-based half-open [s, e) on the rotated contig -> on the input; the
    end may pass the length (the interval crosses the origin)."""
    s2 = (s + offset) % length
    return s2, s2 + (e - s)


def unrotate_gff(path: Path, offsets: dict[str, int], lengths: dict[str, int]) -> None:
    """Map a GFF3 written on rotated contigs back onto the input contigs."""
    lines = Path(path).read_text().splitlines()
    parent: dict[str, str] = {}
    for line in lines:
        f = line.split("\t")
        if line.startswith("#") or len(f) < 9:
            continue
        attrs = dict(ATTR_RE.findall(f[8]))
        if "ID" in attrs and "Parent" in attrs:
            parent.setdefault(attrs["ID"], attrs["Parent"].split(",")[0])

    def root(fid: str) -> str:
        seen = set()
        while fid in parent and fid not in seen:
            seen.add(fid)
            fid = parent[fid]
        return fid

    # where each top-level feature starts on the input, per contig
    root_start: dict[tuple[str, str], int] = {}
    for line in lines:
        f = line.split("\t")
        if line.startswith("#") or len(f) < 9 or f[0] not in offsets or f[2] == "region":
            continue
        attrs = dict(ATTR_RE.findall(f[8]))
        if "Parent" not in attrs and "ID" in attrs:
            key = (f[0], attrs["ID"])
            s = interval(int(f[3]) - 1, int(f[4]), offsets[f[0]], lengths[f[0]])[0]
            root_start[key] = min(root_start.get(key, s), s)

    head, body = [], []
    for line in lines:
        f = line.split("\t")
        if line.startswith("#") or len(f) < 9:
            head.append(line)
            continue
        if f[0] in offsets:
            if f[2] == "region":
                if "Is_circular=" not in f[8]:
                    f[8] += ";Is_circular=true"
            else:
                n = lengths[f[0]]
                s, e = interval(int(f[3]) - 1, int(f[4]), offsets[f[0]], n)
                attrs = dict(ATTR_RE.findall(f[8]))
                rid = root(attrs.get("Parent", "").split(",")[0] or attrs.get("ID", ""))
                rs = root_start.get((f[0], rid))
                if rs is not None and s < rs:
                    s, e = s + n, e + n  # past the origin, inside its gene
                f[3], f[4] = str(s + 1), str(e)
        body.append(f)
    order = {"region": 0}
    body.sort(key=lambda f: (f[0], int(f[3]), order.get(f[2], 1)))
    Path(path).write_text("\n".join(head + ["\t".join(f) for f in body]) + "\n")


def unrotate_table(
    path: Path,
    organelle: str,
    offsets: dict[str, int],
    lengths: dict[str, int],
    exons: pd.DataFrame | None = None,
) -> None:
    """Map one output table's coordinates back onto the input contigs, for
    the rows of `organelle` (mito and plastid contigs can share names).
    reconstructed_edits.tsv has no contig column: its edits take their
    exon's contig (gene, slot) from `exons`, reconstructed_exons as it was
    before being mapped back."""
    path = Path(path)
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        return
    if path.name == "reconstructed_edits.tsv":
        if exons is None or exons.empty:
            return
        where = exons.drop_duplicates(["species", "organelle", "gene", "slot"]).set_index(
            ["species", "organelle", "gene", "slot"]
        )["contig"]
        contig = [where.get(k) for k in zip(df.species, df.organelle, df.gene, df.slot)]
        df = df.assign(_contig=contig)
        spec = ("_contig", [], ["genomic_pos"])
    elif path.name in TABLE_COORDS:
        spec = TABLE_COORDS[path.name]
    else:
        return
    ccol, intervals, points = spec
    for i, row in df.iterrows():
        c = row[ccol]
        if row["organelle"] != organelle or c not in offsets:
            continue
        k, n = offsets[c], lengths[c]
        for s_col, e_col in intervals:
            if pd.notna(row[s_col]):
                df.at[i, s_col], df.at[i, e_col] = interval(int(row[s_col]), int(row[e_col]), k, n)
        for p_col in points:
            if pd.notna(row[p_col]):
                df.at[i, p_col] = point(int(row[p_col]), k, n)
    df = df.drop(columns=[c for c in ("_contig",) if c in df.columns])
    df.to_csv(path, sep="\t", index=False)
