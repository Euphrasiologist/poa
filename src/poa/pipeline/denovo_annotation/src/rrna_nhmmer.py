"""rRNA calling: nhmmer against poa's bundled plant rRNA HMMs -> GFF3.

Python port of barrnap 0.9 run as `barrnap --kingdom plant`, where `plant`
was a database added to one cluster's barrnap install (stock barrnap has
no such kingdom): barrnap's own bacterial 16S/23S/5S models plus Rfam's
RF00001 5S model, now bundled as poa/data/rrna/plant.hmm (see NOTICE
there). Only nhmmer (HMMER) is needed. Checked byte-identical against that
barrnap on the 1250-species dataset (tests/test_rrna_nhmmer.py keeps a
fixture).

barrnap's logic, as used here:
  - nhmmer --cpu N -E <evalue> --w_length 3878 --tblout (MAXLEN is 1.2 x
    the longest model length barrnap knows);
  - a hit shorter than int(0.25 x model length) is rejected; shorter than
    int(0.8 x model length) is kept with product "... (partial)" and a
    note `aligned only <int percent> percent of the <product>`;
  - GFF score column is nhmmer's E-value string as printed (barrnap's
    choice, kept);
  - sorted by seqid then start, source `barrnap:0.9`.
"""
from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

PLANT_HMM = Path(__file__).resolve().parents[3] / "data" / "rrna" / "plant.hmm"
LENG = {"5S_rRNA": 119, "5S_rRNA_all_life": 120, "16S_rRNA": 1585, "23S_rRNA": 3232,
        "5_8S_rRNA": 156, "18S_rRNA": 1869, "28S_rRNA": 2912, "12S_rRNA": 954}
MAXLEN = int(1.2 * max(LENG.values()))
LENCUTOFF, REJECT = 0.8, 0.25
SOURCE = "barrnap:0.9"


def nhmmer_exe() -> str:
    exe = shutil.which("nhmmer") or "/software/team301/hmmer-3.4/src/nhmmer"
    if not Path(exe).exists():
        raise FileNotFoundError("nhmmer (HMMER) not found on PATH")
    return exe


def tblout_to_gff(tblout_text: str) -> str:
    feats = []
    for line in tblout_text.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        x = line.split()
        if len(x) < 8 or not x[6].isdigit():
            raise ValueError(f"bad line in nhmmer output - {line}")
        a, b = int(x[6]), int(x[7])
        begin, end, strand = (a, b, "+") if a < b else (b, a, "-")
        seqid, gene = x[0], x[2]
        score = x[12] if len(x) > 12 else "."
        if gene not in LENG:
            raise ValueError(f"Detected unknown gene '{gene}' in scan")
        prod = gene.replace("_r", " ribosomal ", 1).replace("5_8", "5.8", 1)
        length = end - begin + 1
        if length < int(REJECT * LENG[gene]):
            continue
        tags = f"Name={gene};product={prod}"
        if length < int(LENCUTOFF * LENG[gene]):
            note = f"aligned only {int(100 * length / LENG[gene])} percent of the {prod}"
            tags = f"Name={gene};product={prod} (partial);note={note}"
        feats.append((seqid, begin, "\t".join([seqid, SOURCE, "rRNA", str(begin), str(end), score, strand, ".", tags])))
    feats.sort(key=lambda f: (f[0], f[1]))  # stable, like barrnap's sort
    return "##gff-version 3\n" + "".join(f[2] + "\n" for f in feats)


def call_rrna(fasta: str | Path, evalue: str = "1e-6", threads: str | int = 1, timeout: int = 600) -> str:
    """GFF3 text of the rRNA genes in `fasta`."""
    with tempfile.NamedTemporaryFile(suffix=".tblout") as tbl:
        cmd = [nhmmer_exe(), "--cpu", str(threads), "-E", str(evalue), "--w_length", str(MAXLEN),
               "-o", "/dev/null", "--tblout", tbl.name, str(PLANT_HMM), str(fasta)]
        subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=timeout)
        return tblout_to_gff(Path(tbl.name).read_text())
