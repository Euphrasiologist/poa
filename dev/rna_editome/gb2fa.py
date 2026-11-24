"""gb2fa.py LIST GBDIR OUTDIR : NAME ACCESSION lines -> OUTDIR/NAME.fasta from GBDIR/ACCESSION.gb (run with python -I)"""
import sys
from pathlib import Path
from Bio import SeqIO
for l in open(sys.argv[1]):
    n, a = l.split(); rec = SeqIO.read(Path(sys.argv[2]) / f"{a}.gb", "genbank")
    Path(sys.argv[3], f"{n}.fasta").write_text(f">{rec.id}\n{str(rec.seq).upper()}\n"); print(n, a, len(rec.seq), sum(f.type == "CDS" for f in rec.features))
