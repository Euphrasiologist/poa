#!/bin/bash
# map.sh SPECIES TAG R1 R2 : map RNA-seq to the species' mito + plastid contigs, pile up the mito contigs.
# Writes TAG.pile.tsv.gz (contig pos ref depth A C G T, counts of bases with baseQ>=20 from reads with MAPQ>=20).
set -euo pipefail
source /etc/profile.d/modules.sh
module load bwa-mem2/2.2.1--hd03093a_2 samtools/1.20--h50ea8bc_0
S=$1; T=$2; R1=$3; R2=${4:-}
D=/lustre/scratch122/tol/teams/blaxter/users/mb39/plant_organellar_database/data
W=$(dirname $0)/$T.work2; mkdir -p $W
sed 's/^>/>mito_/' $D/mito/$S/*.mito.ctg.fasta > $W/ref.fa
sed 's/^>/>pltd_/' $D/plastid/$S/*.pltd.ctg.fasta >> $W/ref.fa
bwa-mem2 index $W/ref.fa 2> /dev/null
samtools faidx $W/ref.fa
# only mapped reads are kept; the nuclear majority is dropped on the way
if [[ $R1 == irods:* ]]; then
  # an unaligned CRAM on iRODS: fetch, then stream it as interleaved pairs
  module load ISG/IRODS/1.0
  for try in 1 2 3; do rm -f $W/reads.cram; timeout 40m iget -K "${R1#irods:}" $W/reads.cram && break; done
  [ -s $W/reads.cram ]
  READS="samtools fastq -@ 2 -n $W/reads.cram"; P=-p
else
  READS="cat /dev/null"; P=
fi
{ if [[ -n $P ]]; then $READS 2> /dev/null | bwa-mem2 mem -t ${NT:-4} -p $W/ref.fa - ; else bwa-mem2 mem -t ${NT:-4} $W/ref.fa $R1 $R2; fi; } 2> $W/bwa.log | samtools view -b -F 4 - | samtools sort -@ 2 -m 1G -o $W/map.bam -
samtools index $W/map.bam
# keep the reads on the mito contigs (with their split alignments): junction evidence for spliced genes
samtools view -b -o $(dirname $0)/$T.mito.bam $W/map.bam $(grep '^>mito_' $W/ref.fa | tr -d '>' | awk '{print $1}')
samtools index $(dirname $0)/$T.mito.bam
samtools idxstats $W/map.bam > $(dirname $0)/$T.idxstats
grep '^>mito_' $W/ref.fa | tr -d '>' | awk '{print $1}' > $W/mito.txt
for c in $(cat $W/mito.txt); do
  samtools mpileup -r $c -f $W/ref.fa -q ${MQ:-0} -Q 20 -d 0 -a -B --ff UNMAP,SECONDARY,QCFAIL,SUPPLEMENTARY $W/map.bam 2> /dev/null
done | python3 -c '
import sys, re
for line in sys.stdin:
    c, p, ref, n, bases, _ = line.rstrip("\n").split("\t")
    bases = re.sub(r"\^.", "", bases)
    out = []; i = 0
    while i < len(bases):
        ch = bases[i]
        if ch in "+-":
            j = i + 1
            while bases[j].isdigit(): j += 1
            i = j + int(bases[i + 1:j]); continue
        out.append(ch); i += 1
    s = "".join(out).upper(); r = ref.upper()
    cnt = {b: s.count(b) for b in "ACGT"}
    cnt[r] = cnt.get(r, 0) + s.count(".") + s.count(",")
    print(c.replace("mito_", ""), p, r, sum(cnt[b] for b in "ACGT"), cnt["A"], cnt["C"], cnt["G"], cnt["T"], sep="\t")
' | gzip > $W/pile.tsv.gz
[ -s $(dirname $0)/$T.pile.tsv.gz ] || mv $W/pile.tsv.gz $(dirname $0)/$T.pile.tsv.gz
rm -rf $W
