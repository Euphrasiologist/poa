#!/bin/bash
# map4.sh SPECIES RUN URL1 MD5_1 URL2 MD5_2 : a public paired RNA-seq run (ENA fastq) mapped to the species' own mito +
# plastid contigs. Writes pub.SPECIES.RUN.pile.tsv.gz, .pltd.pile.tsv.gz (contig pos ref depth A C G T) and .mito.bam.
set -euo pipefail
source /etc/profile.d/modules.sh
module load bwa-mem2/2.2.1--hd03093a_2 samtools/1.20--h50ea8bc_0
S=$1; RUN=$2; H=$(dirname $0); T=pub.$S.$RUN
D=/lustre/scratch122/tol/teams/blaxter/users/mb39/plant_organellar_database/data
W=$H/$T.work; mkdir -p $W
get() {  # get URL MD5 OUT
  for try in 1 2 3 4 5 6 7 8 9 10 11 12; do
    # resume a partial file; if the server will not resume (curl 33), start that file again
    rc=0; curl -sS -A curl -L -C - --retry 5 --retry-delay 30 -o $3 "https://$1" || rc=$?
    [ "$(md5sum < $3 | cut -d' ' -f1)" == "$2" ] && return 0
    # a file curl thinks is whole (or that the server will not resume: 33, 22, 36) but with the wrong checksum: start it again
    case $rc in 0|22|33|36) rm -f $3;; esac
    sleep 120
  done
  echo "download failed: $1" >&2; return 1
}
get $3 $4 $W/r1.fq.gz; get $5 $6 $W/r2.fq.gz
sed 's/^>/>mito_/' $(ls $D/mito/$S/*.mito.ctg.fasta | head -1) > $W/ref.fa
sed 's/^>/>pltd_/' $(ls $D/plastid/$S/*.pltd.ctg.fasta | head -1) >> $W/ref.fa
bwa-mem2 index $W/ref.fa 2> /dev/null
samtools faidx $W/ref.fa
bwa-mem2 mem -t ${NT:-2} $W/ref.fa $W/r1.fq.gz $W/r2.fq.gz 2> $W/bwa.log | samtools view -b -F 4 - | samtools sort -@ 1 -m 1G -T /tmp/$T.$$ -o $W/map.bam -
rm -f $W/r1.fq.gz $W/r2.fq.gz
samtools index $W/map.bam
samtools view -b -o $H/$T.mito.bam $W/map.bam $(grep '^>mito_' $W/ref.fa | tr -d '>' | awk '{print $1}'); samtools index $H/$T.mito.bam
samtools idxstats $W/map.bam > $H/$T.idxstats
cat > $W/pile.py <<'PY'
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
    if r in cnt: cnt[r] += s.count(".") + s.count(",")
    print(c.split("_", 1)[1], p, r, sum(cnt.values()), cnt["A"], cnt["C"], cnt["G"], cnt["T"], sep="\t")
PY
pile() {
  for c in $(grep "^>$1_" $W/ref.fa | tr -d '>' | awk '{print $1}'); do
    samtools mpileup -A -r $c -f $W/ref.fa -q 0 -Q 20 -d 0 -a -B --ff UNMAP,SECONDARY,QCFAIL,SUPPLEMENTARY $W/map.bam 2> /dev/null
  done | python3 $W/pile.py | gzip > $W/pile.$1.tsv.gz
  mv $W/pile.$1.tsv.gz $2
}
pile mito $H/$T.pile.tsv.gz
pile pltd $H/$T.pltd.pile.tsv.gz
rm -rf $W
