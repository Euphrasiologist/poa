#!/bin/bash
# prep_arr.sh : array element $LSB_JOBINDEX of ends/genomes.tsv (NAME ORG FASTA ACC): a poa run with the current tools, keeping the
# work directory (orfedit's manifests and the contig FASTA they name) in ends/runs/NAME
set -uo pipefail
R=/lustre/scratch122/tol/teams/blaxter/users/mb39/plant_organelle_annotator; H=$R/build/rnaedit
read S ORG FA ACC < <(sed -n "${LSB_JOBINDEX}p" $H/ends/genomes.tsv)
O=$H/ends/runs/$S; rm -rf $O; mkdir -p $(dirname $O)
E="env -i HOME=$HOME PYTHONNOUSERSITE=1 TMPDIR=/tmp PYTHONPATH=$R/src PATH=$H/final3/bin:$R/build/pkgenv/bin:/usr/bin:/bin"
$E poa run --$ORG $FA -o $O -n $S --lineage auto -t 2 > $O.log 2>&1; rc=$?
# keep the editing manifests and the data tree they point into; drop the rest of the work directory
find $O/work -type f | grep -v '/editing/work/manifest\|\.fasta$\|\.fa$\|\.fai$' | xargs -r rm -f; find $O/work -type d -empty -delete
exit $rc
