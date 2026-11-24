#!/bin/bash
# site_score.sh KIND TAG... : score held-out runs. KIND ds: dataset angiosperms (pilot/ang_rna_p.list) against RNA-seq, both
# organelles, and every predicted edit against the reads; gb: the 19 RefSeq genomes against their annotation
set -uo pipefail
R=/lustre/scratch122/tol/teams/blaxter/users/mb39/plant_organelle_annotator; H=$R/build/rnaedit; cd $H/pilot
E="env -i HOME=$HOME PYTHONNOUSERSITE=1 TMPDIR=/tmp PYTHONPATH=$R/src PATH=$H/final10/bin:$R/build/pkgenv/bin:/usr/bin:/bin USE_DATASET=1"
K=$1; shift; SP=$(cat ${LIST:-ang_rna_p.list} | tr "\n" " ")
for T in "$@"; do
  if [ "$K" == "ds" ]; then
    $E python assess.py $T $SP > assess.$T.log 2>&1 || { echo "$T mito failed"; tail -3 assess.$T.log; }
    $E python assess_pltd.py $T $SP > assess_pltd.$T.log 2>&1 || { echo "$T pltd failed"; tail -3 assess_pltd.$T.log; }
    (cd $H; for o in mito pltd; do $E python pilot/call_frac.py $T $o $SP | head -1; done)
  else
    (cd $H; echo "== $T"; $E python pilot/score_models.py $T poa:$T $(cat pilot/ang_all.list) 2>&1 | grep 'CDS with every exon'
     $E python pilot/prot_score.py $T poa:$T $(cat pilot/ang_edset.list) 2>&1 | grep TOTAL; $E python pilot/edit_score.py $T poa:$T $(cat pilot/ang_edset.list) 2>&1 | grep TOTAL)
  fi
done
