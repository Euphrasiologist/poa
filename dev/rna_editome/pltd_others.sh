#!/bin/bash
# pltd_others.sh SEED : RNA-edited plastid references of fern, horsetail and clubmoss again, from poa runs made with the
# current tools (pilot/SEED/<species>: whole rbcL, genes ending at their stops), installed and built
set -uo pipefail
R=/lustre/scratch122/tol/teams/blaxter/users/mb39/plant_organelle_annotator; H=$R/build/rnaedit; cd $H; SEED=$1
E="env -i HOME=$HOME PYTHONNOUSERSITE=1 TMPDIR=/tmp PYTHONPATH=$R/src PATH=$H/final10/bin:$R/build/pkgenv/bin:/usr/bin:/bin"
F="Asplenium_scolopendrium Asplenium_trichomanes Asplenium_viride Azolla_filiculoides Dryopteris_borreri Dryopteris_dilatata Dryopteris_filix_mas Onoclea_sensibilis Phegopteris_connectilis Polystichum_setiferum Woodsia_ilvensis"
HS="Equisetum_arvense Equisetum_fluviatile Equisetum_hyemale Equisetum_palustre Equisetum_pratense Equisetum_sylvaticum Equisetum_telmateia Equisetum_variegatum"
LY="Huperzia_selago Diphasiastrum_alpinum"
have() { for s in "$@"; do [ -f pilot/$SEED/$s/tables/edited_gene_calls.tsv ] && printf "%s " $s; done; }
for x in "fern $(have $F)" "horsetail $(have $HS)" "lycophyte $(have $LY)"; do
  set -- $x; L=$1; shift
  echo "== $L: $# species"
  $E python pltd_refs.py $SEED stage/${L}_pltd3 "$@" > pilot/pltd_refs3.$L.log 2>&1 || { echo "$L pltd_refs failed"; tail -3 pilot/pltd_refs3.$L.log; continue; }
  grep -h 'single-exon genes' pilot/pltd_refs3.$L.log | cut -c1-200
  bash pilot/build_pltd3.sh $L ${L}_pltd3 | cut -c1-200
  echo "rbcL columns: $(($(wc -l < $R/src/poa/pipeline/editing/reference/lineages/$L/pltd/rbcL.pssm) - 1))"
done
