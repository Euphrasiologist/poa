#!/bin/bash
# final_angio.sh SEED STAGE_SUFFIX : the angiosperm references from every dataset angiosperm mapped so far (pilot/ang_rna.list
# with pileups of both organelles and a seed run pilot/SEED/<species>): RNA-edited proteins and loci for both organelles
# (rna_refs.py on the seed runs' gene models, C-to-U only), installed and built into the package's templates and profiles.
# Writes pilot/ang_rna_final.list (the species used).
set -uo pipefail
R=/lustre/scratch122/tol/teams/blaxter/users/mb39/plant_organelle_annotator; H=$R/build/rnaedit; cd $H; SEED=$1; X=$2
sp=""; for s in $(cat pilot/ang_rna.list); do /bin/ls pub.$s.*[0-9].pile.tsv.gz > /dev/null 2>&1 && /bin/ls pub.$s.*[0-9].pltd.pile.tsv.gz > /dev/null 2>&1 && [ -f pilot/$SEED/$s/tables/edited_gene_calls.tsv ] && sp="$sp $s"; done
echo $sp | tr ' ' '\n' | sort > pilot/ang_rna_final.list; echo "species: $(wc -l < pilot/ang_rna_final.list)"
bash pilot/rna_refs_run3.sh mito $SEED stage/angio_mito_$X $sp | grep 'species:\|single-exon genes'
bash pilot/rna_refs_run3.sh pltd $SEED stage/angio_pltd_$X $sp | grep 'species:\|single-exon genes'
bash pilot/build_angio_rna3.sh angio_mito_$X angio_pltd_$X
