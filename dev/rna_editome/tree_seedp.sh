#!/bin/bash
# tree_seedp.sh : a copy of the package (build/rnaedit/tree_seedp/src) in which fern, horsetail and clubmoss have no plastid
# references of their own, so that their plastid genes are called on the default ones - as they were when those lineages'
# plastid references were first read from RNA-seq. Seed runs for pltd_refs.py are made with it (run2_seedp.sh): a lineage's
# own profile cannot be the seed for itself where it is too short (rbcL, cut at 341 codons by the old window).
set -uo pipefail
R=/lustre/scratch122/tol/teams/blaxter/users/mb39/plant_organelle_annotator; H=$R/build/rnaedit; T=$H/tree_seedp
rm -rf $T; mkdir -p $T; cp -r $R/src $T/; find $T -name __pycache__ -type d -prune -exec rm -rf {} +
P=$T/src/poa/pipeline
sed -i 's|^LINEAGE_PLASTID = ("fern", "horsetail", "lycophyte", "moss", "liverwort")|LINEAGE_PLASTID = ("moss", "liverwort")|' $P/trans_splicing/reference_genes.py; grep -c '^LINEAGE_PLASTID = ("moss", "liverwort")' $P/trans_splicing/reference_genes.py
for L in fern horsetail lycophyte; do rm -rf $P/editing/reference/lineages/$L/pltd $P/trans_splicing/reference/lineages/$L/plastid; done
