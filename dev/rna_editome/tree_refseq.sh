#!/bin/bash
# tree_refseq.sh : a copy of the package (build/rnaedit/tree_refseq/src) whose angiosperm references are the edited RefSeq
# mitogenomes alone, as they were before the RNA-derived genes went in: no raw_rna loci, no rna_proteins, no angiosperm plastid
# profiles of their own (the dataset-built ones serve). Seed runs for rna_refs.py are made with it (run2_seed.sh).
set -uo pipefail
R=/lustre/scratch122/tol/teams/blaxter/users/mb39/plant_organelle_annotator; H=$R/build/rnaedit; T=$H/tree_refseq
rm -rf $T; mkdir -p $T; cp -r $R/src $T/; find $T -name __pycache__ -type d -prune -exec rm -rf {} +
P=$T/src/poa/pipeline; ED=$P/editing/reference/lineages/angiosperm
rm -rf $P/trans_splicing/reference/raw_rna $P/trans_splicing/reference/plastid/raw_rna $ED/pltd; rm -f $ED/rna_proteins.fasta $ED/rna_proteins.pltd.fasta
# and no editing sites, so one cost of an edit everywhere, as then
sed -i 's|^SITE_LINEAGES = (DEFAULT_LINEAGE,)|SITE_LINEAGES = ()|' $P/common/lineage.py; grep -c '^SITE_LINEAGES = ()' $P/common/lineage.py
E="env -i HOME=$HOME PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 TMPDIR=/tmp PYTHONPATH=$T/src PATH=$H/final10/bin:$R/build/pkgenv/bin:/usr/bin:/bin"
W=$H/pilot/build.treerefseq.$$; mkdir -p $W
$E PLANT_ORGANELLE_DATA_ROOT=$W python $P/trans_splicing/src/00_build_exon_profiles.py --lineage angiosperm > $H/pilot/treerefseq.tpl.log 2>&1 || echo "mito template build failed"
$E PLANT_ORGANELLE_DATA_ROOT=$W python $P/editing/src/00b_build_lineage_profiles.py --lineage angiosperm > $H/pilot/treerefseq.ed.log 2>&1 || echo "mito profile build failed"
$E PLANT_ORGANELLE_DATA_ROOT=$W python $P/trans_splicing/src/00_build_exon_profiles.py --lineage angiosperm --organelle pltd > $H/pilot/treerefseq.ptpl.log 2>&1 || echo "plastid template build failed"
rm -rf $W
echo "mito profiles: $(/bin/ls $P/editing/reference/profiles_refseq/mito | wc -l); with an edit column: $(grep -l 'edit$' $P/editing/reference/profiles_refseq/mito/*.pssm 2>/dev/null | wc -l)"
echo "plastid templates: $(/bin/ls $P/trans_splicing/reference/plastid/profiles | tr '\n' ' ')"
$E python -c "import poa, sys; print(poa.__file__)"
grep -h 'canonical exon count' $H/pilot/treerefseq.tpl.log | head -3 | cut -c1-150
