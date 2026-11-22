#!/bin/bash
# Rebuild the bundled fern editing profiles and exon templates from poa's own builders and validate on the DToL ferns.
set -uo pipefail
R=/lustre/scratch122/tol/teams/blaxter/users/mb39/plant_organelle_annotator; I=$R/build/rnaedit/wire; rm -rf $I; mkdir -p $I; cd $I
NEW=$R/build/rev_edit/bin2
E="env -i HOME=$HOME PYTHONNOUSERSITE=1 TMPDIR=/tmp PYTHONPATH=$R/src PATH=$NEW:$R/build/pkgenv/bin:/usr/bin:/bin"
FERN='Asplenium|Dryopteris|Cystopteris|Polystichum|Onoclea|Phegopteris|Polypodium|Struthiopteris|Woodsia|Gymnocarpium|Adiantum|Thelypteris|Pteridium|Azolla|Pilularia|Vandenboschia'
P=$R/src/poa/pipeline/editing/reference/lineages/fern
ts() { $NEW/transsplice scan-batch --manifest $1 --genes-out $2.genes.tsv --exons-out $2.exons.tsv --junctions-out $2.junctions.tsv --edits-out $2.edits.tsv --frameshifts-out $2.fs.tsv --threads 4 --stop-edit-penalty 2.5 --reverse-edit-penalty 5 2>/dev/null; }
oe() { $NEW/orfedit scan-batch --manifest $1 --gene-calls-out $2.calls.tsv --edits-out $2.edits.tsv --frameshifts-out $2.fs.tsv --frameshift-penalty 20 --frameshift-flank 0 --stop-edit-penalty 2.5 --reverse-edit-penalty 5 --threads 4 2>/dev/null; }
# 1. editing profiles
mkdir w; $E PLANT_ORGANELLE_DATA_ROOT=$I/w python $R/src/poa/pipeline/editing/src/00b_build_lineage_profiles.py --lineage fern > build.editing.log 2>&1 || echo "editing build failed"
python3 - <<PY
import pandas as pd, os
m=pd.read_csv("$R/build/rev_edit/fern.fs_val_manifest.tsv",sep="\t")
m=m[m.species.str.contains("$FERN")&(m.organelle=="mito")]
d="$P/mito"; m=m[[os.path.exists(f"{d}/{g.split('|')[0]}.pssm") for g in m.gene]].copy(); m["profile"]=[f"{d}/{g.split('|')[0]}.pssm" for g in m.gene]; m.to_csv("single.tsv",sep="\t",index=False); print("single-exon rows",len(m))
PY
oe single.tsv single
# 2. exon templates: unmasked, mask from the ferns' reconstructions, masked
rm -f $R/src/poa/pipeline/trans_splicing/reference/lineages/fern/whole_gene_mask.tsv
$E PLANT_ORGANELLE_DATA_ROOT=$I/w python $R/src/poa/pipeline/trans_splicing/src/00_build_exon_profiles.py --lineage fern --no-mask > build.tpl.unmasked.log 2>&1 || echo "template build failed"
ts $R/build/fern_tpl/m.fern.tsv spliced.unmasked
$E python -m poa.pipeline.trans_splicing.column_mask spliced.unmasked.genes.tsv --lineage fern --species "$FERN" > mask.log 2>&1
$E PLANT_ORGANELLE_DATA_ROOT=$I/w python $R/src/poa/pipeline/trans_splicing/src/00_build_exon_profiles.py --lineage fern > build.tpl.log 2>&1 || echo "template build failed"
rm -rf w
ts $R/build/fern_tpl/m.fern.tsv spliced
# 3. planted-indel recall
python3 $R/build/recall/plant_single.py single.tsv single.calls.tsv rcs 40 5
oe rcs.manifest.tsv rcs.new; python3 $R/build/recall/score.py single rcs.truth.tsv rcs.new.fs.tsv "single-exon, new fern profiles" > recall.txt 2>&1
python3 $R/build/recall/plant_spliced.py $R/build/fern_tpl/m.fern.tsv spliced rcp 40 6 --only "$FERN"
ts rcp.manifest.tsv rcp.new; python3 $R/build/recall/score.py spliced rcp.truth.tsv rcp.new.fs.tsv "spliced, new fern templates" >> recall.txt 2>&1
rm -rf rcs.genomes rcp.genomes rcs.manifest*.tsv rcp.manifest*.tsv
# 4. end to end, unit tests, smoke test
D=/lustre/scratch122/tol/teams/blaxter/users/mb39/plant_organellar_database
mkdir e2e; $E poa run --mito $(ls $D/data/mito/Asplenium_ruta_muraria/*.mito.ctg.fasta) --pltd $(ls $D/data/plastid/Asplenium_ruta_muraria/*.pltd.ctg.fasta) -o e2e/fern -n Asplenium_ruta_muraria --lineage fern -t 4 > e2e/fern.log 2>&1 || echo "poa run failed"
(cd $R; for t in tests/test_*.py; do m=tests.$(basename $t .py); $E python -m unittest $m 2>&1 | tail -3 | tr "\n" " "; echo " <- $m"; done) > unittest.log 2>&1
sed "s|build/fern_edit|build/rnaedit/wire|g" $R/build/fern_edit/smoke.sh > smoke.sh; bash smoke.sh; tail -2 smoke.log
