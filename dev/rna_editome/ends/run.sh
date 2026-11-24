#!/bin/bash
# run.sh TAG BINDIR [orfedit options...] : orfedit scan-batch over the pooled manifests with poa's options for angiosperm
# profiles (2.5; 1.5/10 for the species that follow the editing sites) plus any given, into ends/TAG.calls.tsv and TAG.fs.tsv
set -uo pipefail
H=/lustre/scratch122/tol/teams/blaxter/users/mb39/plant_organelle_annotator/build/rnaedit/ends; TAG=$1; BIN=$2; shift 2
W=$H/tmp.$TAG; mkdir -p $W
$BIN/orfedit scan-batch --manifest $H/sites.manifest.tsv --gene-calls-out $W/s.calls --edits-out $W/s.edits --frameshifts-out $W/s.fs --threads ${THREADS:-4} --flank 300 --edit-penalty 2.5 --frameshift-penalty 20 --known-edit-penalty 1.5 --novel-edit-penalty 10 "$@" 2> $W/s.log || { echo "$TAG sites failed"; tail -3 $W/s.log; }
$BIN/orfedit scan-batch --manifest $H/flat.manifest.tsv --gene-calls-out $W/f.calls --edits-out $W/f.edits --frameshifts-out $W/f.fs --threads ${THREADS:-4} --flank 300 --edit-penalty 2.5 --frameshift-penalty 20 "$@" 2> $W/f.log || { echo "$TAG flat failed"; tail -3 $W/f.log; }
(cat $W/s.calls; tail -n +2 $W/f.calls) > $H/$TAG.calls.tsv; (cat $W/s.fs; tail -n +2 $W/f.fs) > $H/$TAG.fs.tsv; (cat $W/s.edits; tail -n +2 $W/f.edits) > $H/$TAG.edits.tsv
rm -rf $W; echo "$TAG: $(($(wc -l < $H/$TAG.calls.tsv) - 1)) calls, $(($(wc -l < $H/$TAG.fs.tsv) - 1)) frame changes, $(($(wc -l < $H/$TAG.edits.tsv) - 1)) edits"
