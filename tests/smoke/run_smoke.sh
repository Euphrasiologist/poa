#!/usr/bin/env bash
# End-to-end smoke test: Arabidopsis thaliana mito + plastid (tests/data/)
# through QC -> denovo_annotation -> editing -> trans_splicing ->
# unitig_coords -> gff_export, in a fresh data root separate from this
# checkout, then the mito GFF checked against RefSeq NC_037304.1
# (check_arabidopsis.py). Exits non-zero if any stage or the check fails.
#
# Usage: tests/smoke/run_smoke.sh [DATA_ROOT] [THREADS]
#   DATA_ROOT  where to build the test data root (default: a new temp dir);
#              must not already contain data/
#   THREADS    default 4
# Needs the tools from environment.yml + install.sh on PATH (or this
# cluster's fallbacks in analysis/common/tool_paths.sh), and oatkDB's
# .fam databases via OATKDB_MITO_FAM / OATKDB_PLTD_FAM - see README.md.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
DATA_ROOT="${1:-$(mktemp -d "${TMPDIR:-/tmp}/annotator_smoke.XXXXXX")}"
THREADS="${2:-4}"
SP=Arabidopsis_thaliana

if [[ -e "${DATA_ROOT}/data" ]]; then
  echo "[err] ${DATA_ROOT}/data already exists - pass a fresh directory" >&2
  exit 1
fi
log() { echo "[$(date '+%H:%M:%S')] $*" >&2; }

for organelle in mito plastid; do
  mkdir -p "${DATA_ROOT}/data/${organelle}/${SP}"
  for gz in "${REPO}/tests/data/${organelle}/${SP}"/*.gz; do
    gunzip -c "${gz}" > "${DATA_ROOT}/data/${organelle}/${SP}/$(basename "${gz}" .gz)"
  done
done
export PLANT_ORGANELLE_DATA_ROOT="${DATA_ROOT}"
cd "${DATA_ROOT}"
A="${REPO}/analysis"
GENE_CALLS="${DATA_ROOT}/analysis/denovo_annotation/results/gene_calls.tsv"
log "data root: ${DATA_ROOT}"

log "=== qc_basic_stats ==="
for s in 01_gfa_stats 02_contig_stats 03_gene_matrix 04_core_gene_list 05_qc_summary; do
  python3 "${A}/qc_basic_stats/src/${s}.py" --organelle both
done

log "=== denovo_annotation ==="
python3 "${A}/denovo_annotation/src/00_run_nhmmscan.py" --organelle both --jobs 2 --cpu "$(( (THREADS + 1) / 2 ))"
python3 "${A}/denovo_annotation/src/01_run_trnascan.py" --organelle both --jobs 2 --thread 1
python3 "${A}/denovo_annotation/src/02_run_barrnap.py" --organelle both --jobs 2 --threads 1
python3 "${A}/denovo_annotation/src/03_build_combined_gff.py" --organelle both
python3 "${A}/denovo_annotation/src/04_build_gene_calls.py" --organelle both

log "=== editing (bundled profiles) ==="
python3 "${A}/editing/src/01_scan_editing.py" --organelle both --threads "${THREADS}" --gene-calls "${GENE_CALLS}"

log "=== trans_splicing ==="
python3 "${A}/trans_splicing/src/01_reconstruct.py" --threads "${THREADS}" --gene-calls "${GENE_CALLS}"

log "=== unitig_coords ==="
python3 "${A}/unitig_coords/src/00_build_unitig_map.py" --organelle both --jobs 2
python3 "${A}/unitig_coords/src/01_linearization_qc.py" --organelle both
python3 "${A}/qc_basic_stats/src/06_low_depth_paths.py"

log "=== gff_export ==="
python3 "${A}/gff_export/src/01_build_gff.py" --organelle both --gene-calls "${GENE_CALLS}"

OUT="${DATA_ROOT}/analysis/gff_export/results"
for o in mito pltd; do
  [[ -s "${OUT}/${SP}.${o}.gff" ]] || { echo "[err] no ${o} GFF: ${OUT}/${SP}.${o}.gff" >&2; exit 1; }
done
log "plastid GFF: $(grep -c $'\tgene\t' "${OUT}/${SP}.pltd.gff") genes"

log "=== mito GFF vs RefSeq NC_037304.1 ==="
python3 "${REPO}/tests/smoke/check_arabidopsis.py" \
  --gff "${OUT}/${SP}.mito.gff" \
  --fasta "${DATA_ROOT}/data/mito/${SP}/${SP}.mito.ctg.fasta" \
  --out "${DATA_ROOT}/check_arabidopsis.tsv"
