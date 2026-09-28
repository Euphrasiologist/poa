#!/usr/bin/env bash
# Per-species worker for the gfatk-resolve fallback tier: recruit organelle-relevant
# reads via minimap2, align them to the graph with GraphAligner, then let gfatk
# resolve use that read-path evidence to pick a circuit. Not meant to be run
# directly at scale - 03_recruit_and_resolve.sh submits one bsub job per species
# x organelle that runs this.
#
# Why recruit first: whole-genome raw HiFi read sets for this dataset total
# ~24TB compressed (verified via meta/file.txt) - GraphAligner against the raw
# set would scan all of that per species just to align the handful of MB that
# actually belong to one small organelle graph. minimap2 -x map-hifi against
# the graph's own segment sequences (gfatk fasta, works even when oatk's own
# Pathfinder never produced a ctg.fasta - that's exactly this tier's target
# population) cheaply shrinks that down first.
#
# Recruitment + extraction in ONE pass: minimap2 -a (SAM, which embeds each
# read's sequence) piped straight into `samtools fasta -F 0x904` (drop
# unmapped/secondary/supplementary, keep one primary record per read).
# Measured directly on a real species (Carlina_vulgaris, ~30GB compressed
# reads): a first version that ran minimap2 -x map-hifi for a PAF then a
# second full zcat+seqkit-grep pass to pull matching reads cost 33.5min +
# 31.7min = ~65min, almost entirely spent re-decompressing the same reads
# twice (minimap2 itself only got ~2.7x parallelism out of 16 threads,
# consistent with gzip decompression - not alignment - being the real
# bottleneck). Piping avoids the second pass entirely.
set -euo pipefail

usage() {
  echo "Usage: $0 <species> <organelle> <gfa> <meta_file> <work_dir> <results_dir> <threads>" >&2
  exit 1
}
[[ $# -eq 7 ]] || usage

SPECIES="$1"; ORGANELLE="$2"; GFA="$3"; META="$4"; WORK_DIR="$5"; RESULTS_DIR="$6"; THREADS="$7"

SRC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SRC_DIR}/../../common/tool_paths.sh"

mkdir -p "${WORK_DIR}" "${RESULTS_DIR}/${ORGANELLE}"
PREFIX="${WORK_DIR}/${SPECIES}.${ORGANELLE}"

echo "[info] ${SPECIES} (${ORGANELLE}): building segment reference from GFA"
"${GFATK}" fasta "${GFA}" > "${PREFIX}.segs.fasta"

mapfile -t READ_FILES < <(awk -F'\t' -v s="${SPECIES}" '$1==s {print $2}' "${META}")
if (( ${#READ_FILES[@]} == 0 )); then
  echo "[err] ${SPECIES} (${ORGANELLE}): no raw read files in ${META}" >&2
  exit 1
fi
echo "[info] ${SPECIES} (${ORGANELLE}): recruiting from ${#READ_FILES[@]} read file(s) (one-pass align+extract)"

"${MINIMAP2}" -a -x map-hifi -t "${THREADS}" "${PREFIX}.segs.fasta" "${READ_FILES[@]}" 2> "${PREFIX}.recruit.log" \
  | "${SAMTOOLS}" fasta -F 0x904 - 2> "${PREFIX}.extract.log" \
  | gzip > "${PREFIX}.recruited.fasta.gz"

N_RECRUITED=$(zcat "${PREFIX}.recruited.fasta.gz" | grep -c "^>" || true)
echo "[info] ${SPECIES} (${ORGANELLE}): recruited ${N_RECRUITED} reads"

if (( N_RECRUITED == 0 )); then
  echo "[warn] ${SPECIES} (${ORGANELLE}): zero reads recruited - nothing maps to this graph, skipping GraphAligner/resolve" >&2
  exit 0
fi

echo "[info] ${SPECIES} (${ORGANELLE}): running GraphAligner"
"${GRAPHALIGNER}" -g "${GFA}" -f "${PREFIX}.recruited.fasta.gz" -a "${PREFIX}.gaf" -x vg -t "${THREADS}" \
  > "${PREFIX}.ga.out" 2> "${PREFIX}.ga.err"

N_ALIGNED=$(wc -l < "${PREFIX}.gaf" || echo 0)
echo "[info] ${SPECIES} (${ORGANELLE}): ${N_ALIGNED} read alignments in GAF"

echo "[info] ${SPECIES} (${ORGANELLE}): running gfatk resolve"
"${GFATK}" resolve --gaf "${PREFIX}.gaf" "${GFA}" \
  > "${RESULTS_DIR}/${ORGANELLE}/${SPECIES}.resolve.fasta" 2> "${PREFIX}.resolve.log"

echo "[ok] ${SPECIES} (${ORGANELLE}): done -> ${RESULTS_DIR}/${ORGANELLE}/${SPECIES}.resolve.fasta"
