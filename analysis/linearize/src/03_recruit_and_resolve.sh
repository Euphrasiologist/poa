#!/usr/bin/env bash
# Submits one bsub job per species x organelle row in a resolve_targets.tsv
# (from 02_select_resolve_targets.py) to run the gfatk-resolve fallback tier:
# minimap2 read recruitment -> GraphAligner -> gfatk resolve (see
# _resolve_one.sh for what each job actually does, and analysis/linearize/README.md
# for why this two-step recruit-then-align design exists instead of running
# GraphAligner against the raw whole-genome read sets directly).
set -euo pipefail

SRC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANALYSIS_DIR="$(cd "${SRC_DIR}/../.." && pwd)"
ROOT_DIR="$(cd "${ANALYSIS_DIR}/.." && pwd)"

TARGETS="${ANALYSIS_DIR}/linearize/work/resolve_targets.tsv"
META="${ROOT_DIR}/meta/file.txt"
WORK_DIR="${ANALYSIS_DIR}/linearize/work"
RESULTS_DIR="${ANALYSIS_DIR}/linearize/results/resolve"
LOGDIR="${ANALYSIS_DIR}/linearize/logs"
THREADS=16
MEM_MB=32000
QUEUE="normal"

usage() {
  cat <<EOF
Usage:
  $0 [options]

Options:
  --targets <tsv>   (default: ${TARGETS})
  --meta <file>     (default: ${META})
  --threads <int>   (default: ${THREADS})
  --mem-mb <int>    (default: ${MEM_MB})
  --queue <name>    (default: ${QUEUE})
  --limit <int>     only submit the first N rows (for piloting)
EOF
}

LIMIT=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --targets) TARGETS="$2"; shift 2;;
    --meta) META="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --mem-mb) MEM_MB="$2"; shift 2;;
    --queue) QUEUE="$2"; shift 2;;
    --limit) LIMIT="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "[ERROR] Unknown arg: $1" >&2; usage; exit 1;;
  esac
done

[[ -f "${TARGETS}" ]] || { echo "[ERROR] targets not found: ${TARGETS} (run 02_select_resolve_targets.py first)" >&2; exit 1; }
[[ -f "${META}" ]] || { echo "[ERROR] meta not found: ${META}" >&2; exit 1; }

mkdir -p "${WORK_DIR}" "${RESULTS_DIR}" "${LOGDIR}"

submitted=0
skipped=0

# TARGETS has a header line: species<TAB>organelle<TAB>gfa
# Process substitution (not a pipe) so submitted/skipped survive the loop -
# `... | while read` would run the loop in a subshell and silently reset
# both counters to 0 on exit.
while IFS=$'\t' read -r species organelle gfa; do
  [[ -n "${species}" ]] || continue

  if (( LIMIT > 0 && submitted >= LIMIT )); then
    break
  fi

  out_fasta="${RESULTS_DIR}/${organelle}/${species}.resolve.fasta"
  if [[ -f "${out_fasta}" && "${out_fasta}" -nt "${gfa}" ]]; then
    ((skipped+=1))
    continue
  fi

  # Named so a rerun of this script (e.g. after --targets grows with newly
  # QC-caught species) can't silently double-submit a species whose job is
  # still queued/running - checking the output file alone isn't enough,
  # since an in-flight job hasn't written it yet. Caught this exact bug
  # directly: an early run submitted 35 real duplicate jobs this way.
  JOB_NAME="linearize_resolve.${species}.${organelle}"
  if bjobs -J "${JOB_NAME}" >/dev/null 2>&1; then
    echo "[SKIP] ${species} (${organelle}): already queued/running as ${JOB_NAME}"
    ((skipped+=1))
    continue
  fi

  bsub -J "${JOB_NAME}" -n "${THREADS}" -q "${QUEUE}" \
    -R"span[hosts=1] select[mem>${MEM_MB}] rusage[mem=${MEM_MB}]" \
    -M"${MEM_MB}" \
    -o "${LOGDIR}/${species}.${organelle}.out" \
    -e "${LOGDIR}/${species}.${organelle}.err" \
    "${SRC_DIR}/_resolve_one.sh" "${species}" "${organelle}" "${gfa}" "${META}" "${WORK_DIR}" "${RESULTS_DIR}" "${THREADS}"

  echo "[SUBMIT] ${species} (${organelle})"
  ((submitted+=1))
done < <(tail -n +2 "${TARGETS}")

echo "[DONE] submitted=${submitted} skipped=${skipped}"
