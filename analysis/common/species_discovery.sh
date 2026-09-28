#!/usr/bin/env bash
# Cached wrapper around species_discovery.py --print-table.
# Usage: species_discovery.sh <mito|pltd> [species_list_file]
# Prints a TSV to stdout: species organelle gfa ctg_bed ctg_fasta bed annot run_prefix n_candidates status
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"

ORGANELLE="${1:?usage: species_discovery.sh <mito|pltd> [species_list_file]}"
SPECIES_LIST="${2:-}"

DATA_ROOT="${ROOT_DIR}/data/$([[ "${ORGANELLE}" == mito ]] && echo mito || echo plastid)"
CACHE_DIR="${SCRIPT_DIR}/.cache"
mkdir -p "${CACHE_DIR}"
CACHE_FILE="${CACHE_DIR}/species_${ORGANELLE}.tsv"

# species_list-filtered calls are cheap enough (and filter-dependent) to never cache
if [[ -n "${SPECIES_LIST}" ]]; then
  python3 "${SCRIPT_DIR}/species_discovery.py" --organelle "${ORGANELLE}" \
    --root-dir "${ROOT_DIR}" --species-list "${SPECIES_LIST}" --print-table
  exit 0
fi

if [[ -f "${CACHE_FILE}" ]] && [[ -z "$(find "${DATA_ROOT}" -newer "${CACHE_FILE}" -print -quit 2>/dev/null)" ]]; then
  cat "${CACHE_FILE}"
else
  python3 "${SCRIPT_DIR}/species_discovery.py" --organelle "${ORGANELLE}" \
    --root-dir "${ROOT_DIR}" --print-table > "${CACHE_FILE}.tmp"
  mv "${CACHE_FILE}.tmp" "${CACHE_FILE}"
  cat "${CACHE_FILE}"
fi
