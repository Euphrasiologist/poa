#!/usr/bin/env bash
# Shared CLI flag parsing for bash orchestration scripts in analysis/.
#
# Usage:
#   source "$(dirname "${BASH_SOURCE[0]}")/../../common/args.sh"
#   parse_common_args "$@"
#   # now use: $SPECIES_LIST $ORGANELLE $FORCE $QC_STATUS
#   # any unrecognised args are left in the REMAINING_ARGS array
#
# Flags:
#   --species-list FILE   plain text file, one species name per line (default: all species)
#   --organelle X          mito | pltd | both (default: both)
#   --force                 recompute everything, ignore existing results
#   --qc-status LIST       comma-separated subset of pass,flag,fail (default: pass)
#   --data-root DIR         override repo root (default: auto-detected from this script's location)

parse_common_args() {
  SPECIES_LIST=""
  ORGANELLE="both"
  FORCE=0
  QC_STATUS="pass"
  REMAINING_ARGS=()

  while [[ $# -gt 0 ]]; do
    case "$1" in
      --species-list) SPECIES_LIST="$2"; shift 2 ;;
      --organelle) ORGANELLE="$2"; shift 2 ;;
      --force) FORCE=1; shift ;;
      --qc-status) QC_STATUS="$2"; shift 2 ;;
      *) REMAINING_ARGS+=("$1"); shift ;;
    esac
  done

  if [[ -n "${SPECIES_LIST}" && ! -f "${SPECIES_LIST}" ]]; then
    echo "[err] --species-list file not found: ${SPECIES_LIST}" >&2
    exit 1
  fi
  case "${ORGANELLE}" in
    mito|pltd|both) ;;
    *) echo "[err] --organelle must be mito|pltd|both, got: ${ORGANELLE}" >&2; exit 1 ;;
  esac

  export SPECIES_LIST ORGANELLE FORCE QC_STATUS
}

# Resolve the analysis/ dir and repo ROOT_DIR from any script under analysis/<module>/src/.
resolve_analysis_root() {
  local script_dir
  script_dir="$(cd "$(dirname "${BASH_SOURCE[1]}")" && pwd)"
  ANALYSIS_DIR="$(cd "${script_dir}/../.." && pwd)"
  ROOT_DIR="$(cd "${ANALYSIS_DIR}/.." && pwd)"
  export ANALYSIS_DIR ROOT_DIR
}
