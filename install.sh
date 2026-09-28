#!/usr/bin/env bash
# Installs this project's own Rust command-line tools into an already
# -activated conda/mamba environment (see environment.yml, which covers
# every third-party bioinformatics tool this pipeline needs).
#
# Usage:
#   mamba env create -f environment.yml
#   conda activate plant_organellar_database
#   ./install.sh
#
# Reference databases (oatkDB's mito/plastid gene-family HMMs, and the
# editing/trans_splicing reference profiles) are NOT handled here yet -
# see README.md's "Reference data" section for the current status.
set -euo pipefail

if ! command -v cargo >/dev/null 2>&1; then
  echo "[err] cargo not found on PATH." >&2
  echo "      Run 'mamba env create -f environment.yml && conda activate plant_organellar_database' first" >&2
  echo "      (that env includes the rust/cargo toolchain)." >&2
  exit 1
fi

echo "[info] installing this project's own Rust tools via cargo..."

# Pinned to the versions this pipeline was built and validated against.
# gfatk, hmm_to_gff, filter_tblout are already on crates.io; orfedit and
# transsplice will be once first published (see each repo's README in the
# meantime - ARU-life-sciences/orfedit, ARU-life-sciences/transsplice).
cargo install --locked gfatk@0.4.0
cargo install --locked hmm_to_gff@0.1.1
cargo install --locked filter_tblout@0.1.1
cargo install --locked orfedit@0.1.0
cargo install --locked transsplice@0.1.0

echo "[info] done. Verifying tool_paths.sh resolves everything cleanly:"
bash -c 'source analysis/common/tool_paths.sh'
echo "[info] if no [warn] lines appeared above, every tool this pipeline needs is on PATH."
