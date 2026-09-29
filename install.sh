#!/usr/bin/env bash
# Installs this project's own Rust command-line tools into an already
# -activated conda/mamba environment (see environment.yml, which covers
# every third-party bioinformatics tool this pipeline needs).
#
# Usage:
#   mamba env create -f environment.yml
#   conda activate plant_organelle_annotator
#   ./install.sh
#
# Reference data: trans_splicing's exon templates and editing's profiles
# are bundled in this repo. oatkDB's gene-family .fam databases are not -
# see README.md's "Reference data" section.
set -euo pipefail

if ! command -v cargo >/dev/null 2>&1; then
  echo "[err] cargo not found on PATH." >&2
  echo "      Run 'mamba env create -f environment.yml && conda activate plant_organelle_annotator' first" >&2
  echo "      (that env includes the rust/cargo toolchain)." >&2
  exit 1
fi

echo "[info] installing this project's own Rust tools via cargo..."

# Pinned to the versions this pipeline was validated against (tags in each
# tool's repo). gfatk, orfedit and transsplice come from crates.io;
# hmm_to_gff and filter_tblout aren't published there, so they're built
# from their GitHub release tags.
cargo install --locked gfatk@0.4.0
cargo install --locked orfedit@0.1.0
cargo install --locked transsplice@0.2.1
cargo install --locked --git https://github.com/ARU-life-sciences/hmm_to_gff --tag 0.1.1 hmm_to_gff
cargo install --locked --git https://github.com/ARU-life-sciences/filter_tblout --tag 0.1.1 filter_tblout

echo "[info] done. Verifying tool_paths.sh resolves everything cleanly:"
bash -c 'source analysis/common/tool_paths.sh'
echo "[info] if no [warn] lines appeared above, every tool this pipeline needs is on PATH."
