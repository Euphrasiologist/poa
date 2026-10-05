#!/usr/bin/env bash
# Optional: builds this project's own Rust tools from source with cargo
# instead of bioconda, at the versions environment.yml pins. Needs a Rust
# toolchain (rustup or `mamba install rust`).
#
#   ./install.sh
set -euo pipefail

if ! command -v cargo >/dev/null 2>&1; then
  echo "[err] cargo not found on PATH - install Rust (https://rustup.rs) or 'mamba install rust'." >&2
  exit 1
fi

echo "[info] installing gfatk, orfedit, transsplice via cargo..."
cargo install --locked gfatk@0.6.1
cargo install --locked orfedit@0.1.0
cargo install --locked --git https://github.com/ARU-life-sciences/transsplice --tag 0.2.1 transsplice

echo "[info] done. 'poa check' (after 'pip install .') confirms every tool is found."
