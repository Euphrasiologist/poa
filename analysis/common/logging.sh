#!/usr/bin/env bash
# Shared logging helpers, matching the [info]/[warn]/[err] style used in src/*.sh.

log_info() { echo "[info] $*" >&2; }
log_warn() { echo "[warn] $*" >&2; }
log_err()  { echo "[err] $*" >&2; }
