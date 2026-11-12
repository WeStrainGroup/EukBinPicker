#!/usr/bin/env bash
# Install from this checkout, or fetch the matching release when piped to bash.
set -euo pipefail
VERSION=0.2
CONDA_BIN="${CONDA_EXE:-$(command -v conda || true)}"
[[ -n "$CONDA_BIN" ]] || { echo 'Conda is required.' >&2; exit 2; }
BASE="$("$CONDA_BIN" info --base)"
SOURCE="$(cd "$(dirname "${BASH_SOURCE[0]:-/dev/stdin}")" && pwd)"
STAGE=""
trap '[[ -z "$STAGE" ]] || rm -rf -- "$STAGE"' EXIT
if [[ ! -f "$SOURCE/scripts/install.py" ]]; then
    STAGE="$(mktemp -d)"
    curl -fsSL --retry 3 "https://github.com/WeStrainGroup/EukBinPicker/archive/refs/tags/v${VERSION}.tar.gz" -o "$STAGE/source.tar.gz"
    tar -xzf "$STAGE/source.tar.gz" -C "$STAGE"
    SOURCE="$STAGE/EukBinPicker-${VERSION}"
fi
"$BASE/bin/python" "$SOURCE/scripts/install.py" --conda "$CONDA_BIN" "$@"
