#!/usr/bin/env bash
set -euo pipefail
prefix="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
runtime="$prefix/share/eukbinpicker/runtime/0.1"
if [[ ! -f "$runtime/share/eukbinpicker/installation.json" ]]; then
    echo 'EukBinPicker dependency installation is incomplete; repeat the installation command.' >&2
    exit 1
fi
export PATH="$runtime/bin:$PATH"
exec "$runtime/bin/python" "$runtime/bin/eukbinpicker" "$@"
