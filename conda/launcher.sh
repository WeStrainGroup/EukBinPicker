#!/usr/bin/env bash
set -euo pipefail
prefix="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
runtime="$prefix/share/eukbinpicker/runtime/0.4"
if [[ ! -f "$runtime/share/eukbinpicker/installation.json" ]]; then
    echo 'EukBinPicker dependency installation is incomplete; repeat the installation command.' >&2
    exit 1
fi
# Activate the runtime so nested conda run selects the private tool interpreter.
exec "$prefix/bin/conda" run --no-capture-output -p "$runtime" \
    "$runtime/bin/python" "$runtime/bin/eukbinpicker" "$@"
