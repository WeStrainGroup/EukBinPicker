#!/usr/bin/env bash
set -euo pipefail
prefix="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
runtime="$prefix/share/eukbinminer/runtime/0.40"
if [[ ! -f "$runtime/share/eukbinminer/installation.json" ]]; then
    echo 'EukBinMiner dependency installation is incomplete; repeat the installation command.' >&2
    exit 1
fi
export PATH="$runtime/bin:$PATH"
exec "$runtime/bin/python" "$runtime/bin/eukbinminer" "$@"
