#!/usr/bin/env bash
set -euo pipefail
"$PREFIX/bin/python" "$PREFIX/share/eukbinminer/source/scripts/install.py" \
    --prefix "$PREFIX/share/eukbinminer/runtime/0.40" \
    --conda "$PREFIX/bin/conda"
