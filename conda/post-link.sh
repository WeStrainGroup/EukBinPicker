#!/usr/bin/env bash
set -euo pipefail
"$PREFIX/bin/python" "$PREFIX/share/eukbinpicker/source/scripts/install.py" \
    --prefix "$PREFIX/share/eukbinpicker/runtime/0.4" \
    --conda "$PREFIX/bin/conda"
