#!/usr/bin/env bash
set -euo pipefail
source_dir="$PREFIX/share/eukbinpicker/source"
mkdir -p "$source_dir/bin" "$PREFIX/bin"
cp -R scripts envs docs "$source_dir/"
cp README.md LICENSE CITATION.cff "$source_dir/"
cp bin/eukbinpicker "$source_dir/bin/"
cp conda/launcher.sh "$PREFIX/bin/eukbinpicker"
chmod 755 "$PREFIX/bin/eukbinpicker"
