#!/usr/bin/env bash
set -euo pipefail
source_dir="$PREFIX/share/eukbinminer/source"
mkdir -p "$source_dir/bin" "$PREFIX/bin"
cp -R scripts envs "$source_dir/"
cp README.md LICENSE CITATION.cff "$source_dir/"
cp bin/eukbinminer "$source_dir/bin/"
cp conda/launcher.sh "$PREFIX/bin/eukbinminer"
chmod 755 "$PREFIX/bin/eukbinminer"
