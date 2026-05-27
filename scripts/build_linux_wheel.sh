#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "$script_dir/.." && pwd)"
stage="$root/build/wheel-linux"
out_dir="$root/dist/linux"

case "$stage" in
  "$root"/build/*) ;;
  *) echo "Refusing to remove unexpected path: $stage" >&2; exit 1 ;;
esac

rm -rf "$stage"
mkdir -p "$stage"
mkdir -p "$out_dir"

cp -a \
  "$root/bfee_docking" \
  "$root/pyproject.toml" \
  "$root/setup.py" \
  "$root/MANIFEST.in" \
  "$root/README.md" \
  "$root/LICENSE" \
  "$stage/"

rm -f "$stage/bfee_docking/third_party/vina/vina.exe"
rm -f "$stage/bfee_docking/third_party/smina/smina.exe"
find "$stage/bfee_docking/third_party" -type f \( -name "*.dll" -o -name "*.bat" -o -name "*.jar" \) -delete
rm -rf "$stage/bfee_docking/third_party/obabel"

export BFEE_DOCKING_PLATFORM_TAG=manylinux2014_x86_64
python -m build --wheel --no-isolation --outdir "$out_dir" "$stage"
