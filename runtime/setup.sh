#!/usr/bin/env bash
# Only runs when explicitly invoked on the destination Linux host.
set -euo pipefail
austin_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
python_bin="${AUSTIN_PYTHON:-python3}"
for program in "$python_bin" gcc make curl; do
    command -v "$program" >/dev/null || { printf 'Missing program: %s\n' "$program" >&2; exit 1; }
done
"$python_bin" -c 'import sys; assert sys.version_info >= (3, 11), "Python >= 3.11 required"; assert sys.platform == "linux", "Linux required"'
"$python_bin" -m venv "$austin_dir/.venv-runtime"
"$austin_dir/.venv-runtime/bin/python" -m pip install --upgrade pip
"$austin_dir/.venv-runtime/bin/python" -m pip install -e "$austin_dir/runtime"
mkdir -p "$austin_dir/runtime/build/tap-b"
cp -a "$austin_dir/vendor/tap-b/." "$austin_dir/runtime/build/tap-b/"
make -C "$austin_dir/runtime/build/tap-b" clean
make -C "$austin_dir/runtime/build/tap-b" -j"${AUSTIN_BUILD_JOBS:-4}" parallel
"$austin_dir/.venv-runtime/bin/python" -m pip freeze > "$austin_dir/runtime/build/environment.freeze.txt"
printf 'Installed and compiled. No experiments started. Follow runtime/README.md for data preparation and acceptance.\n'
