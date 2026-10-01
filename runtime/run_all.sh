#!/usr/bin/env bash
set -euo pipefail
austin_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
profile="${1:-research}"
if [[ "$profile" != smoke && "$profile" != research ]]; then
    printf 'Usage: bash runtime/run_all.sh [smoke|research] [workers]\n' >&2
    exit 2
fi
cd "$austin_dir"
if [[ ! -x .venv-runtime/bin/python ]]; then
    printf 'Run bash runtime/bootstrap.sh on this host first.\n' >&2
    exit 1
fi
extra=()
if [[ -n "${2:-}" ]]; then extra=(--workers "$2"); fi
exec .venv-runtime/bin/python runtime/run.py --config "runtime/configs/$profile.toml" run --stage all "${extra[@]}"
