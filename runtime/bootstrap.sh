#!/usr/bin/env bash
# Destination-host preparation only. This file is inert until explicitly invoked.
set -euo pipefail
austin_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$austin_dir"
bash runtime/setup.sh
.venv-runtime/bin/python scripts/fetch_sources.py --download-missing
if [[ ! -f data/processed/coupling/signal_road_power_provisional.csv ]]; then
    printf 'Frozen processed dataset is missing. Include data/processed/ in the clone; see runtime/README.md.\n' >&2
    exit 1
fi
.venv-runtime/bin/python -m unittest discover -s runtime/tests -v
.venv-runtime/bin/python runtime/run.py prepare
.venv-runtime/bin/python runtime/run.py doctor
.venv-runtime/bin/python runtime/run.py plan
printf 'Preparation finished. Explicitly run: bash runtime/run_all.sh smoke\n'
