PY := .venv/bin/python
.PHONY: setup sources build preview validate-road validate-power verify package validate-signal-scenario
setup:
	python3 -m venv .venv
	$(PY) -m pip install -r requirements.txt
sources:
	$(PY) scripts/fetch_sources.py
build: sources
	$(PY) scripts/build_dataset.py
preview:
	$(PY) scripts/make_preview.py
validate-road:
	$(PY) scripts/validate_solvers.py --tapb vendor/tap-b/bin/tap
validate-signal-scenario:
	$(PY) scripts/validate_signal_scenario.py
validate-power:
	$(PY) scripts/validate_solvers.py --opendss --all-feeders --allow-known-source-exceptions
	$(PY) scripts/check_power_exceptions.py
verify:
	$(PY) scripts/verify_release.py
package: verify
	$(PY) scripts/package_release.py
