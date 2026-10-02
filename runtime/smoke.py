"""Explicit physical acceptance only; never launches the A/B/C experiment matrix.

This entry point stays outside austin_runtime so routing-only changes do not
invalidate the numerical fingerprint or existing physical state caches.
"""
from __future__ import annotations

import argparse
import importlib.metadata
import json
import os
import sys
import time

from austin_runtime.common import (
    RUNTIME, atomic_json, config, cpu_budget, file_lock, fingerprint, local_path, read_json,
)
from austin_runtime.progress import phase, report
from austin_runtime.validation import validate, verify_inputs


SCOPE = "Six-region AC and TAP-B acceptance only; no A/B/C experiments or analysis"


def smoke_plan(cfg):
    return dict(scope=SCOPE, experiment_jobs=0, construction_permutations=0,
                full_research_validated=False, output=str(local_path(cfg, "output")),
                checks=["healthy baseline", "one full substation fault",
                        "fresh restoration round trip", "partial transformer derating"],
                note="A cold validation still uses complete circuits; no elapsed-time guarantee.")


def run_smoke(cfg, requested=None):
    catalog_path = local_path(cfg, "prepared") / "catalog.json"
    if not catalog_path.exists():
        raise FileNotFoundError("Run the explicit prepare step first")
    with phase("smoke: verify prepared data and model checksums"):
        catalog = read_json(catalog_path)
        verify_inputs(cfg, catalog)
    fp = fingerprint(cfg, catalog)
    budget = cpu_budget(cfg, requested)
    output = local_path(cfg, "output")
    output.mkdir(parents=True, exist_ok=True)
    with file_lock(output / "run.lock", blocking=False):
        manifest_path = output / "run_manifest.json"
        gate = output / "validation.json"
        saved = read_json(gate) if gate.exists() else None
        for path, previous in (
            (manifest_path, read_json(manifest_path) if manifest_path.exists() else None),
            (gate, saved),
        ):
            if previous is not None and previous.get("fingerprint") != fp:
                raise ValueError(f"Different numerical fingerprint in {path}; use a new output directory")
        if not manifest_path.exists():
            versions = {k: importlib.metadata.version(k)
                        for k in ("numpy", "scipy", "networkx", "OpenDSSDirect.py")}
            atomic_json(manifest_path, dict(fingerprint=fp, config=cfg, budget=budget,
                        python=sys.version, packages=versions, validation_required=True,
                        source="Austin physical smoke; A/B/C experiment jobs are not scheduled"))
        summary = dict(smoke_plan(cfg), fingerprint=fp, status="running", passed=False,
                       validation_reused=False, budget=budget)
        summary_path = output / "smoke_summary.json"
        atomic_json(summary_path, summary)
        started = time.monotonic()
        report(f"smoke: 0 experiment jobs; AC validation uses {budget['validation_power_workers']} processes")
        try:
            if saved is not None and saved.get("passed") is True:
                summary["validation_reused"] = True
                report("smoke: matching physical validation already passed; reusing it")
            else:
                saved = validate(cfg, fp, budget)
                if saved.get("passed") is not True:
                    raise RuntimeError("Physical acceptance did not pass")
            summary.update(status="passed", passed=True,
                           validation_report=str(gate),
                           checks_passed=saved.get("checks", []))
            report("smoke: physical acceptance passed; finished (no construct, optimization or sensitivity jobs)")
            return summary
        except BaseException as exc:
            summary.update(status="failed", error=f"{type(exc).__name__}: {exc}")
            raise
        finally:
            summary["invocation_wall_seconds"] = time.monotonic() - started
            atomic_json(summary_path, summary)


def main(argv=None):
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[name] = "1"
    parser = argparse.ArgumentParser(description=SCOPE)
    parser.add_argument("--config", default=str(RUNTIME / "configs/smoke.toml"))
    parser.add_argument("--output", help="Absolute path or path relative to runtime/")
    parser.add_argument("--workers", type=int,
                        help="Existing resource-budget argument; AC processes use validation_power_workers")
    parser.add_argument("--plan", action="store_true", help="Describe acceptance scope without native solvers")
    args = parser.parse_args(argv)
    cfg = config(args.config)
    if args.output:
        cfg["runtime"]["output"] = args.output
    result = smoke_plan(cfg) if args.plan else run_smoke(cfg, args.workers)
    print(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
