"""Destination-host checks for smoke routing and reuse; no native solver calls."""
from __future__ import annotations

import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from austin_runtime.common import atomic_json, config, read_json


spec = importlib.util.spec_from_file_location("austin_smoke_entry", Path(__file__).resolve().parents[1] / "smoke.py")
smoke = importlib.util.module_from_spec(spec)
spec.loader.exec_module(smoke)


class SmokeAcceptanceOnly(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.output = Path(temporary.name) / "output"
        prepared = Path(temporary.name) / "prepared"
        atomic_json(prepared / "catalog.json", {})
        self.cfg = config()
        self.cfg["runtime"].update(output=str(self.output), prepared=str(prepared))
        self.budget = dict(workers=2, validation_power_workers=4)
        for name, value in (("verify_inputs", True), ("fingerprint", "same-physics"), ("cpu_budget", self.budget)):
            patcher = patch.object(smoke, name, return_value=value)
            patcher.start()
            self.addCleanup(patcher.stop)
        atomic_json(self.output / "run_manifest.json", dict(fingerprint="same-physics", existing_run=True))

    def test_reuses_passed_validation_and_preserves_existing_experiment(self):
        atomic_json(self.output / "validation.json", dict(fingerprint="same-physics", passed=True, checks=["physical"]))
        manifest = (self.output / "run_manifest.json").read_bytes()
        with patch.object(smoke, "validate") as validate:
            result = smoke.run_smoke(self.cfg)
        validate.assert_not_called()
        self.assertTrue(result["validation_reused"])
        self.assertEqual(result["experiment_jobs"], 0)
        self.assertFalse(result["full_research_validated"])
        self.assertEqual((self.output / "run_manifest.json").read_bytes(), manifest)
        self.assertFalse((self.output / "results").exists())

    def test_failed_validation_is_retried_then_smoke_stops(self):
        atomic_json(self.output / "validation.json", dict(fingerprint="same-physics", passed=False))
        with patch.object(smoke, "validate", return_value=dict(passed=True, checks=["physical"])) as validate:
            result = smoke.run_smoke(self.cfg, 2)
        validate.assert_called_once_with(self.cfg, "same-physics", self.budget)
        self.assertTrue(result["passed"])
        self.assertFalse(result["validation_reused"])
        self.assertFalse((self.output / "results").exists())

    def test_failure_is_not_reported_as_pass(self):
        with patch.object(smoke, "validate", side_effect=RuntimeError("AC failed")):
            with self.assertRaisesRegex(RuntimeError, "AC failed"):
                smoke.run_smoke(self.cfg)
        result = read_json(self.output / "smoke_summary.json")
        self.assertFalse(result["passed"])
        self.assertEqual(result["status"], "failed")

    def test_stale_validation_is_not_overwritten(self):
        atomic_json(self.output / "validation.json", dict(fingerprint="other-physics", passed=True))
        before = (self.output / "validation.json").read_bytes()
        with patch.object(smoke, "validate") as validate:
            with self.assertRaisesRegex(ValueError, "Different numerical fingerprint"):
                smoke.run_smoke(self.cfg)
        validate.assert_not_called()
        self.assertEqual((self.output / "validation.json").read_bytes(), before)
        self.assertFalse((self.output / "smoke_summary.json").exists())


if __name__ == "__main__":
    unittest.main()
