"""Native check of the local power operator on one real TAMU region (opt-in, a few minutes).

Run after `cli --config configs/core-v2.toml prepare`:
    AUSTIN_NATIVE_TESTS=1 python -m unittest tests.test_power_local_native -v
Region P2U (14 substations; five transformers and 21 lines above nameplate in the base case).
"""
from __future__ import annotations
import os, random, tempfile, unittest
from pathlib import Path
from austin_runtime.common import RUNTIME, config, local_path, read_json

CONFIG = RUNTIME/"configs/core-v2.toml"
REGION = os.environ.get("AUSTIN_NATIVE_REGION", "P2U")


@unittest.skipUnless(os.environ.get("AUSTIN_NATIVE_TESTS") == "1", "set AUSTIN_NATIVE_TESTS=1 to start OpenDSS")
class LocalOperatorOnTamuRegion(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from austin_runtime.power import PowerEngine
        cls.folder = tempfile.TemporaryDirectory()
        cfg = config(CONFIG); cfg["power"]["isolate_states"] = False; cfg["runtime"]["cache"] = cls.folder.name
        catalog_path = local_path(cfg, "prepared")/"catalog.json"
        if not catalog_path.exists(): raise unittest.SkipTest(f"prepare {CONFIG.name} first: {catalog_path}")
        cls.catalog = read_json(catalog_path)
        cls.engine = PowerEngine(cfg, cls.catalog, "native-unit-test", Path(cls.folder.name)/"scratch")
        cls.subs = sorted(k for k, a in cls.catalog["assets"].items() if a["kind"] == "power" and a["region"] == REGION)
        cls.healthy = cls.engine._operate(REGION, {})

    @classmethod
    def tearDownClass(cls):
        cls.engine.close(); cls.folder.cleanup()

    def test_healthy_base_case_serves_nominal_load(self):
        h = self.healthy
        self.assertTrue(h["feasible"]); self.assertEqual(h["shed_loads"], 0); self.assertFalse(h["regional_fallback"])
        self.assertAlmostEqual(h["served_kw"]/h["nominal_kw"], 1.0, places=6)

    def test_partial_derating_sheds_only_its_own_customers(self):
        for sid in self.subs[:4]:
            state = self.engine._operate(REGION, {sid: 0.4}); own = sid.split(":", 1)[1]
            for other, served in self.healthy["substation_served_kw"].items():
                if other != own: self.assertAlmostEqual(state["substation_served_kw"].get(other, 0.0), served, delta=1.0, msg=f"{sid} shed {other}")
            self.assertLessEqual(state["substation_served_kw"].get(own, 0.0), self.healthy["substation_served_kw"][own] + 1e-6)

    def test_repairs_never_reduce_served_load(self):
        rng = random.Random(20261007)
        for _ in range(2):
            picked = rng.sample(self.subs, 4)
            remaining = {k: (0.0 if rng.random() < .5 else rng.uniform(.2, .8)) for k in picked}
            previous = self.engine._operate(REGION, remaining)["served_kw"]
            for asset in rng.sample(picked, len(picked)):
                remaining.pop(asset); current = self.engine._operate(REGION, remaining)["served_kw"]
                self.assertGreaterEqual(current, previous - 1e-3, f"repairing {asset} reduced supply")
                previous = current
        self.assertAlmostEqual(previous, self.healthy["served_kw"], delta=1e-3)


if __name__ == "__main__":
    unittest.main()
