"""Regression checks for destination host; no native solver is needed here."""
from __future__ import annotations

from concurrent.futures import Future
from pathlib import Path
from types import SimpleNamespace
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from austin_runtime.common import config, cpu_budget
from austin_runtime.power import PowerEngine
from austin_runtime.traffic import TrafficEngine


# The two-line-per-arc layout is emitted by vendor/tap-b/src/fileio.c.
# Parallel arcs must remain separate even when their endpoint IDs are identical.
NATIVE_FLOWS = "(1,2) 3.000000 2.000000 \n3.000000 \n(1,2) 4.000000 2.000000 \n4.000000 \n"
NATIVE_LOG = "Iteration 4: gap 1.000000e-05\nAggregate TSTT: 14.000000\n"


class NativeOutputContract(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        root = Path(self.temporary.name)
        engine = TrafficEngine.__new__(TrafficEngine)
        engine.cfg = config()
        engine.verbose = False
        engine.scratch = root / "scratch"
        engine.scratch.mkdir()
        engine.net = root / "network.tntp"
        engine.net.write_text("\n".join(["~ test header"] * 8) + "\n")
        engine.trips = root / "trips.tntp"
        engine.binary = root / "tap"
        engine.binary_hash = "fixture"
        engine.supply = np.array([0., 7., -7.])
        arc = dict(from_node="1", to_node="2", capacity_source="100", length_source="1",
                   free_flow_time_source="2", bpr_alpha="0.15", bpr_beta="4",
                   speed_source="0", toll_source="0", link_type="1")
        engine.links = [dict(arc, link_id="1"), dict(arc, link_id="2")]
        self.engine = engine

    def native_call(self, contents=NATIVE_FLOWS, filename="flows.txt", *, timeout=False):
        def emit(args, *, cwd, stdout, **kwargs):
            stdout.write(NATIVE_LOG)
            stdout.flush()
            if contents is not None:
                (Path(cwd) / filename).write_text(contents)
            if timeout:
                raise subprocess.TimeoutExpired(args, kwargs["timeout"])
            return SimpleNamespace(returncode=0)
        return emit

    def solve(self):
        return self.engine._solve(np.ones(2), 9999., "fixture")

    def retained_folder(self):
        folders = list(self.engine.scratch.glob("tapb-*"))
        self.assertEqual(len(folders), 1)
        self.assertTrue((folders[0] / "network.tntp").is_file())
        self.assertEqual((folders[0] / "tapb.log").read_text(), NATIVE_LOG)
        self.assertTrue((folders[0] / "failure.txt").is_file())
        return folders[0]

    def test_reads_upstream_flows_and_preserves_parallel_arcs(self):
        with patch("austin_runtime.traffic.subprocess.run", side_effect=self.native_call()):
            costs, flows, report = self.solve()
        np.testing.assert_array_equal(costs, [2., 2.])
        np.testing.assert_array_equal(flows, [3., 4.])
        self.assertEqual(report["tstt"], 14.)
        self.assertEqual(report["max_node_imbalance"], 0.)
        self.assertFalse(list(self.engine.scratch.glob("tapb-*")))

    def test_wrong_output_name_reports_error_and_keeps_evidence(self):
        with patch("austin_runtime.traffic.subprocess.run", side_effect=self.native_call(filename="s.txt")):
            with self.assertRaisesRegex(RuntimeError, "did not write flows.txt.*diagnostic files retained"):
                self.solve()
        self.assertTrue((self.retained_folder() / "s.txt").is_file())

    def test_arc_identity_failure_keeps_actual_output(self):
        malformed = NATIVE_FLOWS.replace("(1,2)", "(2,1)", 1)
        with patch("austin_runtime.traffic.subprocess.run", side_effect=self.native_call(malformed)):
            with self.assertRaisesRegex(RuntimeError, "changed arc order"):
                self.solve()
        self.assertEqual((self.retained_folder() / "flows.txt").read_text(), malformed)

    def test_timeout_keeps_network_and_partial_log(self):
        with patch("austin_runtime.traffic.subprocess.run", side_effect=self.native_call(None, timeout=True)):
            with self.assertRaisesRegex(RuntimeError, "timed out.*diagnostic files retained"):
                self.solve()
        self.retained_folder()


class ImmediateExecutor:
    """Exercise submission/merge/error handling without starting processes."""
    def __init__(self, **kwargs):
        self.options = kwargs
        self.calls = []
        self.shutdown_calls = []

    def submit(self, function, *args):
        self.calls.append(args)
        future = Future()
        try:
            future.set_result(function(*args))
        except BaseException as exc:
            future.set_exception(exc)
        return future

    def shutdown(self, **kwargs):
        self.shutdown_calls.append(kwargs)


class RegionalScheduling(unittest.TestCase):
    def setUp(self):
        self.states = {
            "R1": dict(region="R1", served_kw=10., nominal_kw=20., feasible=True,
                       zone_served_kw={"1": 10.}, substation_served_kw={"s1": 10.}, signal_powered={"a": True}),
            "R2": dict(region="R2", served_kw=30., nominal_kw=40., feasible=True,
                       zone_served_kw={"1": 30.}, substation_served_kw={"s2": 30.}, signal_powered={"b": False}),
        }

    def engine(self, workers):
        engine = PowerEngine.__new__(PowerEngine)
        engine.catalog = dict(regions=["R1", "R2"], assets={
            "power:s1": {"region": "R1"}, "power:s2": {"region": "R2"}})
        engine.region_workers = workers
        engine.verbose = False
        engine._pool = None
        engine._pool_args = ()
        engine._regional_state = lambda region, damage: self.states[region]
        return engine

    def test_parallel_merge_matches_serial_despite_reversed_completion(self):
        damage = {"power:s1": 0., "power:s2": .5, "road:1-2": 0.}
        serial = self.engine(1).evaluate(damage)
        engine = self.engine(2)
        worker = SimpleNamespace(_regional_state=lambda region, state: self.states[region])
        with patch("austin_runtime.power._regional_engine", worker), \
             patch("austin_runtime.power.futures.ProcessPoolExecutor", ImmediateExecutor), \
             patch("austin_runtime.power.futures.as_completed", side_effect=lambda pending: reversed(list(pending))):
            parallel = engine.evaluate(damage)
            pool = engine._pool
            self.assertEqual(pool.options["max_workers"], 2)
            self.assertEqual(pool.options["mp_context"].get_start_method(), "spawn")
            self.assertEqual(pool.calls, [("R1", {"power:s1": 0.}), ("R2", {"power:s2": .5})])
            engine.close()
        self.assertEqual(parallel, serial)
        self.assertEqual(parallel["zone_served_kw"], {"1": 40.})
        self.assertEqual([r["region"] for r in parallel["regions"]], ["R1", "R2"])
        self.assertEqual(pool.shutdown_calls, [dict(wait=True, cancel_futures=False)])

    def test_fresh_restoration_bypasses_cache_in_existing_worker(self):
        engine = self.engine(2)
        pool = engine._pool = ImmediateExecutor()
        called = []
        def fresh(region, damage):
            called.append((region, damage))
            return self.states[region]
        worker = SimpleNamespace(_operate=fresh,
                                 _regional_state=lambda *args: self.fail("fresh restoration used cache"))
        with patch("austin_runtime.power._regional_engine", worker):
            self.assertEqual(engine.fresh_region("R1", {}), self.states["R1"])
        self.assertEqual(called, [("R1", {})])
        self.assertEqual(pool.calls, [("R1", {}, True)])
        engine.close()

    def test_region_failure_closes_pool_and_does_not_return_partial_power(self):
        engine = self.engine(2)
        pool = engine._pool = ImmediateExecutor()
        def fail(region, damage):
            if region == "R2": raise ValueError("regional infeasibility")
            return self.states[region]
        with patch("austin_runtime.power._regional_engine", SimpleNamespace(_regional_state=fail)):
            with self.assertRaisesRegex(ValueError, "regional infeasibility"):
                engine.evaluate({})
        self.assertIsNone(engine._pool)
        self.assertEqual(pool.shutdown_calls, [dict(wait=True, cancel_futures=True)])

    def test_validation_pool_uses_ram_budget_independently_of_smoke_workers(self):
        cfg = config()
        with patch("austin_runtime.common.os.sched_getaffinity", return_value=set(range(32))), \
             patch("austin_runtime.common.Path.exists", lambda path: str(path) == "/proc/meminfo"), \
             patch("austin_runtime.common.Path.read_text", return_value=f"MemAvailable: {95 * 2**20} kB\n"):
            budget = cpu_budget(cfg, requested=2)
            self.assertEqual(budget["workers"], 2)
            self.assertEqual(budget["validation_power_workers"], 5)
            cfg["runtime"]["validation_power_workers"] = 6
            with self.assertRaisesRegex(ValueError, "validation_power_workers"):
                cpu_budget(cfg, requested=2)


if __name__ == "__main__":
    unittest.main()
