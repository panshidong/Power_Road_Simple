from __future__ import annotations

"""Probe for Task A hypothesis 2 (does an interdependency-aware joint ranking beat
separate rankings?). Under the current Task A design the two rankings are executed
by two specialized crews working in parallel, so the joint list's cross-system
interleaving never affects execution. This probe re-executes the SAME rankings with
ONE multifunction crew (as in the proposal's S0 power-first / S1 road-first / S2
joint design and in Pan & Boyles 2025), where interleaving does matter.

Run with cwd = an isolated worker copy (fixed scratch files). Output: CSV + summary.
"""

import csv
import os
import statistics as st
import sys

from disaster import generate_scenarios
from resilience_measurement import run_model_multi
from task_a_criticality import TaskACriticalityConfig, build_task_a_strategies

N = int(sys.argv[1]) if len(sys.argv) > 1 else 20
OUT_DIR = sys.argv[2] if len(sys.argv) > 2 else "/home/workenv/TaskB_powerfix/results/taskA_single_crew_probe"
os.makedirs(OUT_DIR, exist_ok=True)

scenarios = generate_scenarios(n_scenarios=N, seed0=20260701, bus_count_range=(8, 15), link_count_range=(3, 13),
                               link_drop_range=(0.5, 1.0), out_json=os.path.join(OUT_DIR, "random_disasters.json"))
cfg = TaskACriticalityConfig()


def run(seq, sid, tag, crew_mode, preserve, link_factors):
    kw = dict(result_root=OUT_DIR, run_dir=os.path.join(OUT_DIR, "runs", f"{sid}_{tag}"), message=tag, Scenario=f"{sid}_{tag}",
              strict=True, save_artifacts=False, preserve_sequence_order=preserve, objective="triangle",
              broken_link_factors=dict(link_factors))
    if crew_mode == "multifunction":
        r = run_model_multi(seq, crew_mode="multifunction", multifunction_crews=1, power_crews=0, road_crews=0, **kw)
    else:
        r = run_model_multi(seq, crew_mode="specialized", power_crews=1, road_crews=1, **kw)
    return float(r["triangle_area"])


rows = []
for i, sc in enumerate(scenarios, 1):
    print(f"[probe] {sc.scenario_id} ({i}/{N})", flush=True)
    strats = {s.strategy_id: s for s in build_task_a_strategies(sc, cfg=cfg)}
    s0, s1, s2 = strats["S0_centrality"], strats["S1_separate_shapley"], strats["S2_integrated_shapley"]
    row = {"scenario_id": sc.scenario_id, "n_bus": len(sc.broken_buses), "n_link": len(sc.broken_links)}
    # specialized parallel crews (current Task A execution)
    row["spec_S1"] = run(s1.sequence, sc.scenario_id, "spec_S1", "specialized", True, sc.link_capacity_factors)
    row["spec_S2"] = run(s2.sequence, sc.scenario_id, "spec_S2", "specialized", True, sc.link_capacity_factors)
    # one multifunction crew: separate Shapley stacked power-first / road-first, vs joint Shapley interleaved
    row["one_S1_power_first"] = run(list(s1.power_sequence) + list(s1.road_sequence), sc.scenario_id, "one_S1_pf", "multifunction", True, sc.link_capacity_factors)
    row["one_S1_road_first"] = run(list(s1.road_sequence) + list(s1.power_sequence), sc.scenario_id, "one_S1_rf", "multifunction", True, sc.link_capacity_factors)
    row["one_S2_joint"] = run(list(s2.sequence), sc.scenario_id, "one_S2", "multifunction", True, sc.link_capacity_factors)
    row["one_S0_power_first"] = run(list(s0.power_sequence) + list(s0.road_sequence), sc.scenario_id, "one_S0_pf", "multifunction", True, sc.link_capacity_factors)
    row["one_S0_road_first"] = run(list(s0.road_sequence) + list(s0.power_sequence), sc.scenario_id, "one_S0_rf", "multifunction", True, sc.link_capacity_factors)
    rows.append(row)
    with open(os.path.join(OUT_DIR, "probe_rows.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

keys = [k for k in rows[0] if k not in ("scenario_id", "n_bus", "n_link")]
lines = ["# Single-crew probe (n=%d scenarios, seed0 20260701)\n" % len(rows), "| execution / ranking | mean triangle ± 95% CI |", "|---|---|"]
for k in keys:
    xs = [r[k] for r in rows]
    lines.append(f"| {k} | {st.mean(xs):.1f} ± {1.96*st.stdev(xs)/len(xs)**0.5:.1f} |")


def paired(a, b):
    d = [r[a] - r[b] for r in rows]
    return f"{st.mean(d):+.1f} ± {1.96*st.stdev(d)/len(d)**0.5:.1f} ({a} better in {sum(1 for x in d if x < 0)}/{len(d)})"


lines += ["", "paired differences:",
          "- spec: S2 − S1 = " + paired("spec_S2", "spec_S1"),
          "- one crew: S2 joint − S1 power-first = " + paired("one_S2_joint", "one_S1_power_first"),
          "- one crew: S2 joint − S1 road-first  = " + paired("one_S2_joint", "one_S1_road_first"),
          "- one crew: S2 joint − best of the two S1 stackings = " + f"{st.mean([r['one_S2_joint'] - min(r['one_S1_power_first'], r['one_S1_road_first']) for r in rows]):+.1f}",
          "- one crew: S2 joint − S0 road-first = " + paired("one_S2_joint", "one_S0_road_first")]
open(os.path.join(OUT_DIR, "probe_summary.md"), "w").write("\n".join(lines))
print("\n".join(lines))
