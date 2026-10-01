from __future__ import annotations

"""Task C (O-D road criticality) evaluated against the pre-event Task A tables.

Run with cwd = a private worker directory that is a copy of TaskA_v4/template
(physical-road disruption sampling, TAP-B gap-relaxation fallback, corrected
power rule) plus task_c_od.py / task_c_weighted_triangle.py from /home/workenv/OD.

Strategies per scenario (all are pre-event tables; specialized crews, so only the
within-trade order matters):
  CEN, JSH, IJSH          : Task A v4 tables applied exactly as evaluation_worker.py
  OD_CEN, OD_JSH, OD_IJSH : road block ordered by the static critical-O-D score
                            (task_c_od.build_static_road_importance), power block
                            ordered by the named table's power scores.
Each OD_X vs X pair therefore isolates the road-ranking change with the power
order held fixed.
"""

import argparse
import csv
import json
import os
import re

from disaster import generate_scenarios
from resilience_measurement import run_model_multi
from task_a_criticality import TaskACriticalityConfig, _normalize_link, _sorted_by_score
from task_c_od import TaskCODConfig, build_static_road_importance
from task_c_weighted_triangle import CriticalWeightConfig, compute_weighted_triangle

EVALUATION_SEED0 = 20264001  # identical to TaskA_v4/evaluation_worker.py
# Shifted-damage ensembles: identical definitions and seeds to TaskA_v4/robustness_worker.py
VARIANTS = {
    "small": ((3, 7), (1, 5), (0.5, 1.0), False, 20265001),
    "large": ((15, 20), (12, 18), (0.5, 1.0), False, 20266001),
    "light": ((8, 15), (3, 13), (0.2, 0.6), False, 20267001),
    "clustered": ((8, 15), (3, 13), (0.5, 1.0), True, 20268001),
}


def _road_adjacency(roads):
    adjacency = {}
    for u, v in roads:
        adjacency.setdefault(int(u), set()).add(int(v))
        adjacency.setdefault(int(v), set()).add(int(u))
    return adjacency


def _bus_anchor_nodes(path="new_bus_to_link.json"):
    raw = json.load(open(path, encoding="utf-8"))
    return {int(bus): {int(link[0]), int(link[1])} for bus, link in raw.items()}


def _clustered_sample(rng, candidates, anchors, n_power, n_road):
    from collections import deque
    adjacency = _road_adjacency(candidates["links"])
    center = rng.choice(sorted(adjacency))
    ball = {center}
    frontier = deque([center])
    while frontier:
        road_pool = [road for road in candidates["links"] if road[0] in ball and road[1] in ball]
        power_pool = [bus for bus in candidates["buses"] if anchors.get(bus, set()) & ball]
        if len(road_pool) >= n_road and len(power_pool) >= n_power:
            return rng.sample(power_pool, n_power), rng.sample(road_pool, n_road)
        node = frontier.popleft()
        for neighbor in sorted(adjacency.get(node, ())):
            if neighbor not in ball:
                ball.add(neighbor)
                frontier.append(neighbor)
    raise RuntimeError("Could not construct clustered disruption with requested counts")


def generate_variant(variant, count):
    import random
    from disaster import DisasterScenario, load_candidates
    bus_range, road_range, drop_range, clustered, seed0 = VARIANTS[variant]
    candidates = load_candidates()
    anchors = _bus_anchor_nodes()
    out = []
    for index in range(count):
        seed = seed0 + index
        rng = random.Random(seed)
        n_power = rng.randint(*bus_range)
        n_road = rng.randint(*road_range)
        if clustered:
            buses, roads = _clustered_sample(rng, candidates, anchors, n_power, n_road)
        else:
            buses = rng.sample(candidates["buses"], n_power)
            roads = rng.sample(candidates["links"], n_road)
        factors = {road: max(0.0, 1.0 - rng.uniform(*drop_range)) for road in roads}
        out.append(DisasterScenario(scenario_id=f"{variant}_{index + 1:03d}", seed=seed, broken_buses=buses,
                                    broken_links=roads, link_capacity_factors=factors))
    return out
POWER_RULES = ("CEN", "JSH", "IJSH")
DEFAULT_WEIGHTS = (1.0, 0.6, 0.35)


def asset_key(asset):
    if isinstance(asset, tuple):
        u, v = sorted((int(asset[0]), int(asset[1])))
        return f"road:{u}-{v}"
    return f"power:{int(asset)}"


def metric(run, key):
    return float(run["metric_catalog"].get(key, float("nan")))


def od_lookup(table, link):
    return float(table.get(link, table.get((link[1], link[0]), 0.0)))


def fallback_lines():
    p = os.path.join(os.getcwd(), "tapb_fallback.log")
    if not os.path.exists(p):
        return 0
    with open(p, encoding="utf-8") as f:
        return sum(1 for _ in f)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--worker-id", type=int, required=True)
    ap.add_argument("--n-workers", type=int, required=True)
    ap.add_argument("--n-scenarios", type=int, default=300)
    ap.add_argument("--result-root", required=True)
    ap.add_argument("--tables", required=True)
    ap.add_argument("--checkpoint", type=int, default=400)
    ap.add_argument("--k", type=int, default=3)
    ap.add_argument("--strategies", default="CEN,JSH,IJSH,OD_CEN,OD_JSH,OD_IJSH")
    ap.add_argument("--critical-buses", default="8,17")
    ap.add_argument("--shelter-nodes", default="24")
    ap.add_argument("--weight-buses", default=None, help="buses upweighted in the supplemental metric (default: the critical buses)")
    ap.add_argument("--pair-mode", default="depot_to_critical", choices=("depot_to_critical", "all_pairs"))
    ap.add_argument("--variant", default=None, help="shifted-damage ensemble name (small/large/light/clustered)")
    ap.add_argument("--reference-manifest", default="/home/workenv/TaskA_v4/results/evaluation/worker_00/evaluation_disasters.json")
    args = ap.parse_args()

    os.makedirs(args.result_root, exist_ok=True)
    tables = json.load(open(args.tables, encoding="utf-8"))["checkpoints"][str(args.checkpoint)]["tables"]
    wanted = [s.strip() for s in args.strategies.split(",") if s.strip()]

    if args.variant:
        scenarios = generate_variant(args.variant, args.n_scenarios)
    else:
        scenarios = generate_scenarios(
            n_scenarios=args.n_scenarios,
            seed0=EVALUATION_SEED0,
            bus_count_range=(8, 15),
            link_count_range=(3, 13),
            link_drop_range=(0.5, 1.0),
            out_json=os.path.join(args.result_root, "evaluation_disasters.json"),
        )
    # Cross-check against the manifest the Task A evaluation actually used.
    if not args.variant and args.reference_manifest and os.path.exists(args.reference_manifest):
        ref = {s["scenario_id"]: s for s in json.load(open(args.reference_manifest, encoding="utf-8"))["scenarios"]}
        for sc in scenarios:
            r = ref.get(sc.scenario_id)
            if r is None:
                continue
            assert sorted(map(int, sc.broken_buses)) == sorted(map(int, r["broken_buses"])), sc.scenario_id
            ref_links = {(int(l["u"]), int(l["v"])): float(l["remaining_capacity_factor"]) for l in r["broken_links"]}
            mine = {(int(u), int(v)): float(f) for (u, v), f in sc.link_capacity_factors.items()}
            assert mine == ref_links, (sc.scenario_id, mine, ref_links)

    cb = tuple(int(x) for x in args.critical_buses.split(",") if x.strip()); sh = tuple(int(x) for x in args.shelter_nodes.split(","))
    od_cfg = TaskCODConfig(k_paths=args.k, path_rank_weights=tuple(DEFAULT_WEIGHTS[: args.k]) or (1.0,), pair_mode=args.pair_mode, critical_buses=cb, shelter_nodes=sh)
    static_od = build_static_road_importance(od_cfg, TaskACriticalityConfig())
    wb = tuple(int(x) for x in args.weight_buses.split(",") if x.strip()) if args.weight_buses else cb
    weight_cfg = CriticalWeightConfig(critical_buses=wb, shelter_node=sh[0])

    with open(os.path.join(args.result_root, f"static_od_table_k{args.k}_{args.pair_mode}.json"), "w", encoding="utf-8") as f:
        json.dump({f"{u}-{v}": val for (u, v), val in sorted(static_od.items())}, f, indent=1)

    for index, scenario in enumerate(scenarios, start=1):
        if (index - 1) % args.n_workers != args.worker_id:
            continue
        out_path = os.path.join(args.result_root, f"rows_{scenario.scenario_id}.csv")
        if os.path.exists(out_path):
            continue
        print(f"[w{args.worker_id:02d}] {scenario.scenario_id} ({index}/{len(scenarios)}) k={args.k}", flush=True)

        power_assets = [int(b) for b in scenario.broken_buses]
        road_assets = [_normalize_link(l) for l in scenario.broken_links]
        assets = power_assets + road_assets
        rows = []
        try:
          for strategy in wanted:
              if strategy in POWER_RULES:
                  scores = {a: float(tables[strategy][asset_key(a)]) for a in assets}
                  sequence = _sorted_by_score(assets, scores)
              elif strategy.startswith("MIX_"):
                  # MIX_<roadrule>road_<powerrule>power: decomposition run (road order from one
                  # Task A table, power order from another), used to attribute IJSH vs JSH.
                  m = re.fullmatch(r"MIX_(\w+)road_(\w+)power", strategy)
                  r_rule, p_rule = m.group(1), m.group(2)
                  p_scores = {a: float(tables[p_rule][asset_key(a)]) for a in power_assets}
                  r_scores = {a: float(tables[r_rule][asset_key(a)]) for a in road_assets}
                  sequence = _sorted_by_score(power_assets, p_scores) + _sorted_by_score(road_assets, r_scores)
                  scores = {**p_scores, **r_scores}
              else:
                  rule = strategy.split("_", 1)[1]
                  p_scores = {a: float(tables[rule][asset_key(a)]) for a in power_assets}
                  r_scores = {a: od_lookup(static_od, a) for a in road_assets}
                  sequence = _sorted_by_score(power_assets, p_scores) + _sorted_by_score(road_assets, r_scores)
                  scores = {**p_scores, **r_scores}
              fb0 = fallback_lines()
              run = run_model_multi(
                  sequence,
                  result_root=args.result_root,
                  run_dir=os.path.join(args.result_root, "runs", f"{scenario.scenario_id}_{strategy}"),
                  message=f"Task C offline {strategy}",
                  Scenario=f"{scenario.scenario_id}_{strategy}",
                  strict=True,
                  save_artifacts=False,
                  crew_mode="specialized",
                  power_crews=1,
                  road_crews=1,
                  preserve_sequence_order=True,
                  broken_link_factors=dict(scenario.link_capacity_factors),
                  objective="triangle",
              )
              weighted = compute_weighted_triangle(run, cfg=weight_cfg)
              rows.append(
                  {
                      "scenario_id": scenario.scenario_id,
                      "scenario_seed": scenario.seed,
                      "n_power": len(power_assets),
                      "n_road": len(road_assets),
                      "k": args.k,
                      "strategy_id": strategy,
                      "triangle_area": float(run["triangle_area"]),
                      "weighted_triangle_area": float(weighted["weighted_triangle_area"]),
                      "shelter_access_initial": float(weighted["shelter_access_initial"]),
                      "shelter_access_final": float(weighted["shelter_access_final"]),
                      "gini_restore": metric(run, "equity:gini_restore"),
                      "min_time_avg_cri": metric(run, "equity:min_time_avg_cri"),
                      "p90_access_restore": metric(run, "critical_access:p90_access_restore"),
                      "n_road_nonzero_od": sum(1 for a in road_assets if od_lookup(static_od, a) > 0.0),
                      "tapb_fallbacks": fallback_lines() - fb0,
                      "power_sequence": repr([a for a in run["sequence"] if not isinstance(a, tuple)]),
                      "road_sequence": repr([a for a in run["sequence"] if isinstance(a, tuple)]),
                  }
              )
        except Exception as exc:  # heavily damaged variant that TAP-B cannot solve: record and skip, as robustness_worker.py does
            open(os.path.join(args.result_root, f"FAILED_{scenario.scenario_id}.txt"), "w").write(f"{type(exc).__name__}: {exc}\n")
            print(f"[w{args.worker_id:02d}] SKIP {scenario.scenario_id}: {type(exc).__name__}: {str(exc)[:200]}", flush=True)
            continue
        tmp = out_path + ".tmp"
        with open(tmp, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        os.replace(tmp, out_path)
    print(f"[w{args.worker_id:02d}] DONE", flush=True)


if __name__ == "__main__":
    main()
