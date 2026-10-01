from __future__ import annotations

"""Evaluate frozen v4 tables under shifted disruption distributions."""

import argparse
import csv
import json
import os
import random
from collections import deque

from disaster import DisasterScenario, load_candidates
from resilience_measurement import run_model_multi
from task_a_criticality import _sorted_by_score

VARIANTS = {
    "small": ((3, 7), (1, 5), (0.5, 1.0), False, 20265001),
    "large": ((15, 20), (12, 18), (0.5, 1.0), False, 20266001),
    "light": ((8, 15), (3, 13), (0.2, 0.6), False, 20267001),
    "clustered": ((8, 15), (3, 13), (0.5, 1.0), True, 20268001),
}


def asset_key(asset):
    if isinstance(asset, tuple):
        u, v = sorted((int(asset[0]), int(asset[1])))
        return f"road:{u}-{v}"
    return f"power:{int(asset)}"


def road_adjacency(roads):
    adjacency = {}
    for u, v in roads:
        adjacency.setdefault(int(u), set()).add(int(v))
        adjacency.setdefault(int(v), set()).add(int(u))
    return adjacency


def bus_anchor_nodes(path="new_bus_to_link.json"):
    raw = json.load(open(path, encoding="utf-8"))
    return {int(bus): {int(link[0]), int(link[1])} for bus, link in raw.items()}


def clustered_sample(rng, candidates, anchors, n_power, n_road):
    adjacency = road_adjacency(candidates["links"])
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


def generate_variant(variant: str, count: int):
    bus_range, road_range, drop_range, clustered, seed0 = VARIANTS[variant]
    candidates = load_candidates()
    anchors = bus_anchor_nodes()
    scenarios = []
    for index in range(count):
        seed = seed0 + index
        rng = random.Random(seed)
        n_power = rng.randint(*bus_range)
        n_road = rng.randint(*road_range)
        if clustered:
            buses, roads = clustered_sample(rng, candidates, anchors, n_power, n_road)
        else:
            buses = rng.sample(candidates["buses"], n_power)
            roads = rng.sample(candidates["links"], n_road)
        factors = {
            road: max(0.0, 1.0 - rng.uniform(*drop_range))
            for road in roads
        }
        scenarios.append(
            DisasterScenario(
                scenario_id=f"{variant}_{index + 1:03d}",
                seed=seed,
                broken_buses=buses,
                broken_links=roads,
                link_capacity_factors=factors,
            )
        )
    return scenarios


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--worker-id", type=int, required=True)
    ap.add_argument("--n-workers", type=int, required=True)
    ap.add_argument("--n-scenarios", type=int, default=100, help="scenarios per variant")
    ap.add_argument("--result-root", required=True)
    ap.add_argument("--tables", required=True)
    ap.add_argument("--checkpoint", type=int, default=400)
    args = ap.parse_args()

    os.makedirs(args.result_root, exist_ok=True)
    tables = json.load(open(args.tables, encoding="utf-8"))["checkpoints"][str(args.checkpoint)]["tables"]
    jobs = [(variant, scenario) for variant in VARIANTS for scenario in generate_variant(variant, args.n_scenarios)]

    for job_index, (variant, scenario) in enumerate(jobs):
        if job_index % args.n_workers != args.worker_id:
            continue
        out_path = os.path.join(args.result_root, f"rows_{scenario.scenario_id}.csv")
        failed_path = os.path.join(args.result_root, f"FAILED_{scenario.scenario_id}.txt")
        if os.path.exists(out_path):
            continue
        print(f"[robustness w{args.worker_id:02d}] {scenario.scenario_id}", flush=True)
        assets = list(map(int, scenario.broken_buses)) + [tuple(map(int, road)) for road in scenario.broken_links]
        rows = []
        try:
            for strategy in ("CEN", "JSH", "IJSH"):
                scores = {asset: float(tables[strategy][asset_key(asset)]) for asset in assets}
                sequence = _sorted_by_score(assets, scores)
                run = run_model_multi(
                    sequence,
                    result_root=args.result_root,
                    run_dir=os.path.join(args.result_root, "runs", f"{scenario.scenario_id}_{strategy}"),
                    message=f"Task A v4 robustness {strategy}",
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
                rows.append(
                    {
                        "variant": variant,
                        "scenario_id": scenario.scenario_id,
                        "scenario_seed": scenario.seed,
                        "n_power": len(scenario.broken_buses),
                        "n_road": len(scenario.broken_links),
                        "strategy_id": strategy,
                        "triangle_area": float(run["triangle_area"]),
                    }
                )
        except Exception as exc:
            with open(failed_path, "w", encoding="utf-8") as handle:
                handle.write(f"{type(exc).__name__}: {exc}\n")
            print(f"[robustness w{args.worker_id:02d}] FAILED {scenario.scenario_id}: {exc}", flush=True)
            continue

        tmp_path = out_path + ".tmp"
        with open(tmp_path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        os.replace(tmp_path, out_path)
        if os.path.exists(failed_path):
            os.remove(failed_path)

    print(f"[robustness w{args.worker_id:02d}] DONE", flush=True)


if __name__ == "__main__":
    main()
