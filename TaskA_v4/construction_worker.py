from __future__ import annotations

"""Build scenario-level joint Shapley scores for the table-construction ensemble.

Run from a private worker directory containing the simulator template.  The
driver assigns disjoint scenarios to workers, while every worker generates the
same deterministic ensemble manifest.
"""

import argparse
import json
import os

from access_value import AccessAwareValue, AccessConfig
from disaster import generate_scenarios
from task_a_criticality import (
    FullFunctionalityValue,
    TaskACriticalityConfig,
    _damaged_link_factors,
    _normalize_link,
    centrality_scores,
    read_road_network,
    sampled_shapley,
)

CONSTRUCTION_SEED0 = 20263001


def asset_key(asset):
    if isinstance(asset, tuple):
        u, v = sorted((int(asset[0]), int(asset[1])))
        return f"road:{u}-{v}"
    return f"power:{int(asset)}"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--worker-id", type=int, required=True)
    ap.add_argument("--n-workers", type=int, required=True)
    ap.add_argument("--n-scenarios", type=int, default=400)
    ap.add_argument("--result-root", required=True)
    args = ap.parse_args()

    os.makedirs(args.result_root, exist_ok=True)
    scenarios = generate_scenarios(
        n_scenarios=args.n_scenarios,
        seed0=CONSTRUCTION_SEED0,
        bus_count_range=(8, 15),
        link_count_range=(3, 13),
        link_drop_range=(0.5, 1.0),
        out_json=os.path.join(args.result_root, "construction_disasters.json"),
    )
    cfg = TaskACriticalityConfig(shapley_samples=120)

    for index, scenario in enumerate(scenarios, start=1):
        if (index - 1) % args.n_workers != args.worker_id:
            continue
        out_path = os.path.join(args.result_root, f"scores_{scenario.scenario_id}.json")
        if os.path.exists(out_path):
            continue
        print(f"[construction w{args.worker_id:02d}] {scenario.scenario_id} ({index}/{len(scenarios)})", flush=True)

        power_assets = [int(bus) for bus in scenario.broken_buses]
        road_assets = [_normalize_link(link) for link in scenario.broken_links]
        all_assets = list(power_assets) + list(road_assets)
        link_factors = _damaged_link_factors(scenario)

        power_cent, road_cent = centrality_scores(read_road_network(cfg.road_net_path))
        cen_scores = {
            **{bus: float(power_cent.get(bus, 0.0)) for bus in power_assets},
            **{
                link: max(
                    float(road_cent.get(link, 0.0)),
                    float(road_cent.get((link[1], link[0]), 0.0)),
                )
                for link in road_assets
            },
        }
        joint = FullFunctionalityValue(
            power_assets=power_assets,
            road_assets=road_assets,
            link_factors=link_factors,
            cfg=cfg,
        )
        jsh_scores = sampled_shapley(
            all_assets,
            joint.integrated_value,
            samples=cfg.shapley_samples,
            seed=cfg.shapley_seed + scenario.seed + 303,
        )
        aware = AccessAwareValue(
            power_assets=power_assets,
            road_assets=road_assets,
            link_factors=link_factors,
            cfg=cfg,
            acc=AccessConfig(alpha=1.0),
        )
        ijsh_scores = sampled_shapley(
            all_assets,
            aware.joint_value,
            samples=cfg.shapley_samples,
            seed=cfg.shapley_seed + scenario.seed + 303,
        )
        selected = {"CEN": cen_scores, "JSH": jsh_scores, "IJSH": ijsh_scores}
        record = {
            "scenario_id": scenario.scenario_id,
            "scenario_seed": scenario.seed,
            "broken_buses": list(scenario.broken_buses),
            "broken_roads": [list(link) for link in scenario.broken_links],
            "scores": {
                name: {asset_key(asset): float(value) for asset, value in scores.items()}
                for name, scores in selected.items()
            },
        }
        tmp_path = out_path + ".tmp"
        with open(tmp_path, "w", encoding="utf-8") as handle:
            json.dump(record, handle, indent=1)
        os.replace(tmp_path, out_path)

    print(f"[construction w{args.worker_id:02d}] DONE", flush=True)


if __name__ == "__main__":
    main()
