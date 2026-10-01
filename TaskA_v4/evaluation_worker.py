from __future__ import annotations

"""Apply frozen offline tables to an independent disruption ensemble."""

import argparse
import csv
import json
import os

from disaster import generate_scenarios
from resilience_measurement import run_model_multi
from task_a_criticality import _sorted_by_score, _split_assets

EVALUATION_SEED0 = 20264001


def asset_key(asset):
    if isinstance(asset, tuple):
        u, v = sorted((int(asset[0]), int(asset[1])))
        return f"road:{u}-{v}"
    return f"power:{int(asset)}"


def metric(run, key):
    return float(run["metric_catalog"].get(key, float("nan")))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--worker-id", type=int, required=True)
    ap.add_argument("--n-workers", type=int, required=True)
    ap.add_argument("--n-scenarios", type=int, default=300)
    ap.add_argument("--result-root", required=True)
    ap.add_argument("--tables", required=True)
    ap.add_argument("--checkpoint", type=int, default=400)
    args = ap.parse_args()

    os.makedirs(args.result_root, exist_ok=True)
    table_file = json.load(open(args.tables, encoding="utf-8"))
    checkpoint = table_file["checkpoints"][str(args.checkpoint)]["tables"]
    scenarios = generate_scenarios(
        n_scenarios=args.n_scenarios,
        seed0=EVALUATION_SEED0,
        bus_count_range=(8, 15),
        link_count_range=(3, 13),
        link_drop_range=(0.5, 1.0),
        out_json=os.path.join(args.result_root, "evaluation_disasters.json"),
    )

    for index, scenario in enumerate(scenarios, start=1):
        if (index - 1) % args.n_workers != args.worker_id:
            continue
        out_path = os.path.join(args.result_root, f"rows_{scenario.scenario_id}.csv")
        if os.path.exists(out_path):
            continue
        print(f"[evaluation w{args.worker_id:02d}] {scenario.scenario_id} ({index}/{len(scenarios)})", flush=True)

        assets = list(map(int, scenario.broken_buses)) + [tuple(map(int, link)) for link in scenario.broken_links]
        rows = []
        for strategy in ("CEN", "JSH", "IJSH"):
            missing = [asset_key(asset) for asset in assets if asset_key(asset) not in checkpoint[strategy]]
            if missing:
                raise RuntimeError(f"{strategy} table misses {missing}")
            scores = {asset: float(checkpoint[strategy][asset_key(asset)]) for asset in assets}
            sequence = _sorted_by_score(assets, scores)
            run = run_model_multi(
                sequence,
                result_root=args.result_root,
                run_dir=os.path.join(args.result_root, "runs", f"{scenario.scenario_id}_{strategy}"),
                message=f"Task A v4 {strategy}",
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
            power_sequence, road_sequence = _split_assets(sequence)
            rows.append(
                {
                    "scenario_id": scenario.scenario_id,
                    "scenario_seed": scenario.seed,
                    "n_power": len(scenario.broken_buses),
                    "n_road": len(scenario.broken_links),
                    "strategy_id": strategy,
                    "triangle_area": float(run["triangle_area"]),
                    "gini_restore": metric(run, "equity:gini_restore"),
                    "min_time_avg_cri": metric(run, "equity:min_time_avg_cri"),
                    "p90_access_restore": metric(run, "critical_access:p90_access_restore"),
                    "power_sequence": repr(power_sequence),
                    "road_sequence": repr(road_sequence),
                }
            )

        tmp_path = out_path + ".tmp"
        with open(tmp_path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        os.replace(tmp_path, out_path)

    print(f"[evaluation w{args.worker_id:02d}] DONE", flush=True)


if __name__ == "__main__":
    main()
