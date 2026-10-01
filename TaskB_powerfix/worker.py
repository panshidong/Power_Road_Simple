from __future__ import annotations

"""
One worker of the parallel Task B (equity / trade-off) rerun.

Run with cwd = its own private copy of the code+data template, because the
simulator uses fixed relative scratch paths (work/SiouxFalls_net*.txt, s.txt,
bus_location.json, taskB_*.json). Each worker processes the scenarios whose
1-based index satisfies (idx - 1) % n_workers == worker_id, using the exact
fast100 configuration that produced the equity paper's `newruns` batch
(runner_fast100.build_fast100_config: 7 experiments + baseline, SA max_iter=80,
T0=1.4, alpha=0.97, seed0=20260325). Resume-safe: a scenario whose
rows_<scenario_id>.csv already exists is skipped.
"""

import argparse
import csv
import os
from dataclasses import replace

from batch_tradeoff_runner import _scenario_row
from disaster import generate_scenarios, scenario_sequence
from runner_fast100 import build_fast100_config
from tradeoff_runner import run_tradeoff_study


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--worker-id", type=int, required=True)
    ap.add_argument("--n-workers", type=int, required=True)
    ap.add_argument("--result-root", required=True, help="absolute path; per-worker")
    args = ap.parse_args()

    cfg = build_fast100_config()
    os.makedirs(args.result_root, exist_ok=True)
    scenario_run_root = os.path.join(args.result_root, "scenario_runs")
    os.makedirs(scenario_run_root, exist_ok=True)

    scenarios = generate_scenarios(
        n_scenarios=cfg.n_scenarios,
        seed0=cfg.seed0,
        bus_count_range=cfg.bus_count_range,
        link_count_range=cfg.link_count_range,
        link_drop_range=cfg.link_drop_range,
        bus_candidate_path=cfg.bus_candidate_path,
        base_net_path=cfg.base_net_path,
        out_json=os.path.join(args.result_root, "random_disasters.json"),
    )

    mine = [(idx, s) for idx, s in enumerate(scenarios, start=1) if (idx - 1) % args.n_workers == args.worker_id]
    print(f"[worker {args.worker_id}] {len(mine)} scenarios: {[s.scenario_id for _, s in mine]}", flush=True)

    for idx, scenario in mine:
        rows_csv = os.path.join(args.result_root, f"rows_{scenario.scenario_id}.csv")
        if os.path.exists(rows_csv):
            print(f"[worker {args.worker_id}] skip {scenario.scenario_id} (done)", flush=True)
            continue
        print(f"[worker {args.worker_id}] running {scenario.scenario_id} (global {idx}/{len(scenarios)}) ...", flush=True)
        scenario_cfg = replace(
            cfg.study_cfg,
            result_root=scenario_run_root,
            base_sequence=scenario_sequence(scenario),
            broken_link_factors=dict(scenario.link_capacity_factors),
        )
        info = run_tradeoff_study(scenario_cfg)
        out_rows = []
        with open(info["csv_path"], newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                if row.get("row_kind") != "optimized":
                    continue
                out_rows.append(_scenario_row(row, scenario=scenario, study_result_dir=info["result_dir"]))
        tmp = rows_csv + ".tmp"
        with open(tmp, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
            w.writeheader()
            w.writerows(out_rows)
        os.replace(tmp, rows_csv)
        print(f"[worker {args.worker_id}] finished {scenario.scenario_id}: {len(out_rows)} rows", flush=True)

    print(f"[worker {args.worker_id}] ALL DONE", flush=True)


if __name__ == "__main__":
    main()
