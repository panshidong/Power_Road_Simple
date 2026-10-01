from __future__ import annotations

"""
Paired-comparison analysis of the Task C main study (100 scenarios, S0-S3, v2
redesign). Read-only against results/task_c_main_100/scenario_strategy_rows.csv.
"""

import csv
import math
from collections import defaultdict
from typing import Dict, List

METRICS = ["triangle_area", "weighted_triangle_area"]
STRATEGIES = ["S0_centrality", "S1_separate_shapley", "S2_integrated_shapley", "S3_od_representation"]


def load_rows(path: str) -> List[Dict[str, str]]:
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def mean(xs: List[float]) -> float:
    return sum(xs) / len(xs) if xs else float("nan")


def std(xs: List[float]) -> float:
    if len(xs) <= 1:
        return 0.0
    mu = mean(xs)
    return math.sqrt(sum((x - mu) ** 2 for x in xs) / (len(xs) - 1))


def ci95(xs: List[float]) -> float:
    if len(xs) <= 1:
        return 0.0
    return 1.96 * std(xs) / math.sqrt(len(xs))


def main() -> None:
    rows = load_rows("results/task_c_main_100/scenario_strategy_rows.csv")
    by_scenario_strategy: Dict[str, Dict[str, Dict[str, float]]] = defaultdict(dict)
    for row in rows:
        sid = row["scenario_id"]
        strat = row["strategy_id"]
        by_scenario_strategy[sid][strat] = {m: float(row[m]) for m in METRICS}

    scenario_ids = sorted(by_scenario_strategy.keys())
    complete_scenarios = [
        sid for sid in scenario_ids if all(s in by_scenario_strategy[sid] for s in STRATEGIES)
    ]
    print(f"n scenarios with all 4 strategies: {len(complete_scenarios)} / {len(scenario_ids)}")

    print("\n=== Marginal means ===")
    for strat in STRATEGIES:
        for metric in METRICS:
            xs = [by_scenario_strategy[sid][strat][metric] for sid in complete_scenarios]
            print(f"{strat:24s} {metric:24s} mean={mean(xs):12.3f} ci95={ci95(xs):10.3f}")
        print()

    print("\n=== Paired differences: S3 - baseline (per-scenario) ===")
    for baseline in ["S0_centrality", "S1_separate_shapley", "S2_integrated_shapley"]:
        print(f"\n-- S3_od_representation vs {baseline} --")
        for metric in METRICS:
            diffs = [
                by_scenario_strategy[sid]["S3_od_representation"][metric] - by_scenario_strategy[sid][baseline][metric]
                for sid in complete_scenarios
            ]
            n_better = sum(1 for d in diffs if d < 0)
            n_worse = sum(1 for d in diffs if d > 0)
            n_tied = sum(1 for d in diffs if d == 0)
            m, c = mean(diffs), ci95(diffs)
            sig = "SIGNIFICANT" if abs(m) > c else "not significant"
            print(
                f"{metric:24s} mean_diff={m:12.3f} ci95={c:10.3f} [{sig}] "
                f"(S3 better {n_better}, worse {n_worse}, tied {n_tied} of {len(diffs)})"
            )

    print("\n=== Relative change in mean triangle metrics (S3 vs baseline) ===")
    for metric in METRICS:
        s3_mean = mean([by_scenario_strategy[sid]["S3_od_representation"][metric] for sid in complete_scenarios])
        for baseline in ["S0_centrality", "S1_separate_shapley", "S2_integrated_shapley"]:
            base_mean = mean([by_scenario_strategy[sid][baseline][metric] for sid in complete_scenarios])
            pct = (s3_mean - base_mean) / base_mean * 100.0
            print(f"{metric:24s} vs {baseline:24s} base={base_mean:10.3f} S3={s3_mean:10.3f} pct={pct:8.2f}%")

    print("\n=== How often is S3's road sequence identical to S0's? (degeneracy check) ===")
    identical = 0
    for row in rows:
        pass
    road_seq_by = {}
    for row in rows:
        if row["strategy_id"] in ("S0_centrality", "S3_od_representation"):
            road_seq_by.setdefault(row["scenario_id"], {})[row["strategy_id"]] = row["road_sequence"]
    n_identical = sum(
        1
        for sid, d in road_seq_by.items()
        if "S0_centrality" in d and "S3_od_representation" in d and d["S0_centrality"] == d["S3_od_representation"]
    )
    print(f"S3 road sequence identical to S0 in {n_identical} / {len(road_seq_by)} scenarios")


if __name__ == "__main__":
    main()
