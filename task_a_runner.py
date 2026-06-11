from __future__ import annotations

import csv
import math
import os
from dataclasses import dataclass, field
from datetime import datetime
from typing import Any, Dict, List, Tuple

from disaster import DisasterScenario, generate_scenarios, scenario_sequence
from resilience_measurement import optimize_sequence_sa, run_model_multi
from simulated_annealing import SAConfig
from task_a_criticality import (
    TaskACriticalityConfig,
    TaskAStrategy,
    build_task_a_strategies,
    strategy_ranking_rows,
)


def _ts() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def _sequence_key(seq: List[Any]) -> str:
    return repr(list(seq))


def _safe_float(value: Any) -> float:
    try:
        return float(value)
    except Exception:
        return float("nan")


def _mean(values: List[float]) -> float:
    xs = [float(v) for v in values if not math.isnan(float(v))]
    return sum(xs) / len(xs) if xs else float("nan")


def _std(values: List[float]) -> float:
    xs = [float(v) for v in values if not math.isnan(float(v))]
    if len(xs) <= 1:
        return 0.0
    mu = _mean(xs)
    return math.sqrt(sum((x - mu) ** 2 for x in xs) / (len(xs) - 1))


def _ci95(values: List[float]) -> float:
    xs = [float(v) for v in values if not math.isnan(float(v))]
    if len(xs) <= 1:
        return 0.0
    return 1.96 * _std(xs) / math.sqrt(len(xs))


def _metric_or_nan(run: Dict[str, Any], key: str) -> float:
    return float(run["metric_catalog"].get(key, float("nan")))


def _tapb_snapshot_count(run: Dict[str, Any]) -> int:
    paths = []
    for event in run.get("event_log", []):
        state = event.get("state", {}) if isinstance(event, dict) else {}
        path = state.get("s_txt_path")
        if path:
            paths.append(path)
    return len(set(paths))


@dataclass
class TaskAConfig:
    n_scenarios: int = 10
    seed0: int = 20260402
    result_root: str = "results"
    bus_count_range: Tuple[int, int] = (8, 15)
    link_count_range: Tuple[int, int] = (3, 13)
    link_drop_range: Tuple[float, float] = (0.5, 1.0)
    bus_candidate_path: str = "new_bus_to_link.json"
    base_net_path: str = "tap-b/net/SiouxFalls_net.txt"
    include_heuristic: bool = True
    save_strategy_artifacts: bool = False
    save_heuristic_artifacts: bool = False
    strict: bool = True
    crew_mode: str = "specialized"
    power_crews: int = 1
    road_crews: int = 1
    criticality: TaskACriticalityConfig = field(default_factory=TaskACriticalityConfig)
    heuristic_sa: SAConfig = field(
        default_factory=lambda: SAConfig(seed=20260402, max_iter=60, T0=2.0, alpha=0.98, neighbor="swap")
    )


def _row_from_run(
    *,
    scenario: DisasterScenario,
    strategy_id: str,
    strategy_label: str,
    score_method: str,
    row_kind: str,
    run: Dict[str, Any],
    preserve_sequence_order: bool,
    run_dir: str,
) -> Dict[str, Any]:
    return {
        "scenario_id": scenario.scenario_id,
        "scenario_seed": scenario.seed,
        "broken_bus_count": len(scenario.broken_buses),
        "broken_link_count": len(scenario.broken_links),
        "strategy_id": strategy_id,
        "strategy_label": strategy_label,
        "score_method": score_method,
        "row_kind": row_kind,
        "triangle_area": float(run["triangle_area"]),
        "objective_value": float(run["objective_value"]),
        "gini_restore": _metric_or_nan(run, "equity:gini_restore"),
        "var_restore": _metric_or_nan(run, "equity:var_restore"),
        "p90_restore": _metric_or_nan(run, "equity:p90_restore"),
        "min_time_avg_cri": _metric_or_nan(run, "equity:min_time_avg_cri"),
        "p90_access_restore": _metric_or_nan(run, "critical_access:p90_access_restore"),
        "share_access_initial": _metric_or_nan(run, "critical_access:share_access_initial"),
        "share_access_final": _metric_or_nan(run, "critical_access:share_access_final"),
        "time_avg_share_access": _metric_or_nan(run, "critical_access:time_avg_share_access"),
        "crew_mode": run.get("crew_mode", ""),
        "power_crews": run.get("power_crews", ""),
        "road_crews": run.get("road_crews", ""),
        "multifunction_crews": run.get("active_multifunction_crews", run.get("multifunction_crews", "")),
        "crew_execution_rule": "specialized crews each take the highest-ranked remaining asset they can repair",
        "preserve_sequence_order": int(bool(preserve_sequence_order)),
        "event_count": len(run.get("event_log", [])),
        "tapb_state_snapshots": _tapb_snapshot_count(run),
        "sequence": _sequence_key(run["sequence"]),
        "power_sequence": _sequence_key(run.get("power_sequence", [])),
        "road_sequence": _sequence_key(run.get("road_sequence", [])),
        "run_dir": run_dir,
    }


def _run_strategy(
    *,
    cfg: TaskAConfig,
    scenario: DisasterScenario,
    strategy: TaskAStrategy,
    scenario_root: str,
) -> Dict[str, Any]:
    run_dir = os.path.join(scenario_root, strategy.strategy_id)
    run = run_model_multi(
        strategy.sequence,
        result_root=scenario_root,
        run_dir=run_dir,
        message=f"Task A {strategy.strategy_label}",
        Scenario=f"{scenario.scenario_id}_{strategy.strategy_id}",
        strict=cfg.strict,
        save_artifacts=cfg.save_strategy_artifacts,
        crew_mode=cfg.crew_mode,
        power_crews=cfg.power_crews,
        road_crews=cfg.road_crews,
        preserve_sequence_order=strategy.preserve_sequence_order,
        broken_link_factors=dict(scenario.link_capacity_factors),
        objective="triangle",
    )
    return _row_from_run(
        scenario=scenario,
        strategy_id=strategy.strategy_id,
        strategy_label=strategy.strategy_label,
        score_method=strategy.score_method,
        row_kind="strategy",
        run=run,
        preserve_sequence_order=strategy.preserve_sequence_order,
        run_dir=run_dir,
    )


def _run_heuristic_benchmark(
    *,
    cfg: TaskAConfig,
    scenario: DisasterScenario,
    strategies: List[TaskAStrategy],
    scenario_root: str,
) -> Dict[str, Any]:
    s2 = next((strategy for strategy in strategies if strategy.strategy_id == "S2_integrated_shapley"), None)
    base_sequence = list(s2.sequence if s2 is not None else scenario_sequence(scenario))
    res = optimize_sequence_sa(
        base_sequence=base_sequence,
        result_root=scenario_root,
        message="Task A heuristic benchmark with specialized crews",
        Scenario=f"{scenario.scenario_id}_heuristic_triangle_sa",
        objective="triangle",
        sa=cfg.heuristic_sa,
        strict=cfg.strict,
        save_baseline=False,
        save_best_artifacts=cfg.save_heuristic_artifacts,
        save_best_debug=False,
        crew_mode=cfg.crew_mode,
        power_crews=cfg.power_crews,
        road_crews=cfg.road_crews,
        preserve_sequence_order=True,
        broken_link_factors=dict(scenario.link_capacity_factors),
    )
    best_run = res["best_run"]
    return _row_from_run(
        scenario=scenario,
        strategy_id="heuristic_triangle_sa",
        strategy_label="Heuristic benchmark: SA minimizes triangle",
        score_method="simulated annealing over repair priority list with specialized crews",
        row_kind="heuristic",
        run=best_run,
        preserve_sequence_order=True,
        run_dir=res["session_dir"],
    )


def _aggregate_rows(rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    metric_names = [
        "triangle_area",
        "gini_restore",
        "var_restore",
        "p90_restore",
        "min_time_avg_cri",
        "p90_access_restore",
        "share_access_initial",
        "share_access_final",
        "time_avg_share_access",
    ]
    groups: Dict[str, List[Dict[str, Any]]] = {}
    for row in rows:
        groups.setdefault(row["strategy_id"], []).append(row)

    out: List[Dict[str, Any]] = []
    for strategy_id, items in sorted(groups.items()):
        sample = items[0]
        agg: Dict[str, Any] = {
            "strategy_id": strategy_id,
            "strategy_label": sample["strategy_label"],
            "score_method": sample["score_method"],
            "row_kind": sample["row_kind"],
            "n_scenarios": len(items),
        }
        for metric in metric_names:
            values = [_safe_float(row[metric]) for row in items]
            agg[f"{metric}_mean"] = _mean(values)
            agg[f"{metric}_std"] = _std(values)
            agg[f"{metric}_ci95"] = _ci95(values)
        out.append(agg)
    out.sort(key=lambda row: float(row["triangle_area_mean"]))
    return out


def _write_csv(path: str, rows: List[Dict[str, Any]]) -> None:
    if not rows:
        return
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def _write_summary(path: str, cfg: TaskAConfig, aggregate_rows: List[Dict[str, Any]]) -> None:
    with open(path, "w", encoding="utf-8") as f:
        f.write("# Task A Criticality Strategy Summary\n\n")
        f.write(f"- Scenario count: `{cfg.n_scenarios}`\n")
        f.write(f"- Scenario seed start: `{cfg.seed0}`\n")
        f.write(f"- Power failures per scenario: `{cfg.bus_count_range}`\n")
        f.write(f"- Road failures per scenario: `{cfg.link_count_range}`\n")
        f.write(f"- Road capacity-drop range: `{cfg.link_drop_range}`\n")
        f.write(f"- Crew model: `{cfg.crew_mode}` with power crews `{cfg.power_crews}` and road crews `{cfg.road_crews}`\n")
        f.write("- Execution rule: strategy rankings are priority lists; specialized power and road crews each take the highest-ranked remaining asset they are eligible to repair.\n")
        f.write("- S0 centrality uses different network-specific methods: power radial downstream service impact and road weighted edge betweenness.\n")
        f.write("- S1 and S2 Shapley values use coalition marginal gains in simulator functionality, not centrality proxy scores.\n")
        f.write("- Roadway damage is interpreted as physical two-way roadway blockage; a damaged link affects both travel directions.\n")
        f.write("- S1 uses pure separate-network Shapley games: power Shapley excludes road failures, and road Shapley excludes power failures.\n")
        f.write("- S2 uses one integrated Shapley ranking, but repairs are still executed by separate specialized crews.\n")
        f.write("- The simulator re-runs TAP-B at each changed recovery state; `tapb_state_snapshots` in the CSV records the unique background network states evaluated per run.\n")
        f.write(f"- Road TSTT baseline follows the existing simulator baseline: `{cfg.criticality.baseline_tstt}`\n")
        f.write(f"- Shapley samples: `{cfg.criticality.shapley_samples}`\n")
        f.write(f"- Full-functionality Shapley: `{cfg.criticality.use_full_functionality_shapley}`\n")
        f.write(f"- Integrated Shapley power weight: `{cfg.criticality.integrated_power_weight}`\n")
        f.write(f"- Heuristic benchmark: `{cfg.include_heuristic}`\n\n")
        f.write("## Aggregate Results\n\n")
        f.write("| Strategy | Kind | Triangle mean +/- CI | Gini mean +/- CI | P90 access mean +/- CI |\n")
        f.write("|---|---|---:|---:|---:|\n")
        for row in aggregate_rows:
            f.write(
                f"| {row['strategy_id']} | {row['row_kind']} | "
                f"{float(row['triangle_area_mean']):.3f} +/- {float(row['triangle_area_ci95']):.3f} | "
                f"{float(row['gini_restore_mean']):.3f} +/- {float(row['gini_restore_ci95']):.3f} | "
                f"{float(row['p90_access_restore_mean']):.3f} +/- {float(row['p90_access_restore_ci95']):.3f} |\n"
            )


def run_task_a_study(cfg: TaskAConfig) -> Dict[str, str]:
    result_dir = os.path.join(cfg.result_root, f"task_a_criticality_{_ts()}")
    os.makedirs(result_dir, exist_ok=True)
    scenario_root = os.path.join(result_dir, "scenario_runs")
    os.makedirs(scenario_root, exist_ok=True)

    scenarios = generate_scenarios(
        n_scenarios=cfg.n_scenarios,
        seed0=cfg.seed0,
        bus_count_range=cfg.bus_count_range,
        link_count_range=cfg.link_count_range,
        link_drop_range=cfg.link_drop_range,
        bus_candidate_path=cfg.bus_candidate_path,
        base_net_path=cfg.base_net_path,
        out_json=os.path.join(result_dir, "random_disasters.json"),
    )

    manifest_rows: List[Dict[str, Any]] = []
    strategy_rows: List[Dict[str, Any]] = []
    ranking_rows: List[Dict[str, Any]] = []

    for idx, scenario in enumerate(scenarios, start=1):
        print(f"[task-a] running {scenario.scenario_id} ({idx}/{len(scenarios)}) ...")
        current_root = os.path.join(scenario_root, scenario.scenario_id)
        os.makedirs(current_root, exist_ok=True)
        strategies = build_task_a_strategies(scenario, cfg=cfg.criticality)

        manifest_rows.append(
            {
                "scenario_id": scenario.scenario_id,
                "scenario_seed": scenario.seed,
                "broken_bus_count": len(scenario.broken_buses),
                "broken_link_count": len(scenario.broken_links),
                "broken_buses": _sequence_key(list(scenario.broken_buses)),
                "broken_links": _sequence_key(list(scenario.broken_links)),
                "scenario_result_dir": current_root,
            }
        )

        for strategy in strategies:
            ranking_rows.extend(strategy_ranking_rows(scenario, strategy))
            strategy_rows.append(_run_strategy(cfg=cfg, scenario=scenario, strategy=strategy, scenario_root=current_root))

        if cfg.include_heuristic:
            strategy_rows.append(
                _run_heuristic_benchmark(
                    cfg=cfg,
                    scenario=scenario,
                    strategies=strategies,
                    scenario_root=current_root,
                )
            )

    aggregate_rows = _aggregate_rows(strategy_rows)

    manifest_csv = os.path.join(result_dir, "scenario_manifest.csv")
    strategy_csv = os.path.join(result_dir, "scenario_strategy_rows.csv")
    ranking_csv = os.path.join(result_dir, "criticality_rankings.csv")
    aggregate_csv = os.path.join(result_dir, "aggregate_task_a_summary.csv")
    summary_md = os.path.join(result_dir, "task_a_summary.md")

    _write_csv(manifest_csv, manifest_rows)
    _write_csv(strategy_csv, strategy_rows)
    _write_csv(ranking_csv, ranking_rows)
    _write_csv(aggregate_csv, aggregate_rows)
    _write_summary(summary_md, cfg, aggregate_rows)

    return {
        "result_dir": result_dir,
        "scenario_manifest_csv": manifest_csv,
        "strategy_rows_csv": strategy_csv,
        "ranking_csv": ranking_csv,
        "aggregate_csv": aggregate_csv,
        "summary_md": summary_md,
    }


def main() -> None:
    info = run_task_a_study(TaskAConfig())
    print("Task A study complete.")
    print("Result dir:", info["result_dir"])
    print("Scenario manifest:", info["scenario_manifest_csv"])
    print("Strategy rows:", info["strategy_rows_csv"])
    print("Criticality rankings:", info["ranking_csv"])
    print("Aggregate CSV:", info["aggregate_csv"])
    print("Summary:", info["summary_md"])


if __name__ == "__main__":
    main()
