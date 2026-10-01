from __future__ import annotations

import csv
import math
import os
from dataclasses import dataclass, field
from datetime import datetime
from typing import Any, Dict, List, Tuple

from disaster import DisasterScenario, generate_scenarios
from resilience_measurement import run_model_multi
from task_a_criticality import TaskACriticalityConfig, TaskAStrategy, build_task_a_strategies, strategy_ranking_rows
from task_c_od import Link, TaskCODConfig, build_static_road_importance, build_task_c_strategies
from task_c_weighted_triangle import CriticalWeightConfig, compute_weighted_triangle

"""
Task C experiment runner: evaluates the S3 (O-D representation) road-repair
strategy from task_c_od.py against Task A's S0/S1/S2 baselines under the same
disruption scenarios and the same resilience-triangle simulator, per proposal
Section 4.4 / RQ3. This module only reads task_a_criticality.py and
resilience_measurement.py; it does not modify them.

Every run's outcome is reported under two metrics:
  - triangle_area: the standard, unweighted resilience triangle already used
    elsewhere in this dissertation (primary efficiency metric).
  - weighted_triangle_area: a supplemental, post-hoc composite triangle that
    upweights a small set of designated critical power buses and adds a
    shelter-accessibility term (task_c_weighted_triangle.py). This does not
    change which sequence gets run; it only changes how the SAME run is
    scored, so both numbers come from the same simulation.
"""


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
class TaskCConfig:
    n_scenarios: int = 10
    seed0: int = 20260701
    result_root: str = "results"
    result_dir: str | None = None
    resume: bool = True
    bus_count_range: Tuple[int, int] = (8, 15)
    link_count_range: Tuple[int, int] = (3, 13)
    link_drop_range: Tuple[float, float] = (0.5, 1.0)
    bus_candidate_path: str = "new_bus_to_link.json"
    base_net_path: str = "tap-b/net/SiouxFalls_net.txt"
    include_baselines: bool = True
    save_strategy_artifacts: bool = False
    strict: bool = True
    crew_mode: str = "specialized"
    power_crews: int = 1
    road_crews: int = 1
    criticality: TaskACriticalityConfig = field(default_factory=TaskACriticalityConfig)
    od: TaskCODConfig = field(default_factory=TaskCODConfig)
    weight: CriticalWeightConfig = field(default_factory=CriticalWeightConfig)
    sensitivity_tag: str = "base"


def _row_from_run(
    *,
    scenario: DisasterScenario,
    strategy_id: str,
    strategy_label: str,
    score_method: str,
    sensitivity_tag: str,
    run: Dict[str, Any],
    weighted: Dict[str, float],
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
        "sensitivity_tag": sensitivity_tag,
        "triangle_area": float(run["triangle_area"]),
        "weighted_triangle_area": float(weighted["weighted_triangle_area"]),
        "shelter_access_initial": float(weighted["shelter_access_initial"]),
        "shelter_access_final": float(weighted["shelter_access_final"]),
        "weighted_power_initial": float(weighted["weighted_power_initial"]),
        "weighted_power_final": float(weighted["weighted_power_final"]),
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
    cfg: TaskCConfig,
    scenario: DisasterScenario,
    strategy: TaskAStrategy,
    scenario_root: str,
) -> Dict[str, Any]:
    run_dir = os.path.join(scenario_root, strategy.strategy_id)
    run = run_model_multi(
        strategy.sequence,
        result_root=scenario_root,
        run_dir=run_dir,
        message=f"Task C {strategy.strategy_label}",
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
    weighted = compute_weighted_triangle(run, cfg=cfg.weight)
    return _row_from_run(
        scenario=scenario,
        strategy_id=strategy.strategy_id,
        strategy_label=strategy.strategy_label,
        score_method=strategy.score_method,
        sensitivity_tag=cfg.sensitivity_tag,
        run=run,
        weighted=weighted,
        preserve_sequence_order=strategy.preserve_sequence_order,
        run_dir=run_dir,
    )


def _aggregate_rows(rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    metric_names = [
        "triangle_area",
        "weighted_triangle_area",
        "shelter_access_initial",
        "shelter_access_final",
        "weighted_power_initial",
        "weighted_power_final",
        "gini_restore",
        "var_restore",
        "p90_restore",
        "min_time_avg_cri",
        "p90_access_restore",
        "share_access_initial",
        "share_access_final",
        "time_avg_share_access",
    ]
    groups: Dict[Tuple[str, str], List[Dict[str, Any]]] = {}
    for row in rows:
        groups.setdefault((row["sensitivity_tag"], row["strategy_id"]), []).append(row)

    out: List[Dict[str, Any]] = []
    for (sensitivity_tag, strategy_id), items in sorted(groups.items()):
        sample = items[0]
        agg: Dict[str, Any] = {
            "sensitivity_tag": sensitivity_tag,
            "strategy_id": strategy_id,
            "strategy_label": sample["strategy_label"],
            "score_method": sample["score_method"],
            "n_scenarios": len(items),
        }
        for metric in metric_names:
            values = [_safe_float(row[metric]) for row in items]
            agg[f"{metric}_mean"] = _mean(values)
            agg[f"{metric}_std"] = _std(values)
            agg[f"{metric}_ci95"] = _ci95(values)
        out.append(agg)
    out.sort(key=lambda row: (row["sensitivity_tag"], float(row["triangle_area_mean"])))
    return out


def _write_csv(path: str, rows: List[Dict[str, Any]]) -> None:
    if not rows:
        return
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def _read_csv(path: str) -> List[Dict[str, Any]]:
    if not os.path.exists(path):
        return []
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _checkpoint_paths(result_dir: str) -> Dict[str, str]:
    return {
        "scenario_manifest_csv": os.path.join(result_dir, "scenario_manifest.csv"),
        "strategy_rows_csv": os.path.join(result_dir, "scenario_strategy_rows.csv"),
        "ranking_csv": os.path.join(result_dir, "criticality_rankings.csv"),
        "aggregate_csv": os.path.join(result_dir, "aggregate_task_c_summary.csv"),
        "summary_md": os.path.join(result_dir, "task_c_summary.md"),
    }


def _save_checkpoint(
    *,
    paths: Dict[str, str],
    cfg: TaskCConfig,
    manifest_rows: List[Dict[str, Any]],
    strategy_rows: List[Dict[str, Any]],
    ranking_rows: List[Dict[str, Any]],
) -> None:
    aggregate_rows = _aggregate_rows(strategy_rows)
    _write_csv(paths["scenario_manifest_csv"], manifest_rows)
    _write_csv(paths["strategy_rows_csv"], strategy_rows)
    _write_csv(paths["ranking_csv"], ranking_rows)
    _write_csv(paths["aggregate_csv"], aggregate_rows)
    _write_summary(paths["summary_md"], cfg, aggregate_rows)


def _write_summary(path: str, cfg: TaskCConfig, aggregate_rows: List[Dict[str, Any]]) -> None:
    with open(path, "w", encoding="utf-8") as f:
        f.write("# Task C O-D Representation Strategy Summary\n\n")
        f.write(f"- Scenario count: `{cfg.n_scenarios}`\n")
        f.write(f"- Scenario seed start: `{cfg.seed0}`\n")
        f.write(f"- Power failures per scenario: `{cfg.bus_count_range}`\n")
        f.write(f"- Road failures per scenario: `{cfg.link_count_range}`\n")
        f.write(f"- Road capacity-drop range: `{cfg.link_drop_range}`\n")
        f.write(f"- Crew model: `{cfg.crew_mode}` with power crews `{cfg.power_crews}` and road crews `{cfg.road_crews}`\n")
        f.write(
            "- S0/S1/S2 are Task A's baselines (imported unmodified); S3 is the O-D "
            "representation strategy from task_c_od.py.\n"
        )
        f.write(
            "- S3 pairs an O-D-representation road score with Task A's unmodified "
            "power radial-service-impact score, isolating the road-side method change.\n"
        )
        f.write(f"- O-D depot node(s): `{list(cfg.od.depot_nodes)}`\n")
        f.write(f"- O-D critical destination set: buses `{list(cfg.od.critical_buses)}`, shelter nodes `{list(cfg.od.shelter_nodes)}` (fixed, damage-independent)\n")
        f.write(f"- O-D candidate paths (k): `{cfg.od.k_paths}`\n")
        f.write(f"- O-D path-rank weights: `{list(cfg.od.path_rank_weights)}`\n")
        f.write(f"- O-D allocation method: `{cfg.od.allocation_method}`\n")
        f.write(
            f"- Weighted-triangle config: critical buses `{list(cfg.weight.critical_buses)}` "
            f"weight x`{cfg.weight.bus_weight_multiplier}`, shelter node `{cfg.weight.shelter_node}`, "
            f"term weights road/power/shelter = `{cfg.weight.weight_road}`/`{cfg.weight.weight_power}`/`{cfg.weight.weight_shelter}`\n\n"
        )
        f.write("## Aggregate Results\n\n")
        f.write("| Sensitivity tag | Strategy | Triangle mean +/- CI | Weighted triangle mean +/- CI | Gini mean +/- CI |\n")
        f.write("|---|---|---:|---:|---:|\n")
        for row in aggregate_rows:
            f.write(
                f"| {row['sensitivity_tag']} | {row['strategy_id']} | "
                f"{float(row['triangle_area_mean']):.3f} +/- {float(row['triangle_area_ci95']):.3f} | "
                f"{float(row['weighted_triangle_area_mean']):.3f} +/- {float(row['weighted_triangle_area_ci95']):.3f} | "
                f"{float(row['gini_restore_mean']):.3f} +/- {float(row['gini_restore_ci95']):.3f} |\n"
            )


def run_task_c_study(cfg: TaskCConfig) -> Dict[str, str]:
    result_dir = cfg.result_dir or os.path.join(cfg.result_root, f"task_c_od_{_ts()}")
    os.makedirs(result_dir, exist_ok=True)
    scenario_root = os.path.join(result_dir, "scenario_runs")
    os.makedirs(scenario_root, exist_ok=True)
    paths = _checkpoint_paths(result_dir)

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

    # The O-D set (and therefore the road-link importance table) is fixed and
    # damage-independent, so it is computed once here rather than per scenario.
    static_road_importance: Dict[Link, float] = build_static_road_importance(cfg.od, cfg.criticality)

    manifest_rows: List[Dict[str, Any]] = _read_csv(paths["scenario_manifest_csv"]) if cfg.resume else []
    strategy_rows: List[Dict[str, Any]] = _read_csv(paths["strategy_rows_csv"]) if cfg.resume else []
    ranking_rows: List[Dict[str, Any]] = _read_csv(paths["ranking_csv"]) if cfg.resume else []

    expected_strategy_ids = {"S3_od_representation"}
    if cfg.include_baselines:
        expected_strategy_ids |= {"S0_centrality", "S1_separate_shapley", "S2_integrated_shapley"}
    completed_scenarios = {
        scenario_id
        for scenario_id in {str(row.get("scenario_id", "")) for row in strategy_rows}
        if {
            str(row.get("strategy_id", ""))
            for row in strategy_rows
            if str(row.get("scenario_id", "")) == scenario_id and str(row.get("sensitivity_tag", "")) == cfg.sensitivity_tag
        }
        >= expected_strategy_ids
    }

    for idx, scenario in enumerate(scenarios, start=1):
        if scenario.scenario_id in completed_scenarios:
            print(f"[task-c] skipping completed {scenario.scenario_id} ({idx}/{len(scenarios)}, tag={cfg.sensitivity_tag})")
            continue
        print(f"[task-c] running {scenario.scenario_id} ({idx}/{len(scenarios)}, tag={cfg.sensitivity_tag}) ...")
        current_root = os.path.join(scenario_root, scenario.scenario_id)
        os.makedirs(current_root, exist_ok=True)

        strategies: List[TaskAStrategy] = list(
            build_task_c_strategies(
                scenario,
                od_cfg=cfg.od,
                task_a_cfg=cfg.criticality,
                static_road_importance=static_road_importance,
            )
        )
        if cfg.include_baselines:
            strategies = list(build_task_a_strategies(scenario, cfg=cfg.criticality)) + strategies

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

        _save_checkpoint(
            paths=paths,
            cfg=cfg,
            manifest_rows=manifest_rows,
            strategy_rows=strategy_rows,
            ranking_rows=ranking_rows,
        )

    return {"result_dir": result_dir, **paths}


def run_task_c_sensitivity(
    base_cfg: TaskCConfig,
    *,
    k_values: Tuple[int, ...] = (1, 2, 3, 5),
) -> Dict[str, str]:
    """Sweep k (candidate path count) for the fixed critical-O-D set. Reuses
    run_task_c_study for each value and concatenates results into one combined
    CSV."""
    result_dir = base_cfg.result_dir or os.path.join(base_cfg.result_root, f"task_c_od_sensitivity_{_ts()}")
    os.makedirs(result_dir, exist_ok=True)
    combined_rows: List[Dict[str, Any]] = []

    for k in k_values:
        tag = f"k{k}"
        variant_cfg = TaskCConfig(
            n_scenarios=base_cfg.n_scenarios,
            seed0=base_cfg.seed0,
            result_root=base_cfg.result_root,
            result_dir=os.path.join(result_dir, tag),
            resume=base_cfg.resume,
            bus_count_range=base_cfg.bus_count_range,
            link_count_range=base_cfg.link_count_range,
            link_drop_range=base_cfg.link_drop_range,
            bus_candidate_path=base_cfg.bus_candidate_path,
            base_net_path=base_cfg.base_net_path,
            include_baselines=False,
            save_strategy_artifacts=False,
            strict=base_cfg.strict,
            crew_mode=base_cfg.crew_mode,
            power_crews=base_cfg.power_crews,
            road_crews=base_cfg.road_crews,
            criticality=base_cfg.criticality,
            od=TaskCODConfig(
                road_net_path=base_cfg.od.road_net_path,
                bus_location_path=base_cfg.od.bus_location_path,
                depot_nodes=base_cfg.od.depot_nodes,
                critical_buses=base_cfg.od.critical_buses,
                shelter_nodes=base_cfg.od.shelter_nodes,
                k_paths=k,
                path_rank_weights=tuple(base_cfg.od.path_rank_weights[:k]) or (1.0,),
                allocation_method=base_cfg.od.allocation_method,
            ),
            weight=base_cfg.weight,
            sensitivity_tag=tag,
        )
        info = run_task_c_study(variant_cfg)
        combined_rows.extend(_read_csv(info["strategy_rows_csv"]))

    combined_csv = os.path.join(result_dir, "taskC_od_sensitivity_summary.csv")
    _write_csv(combined_csv, combined_rows)
    aggregate_csv = os.path.join(result_dir, "taskC_od_sensitivity_aggregate.csv")
    _write_csv(aggregate_csv, _aggregate_rows(combined_rows))

    return {"result_dir": result_dir, "combined_csv": combined_csv, "aggregate_csv": aggregate_csv}


def main() -> None:
    info = run_task_c_study(TaskCConfig())
    print("Task C study complete.")
    print("Result dir:", info["result_dir"])
    print("Scenario manifest:", info["scenario_manifest_csv"])
    print("Strategy rows:", info["strategy_rows_csv"])
    print("Criticality rankings:", info["ranking_csv"])
    print("Aggregate CSV:", info["aggregate_csv"])
    print("Summary:", info["summary_md"])


if __name__ == "__main__":
    main()
