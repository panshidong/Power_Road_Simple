from __future__ import annotations

import csv
import math
import os
import textwrap
from dataclasses import dataclass, field, replace
from datetime import datetime
from typing import Any, Dict, List, Tuple

from disaster import DisasterScenario, generate_scenarios, scenario_sequence
from simulated_annealing import SAConfig
from tradeoff_runner import TradeoffConfig, run_tradeoff_study, _display_experiment_label, _display_rule_type


def _ts() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def _wrap_title(title: str, width: int = 78) -> str:
    lines: List[str] = []
    for chunk in str(title).split("\n"):
        lines.extend(textwrap.wrap(chunk, width=width) or [""])
    return "\n".join(lines)


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
    n = len(xs)
    if n <= 1:
        return 0.0
    return 1.96 * _std(xs) / math.sqrt(n)


def _label_for_metric(metric: str) -> str:
    labels = {
        "triangle_area": "Triangle area",
        "gini_restore": "Gini restore",
        "p90_restore": "P90 restore",
        "min_time_avg_cri": "Min time-avg CRI",
        "p90_access_restore": "P90 access restore",
        "time_avg_share_access": "Time-avg share access",
    }
    return labels.get(metric, metric)


def _plot_errorbar_tradeoff(
    rows: List[Dict[str, Any]],
    *,
    out_png: str,
    x_metric: str,
    y_metric: str,
) -> str:
    try:
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        return ""

    palette = {
        "baseline": "#1f77b4",
        "single": "#ff7f0e",
        "weighted_sum": "#2ca02c",
        "guardrail": "#d62728",
    }
    markers = {
        "baseline": "s",
        "single": "o",
        "weighted_sum": "^",
        "guardrail": "D",
    }

    fig, ax = plt.subplots(figsize=(18, 11), constrained_layout=True)
    grouped: Dict[str, List[Dict[str, Any]]] = {}
    for row in rows:
        grouped.setdefault(row["rule_type"], []).append(row)

    for rule_type, items in grouped.items():
        items = sorted(items, key=lambda row: float(row[f"{x_metric}_mean"]))
        seen_label = False
        xs = [float(row[f"{x_metric}_mean"]) for row in items]
        ys = [float(row[f"{y_metric}_mean"]) for row in items]
        if len(items) >= 2:
            ax.plot(xs, ys, color=palette.get(rule_type, "#333333"), alpha=0.35, linewidth=1.8)
        for row in items:
            x = float(row[f"{x_metric}_mean"])
            y = float(row[f"{y_metric}_mean"])
            xerr = float(row[f"{x_metric}_ci95"])
            yerr = float(row[f"{y_metric}_ci95"])
            ax.errorbar(
                x,
                y,
                xerr=xerr,
                yerr=yerr,
                fmt=markers.get(rule_type, "o"),
                color=palette.get(rule_type, "#333333"),
                ecolor=palette.get(rule_type, "#333333"),
                capsize=5,
                elinewidth=1.8,
                markersize=10,
                alpha=0.9,
                label=(_display_rule_type(rule_type) if not seen_label else None),
            )
            seen_label = True
            ax.annotate(
                _display_experiment_label(row["experiment_id"]),
                (x, y),
                xytext=(9, 6),
                textcoords="offset points",
                fontsize=12,
                alpha=0.9,
            )

    ax.set_xlabel(f"{_label_for_metric(x_metric)} (mean ± 95% CI)", fontsize=18, labelpad=10)
    ax.set_ylabel(f"{_label_for_metric(y_metric)} (mean ± 95% CI)", fontsize=18, labelpad=10)
    ax.set_title(_wrap_title("Trade-off Curve Across 100 Disaster Scenarios with 95% Confidence Intervals"), fontsize=22, pad=20)
    ax.tick_params(axis="both", labelsize=14)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=4, frameon=True, fontsize=14)
    fig.savefig(out_png, dpi=600, bbox_inches="tight", pad_inches=0.35)
    plt.close(fig)
    return out_png


@dataclass
class BatchTradeoffConfig:
    n_scenarios: int = 100
    seed0: int = 20260325
    result_root: str = "results"
    bus_count_range: Tuple[int, int] = (8, 15)
    link_count_range: Tuple[int, int] = (3, 13)
    link_drop_range: Tuple[float, float] = (0.5, 1.0)
    bus_candidate_path: str = "new_bus_to_link.json"
    base_net_path: str = "tap-b/net/SiouxFalls_net.txt"
    aggregate_x_metric: str = "triangle_area"
    aggregate_y_metric: str = "gini_restore"
    study_cfg: TradeoffConfig = field(
        default_factory=lambda: TradeoffConfig(
            crew_mode="specialized",
            power_crews=1,
            road_crews=1,
            multifunction_crews=1,
            cri_weight_pairs=[(0.133, 0.867)],
            run_sensitivity=False,
            save_reference_artifacts=False,
            save_baseline=False,
            save_best_artifacts=False,
            save_best_debug=False,
            sa=SAConfig(seed=0, max_iter=80, T0=2.0, alpha=0.98, neighbor="swap"),
        )
    )


def _scenario_row(row: Dict[str, Any], *, scenario: DisasterScenario, study_result_dir: str) -> Dict[str, Any]:
    out = dict(row)
    out["scenario_id"] = scenario.scenario_id
    out["scenario_seed"] = scenario.seed
    out["broken_bus_count"] = len(scenario.broken_buses)
    out["broken_link_count"] = len(scenario.broken_links)
    factors = list(scenario.link_capacity_factors.values())
    out["mean_remaining_capacity_factor"] = _mean(factors)
    out["mean_capacity_drop_fraction"] = _mean([1.0 - f for f in factors]) if factors else float("nan")
    out["study_result_dir"] = study_result_dir
    return out


def _aggregate_rows(rows: List[Dict[str, Any]], *, x_metric: str, y_metric: str) -> List[Dict[str, Any]]:
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
        groups.setdefault(row["experiment_id"], []).append(row)

    aggregated: List[Dict[str, Any]] = []
    for experiment_id, items in sorted(groups.items()):
        sample = items[0]
        agg: Dict[str, Any] = {
            "experiment_id": experiment_id,
            "rule_type": sample["rule_type"],
            "rule_label": sample["rule_label"],
            "description": sample["description"],
            "n_scenarios": len(items),
        }
        for metric in metric_names:
            values = [_safe_float(row[metric]) for row in items]
            agg[f"{metric}_mean"] = _mean(values)
            agg[f"{metric}_std"] = _std(values)
            agg[f"{metric}_ci95"] = _ci95(values)
        aggregated.append(agg)

    aggregated.sort(key=lambda row: float(row[f"{x_metric}_mean"]))
    return aggregated


def run_batch_tradeoff_study(cfg: BatchTradeoffConfig) -> Dict[str, str]:
    result_dir = os.path.join(cfg.result_root, f"batch_tradeoff_{_ts()}")
    os.makedirs(result_dir, exist_ok=True)
    scenario_run_root = os.path.join(result_dir, "scenario_runs")
    os.makedirs(scenario_run_root, exist_ok=True)

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

    scenario_rows: List[Dict[str, Any]] = []
    manifest_rows: List[Dict[str, Any]] = []

    for idx, scenario in enumerate(scenarios, start=1):
        scenario_cfg = replace(
            cfg.study_cfg,
            result_root=scenario_run_root,
            base_sequence=scenario_sequence(scenario),
            broken_link_factors=dict(scenario.link_capacity_factors),
        )
        print(f"[batch] running {scenario.scenario_id} ({idx}/{len(scenarios)}) ...")
        study_info = run_tradeoff_study(scenario_cfg)
        manifest_rows.append(
            {
                "scenario_id": scenario.scenario_id,
                "scenario_seed": scenario.seed,
                "broken_bus_count": len(scenario.broken_buses),
                "broken_link_count": len(scenario.broken_links),
                "study_result_dir": study_info["result_dir"],
                "csv_path": study_info["csv_path"],
            }
        )
        with open(study_info["csv_path"], newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                if row.get("row_kind") != "optimized":
                    continue
                scenario_rows.append(_scenario_row(row, scenario=scenario, study_result_dir=study_info["result_dir"]))

    manifest_csv = os.path.join(result_dir, "scenario_manifest.csv")
    with open(manifest_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(manifest_rows[0].keys()))
        writer.writeheader()
        writer.writerows(manifest_rows)

    scenario_csv = os.path.join(result_dir, "scenario_tradeoff_rows.csv")
    with open(scenario_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(scenario_rows[0].keys()))
        writer.writeheader()
        writer.writerows(scenario_rows)

    aggregate_rows = _aggregate_rows(
        [row for row in scenario_rows if row["row_kind"] == "optimized"],
        x_metric=cfg.aggregate_x_metric,
        y_metric=cfg.aggregate_y_metric,
    )
    aggregate_csv = os.path.join(result_dir, "aggregate_tradeoff_summary.csv")
    with open(aggregate_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(aggregate_rows[0].keys()))
        writer.writeheader()
        writer.writerows(aggregate_rows)

    errorbar_png = os.path.join(result_dir, "aggregate_tradeoff_errorbars.png")
    plot_path = _plot_errorbar_tradeoff(
        aggregate_rows,
        out_png=errorbar_png,
        x_metric=cfg.aggregate_x_metric,
        y_metric=cfg.aggregate_y_metric,
    )

    summary_md = os.path.join(result_dir, "batch_summary.md")
    with open(summary_md, "w", encoding="utf-8") as f:
        f.write("# Batch Trade-off Summary\n\n")
        f.write(f"- Scenario count: `{cfg.n_scenarios}`\n")
        f.write(f"- Scenario seed start: `{cfg.seed0}`\n")
        f.write(f"- Power failures per disaster: `{cfg.bus_count_range}`\n")
        f.write(f"- Road failures per disaster: `{cfg.link_count_range}`\n")
        f.write(f"- Road capacity-drop fraction range: `{cfg.link_drop_range}`\n")
        f.write(f"- Crew mode: `{cfg.study_cfg.crew_mode}` (power=`{cfg.study_cfg.power_crews}`, road=`{cfg.study_cfg.road_crews}`)\n")
        f.write(f"- CRI weights: `{cfg.study_cfg.cri_weight_pairs}`\n")
        f.write(f"- Scenario manifest: `{os.path.basename(manifest_csv)}`\n")
        f.write(f"- Per-scenario rows: `{os.path.basename(scenario_csv)}`\n")
        f.write(f"- Aggregate summary: `{os.path.basename(aggregate_csv)}`\n")
        f.write(f"- Error-bar trade-off figure: `{os.path.basename(plot_path) if plot_path else ''}`\n")
        f.write("\n## Aggregate Means\n\n")
        f.write("| Experiment | Type | Triangle mean ± CI | Gini mean ± CI | P90 access mean ± CI |\n")
        f.write("|---|---|---:|---:|---:|\n")
        for row in aggregate_rows:
            f.write(
                f"| {row['experiment_id']} | {row['rule_type']} | "
                f"{float(row['triangle_area_mean']):.3f} ± {float(row['triangle_area_ci95']):.3f} | "
                f"{float(row['gini_restore_mean']):.3f} ± {float(row['gini_restore_ci95']):.3f} | "
                f"{float(row['p90_access_restore_mean']):.3f} ± {float(row['p90_access_restore_ci95']):.3f} |\n"
            )

    return {
        "result_dir": result_dir,
        "scenario_manifest_csv": manifest_csv,
        "scenario_rows_csv": scenario_csv,
        "aggregate_csv": aggregate_csv,
        "errorbar_png": plot_path,
        "summary_md": summary_md,
    }


def main() -> None:
    info = run_batch_tradeoff_study(BatchTradeoffConfig())
    print("Batch trade-off study complete.")
    print("Result dir:", info["result_dir"])
    print("Scenario manifest:", info["scenario_manifest_csv"])
    print("Scenario rows:", info["scenario_rows_csv"])
    print("Aggregate CSV:", info["aggregate_csv"])
    print("Error-bar plot:", info["errorbar_png"])
    print("Summary:", info["summary_md"])


if __name__ == "__main__":
    main()
