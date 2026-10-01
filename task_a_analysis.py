from __future__ import annotations

import csv
import json
import math
import os
import statistics
from collections import defaultdict
from typing import Any, Dict, Iterable, List, Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


RESULT_DIR = "results/task_a_criticality_final_100"
STABILITY_DIR = "results/task_a_shapley_stability_10"
ANALYSIS_DIR = os.path.join(RESULT_DIR, "analysis")

STRATEGY_ORDER = [
    "S0_centrality",
    "S1_separate_shapley",
    "S2_integrated_shapley",
    "heuristic_triangle_sa",
]
STRATEGY_LABELS = {
    "S0_centrality": "S0: network-specific\ncentrality",
    "S1_separate_shapley": "S1: sector-specific\nShapley",
    "S2_integrated_shapley": "S2: integrated\nShapley",
    "heuristic_triangle_sa": "Heuristic SA\nbenchmark",
}
STRATEGY_SHORT = {
    "S0_centrality": "S0 centrality",
    "S1_separate_shapley": "S1 separate Shapley",
    "S2_integrated_shapley": "S2 integrated Shapley",
    "heuristic_triangle_sa": "Heuristic SA",
}
COLORS = {
    "S0_centrality": "#8A8F98",
    "S1_separate_shapley": "#4C78A8",
    "S2_integrated_shapley": "#F58518",
    "heuristic_triangle_sa": "#54A24B",
}


def _read_csv(path: str) -> List[Dict[str, str]]:
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _write_csv(path: str, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        return
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def _sample_std(values: Sequence[float]) -> float:
    return statistics.stdev(values) if len(values) > 1 else 0.0


def _t_critical_975(n: int) -> float:
    # Exact values for the sample sizes used here; normal limit is adequate above 120.
    lookup = {
        2: 12.706,
        3: 4.303,
        4: 3.182,
        5: 2.776,
        6: 2.571,
        7: 2.447,
        8: 2.365,
        9: 2.306,
        10: 2.262,
        20: 2.093,
        30: 2.045,
        50: 2.010,
        100: 1.984,
    }
    if n in lookup:
        return lookup[n]
    if n > 120:
        return 1.96
    keys = sorted(key for key in lookup if key <= n)
    return lookup[keys[-1]] if keys else 1.96


def _mean_ci(values: Sequence[float]) -> tuple[float, float]:
    mean = statistics.fmean(values)
    half = _t_critical_975(len(values)) * _sample_std(values) / math.sqrt(len(values))
    return mean, half


def _exact_sign_p(wins: int, losses: int) -> float:
    n = wins + losses
    if n == 0:
        return 1.0
    k = min(wins, losses)
    tail = sum(math.comb(n, i) for i in range(k + 1)) / (2.0**n)
    return min(1.0, 2.0 * tail)


def _strategy_values(rows: Sequence[Mapping[str, str]]) -> Dict[str, Dict[str, float]]:
    output: Dict[str, Dict[str, float]] = defaultdict(dict)
    for row in rows:
        output[str(row["scenario_id"])][str(row["strategy_id"])] = float(row["triangle_area"])
    return dict(output)


def _verify(rows: Sequence[Mapping[str, str]]) -> Dict[str, Any]:
    values = _strategy_values(rows)
    if len(values) != 100:
        raise ValueError(f"Expected 100 scenarios, found {len(values)}")
    required = set(STRATEGY_ORDER)
    incomplete = {
        scenario: sorted(required - set(strategy_values))
        for scenario, strategy_values in values.items()
        if not required.issubset(strategy_values)
    }
    if incomplete:
        raise ValueError(f"Incomplete strategy rows: {incomplete}")
    if len(rows) != 400:
        raise ValueError(f"Expected exactly 400 scenario-strategy rows, found {len(rows)}")
    snapshot_counts = [int(float(row["tapb_state_snapshots"])) for row in rows]
    event_counts = [int(float(row["event_count"])) for row in rows]
    return {
        "n_scenarios": len(values),
        "n_strategy_runs": len(rows),
        "tapb_snapshot_min": min(snapshot_counts),
        "tapb_snapshot_max": max(snapshot_counts),
        "event_count_min": min(event_counts),
        "event_count_max": max(event_counts),
    }


def _descriptive_rows(values: Mapping[str, Mapping[str, float]]) -> List[Dict[str, Any]]:
    output: List[Dict[str, Any]] = []
    for strategy in STRATEGY_ORDER:
        sample = [scenario_values[strategy] for scenario_values in values.values()]
        mean, ci = _mean_ci(sample)
        output.append(
            {
                "strategy_id": strategy,
                "strategy_label": STRATEGY_SHORT[strategy],
                "n": len(sample),
                "triangle_mean": mean,
                "triangle_std": _sample_std(sample),
                "triangle_ci95_halfwidth": ci,
                "triangle_ci95_lower": mean - ci,
                "triangle_ci95_upper": mean + ci,
                "triangle_median": statistics.median(sample),
                "triangle_q1": float(np.quantile(sample, 0.25)),
                "triangle_q3": float(np.quantile(sample, 0.75)),
                "triangle_min": min(sample),
                "triangle_max": max(sample),
            }
        )
    return output


def _paired_row(
    values: Mapping[str, Mapping[str, float]],
    *,
    candidate: str,
    reference: str,
) -> Dict[str, Any]:
    improvements = [
        scenario_values[reference] - scenario_values[candidate]
        for scenario_values in values.values()
    ]
    percentages = [
        100.0 * (scenario_values[reference] - scenario_values[candidate]) / scenario_values[reference]
        for scenario_values in values.values()
        if scenario_values[reference] != 0.0
    ]
    wins = sum(value > 1e-9 for value in improvements)
    losses = sum(value < -1e-9 for value in improvements)
    ties = len(improvements) - wins - losses
    mean, ci = _mean_ci(improvements)
    percent_mean, percent_ci = _mean_ci(percentages)
    return {
        "candidate": candidate,
        "candidate_label": STRATEGY_SHORT[candidate],
        "reference": reference,
        "reference_label": STRATEGY_SHORT[reference],
        "n": len(improvements),
        "mean_triangle_reduction": mean,
        "triangle_reduction_ci95_halfwidth": ci,
        "triangle_reduction_ci95_lower": mean - ci,
        "triangle_reduction_ci95_upper": mean + ci,
        "median_triangle_reduction": statistics.median(improvements),
        "mean_percent_reduction": percent_mean,
        "percent_reduction_ci95_halfwidth": percent_ci,
        "median_percent_reduction": statistics.median(percentages),
        "wins": wins,
        "losses": losses,
        "ties": ties,
        "win_rate": wins / len(improvements),
        "exact_sign_test_p": _exact_sign_p(wins, losses),
    }


def _paired_rows(values: Mapping[str, Mapping[str, float]]) -> List[Dict[str, Any]]:
    comparisons = [
        ("S1_separate_shapley", "S0_centrality"),
        ("S2_integrated_shapley", "S0_centrality"),
        ("S2_integrated_shapley", "S1_separate_shapley"),
        ("heuristic_triangle_sa", "S2_integrated_shapley"),
    ]
    return [_paired_row(values, candidate=candidate, reference=reference) for candidate, reference in comparisons]


def _scenario_rows(
    values: Mapping[str, Mapping[str, float]],
    manifest: Sequence[Mapping[str, str]],
) -> List[Dict[str, Any]]:
    damage = {row["scenario_id"]: row for row in manifest}
    output: List[Dict[str, Any]] = []
    for scenario_id in sorted(values):
        row = damage[scenario_id]
        item: Dict[str, Any] = {
            "scenario_id": scenario_id,
            "broken_bus_count": int(row["broken_bus_count"]),
            "broken_link_count": int(row["broken_link_count"]),
        }
        item.update(values[scenario_id])
        item["S1_vs_S0_percent"] = 100.0 * (
            values[scenario_id]["S0_centrality"] - values[scenario_id]["S1_separate_shapley"]
        ) / values[scenario_id]["S0_centrality"]
        item["S2_vs_S0_percent"] = 100.0 * (
            values[scenario_id]["S0_centrality"] - values[scenario_id]["S2_integrated_shapley"]
        ) / values[scenario_id]["S0_centrality"]
        item["S2_vs_S1_percent"] = 100.0 * (
            values[scenario_id]["S1_separate_shapley"] - values[scenario_id]["S2_integrated_shapley"]
        ) / values[scenario_id]["S1_separate_shapley"]
        output.append(item)
    return output


def _style_axes(ax: Any) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.7, alpha=0.7)
    ax.set_axisbelow(True)


def _figure_means(descriptive: Sequence[Mapping[str, Any]], output_path: str) -> None:
    fig, ax = plt.subplots(figsize=(8.4, 5.0))
    x = np.arange(len(descriptive))
    means = [float(row["triangle_mean"]) for row in descriptive]
    errors = [float(row["triangle_ci95_halfwidth"]) for row in descriptive]
    colors = [COLORS[str(row["strategy_id"])] for row in descriptive]
    ax.bar(x, means, yerr=errors, capsize=5, color=colors, width=0.68, edgecolor="white")
    ax.set_xticks(x, [STRATEGY_LABELS[str(row["strategy_id"])] for row in descriptive])
    ax.set_ylabel("Mean resilience-loss area")
    ax.set_title("Recovery performance across 100 disruption scenarios")
    for index, (mean, error) in enumerate(zip(means, errors)):
        ax.text(index, mean + error + max(means) * 0.02, f"{mean:,.0f}", ha="center", va="bottom", fontsize=9)
    ax.margins(y=0.16)
    _style_axes(ax)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _figure_paired(scenario_rows: Sequence[Mapping[str, Any]], output_path: str) -> None:
    columns = ["S1_vs_S0_percent", "S2_vs_S0_percent", "S2_vs_S1_percent"]
    labels = ["S1 relative to S0", "S2 relative to S0", "S2 relative to S1"]
    data = [[float(row[column]) for row in scenario_rows] for column in columns]
    fig, ax = plt.subplots(figsize=(8.4, 5.0))
    box = ax.boxplot(data, patch_artist=True, showfliers=True, whis=(5, 95), widths=0.55)
    for patch, color in zip(box["boxes"], [COLORS["S1_separate_shapley"], COLORS["S2_integrated_shapley"], "#ECA82C"]):
        patch.set_facecolor(color)
        patch.set_alpha(0.78)
    for index, values_at_comparison in enumerate(data, start=1):
        ax.scatter(index, statistics.fmean(values_at_comparison), marker="D", s=34, color="black", zorder=3)
    ax.axhline(0.0, color="#333333", linewidth=1.0)
    ax.set_xticks(range(1, len(labels) + 1), labels)
    ax.set_ylabel("Reduction in resilience-loss area (%)")
    ax.set_title("Paired scenario-level changes (positive values favor the named strategy)")
    _style_axes(ax)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _figure_stability(stability: Sequence[Mapping[str, str]], output_path: str) -> None:
    fig, ax = plt.subplots(figsize=(8.4, 5.0))
    styles = {
        ("S1_separate_power", "power"): ("S1 power", "#4C78A8", "o"),
        ("S1_separate_road", "road"): ("S1 road", "#72A0C1", "s"),
        ("S2_integrated", "power"): ("S2 power", "#F58518", "o"),
        ("S2_integrated", "road"): ("S2 road", "#ECA82C", "s"),
    }
    for key, (label, color, marker) in styles.items():
        rows = sorted(
            [row for row in stability if (row["method"], row["sector"]) == key],
            key=lambda row: int(row["sample_size"]),
        )
        ax.errorbar(
            [int(row["sample_size"]) for row in rows],
            [float(row["spearman_rho_mean"]) for row in rows],
            yerr=[float(row["spearman_rho_ci95"]) for row in rows],
            label=label,
            color=color,
            marker=marker,
            linewidth=2,
            capsize=4,
        )
    ax.set_xticks([30, 60, 120, 240])
    ax.set_ylim(0.78, 1.04)
    ax.set_xlabel("Sampled permutations")
    ax.set_ylabel("Mean Spearman correlation with 240-sample ranking")
    ax.set_title("Shapley ranking stability across ten scenarios")
    ax.legend(frameon=False, ncol=2)
    _style_axes(ax)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    os.makedirs(ANALYSIS_DIR, exist_ok=True)
    rows = _read_csv(os.path.join(RESULT_DIR, "scenario_strategy_rows.csv"))
    manifest = _read_csv(os.path.join(RESULT_DIR, "scenario_manifest.csv"))
    verification = _verify(rows)
    values = _strategy_values(rows)
    descriptive = _descriptive_rows(values)
    paired = _paired_rows(values)
    scenarios = _scenario_rows(values, manifest)

    stability_path = os.path.join(STABILITY_DIR, "stability_aggregate.csv")
    stability = _read_csv(stability_path)
    if not stability:
        raise ValueError(f"Missing stability results: {stability_path}")

    _write_csv(os.path.join(ANALYSIS_DIR, "descriptive_statistics.csv"), descriptive)
    _write_csv(os.path.join(ANALYSIS_DIR, "paired_comparisons.csv"), paired)
    _write_csv(os.path.join(ANALYSIS_DIR, "scenario_comparisons.csv"), scenarios)
    _figure_means(descriptive, os.path.join(ANALYSIS_DIR, "figure_1_mean_triangle.png"))
    _figure_paired(scenarios, os.path.join(ANALYSIS_DIR, "figure_2_paired_changes.png"))
    _figure_stability(stability, os.path.join(ANALYSIS_DIR, "figure_3_shapley_stability.png"))

    manifest_bus = [int(row["broken_bus_count"]) for row in manifest]
    manifest_link = [int(row["broken_link_count"]) for row in manifest]
    payload = {
        "verification": verification,
        "scenario_design": {
            "seed0": 20260402,
            "n_scenarios": len(manifest),
            "power_failure_range": [min(manifest_bus), max(manifest_bus)],
            "road_failure_range": [min(manifest_link), max(manifest_link)],
            "power_failure_mean": statistics.fmean(manifest_bus),
            "road_failure_mean": statistics.fmean(manifest_link),
            "road_capacity_drop_range": [0.5, 1.0],
            "shapley_production_samples": 120,
            "stability_samples": [30, 60, 120, 240],
            "stability_scenarios": 10,
            "heuristic_iterations": 60,
        },
        "descriptive": descriptive,
        "paired": paired,
        "stability": stability,
    }
    with open(os.path.join(ANALYSIS_DIR, "results_summary.json"), "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
