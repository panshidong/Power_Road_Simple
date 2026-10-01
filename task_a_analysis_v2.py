from __future__ import annotations

import csv
import os
import statistics
from collections import defaultdict
from typing import Any, Dict, List, Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


RESULT_DIR = "results/task_a_criticality_final_100"
ANALYSIS_DIR = os.path.join(RESULT_DIR, "analysis_v2")

STRATEGIES = [
    "S0_centrality",
    "S1_separate_shapley",
    "S2_integrated_shapley",
    "heuristic_triangle_sa",
]
LABELS = {
    "S0_centrality": "S0 centrality",
    "S1_separate_shapley": "S1 separate Shapley",
    "S2_integrated_shapley": "S2 integrated Shapley",
    "heuristic_triangle_sa": "Heuristic SA",
}
COLORS = {
    "S0_centrality": "#7F8792",
    "S1_separate_shapley": "#4C78A8",
    "S2_integrated_shapley": "#F58518",
    "heuristic_triangle_sa": "#54A24B",
}


def read_csv(path: str) -> List[Dict[str, str]]:
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: str, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        return
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def style_axes(ax: Any) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(color="#D9D9D9", linewidth=0.7, alpha=0.7)
    ax.set_axisbelow(True)


def load_values() -> tuple[Dict[str, Dict[str, float]], Dict[str, Dict[str, int]]]:
    values: Dict[str, Dict[str, float]] = defaultdict(dict)
    for row in read_csv(os.path.join(RESULT_DIR, "scenario_strategy_rows.csv")):
        values[row["scenario_id"]][row["strategy_id"]] = float(row["triangle_area"])
    damage: Dict[str, Dict[str, int]] = {}
    for row in read_csv(os.path.join(RESULT_DIR, "scenario_manifest.csv")):
        damage[row["scenario_id"]] = {
            "power": int(row["broken_bus_count"]),
            "road": int(row["broken_link_count"]),
        }
    if len(values) != 100 or any(set(item) != set(STRATEGIES) for item in values.values()):
        raise ValueError("The completed main-study table must contain four strategies for 100 scenarios")
    return dict(values), damage


def make_scatter(values: Mapping[str, Mapping[str, float]], damage: Mapping[str, Mapping[str, int]]) -> None:
    scenarios = sorted(values)
    x = np.array([values[item]["S0_centrality"] for item in scenarios])
    y = np.array([values[item]["S2_integrated_shapley"] for item in scenarios])
    count = np.array([damage[item]["power"] + damage[item]["road"] for item in scenarios])
    lower = min(float(x.min()), float(y.min())) * 0.85
    upper = max(float(x.max()), float(y.max())) * 1.15

    fig, ax = plt.subplots(figsize=(7.8, 5.8))
    points = ax.scatter(x, y, c=count, cmap="viridis", s=43, alpha=0.82, edgecolors="white", linewidths=0.35)
    ax.plot([lower, upper], [lower, upper], color="#333333", linestyle="--", linewidth=1.2, label="Equal loss")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lower, upper)
    ax.set_ylim(lower, upper)
    ax.set_xlabel("S0 resilience-loss area (log scale)")
    ax.set_ylabel("S2 resilience-loss area (log scale)")
    ax.set_title("Paired S0 and S2 outcomes across 100 scenarios")
    ax.legend(frameon=False, loc="upper left")
    colorbar = fig.colorbar(points, ax=ax, pad=0.02)
    colorbar.set_label("Number of sampled damaged components")
    style_axes(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(ANALYSIS_DIR, "figure_4_s0_s2_scatter.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def make_ecdf(values: Mapping[str, Mapping[str, float]]) -> None:
    fig, ax = plt.subplots(figsize=(7.8, 5.5))
    for strategy in STRATEGIES:
        ordered = np.sort([item[strategy] for item in values.values()])
        probability = np.arange(1, len(ordered) + 1) / len(ordered)
        ax.step(ordered, probability, where="post", label=LABELS[strategy], color=COLORS[strategy], linewidth=2)
    ax.set_xscale("log")
    ax.set_xlabel("Resilience-loss area (log scale; lower is better)")
    ax.set_ylabel("Empirical cumulative probability")
    ax.set_title("Outcome distributions and upper-tail behavior")
    ax.legend(frameon=False, loc="lower right")
    style_axes(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(ANALYSIS_DIR, "figure_5_loss_ecdf.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def damage_strata(values: Mapping[str, Mapping[str, float]], damage: Mapping[str, Mapping[str, int]]) -> List[Dict[str, Any]]:
    definitions = [
        ("Low (13–17)", lambda total: total <= 17),
        ("Middle (18–21)", lambda total: 18 <= total <= 21),
        ("High (22–27)", lambda total: total >= 22),
    ]
    output: List[Dict[str, Any]] = []
    for label, predicate in definitions:
        ids = [item for item in sorted(values) if predicate(damage[item]["power"] + damage[item]["road"])]
        s0 = [values[item]["S0_centrality"] for item in ids]
        s1 = [values[item]["S1_separate_shapley"] for item in ids]
        s2 = [values[item]["S2_integrated_shapley"] for item in ids]
        heuristic = [values[item]["heuristic_triangle_sa"] for item in ids]
        d20 = [a - b for a, b in zip(s0, s2)]
        d21 = [a - b for a, b in zip(s1, s2)]
        output.append({
            "damage_stratum": label,
            "n": len(ids),
            "S0_mean": statistics.fmean(s0),
            "S1_mean": statistics.fmean(s1),
            "S2_mean": statistics.fmean(s2),
            "heuristic_mean": statistics.fmean(heuristic),
            "S2_vs_S0_mean_reduction": statistics.fmean(d20),
            "S2_vs_S0_wins": sum(value > 0 for value in d20),
            "S2_vs_S1_mean_reduction": statistics.fmean(d21),
            "S2_vs_S1_wins": sum(value > 0 for value in d21),
        })
    return output


def extreme_cases(values: Mapping[str, Mapping[str, float]], damage: Mapping[str, Mapping[str, int]]) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for scenario in sorted(values):
        item = values[scenario]
        rows.append({
            "scenario_id": scenario,
            "total_damage": damage[scenario]["power"] + damage[scenario]["road"],
            "S0": item["S0_centrality"],
            "S1": item["S1_separate_shapley"],
            "S2": item["S2_integrated_shapley"],
            "heuristic": item["heuristic_triangle_sa"],
            "S1_minus_S0": item["S1_separate_shapley"] - item["S0_centrality"],
            "S2_minus_S0": item["S2_integrated_shapley"] - item["S0_centrality"],
            "heuristic_minus_S2": item["heuristic_triangle_sa"] - item["S2_integrated_shapley"],
        })
    return sorted(rows, key=lambda row: abs(float(row["S1_minus_S0"])), reverse=True)


def main() -> None:
    os.makedirs(ANALYSIS_DIR, exist_ok=True)
    values, damage = load_values()
    make_scatter(values, damage)
    make_ecdf(values)
    write_csv(os.path.join(ANALYSIS_DIR, "damage_strata.csv"), damage_strata(values, damage))
    write_csv(os.path.join(ANALYSIS_DIR, "extreme_cases.csv"), extreme_cases(values, damage))
    print(f"Saved v2 analysis artifacts to {ANALYSIS_DIR}")


if __name__ == "__main__":
    main()
