from __future__ import annotations

"""Figures for the Task C paper draft, all read-only against
results/task_c_main_100/aggregate_task_c_summary.csv.

  od_efficiency_figure.png : mean standard resilience-triangle loss by strategy (dot + 95% CI)
  od_weighted_figure.png   : mean critical-weighted composite triangle by strategy (dot + 95% CI)
  od_tradeoff_figure.png   : standard triangle vs restoration-time Gini (supplemental)

Palette: validated categorical slots 1-4 (blue, orange, aqua, yellow); identity is
also carried by marker shape and direct labels, never by hue alone.
"""

import csv

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ORDER = ["S0_centrality", "S1_separate_shapley", "S2_integrated_shapley", "S3_od_representation"]
LABELS = {
    "S0_centrality": "S0 topological centrality",
    "S1_separate_shapley": "S1 separate-network Shapley",
    "S2_integrated_shapley": "S2 integrated Shapley",
    "S3_od_representation": "S3 O-D representation",
}
SHORT = {"S0_centrality": "S0", "S1_separate_shapley": "S1", "S2_integrated_shapley": "S2", "S3_od_representation": "S3"}
COLORS = {
    "S0_centrality": "#2a78d6",
    "S1_separate_shapley": "#eb6834",
    "S2_integrated_shapley": "#1baf7a",
    "S3_od_representation": "#eda100",
}
MARKERS = {"S0_centrality": "o", "S1_separate_shapley": "s", "S2_integrated_shapley": "^", "S3_od_representation": "D"}
INK = "#222222"
MUTED = "#6b6b6b"


def load():
    with open("results/task_c_main_100/aggregate_task_c_summary.csv", newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    return {r["strategy_id"]: r for r in rows}


def _style(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.spines["left"].set_color("#cccccc")
    ax.spines["bottom"].set_color("#cccccc")
    ax.tick_params(colors=INK, labelsize=9)
    ax.grid(True, axis="y", linestyle="-", linewidth=0.6, color="#e6e6e6")
    ax.set_axisbelow(True)


def dot_ci(agg, metric, ylabel, title, out):
    fig, ax = plt.subplots(figsize=(6.0, 3.8), dpi=200)
    xs = list(range(len(ORDER)))
    for x, sid in zip(xs, ORDER):
        m = float(agg[sid][f"{metric}_mean"])
        c = float(agg[sid][f"{metric}_ci95"])
        ax.errorbar(x, m, yerr=c, fmt=MARKERS[sid], color=COLORS[sid], ecolor=COLORS[sid],
                    elinewidth=1.5, capsize=4, markersize=8, markeredgecolor="white", markeredgewidth=1.0)
        ax.annotate(f"{m:.0f}", (x, m), textcoords="offset points", xytext=(10, 0), va="center",
                    fontsize=9, color=INK)
    ax.set_xticks(xs)
    tick_labels = {
        "S0_centrality": "S0\ntopological\ncentrality",
        "S1_separate_shapley": "S1\nseparate-network\nShapley",
        "S2_integrated_shapley": "S2\nintegrated\nShapley",
        "S3_od_representation": "S3\nO-D\nrepresentation",
    }
    ax.set_xticklabels([tick_labels[s] for s in ORDER], fontsize=8.5, color=INK)
    ax.set_ylabel(ylabel, fontsize=9.5, color=INK)
    ax.set_title(title, fontsize=10.5, color=INK, loc="left")
    ax.set_ylim(bottom=0)
    _style(ax)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    print("Wrote", out)


def tradeoff(agg, out):
    fig, ax = plt.subplots(figsize=(6.0, 4.6), dpi=200)
    for sid in ORDER:
        x = float(agg[sid]["triangle_area_mean"]); xe = float(agg[sid]["triangle_area_ci95"])
        y = float(agg[sid]["gini_restore_mean"]); ye = float(agg[sid]["gini_restore_ci95"])
        ax.errorbar(x, y, xerr=xe, yerr=ye, fmt=MARKERS[sid], color=COLORS[sid], ecolor=COLORS[sid],
                    elinewidth=1.2, capsize=3, markersize=8, markeredgecolor="white", markeredgewidth=1.0,
                    label=LABELS[sid])
        ax.annotate(SHORT[sid], (x, y), textcoords="offset points", xytext=(8, 8), fontsize=9, color=INK)
    ax.set_xlabel("Mean resilience-triangle loss (lower is better)", fontsize=9.5, color=INK)
    ax.set_ylabel("Mean restoration-time Gini (lower is better)", fontsize=9.5, color=INK)
    ax.set_title("Resilience efficiency vs. restoration-time inequality (100 scenarios, 95% CI)",
                 fontsize=10, color=INK, loc="left")
    ax.legend(loc="upper right", fontsize=8.5, frameon=False)
    ax.grid(True, linestyle="-", linewidth=0.6, color="#e6e6e6")
    _style(ax)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    print("Wrote", out)


def main() -> None:
    agg = load()
    dot_ci(agg, "triangle_area", "Mean resilience-triangle loss (lower is better)",
           "Standard resilience-triangle loss by strategy (100 scenarios, 95% CI)", "od_efficiency_figure.png")
    dot_ci(agg, "weighted_triangle_area", "Mean critical-weighted composite triangle (lower is better)",
           "Critical-weighted composite triangle by strategy (100 scenarios, 95% CI)", "od_weighted_figure.png")
    tradeoff(agg, "od_tradeoff_figure.png")


if __name__ == "__main__":
    main()
