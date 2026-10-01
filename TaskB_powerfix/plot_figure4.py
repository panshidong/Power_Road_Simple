from __future__ import annotations

"""Paper Figure 4 (main-simulation trade-off, mean ± 95% CI) with fonts sized for a
6-inch-wide figure in the manuscript (reviewer comment: labels too small)."""

import csv
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
AGG = os.path.join(HERE, "results", "merged", "aggregate_tradeoff_summary.csv")
OUT = os.path.join(HERE, "results", "merged", "figure4_paper.png")

LABEL = {
    "baseline_reference": "Baseline", "single_triangle": "Triangle", "single_gini_restore": "Gini",
    "single_maximin_time_avg_cri": "Maximin CRI", "weighted_maximin_time_avg_cri_l050": "Weighted maximin",
    "guardrail_gini_restore": "Gini guardrail", "weighted_gini_restore_l050": "Weighted Gini",
    "single_p90_access_restore": "P90 access",
}
TYPE_LABEL = {"baseline": "Baseline", "single": "Single objective", "weighted_sum": "Weighted sum", "guardrail": "Guardrail"}
COLOR = {"baseline": "#2a78d6", "single": "#eb6834", "weighted_sum": "#1baf7a", "guardrail": "#eda100"}
MARK = {"baseline": "s", "single": "o", "weighted_sum": "^", "guardrail": "D"}
OFFSET = {"single_triangle": (12, 6), "single_p90_access_restore": (-34, -52), "weighted_maximin_time_avg_cri_l050": (-40, -26),
          "guardrail_gini_restore": (12, -20), "weighted_gini_restore_l050": (12, 10), "single_maximin_time_avg_cri": (12, 8),
          "baseline_reference": (12, 8), "single_gini_restore": (-60, 12)}


def main():
    rows = list(csv.DictReader(open(AGG, newline="", encoding="utf-8")))
    fig, ax = plt.subplots(figsize=(9, 6.5), dpi=300)
    seen = set()
    for r in rows:
        t = r["rule_type"]
        x, xe = float(r["triangle_area_mean"]), float(r["triangle_area_ci95"])
        y, ye = float(r["gini_restore_mean"]), float(r["gini_restore_ci95"])
        ax.errorbar(x, y, xerr=xe, yerr=ye, fmt=MARK[t], color=COLOR[t], ecolor=COLOR[t], capsize=5, elinewidth=1.6,
                    markersize=11, markeredgecolor="white", label=TYPE_LABEL[t] if t not in seen else None)
        seen.add(t)
        ax.annotate(LABEL[r["experiment_id"]], (x, y), xytext=OFFSET.get(r["experiment_id"], (10, 8)),
                    textcoords="offset points", fontsize=14, color="#222222")
    ax.set_xlabel("Resilience-triangle area (mean ± 95% CI)", fontsize=16)
    ax.set_ylabel("Restoration-time Gini (mean ± 95% CI)", fontsize=16)
    ax.tick_params(labelsize=13)
    ax.grid(True, linestyle="-", linewidth=0.6, color="#e6e6e6")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(loc="upper right", fontsize=13, frameon=False)
    fig.tight_layout()
    fig.savefig(OUT)
    print("Wrote", OUT)


if __name__ == "__main__":
    main()
