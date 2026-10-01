from __future__ import annotations

"""Figures for the OD paper (v4, pre-event-table baselines), read from results/merged/stats.json.

  od4_efficiency_figure.png : mean triangle loss, six strategies (dot + paired-bootstrap 95% CI)
  od4_paired_figure.png     : paired mean reduction of the O-D road order vs each table, holding the
                              power order fixed (dot + 95% CI, symlog scale not needed: all positive)
  od4_weighted_figure.png   : mean critical-weighted composite triangle, six strategies
  od4_tradeoff_figure.png   : triangle loss vs restoration-time Gini (supplemental)
Palette: validated categorical slots (blue, orange, aqua, yellow); identity also carried by marker
shape and direct labels.
"""

import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = "/home/workenv/OD"
ORDER = ["CEN", "JSH", "IJSH", "OD_ALL", "OD_CEN"]
TICK = {"CEN": "Structural", "JSH": "Service\nShapley", "IJSH": "Access-\naugmented\nShapley",
        "OD_CEN": "Origin-\ndestination\n(depot-origin\ntrips)", "OD_ALL": "Origin-\ndestination\n(all trips)", "OD_JSH": "O-D road,\nservice-\nShapley power", "OD_IJSH": "O-D road,\naccess-\naugmented power"}
NAME = {"CEN": "Structural", "JSH": "Service Shapley", "IJSH": "Access-augmented", "OD_CEN": "O-D, depot-origin trips", "OD_ALL": "O-D, all trips",
        "OD_JSH": "O-D road, service-Shapley power", "OD_IJSH": "O-D road, access-augmented power"}
COL = {"CEN": "#2a78d6", "JSH": "#eb6834", "IJSH": "#1baf7a", "OD_CEN": "#eda100", "OD_ALL": "#eda100", "OD_JSH": "#eda100", "OD_IJSH": "#eda100"}
MK = {"CEN": "o", "JSH": "s", "IJSH": "^", "OD_CEN": "D", "OD_ALL": "v", "OD_JSH": "D", "OD_IJSH": "D"}
INK = "#222222"


def _style(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.spines["left"].set_color("#cccccc"); ax.spines["bottom"].set_color("#cccccc")
    ax.tick_params(colors=INK, labelsize=9)
    ax.grid(True, axis="y", linestyle="-", linewidth=0.6, color="#e6e6e6"); ax.set_axisbelow(True)


def dot_ci(desc, ylabel, title, out, fmt="{:.0f}"):
    d = {r["strategy"]: r for r in desc}
    fig, ax = plt.subplots(figsize=(6.6, 3.9), dpi=200)
    for x, s in enumerate(ORDER):
        m, lo, hi = d[s]["mean"], d[s]["ci_lo"], d[s]["ci_hi"]
        ax.errorbar(x, m, yerr=[[m - lo], [hi - m]], fmt=MK[s], color=COL[s], ecolor=COL[s], elinewidth=1.5, capsize=4,
                    markersize=8, markeredgecolor="white", markeredgewidth=1.0)
        ax.annotate(fmt.format(m), (x, m), textcoords="offset points", xytext=(9, 0), va="center", fontsize=8.5, color=INK)
        ax.set_xticks(range(len(ORDER))); ax.set_xticklabels([TICK[s] for s in ORDER], fontsize=8, color=INK)
    ax.set_ylabel(ylabel, fontsize=9.5, color=INK); ax.set_title(title, fontsize=10, color=INK, loc="left")
    ax.set_ylim(bottom=0); _style(ax); fig.tight_layout(); fig.savefig(out); plt.close(fig); print("Wrote", out)


def paired_fig(pairs, out):
    want = [("OD_CEN", "CEN"), ("OD_CEN", "JSH"), ("OD_CEN", "IJSH")]
    fig, ax = plt.subplots(figsize=(6.0, 3.6), dpi=200)
    for x, (c, r) in enumerate(want):
        p = next(q for q in pairs if q["candidate"] == c and q["reference"] == r)
        ax.errorbar(x, p["mean_reduction"], yerr=[[p["mean_reduction"] - p["ci_lo"]], [p["ci_hi"] - p["mean_reduction"]]],
                    fmt="D", color=COL[r], ecolor=COL[r], elinewidth=1.5, capsize=4, markersize=8, markeredgecolor="white")
        ax.plot(x, p["median"], marker="s", color=COL[r], markersize=6, markerfacecolor="white")
        ax.annotate(f"{p['pct']:.1f}%\n{p['wins']}/{p['losses']}/{p['ties']}", (x, p["ci_hi"]), textcoords="offset points",
                    xytext=(0, 6), ha="center", fontsize=8.5, color=INK)
    ax.axhline(0, color="#888888", linewidth=0.8)
    ax.set_xticks(range(3)); ax.set_xticklabels(["vs. structural\npriorities", "vs. service Shapley\npriorities", "vs. access-augmented\nShapley priorities"], fontsize=8.5, color=INK)
    ax.set_ylabel("Service loss avoided\n(positive favours origin-destination priorities)", fontsize=9, color=INK)
    ax.set_title("Origin-destination priorities vs. each baseline rule (300 disruptions)", fontsize=10, color=INK, loc="left")
    ax.set_ylim(top=ax.get_ylim()[1] * 1.18)
    _style(ax); fig.tight_layout(); fig.savefig(out); plt.close(fig); print("Wrote", out)


def tradeoff(st, out):
    t = {r["strategy"]: r for r in st["descriptive"]}; g = {r["strategy"]: r for r in st["descriptive_gini"]}
    fig, ax = plt.subplots(figsize=(6.0, 4.4), dpi=200)
    for s in ORDER:
        x, y = t[s]["mean"], g[s]["mean"]
        ax.errorbar(x, y, xerr=[[x - t[s]["ci_lo"]], [t[s]["ci_hi"] - x]], yerr=[[y - g[s]["ci_lo"]], [g[s]["ci_hi"] - y]], fmt=MK[s], color=COL[s],
                    ecolor=COL[s], elinewidth=1.1, capsize=3, markersize=8, markeredgecolor="white", markeredgewidth=1.0)
        ax.annotate(NAME[s], (x, y), textcoords="offset points", xytext=(7, 6), fontsize=8.5, color=INK)
    ax.set_xlabel("Mean resilience loss (lower is better)", fontsize=9.5, color=INK)
    ax.set_ylabel("Gini of zone restoration times (lower is more even)", fontsize=9.5, color=INK)
    ax.set_title("Resilience loss against evenness of restoration (300 disruptions, 95% intervals)", fontsize=10, color=INK, loc="left")
    ax.grid(True, linestyle="-", linewidth=0.6, color="#e6e6e6"); _style(ax); fig.tight_layout(); fig.savefig(out); plt.close(fig); print("Wrote", out)


import sys
SFX = "_allpairs" if "--allpairs" in sys.argv else ""


def main():
    st = json.load(open(os.path.join(HERE, "results", "merged", f"stats{SFX}.json")))
    ap = json.load(open(os.path.join(HERE, "results", "merged", "stats_allpairs.json")))
    for key in ("descriptive", "descriptive_weighted", "descriptive_gini"):
        row = dict(next(r for r in ap[key] if r["strategy"] == "OD_CEN")); row["strategy"] = "OD_ALL"; st[key].append(row)
    n = st["n"]
    dot_ci(st["descriptive"], "Mean resilience loss (lower is better)", f"Resilience loss by priority rule ({n} disruptions, 95% interval)", os.path.join(OUT, "od4_efficiency_figure.png"))
    paired_fig(st["paired"], os.path.join(OUT, "od4_paired_figure.png"))
    dot_ci(st["descriptive_weighted"], "Mean critical-weighted resilience loss (lower is better)", f"Critical-weighted resilience loss by priority rule ({n} disruptions, 95% interval)", os.path.join(OUT, "od4_weighted_figure.png"))
    tradeoff(st, os.path.join(OUT, "od4_tradeoff_figure.png"))


if __name__ == "__main__":
    main()


def shift_fig(out):
    st = json.load(open(os.path.join(HERE, "results", "merged", f"shift_stats{SFX}.json")))
    main = json.load(open(os.path.join(HERE, "results", "merged", f"stats{SFX}.json")))
    conds = [("main", "Main\n(300)"), ("small", "Smaller\n(100)"), ("large", "Larger\n(100)"), ("light", "Lighter\n(100)"), ("clustered", "Clustered\n(100)")]
    refs = [("CEN", "vs. structural priorities"), ("JSH", "vs. service Shapley priorities"), ("IJSH", "vs. access-augmented Shapley priorities")]
    fig, ax = plt.subplots(figsize=(6.6, 3.9), dpi=200)
    for j, (r, lab) in enumerate(refs):
        xs, ys, lo, hi = [], [], [], []
        for i, (c, _) in enumerate(conds):
            pairs = main["paired"] if c == "main" else st[c]["paired"]
            p = next(q for q in pairs if q["candidate"] == "OD_CEN" and q["reference"] == r)
            xs.append(i + (j - 1) * 0.22); ys.append(p["pct"]); lo.append(p["pct"] - p["pct_lo"]); hi.append(p["pct_hi"] - p["pct"])
        ax.errorbar(xs, ys, yerr=[lo, hi], fmt=MK[r], color=COL[r], ecolor=COL[r], elinewidth=1.3, capsize=3, markersize=7,
                    markeredgecolor="white", linestyle="none", label=lab)
    ax.axhline(0, color="#888888", linewidth=0.8)
    ax.set_xticks(range(len(conds))); ax.set_xticklabels([t for _, t in conds], fontsize=8.5, color=INK)
    ax.set_ylabel("Reduction in mean resilience loss (%)\n(positive favours origin-destination priorities)", fontsize=9, color=INK)
    ax.set_title("Origin-destination priorities vs. each rule under different damage conditions", fontsize=10, color=INK, loc="left")
    ax.legend(loc="upper left", fontsize=8, frameon=False)
    _style(ax); fig.tight_layout(); fig.savefig(out); plt.close(fig); print("Wrote", out)


if __name__ == "__main__":
    shift_fig(os.path.join(OUT, "od4_shift_figure.png"))
