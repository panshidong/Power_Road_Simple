from __future__ import annotations

"""Audit, analyze, and plot the completed v4 experiment."""

import argparse
import ast
import csv
import glob
import json
import math
import os
import statistics
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

STRATEGIES = ("CEN", "JSH", "IJSH")
LABELS = {
    "CEN": "Structural\ncentrality (CEN)",
    "JSH": "Expected joint\nShapley (JSH)",
    "IJSH": "Interdependency-aware\nShapley (IJSH)",
}
COLORS = {"CEN": "#8A8F98", "JSH": "#4C78A8", "IJSH": "#F58518"}


def mean(values):
    return float(statistics.fmean(values))


def sample_sd(values):
    return float(statistics.stdev(values)) if len(values) > 1 else 0.0


def exact_sign_p(reductions):
    nonzero = [value for value in reductions if abs(value) > 1e-12]
    n = len(nonzero)
    wins = sum(value > 0 for value in nonzero)
    if n == 0:
        return 1.0
    tail = min(wins, n - wins)
    return min(1.0, 2.0 * sum(math.comb(n, k) for k in range(tail + 1)) / (2 ** n))


def trimmed_mean(values, fraction_each_tail=0.025):
    ordered = sorted(float(value) for value in values)
    cut = int(len(ordered) * fraction_each_tail)
    kept = ordered[cut:len(ordered) - cut] if cut else ordered
    return mean(kept)


def bootstrap_interval(values, statistic, rng, draws=20000):
    array = np.asarray(values, dtype=float)
    estimates = np.empty(draws, dtype=float)
    for start in range(0, draws, 1000):
        count = min(1000, draws - start)
        samples = array[rng.integers(0, len(array), size=(count, len(array)))]
        if statistic == "mean":
            estimates[start:start + count] = samples.mean(axis=1)
        elif statistic == "trimmed":
            samples.sort(axis=1)
            cut = int(len(array) * 0.025)
            estimates[start:start + count] = samples[:, cut:len(array) - cut].mean(axis=1) if cut else samples.mean(axis=1)
        else:
            raise ValueError(statistic)
    return [float(value) for value in np.quantile(estimates, (0.025, 0.975))]


def read_worker_rows(root, pattern="rows_*.csv"):
    rows = []
    for path in sorted(glob.glob(os.path.join(root, "worker_*", pattern))):
        rows.extend(csv.DictReader(open(path, newline="", encoding="utf-8")))
    return rows


def index_rows(rows, include_variant=False):
    indexed = defaultdict(dict)
    for row in rows:
        scenario = row["scenario_id"] if not include_variant else f"{row['variant']}::{row['scenario_id']}"
        key = (scenario, row["strategy_id"])
        if row["strategy_id"] in indexed[scenario]:
            raise RuntimeError(f"Duplicate result row {key}")
        indexed[scenario][row["strategy_id"]] = row
    return indexed


def paired_record(indexed, scenarios, candidate, reference, rng, comparison, variant="main"):
    reductions = [
        float(indexed[scenario][reference]["triangle_area"])
        - float(indexed[scenario][candidate]["triangle_area"])
        for scenario in scenarios
    ]
    reference_values = [float(indexed[scenario][reference]["triangle_area"]) for scenario in scenarios]
    ci_low, ci_high = bootstrap_interval(reductions, "mean", rng)
    trim_low, trim_high = bootstrap_interval(reductions, "trimmed", rng)
    return {
        "variant": variant,
        "comparison": comparison,
        "candidate": candidate,
        "reference": reference,
        "n": len(reductions),
        "mean_reduction": mean(reductions),
        "reference_mean": mean(reference_values),
        "mean_reduction_ci95_lower": ci_low,
        "mean_reduction_ci95_upper": ci_high,
        "percent_reduction_of_reference_mean": 100.0 * mean(reductions) / mean(reference_values),
        "median_reduction": float(statistics.median(reductions)),
        "trimmed_mean_reduction": trimmed_mean(reductions),
        "trimmed_ci95_lower": trim_low,
        "trimmed_ci95_upper": trim_high,
        "wins": sum(value > 1e-12 for value in reductions),
        "losses": sum(value < -1e-12 for value in reductions),
        "ties": sum(abs(value) <= 1e-12 for value in reductions),
        "exact_sign_p": exact_sign_p(reductions),
        "reductions": reductions,
    }


def average_ranks(values):
    values = np.asarray(values, dtype=float)
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=float)
    position = 0
    while position < len(values):
        end = position + 1
        while end < len(values) and values[order[end]] == values[order[position]]:
            end += 1
        ranks[order[position:end]] = (position + end - 1) / 2.0
        position = end
    return ranks


def spearman(first, second):
    if len(first) < 2:
        return 1.0
    x = average_ranks(first)
    y = average_ranks(second)
    if float(np.std(x)) == 0.0 or float(np.std(y)) == 0.0:
        return 1.0 if np.array_equal(x, y) else 0.0
    return float(np.corrcoef(x, y)[0, 1])


def convergence_rows(table_payload, manifest):
    checkpoints = sorted(map(int, table_payload["checkpoints"]))
    final_tables = table_payload["checkpoints"][str(checkpoints[-1])]["tables"]
    rows = []
    for checkpoint in checkpoints:
        tables = table_payload["checkpoints"][str(checkpoint)]["tables"]
        for strategy in ("JSH", "IJSH"):
            for asset_type in ("power", "road"):
                correlations = []
                top_one = []
                top_three = []
                for scenario in manifest["scenarios"]:
                    if asset_type == "power":
                        assets = [f"power:{int(bus)}" for bus in scenario["broken_buses"]]
                    else:
                        assets = [f"road:{min(int(link['u']), int(link['v']))}-{max(int(link['u']), int(link['v']))}" for link in scenario["broken_links"]]
                    current = tables[strategy]
                    final = final_tables[strategy]
                    # Early checkpoints may not yet have encountered every
                    # physical component.  Convergence is therefore evaluated
                    # only where the complete realized damage set is covered;
                    # deployment itself uses the final, fully covered table.
                    if any(asset not in current or asset not in final for asset in assets):
                        continue
                    correlations.append(spearman([current[a] for a in assets], [final[a] for a in assets]))
                    current_order = sorted(assets, key=lambda asset: (-float(current[asset]), asset))
                    final_order = sorted(assets, key=lambda asset: (-float(final[asset]), asset))
                    top_one.append(current_order[0] == final_order[0])
                    width = min(3, len(assets))
                    top_three.append(len(set(current_order[:width]) & set(final_order[:width])) / width)
                rows.append(
                    {
                        "checkpoint": checkpoint,
                        "strategy": strategy,
                        "asset_type": asset_type,
                        "spearman_mean": mean(correlations),
                        "top1_agreement": mean(top_one),
                        "top3_overlap": mean(top_three),
                        "n_scenarios_covered": len(correlations),
                    }
                )
    return rows


def ranking_change_rows(indexed, scenarios):
    """Summarize whether a richer score changes either trade's executable order."""
    rows = []
    for candidate, reference, comparison in (
        ("JSH", "CEN", "JSH vs CEN"),
        ("IJSH", "JSH", "IJSH vs JSH"),
    ):
        power_changed = 0
        road_changed = 0
        either_changed = 0
        for scenario in scenarios:
            cand_power = ast.literal_eval(indexed[scenario][candidate]["power_sequence"])
            ref_power = ast.literal_eval(indexed[scenario][reference]["power_sequence"])
            cand_road = ast.literal_eval(indexed[scenario][candidate]["road_sequence"])
            ref_road = ast.literal_eval(indexed[scenario][reference]["road_sequence"])
            p_change = cand_power != ref_power
            r_change = cand_road != ref_road
            power_changed += int(p_change)
            road_changed += int(r_change)
            either_changed += int(p_change or r_change)
        rows.append({
            "comparison": comparison,
            "n": len(scenarios),
            "power_order_changed": power_changed,
            "road_order_changed": road_changed,
            "either_order_changed": either_changed,
        })
    return rows


def save_csv(path, rows, exclude=()):
    serial = [{key: value for key, value in row.items() if key not in exclude} for row in rows]
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(serial[0]))
        writer.writeheader()
        writer.writerows(serial)


def workflow_figure(path):
    fig, ax = plt.subplots(figsize=(12.5, 3.4))
    ax.set_xlim(0, 12.75)
    ax.set_ylim(0, 3.4)
    ax.axis("off")
    boxes = [
        (0.25, "Disruption\ndistribution", "#E9EEF5"),
        (2.35, "Monte Carlo\ntable-construction\nensemble", "#DDEAF7"),
        (4.65, "Scenario joint\nShapley values", "#CFE3F5"),
        (6.85, "Expected\npre-event criticality\ntable", "#BDD7EE"),
        (9.05, "Filter to\nrealized damaged\nassets", "#FCE4D6"),
        (11.05, "Coupled recovery\nand system loss", "#F8CBAD"),
    ]
    for x, label, color in boxes:
        patch = FancyBboxPatch((x, 1.05), 1.55, 1.3, boxstyle="round,pad=0.05,rounding_size=0.08", facecolor=color, edgecolor="#455A64", linewidth=1.2)
        ax.add_patch(patch)
        ax.text(x + 0.775, 1.70, label, ha="center", va="center", fontsize=8.7)
    for (x1, _, _), (x2, _, _) in zip(boxes[:-1], boxes[1:]):
        ax.add_patch(FancyArrowPatch((x1 + 1.55, 1.70), (x2, 1.70), arrowstyle="-|>", mutation_scale=13, color="#455A64", linewidth=1.2))
    ax.text(3.9, 2.75, "Pre-event table construction", ha="center", fontsize=11, weight="bold", color="#315B7D")
    ax.text(10.55, 2.75, "Post-event application and evaluation", ha="center", fontsize=11, weight="bold", color="#A85426")
    ax.plot([0.1, 8.65], [2.55, 2.55], color="#6C8EBF", linewidth=1.5)
    ax.plot([8.85, 12.65], [2.55, 2.55], color="#D79B00", linewidth=1.5)
    fig.savefig(path, dpi=240, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def distribution_figure(path, indexed, scenarios):
    data = [[float(indexed[s][strategy]["triangle_area"]) for s in scenarios] for strategy in STRATEGIES]
    fig, ax = plt.subplots(figsize=(8.4, 5.4), constrained_layout=True)
    parts = ax.violinplot(data, positions=range(3), showextrema=False, widths=0.75)
    for body, strategy in zip(parts["bodies"], STRATEGIES):
        body.set_facecolor(COLORS[strategy]); body.set_edgecolor("#444444"); body.set_alpha(0.42)
    rng = np.random.default_rng(73)
    for index, (strategy, values) in enumerate(zip(STRATEGIES, data)):
        jitter = rng.normal(index, 0.055, size=len(values))
        ax.scatter(jitter, values, s=9, alpha=0.24, color=COLORS[strategy], edgecolors="none")
        ax.boxplot(values, positions=[index], widths=0.18, patch_artist=True, showfliers=False,
                   boxprops=dict(facecolor="white", edgecolor="#333333"), medianprops=dict(color="#111111", linewidth=1.8),
                   whiskerprops=dict(color="#333333"), capprops=dict(color="#333333"))
        ax.scatter([index], [mean(values)], marker="D", s=45, color=COLORS[strategy], edgecolor="black", linewidth=0.5, zorder=5)
    ax.set_yscale("log")
    ax.set_xticks(range(3), [LABELS[strategy] for strategy in STRATEGIES])
    ax.set_ylabel("System resilience-loss area (log scale; lower is better)")
    ax.grid(axis="y", alpha=0.25)
    fig.savefig(path, dpi=240, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def paired_figure(path, paired):
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.8), constrained_layout=True)
    rng = np.random.default_rng(81)
    for ax, record, title in zip(axes, paired, ("H1: JSH versus CEN", "H2: IJSH versus JSH")):
        values = np.asarray(record["reductions"])
        ax.axhline(0, color="#555555", linestyle="--", linewidth=1)
        ax.scatter(rng.normal(0, 0.055, size=len(values)), values, s=12, alpha=0.30, color=COLORS[record["candidate"]], edgecolors="none")
        ax.errorbar([0], [record["mean_reduction"]],
                    yerr=[[record["mean_reduction"] - record["mean_reduction_ci95_lower"]], [record["mean_reduction_ci95_upper"] - record["mean_reduction"]]],
                    fmt="D", color="#111111", markerfacecolor=COLORS[record["candidate"]], capsize=5, linewidth=1.8, markersize=7)
        ax.scatter([0.13], [record["median_reduction"]], marker="s", s=45, facecolor="white", edgecolor="#111111", zorder=5)
        ax.set_yscale("symlog", linthresh=100)
        ax.set_xlim(-0.28, 0.30)
        ax.set_xticks([0, 0.13], ["Mean\n(bootstrap CI)", "Median"])
        ax.set_title(title)
        ax.set_ylabel("Paired loss reduction (positive favors candidate)")
        ax.grid(axis="y", alpha=0.22)
    fig.savefig(path, dpi=240, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def convergence_figure(path, rows):
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.3), constrained_layout=False)
    fig.subplots_adjust(left=0.08, right=0.98, top=0.96, bottom=0.25, wspace=0.20)
    for ax, metric, ylabel in zip(axes, ("spearman_mean", "top3_overlap"), ("Mean Spearman correlation", "Mean top-three overlap")):
        for strategy in ("JSH", "IJSH"):
            for asset_type, linestyle in (("power", "-"), ("road", "--")):
                selected = sorted((row for row in rows if row["strategy"] == strategy and row["asset_type"] == asset_type), key=lambda row: row["checkpoint"])
                ax.plot([row["checkpoint"] for row in selected], [row[metric] for row in selected], marker="o", linestyle=linestyle,
                        color=COLORS[strategy], label=f"{strategy} {asset_type}")
        ax.set_xlabel("Construction-ensemble size")
        ax.set_ylabel(ylabel)
        ax.set_ylim(0.75, 1.01)
        ax.grid(alpha=0.25)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4, frameon=True, bbox_to_anchor=(0.5, 0.02))
    fig.savefig(path, dpi=240, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def robustness_figure(path, records):
    variants = ["main", "small", "large", "light", "clustered"]
    comparisons = ["H1: JSH vs CEN", "H2: IJSH vs JSH"]
    fig, ax = plt.subplots(figsize=(9.2, 5.2), constrained_layout=True)
    offsets = {comparisons[0]: -0.12, comparisons[1]: 0.12}
    colors = {comparisons[0]: COLORS["JSH"], comparisons[1]: COLORS["IJSH"]}
    for comparison in comparisons:
        subset = {record["variant"]: record for record in records if record["comparison"] == comparison}
        for index, variant in enumerate(variants):
            if variant not in subset:
                continue
            record = subset[variant]
            base = record["reference_mean"]
            point = record["percent_reduction_of_reference_mean"]
            low = 100.0 * record["mean_reduction_ci95_lower"] / base
            high = 100.0 * record["mean_reduction_ci95_upper"] / base
            ax.errorbar(point, index + offsets[comparison], xerr=[[point - low], [high - point]], fmt="o", color=colors[comparison], capsize=3, label=comparison if index == 0 else None)
    ax.axvline(0, color="#555555", linestyle="--", linewidth=1)
    ax.set_yticks(range(len(variants)), ["Main distribution", "Small damage", "Large damage", "Light damage", "Clustered damage"])
    ax.invert_yaxis()
    ax.set_xlabel("Mean system-loss reduction (%) with bootstrap 95% CI")
    ax.legend(loc="upper right")
    ax.grid(axis="x", alpha=0.25)
    fig.savefig(path, dpi=240, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=os.path.dirname(os.path.abspath(__file__)))
    ap.add_argument("--expected-main", type=int, default=300)
    args = ap.parse_args()
    output = os.path.join(args.root, "analysis")
    os.makedirs(output, exist_ok=True)
    rng = np.random.default_rng(20260904)

    main_rows = read_worker_rows(os.path.join(args.root, "results", "evaluation"), "rows_scenario_*.csv")
    indexed = index_rows(main_rows)
    scenarios = sorted(scenario for scenario, strategies in indexed.items() if all(strategy in strategies for strategy in STRATEGIES))
    if len(scenarios) != args.expected_main or len(main_rows) != args.expected_main * len(STRATEGIES):
        raise RuntimeError(f"Expected {args.expected_main} complete main scenarios; found {len(scenarios)} / {len(main_rows)} rows")

    descriptive = []
    for strategy in STRATEGIES:
        values = [float(indexed[scenario][strategy]["triangle_area"]) for scenario in scenarios]
        ci_low, ci_high = bootstrap_interval(values, "mean", rng)
        descriptive.append({
            "strategy": strategy,
            "n": len(values),
            "mean": mean(values),
            "bootstrap_ci95_lower": ci_low,
            "bootstrap_ci95_upper": ci_high,
            "median": float(statistics.median(values)),
            "sd": sample_sd(values),
            "p90": float(np.quantile(values, 0.90)),
            "p95": float(np.quantile(values, 0.95)),
            "max": max(values),
        })

    paired_main = [
        paired_record(indexed, scenarios, "JSH", "CEN", rng, "H1: JSH vs CEN"),
        paired_record(indexed, scenarios, "IJSH", "JSH", rng, "H2: IJSH vs JSH"),
    ]
    ordered_p = sorted(enumerate(paired_main), key=lambda item: item[1]["exact_sign_p"])
    running = 0.0
    adjusted = [0.0] * len(paired_main)
    for rank, (index, record) in enumerate(ordered_p):
        value = min(1.0, (len(paired_main) - rank) * record["exact_sign_p"])
        running = max(running, value)
        adjusted[index] = running
    for record, value in zip(paired_main, adjusted):
        record["holm_sign_p"] = value

    robustness_records = []
    robustness_rows = read_worker_rows(os.path.join(args.root, "results", "robustness"), "rows_*.csv")
    if robustness_rows:
        robust_indexed = index_rows(robustness_rows, include_variant=True)
        for variant in ("small", "large", "light", "clustered"):
            selected = sorted(key for key, strategies in robust_indexed.items() if key.startswith(variant + "::") and all(strategy in strategies for strategy in STRATEGIES))
            if selected:
                robustness_records.extend([
                    paired_record(robust_indexed, selected, "JSH", "CEN", rng, "H1: JSH vs CEN", variant),
                    paired_record(robust_indexed, selected, "IJSH", "JSH", rng, "H2: IJSH vs JSH", variant),
                ])
    all_effects = paired_main + robustness_records

    tables = json.load(open(os.path.join(args.root, "results", "tables.json"), encoding="utf-8"))
    manifest_path = sorted(glob.glob(os.path.join(args.root, "results", "evaluation", "worker_*", "evaluation_disasters.json")))[0]
    manifest = json.load(open(manifest_path, encoding="utf-8"))
    convergence = convergence_rows(tables, manifest)
    ranking_changes = ranking_change_rows(indexed, scenarios)

    final_checkpoint = str(max(map(int, tables["checkpoints"])))
    final_counts = tables["checkpoints"][final_checkpoint]["counts"]
    table_coverage = []
    for strategy in STRATEGIES:
        for asset_type, prefix in (("power", "power:"), ("road", "road:")):
            counts = [int(value) for asset, value in final_counts[strategy].items() if asset.startswith(prefix)]
            table_coverage.append({
                "strategy": strategy,
                "asset_type": asset_type,
                "components": len(counts),
                "minimum_occurrences": min(counts),
                "median_occurrences": float(statistics.median(counts)),
                "maximum_occurrences": max(counts),
            })

    save_csv(os.path.join(output, "descriptive.csv"), descriptive)
    save_csv(os.path.join(output, "paired_effects.csv"), all_effects, exclude=("reductions",))
    save_csv(os.path.join(output, "table_convergence.csv"), convergence)
    save_csv(os.path.join(output, "ranking_changes.csv"), ranking_changes)
    save_csv(os.path.join(output, "table_coverage.csv"), table_coverage)
    with open(os.path.join(output, "summary.json"), "w", encoding="utf-8") as handle:
        json.dump({
            "descriptive": descriptive,
            "paired": [{k: v for k, v in record.items() if k != "reductions"} for record in all_effects],
            "convergence": convergence,
            "ranking_changes": ranking_changes,
            "table_coverage": table_coverage,
        }, handle, indent=2)

    workflow_figure(os.path.join(output, "figure_1_workflow.png"))
    distribution_figure(os.path.join(output, "figure_2_loss_distributions.png"), indexed, scenarios)
    paired_figure(os.path.join(output, "figure_3_paired_effects.png"), paired_main)
    convergence_figure(os.path.join(output, "figure_4_table_convergence.png"), convergence)
    if robustness_records:
        robustness_figure(os.path.join(output, "figure_5_robustness.png"), all_effects)

    lines = ["# Task A v4 audited results", "", f"Main evaluation ensemble: {len(scenarios)} complete scenarios.", "", "## Descriptive system loss", ""]
    for record in descriptive:
        lines.append(f"- {record['strategy']}: mean {record['mean']:.1f} (bootstrap 95% CI {record['bootstrap_ci95_lower']:.1f} to {record['bootstrap_ci95_upper']:.1f}), median {record['median']:.1f}, p95 {record['p95']:.1f}.")
    lines += ["", "## Primary paired hypotheses", ""]
    for record in paired_main:
        lines.append(
            f"- {record['comparison']}: mean reduction {record['mean_reduction']:.1f} ({record['percent_reduction_of_reference_mean']:.1f}%; "
            f"bootstrap 95% CI {record['mean_reduction_ci95_lower']:.1f} to {record['mean_reduction_ci95_upper']:.1f}); "
            f"median {record['median_reduction']:.1f}; wins/losses/ties {record['wins']}/{record['losses']}/{record['ties']}; "
            f"exact sign p={record['exact_sign_p']:.6g}, Holm p={record['holm_sign_p']:.6g}."
        )
    lines += ["", "## Executable-order changes", ""]
    for record in ranking_changes:
        lines.append(
            f"- {record['comparison']}: power order changed in {record['power_order_changed']}/{record['n']}, "
            f"road order changed in {record['road_order_changed']}/{record['n']}, and at least one changed in "
            f"{record['either_order_changed']}/{record['n']} scenarios."
        )
    if robustness_records:
        lines += ["", "## Shifted-distribution evidence", ""]
        for record in robustness_records:
            lines.append(f"- {record['variant']} {record['comparison']}: {record['percent_reduction_of_reference_mean']:.1f}% mean reduction (bootstrap CI {record['mean_reduction_ci95_lower']:.1f} to {record['mean_reduction_ci95_upper']:.1f} area units), n={record['n']}.")
    with open(os.path.join(output, "RESULTS_SUMMARY.md"), "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
