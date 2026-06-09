from __future__ import annotations

import csv
import os
import textwrap
from dataclasses import dataclass, field
from datetime import datetime
from typing import Any, Dict, List, Tuple

from resilience_measurement import optimize_sequence_sa, run_model_multi
from simulated_annealing import SAConfig
from taskB_setup import generate_taskb_inputs


def _ts() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def _sequence_key(seq: List[Any]) -> str:
    return repr(list(seq))


def _metric_or_nan(run: Dict[str, Any], key: str) -> float:
    return float(run["metric_catalog"].get(key, float("nan")))


def _wrap_title(title: str, width: int = 72) -> str:
    lines: List[str] = []
    for chunk in str(title).split("\n"):
        lines.extend(textwrap.wrap(chunk, width=width) or [""])
    return "\n".join(lines)


def _safe_guardrail_limit(value: float, improvement_factor: float) -> float:
    value = float(value)
    if value <= 0:
        return value
    return value * float(improvement_factor)


def _experiment_description(experiment_id: str, rule_label: str) -> str:
    if experiment_id == "baseline_reference":
        return "Baseline repair sequence used as the reference point."
    if experiment_id.startswith("single_triangle"):
        return "Pure efficiency objective that minimizes resilience-triangle loss."
    if experiment_id.startswith("single_var_restore"):
        return "Equity objective that reduces dispersion in restoration times."
    if experiment_id.startswith("single_gini_restore"):
        return "Equity objective that reduces restoration-time inequality using Gini."
    if experiment_id.startswith("single_p90_restore"):
        return "Equity objective that shortens the time by which 90% of zones recover."
    if experiment_id.startswith("single_maximin_time_avg_cri"):
        return "Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI."
    if experiment_id.startswith("single_p90_access_restore"):
        return "Critical-access objective that shortens the time by which 90% of zones regain access."
    if experiment_id.startswith("weighted_"):
        return f"Weighted trade-off objective. {rule_label.replace('Weighted sum: ', '')}."
    if experiment_id.startswith("guardrail_"):
        return f"Triangle minimization with a fairness guardrail. {rule_label}."
    return rule_label


def _display_rule_type(rule_type: str) -> str:
    labels = {
        "baseline": "Baseline",
        "single": "Single-objective",
        "weighted_sum": "Weighted objective",
        "guardrail": "Guardrail objective",
        "sensitivity": "Sensitivity re-evaluation",
    }
    return labels.get(rule_type, rule_type.replace("_", " ").title())


def _display_experiment_label(experiment_id: str) -> str:
    labels = {
        "baseline_reference": "Baseline sequence",
        "single_triangle": "Efficiency objective",
        "single_var_restore": "Variance equity objective",
        "single_gini_restore": "Gini equity objective",
        "single_maximin_time_avg_cri": "Maximin CRI objective",
        "single_p90_access_restore": "P90 access objective",
        "weighted_gini_restore_l050": "Weighted Gini objective",
        "weighted_maximin_time_avg_cri_l050": "Weighted maximin objective",
        "guardrail_gini_restore": "Gini guardrail objective",
    }
    if experiment_id.startswith("sens::"):
        parts = experiment_id.split("::")
        if len(parts) >= 2:
            return _display_experiment_label(parts[1])
    return labels.get(experiment_id, experiment_id.replace("_", " ").title())


def _pct_change(value: float, baseline: float) -> float:
    baseline = float(baseline)
    value = float(value)
    if baseline == 0:
        return float("nan")
    return (value - baseline) / baseline * 100.0


def _is_dominated(a: Dict[str, Any], b: Dict[str, Any], *, x_key: str, y_key: str) -> bool:
    ax = float(a[x_key]); ay = float(a[y_key])
    bx = float(b[x_key]); by = float(b[y_key])
    return (bx <= ax and by <= ay) and (bx < ax or by < ay)


def _mark_pareto(rows: List[Dict[str, Any]], *, x_key: str, y_key: str) -> None:
    for row in rows:
        dominated = any(_is_dominated(row, other, x_key=x_key, y_key=y_key) for other in rows if other is not row)
        row["is_pareto"] = "1" if not dominated else "0"


def _plot_tradeoff_scatter(rows: List[Dict[str, Any]], out_png: str, *, x_key: str, y_key: str) -> str:
    palette = {
        "baseline": "#1f77b4",
        "single": "#ff7f0e",
        "weighted_sum": "#2ca02c",
        "guardrail": "#d62728",
        "sensitivity": "#7f7f7f",
    }
    markers = {
        "baseline": "s",
        "single": "o",
        "weighted_sum": "^",
        "guardrail": "D",
        "sensitivity": ".",
    }

    try:
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        out_svg = os.path.splitext(out_png)[0] + ".svg"
        _write_tradeoff_svg(rows, out_svg=out_svg, x_key=x_key, y_key=y_key, palette=palette)
        return out_svg

    fig, ax = plt.subplots(figsize=(16, 10), constrained_layout=True)
    groups: Dict[str, List[Dict[str, Any]]] = {}
    for row in rows:
        groups.setdefault(row["rule_type"], []).append(row)

    for rule_type, items in groups.items():
        xs = [float(r[x_key]) for r in items]
        ys = [float(r[y_key]) for r in items]
        ax.scatter(
            xs,
            ys,
            c=palette.get(rule_type, "#333333"),
            marker=markers.get(rule_type, "o"),
            s=130,
            alpha=0.85,
            label=_display_rule_type(rule_type),
        )

    pareto = [r for r in rows if r.get("is_pareto") == "1"]
    pareto_sorted = sorted(pareto, key=lambda r: float(r[x_key]))
    if pareto_sorted:
        ax.plot(
            [float(r[x_key]) for r in pareto_sorted],
            [float(r[y_key]) for r in pareto_sorted],
            color="#111111",
            linewidth=2.4,
            linestyle="--",
            label="Pareto front",
        )

    for row in pareto_sorted[:8]:
        ax.annotate(
            _display_experiment_label(row["experiment_id"]),
            (float(row[x_key]), float(row[y_key])),
            xytext=(8, 6),
            textcoords="offset points",
            fontsize=12,
            alpha=0.9,
        )

    ax.set_xlabel("Triangle area", fontsize=18, labelpad=10)
    ax.set_ylabel("Gini restore", fontsize=18, labelpad=10)
    ax.set_title(_wrap_title("Resilience-Equity Trade-off"), fontsize=22, pad=18)
    ax.tick_params(axis="both", labelsize=14)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3, frameon=True, fontsize=14)
    fig.savefig(out_png, dpi=600, bbox_inches="tight", pad_inches=0.35)
    plt.close(fig)
    return out_png


def _write_tradeoff_svg(
    rows: List[Dict[str, Any]],
    *,
    out_svg: str,
    x_key: str,
    y_key: str,
    palette: Dict[str, str],
) -> None:
    width = 980
    height = 660
    margin_left = 80
    margin_right = 30
    margin_top = 70
    margin_bottom = 120

    xs = [float(r[x_key]) for r in rows]
    ys = [float(r[y_key]) for r in rows]
    xmin = min(xs) if xs else 0.0
    xmax = max(xs) if xs else 1.0
    ymin = min(ys) if ys else 0.0
    ymax = max(ys) if ys else 1.0
    if xmax <= xmin:
        xmax = xmin + 1.0
    if ymax <= ymin:
        ymax = ymin + 1.0
    xmin -= 0.05 * (xmax - xmin)
    xmax += 0.05 * (xmax - xmin)
    ymin = max(0.0, ymin - 0.05 * (ymax - ymin))
    ymax += 0.1 * (ymax - ymin)

    def sx(x: float) -> float:
        return margin_left + (float(x) - xmin) / (xmax - xmin) * (width - margin_left - margin_right)

    def sy(y: float) -> float:
        return height - margin_bottom - (float(y) - ymin) / (ymax - ymin) * (height - margin_top - margin_bottom)

    pareto = sorted([r for r in rows if r.get("is_pareto") == "1"], key=lambda r: float(r[x_key]))
    pareto_line = " ".join(f"{sx(float(r[x_key])):.2f},{sy(float(r[y_key])):.2f}" for r in pareto)

    circles = []
    labels = []
    for row in rows:
        x = sx(float(row[x_key]))
        y = sy(float(row[y_key]))
        color = palette.get(row["rule_type"], "#333333")
        circles.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="5.5" fill="{color}" opacity="0.85"/>')
        if row.get("is_pareto") == "1":
            labels.append(
                f'<text x="{x + 8:.2f}" y="{y - 8:.2f}" font-family="Arial" font-size="11" fill="#222">{row["experiment_id"]}</text>'
            )

    legend_rows = []
    legend_y = height - 126
    for rule_type in ["baseline", "single", "weighted_sum", "guardrail"]:
        color = palette.get(rule_type, "#333333")
        legend_rows.append(
            f'<circle cx="{width-220}" cy="{legend_y}" r="5.5" fill="{color}"/>'
            f'<text x="{width-206}" y="{legend_y + 4}" font-family="Arial" font-size="12">{rule_type}</text>'
        )
        legend_y += 22

    pareto_svg = (
        f'<polyline points="{pareto_line}" fill="none" stroke="#111111" stroke-width="2" stroke-dasharray="6 4"/>'
        if pareto_line
        else ""
    )

    title_lines = _wrap_title("Efficiency-Equity Trade-off", width=60).splitlines() or ["Efficiency-Equity Trade-off"]
    title_svg = "".join(
        f'<tspan x="{width/2:.0f}" dy="{0 if idx == 0 else 20}">{line}</tspan>'
        for idx, line in enumerate(title_lines)
    )

    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
  <rect width="100%" height="100%" fill="white"/>
  <text x="{width/2:.0f}" y="28" text-anchor="middle" font-family="Arial" font-size="18">{title_svg}</text>
  <line x1="{margin_left}" y1="{height-margin_bottom}" x2="{width-margin_right}" y2="{height-margin_bottom}" stroke="#222" stroke-width="2"/>
  <line x1="{margin_left}" y1="{margin_top}" x2="{margin_left}" y2="{height-margin_bottom}" stroke="#222" stroke-width="2"/>
  <text x="{width/2:.0f}" y="{height-18}" text-anchor="middle" font-family="Arial" font-size="14">Triangle area</text>
  <text x="24" y="{height/2:.0f}" text-anchor="middle" font-family="Arial" font-size="14" transform="rotate(-90 24,{height/2:.0f})">Gini restore</text>
  {pareto_svg}
  {''.join(circles)}
  {''.join(labels)}
  <rect x="{width-250}" y="{height-150}" width="210" height="130" fill="white" stroke="#ccc"/>
  {''.join(legend_rows)}
  <line x1="{width-236}" y1="{legend_y + 2}" x2="{width-208}" y2="{legend_y + 2}" stroke="#111111" stroke-width="2" stroke-dasharray="6 4"/>
  <text x="{width-198}" y="{legend_y + 6}" font-family="Arial" font-size="12">Pareto front</text>
</svg>
"""
    with open(out_svg, "w", encoding="utf-8") as f:
        f.write(svg)


@dataclass
class TradeoffConfig:
    base_sequence: List[Any] = field(default_factory=lambda: [(9, 10), 28, 11, 17, 15, 32, (11, 14)])
    result_root: str = "results"
    critical_locations: List[int] = field(default_factory=lambda: [1, 24])
    cri_weight_pairs: List[Tuple[float, float]] = field(default_factory=lambda: [(0.133, 0.867)])
    cri_thresholds: List[float] = field(default_factory=lambda: [0.8, 0.9, 0.95])
    critical_access_thresholds: List[float] = field(default_factory=lambda: [0.8, 0.9, 0.95])
    weighted_lambdas: List[float] = field(default_factory=lambda: [0.25, 0.5, 1.0])
    guardrail_improvement_factor: float = 0.95
    plot_equity_metric: str = "equity:gini_restore"
    sa: SAConfig = field(default_factory=lambda: SAConfig(seed=0, max_iter=80, T0=2.0, alpha=0.98, neighbor="swap"))
    strict: bool = True
    crew_mode: str = "specialized"
    power_crews: int = 1
    road_crews: int = 1
    multifunction_crews: int = 1
    bus_dispatch_mode: str = "link_only"
    bus_location_source: str = "original_bus_location.json"
    bus_to_link_source: str = "new_bus_to_link.json"
    broken_link_factors: Dict[Tuple[int, int], float] = field(default_factory=dict)
    enabled_experiment_ids: List[str] = field(default_factory=list)
    run_sensitivity: bool = True
    save_reference_artifacts: bool = True
    save_baseline: bool = True
    save_best_artifacts: bool = True
    save_best_debug: bool = True


def _build_experiments(reference_metrics: Dict[str, float], cfg: TradeoffConfig) -> List[Dict[str, Any]]:
    specs: List[Dict[str, Any]] = [
        {"experiment_id": "single_triangle", "rule_type": "single", "label": "Single objective: triangle", "objective": "triangle"},
        {"experiment_id": "single_var_restore", "rule_type": "single", "label": "Single objective: var_restore", "objective": "equity:var_restore"},
        {"experiment_id": "single_gini_restore", "rule_type": "single", "label": "Single objective: gini_restore", "objective": "equity:gini_restore"},
        {"experiment_id": "single_p90_restore", "rule_type": "single", "label": "Single objective: p90_restore", "objective": "equity:p90_restore"},
        {
            "experiment_id": "single_maximin_time_avg_cri",
            "rule_type": "single",
            "label": "Single objective: maximize minimum time-avg CRI",
            "objective": "equity:maximin_time_avg_cri_loss",
        },
        {
            "experiment_id": "single_p90_access_restore",
            "rule_type": "single",
            "label": "Single objective: p90_access_restore",
            "objective": "critical_access:p90_access_restore",
        },
    ]

    combo_metrics = [
        ("gini_restore", "equity:gini_restore"),
        ("p90_restore", "equity:p90_restore"),
        ("maximin_time_avg_cri", "equity:maximin_time_avg_cri_loss"),
        ("p90_access_restore", "critical_access:p90_access_restore"),
    ]

    for short_name, metric_key in combo_metrics:
        for lam in cfg.weighted_lambdas:
            specs.append(
                {
                    "experiment_id": f"weighted_{short_name}_l{int(round(lam * 100)):03d}",
                    "rule_type": "weighted_sum",
                    "label": f"Weighted sum: triangle + {lam} * {metric_key}",
                    "objective": "weighted_sum",
                    "objective_weights": {"triangle": 1.0, metric_key: float(lam)},
                    "objective_reference_values": {
                        "triangle": float(reference_metrics["triangle"]),
                        metric_key: float(reference_metrics[metric_key]),
                    },
                }
            )

        specs.append(
            {
                "experiment_id": f"guardrail_{short_name}",
                "rule_type": "guardrail",
                "label": f"Guardrail on {metric_key}",
                "objective": "guardrail",
                "guardrail_primary_metric": "triangle",
                "guardrail_metric": metric_key,
                "guardrail_limit": _safe_guardrail_limit(reference_metrics[metric_key], cfg.guardrail_improvement_factor),
            }
        )

    if not cfg.enabled_experiment_ids:
        return specs

    enabled = set(cfg.enabled_experiment_ids)
    known = {spec["experiment_id"] for spec in specs}
    unknown = sorted(enabled - known)
    if unknown:
        raise ValueError(f"Unknown enabled_experiment_ids: {unknown}")
    return [spec for spec in specs if spec["experiment_id"] in enabled]


def _base_row(run: Dict[str, Any], *, experiment_id: str, rule_type: str, label: str, source_experiment_id: str = "") -> Dict[str, Any]:
    return {
        "row_kind": "optimized",
        "experiment_id": experiment_id,
        "source_experiment_id": source_experiment_id,
        "rule_type": rule_type,
        "rule_label": label,
        "description": _experiment_description(experiment_id, label),
        "objective": run["objective"],
        "objective_value": run["objective_value"],
        "triangle_area": run["triangle_area"],
        "var_restore": _metric_or_nan(run, "equity:var_restore"),
        "gini_restore": _metric_or_nan(run, "equity:gini_restore"),
        "p90_restore": _metric_or_nan(run, "equity:p90_restore"),
        "min_time_avg_cri": _metric_or_nan(run, "equity:min_time_avg_cri"),
        "p90_access_restore": _metric_or_nan(run, "critical_access:p90_access_restore"),
        "share_access_initial": _metric_or_nan(run, "critical_access:share_access_initial"),
        "share_access_final": _metric_or_nan(run, "critical_access:share_access_final"),
        "time_avg_share_access": _metric_or_nan(run, "critical_access:time_avg_share_access"),
        "cri_w_e": run["cri_weights"]["w_e"],
        "cri_w_a": run["cri_weights"]["w_a"],
        "cri_threshold": run["cri_threshold"],
        "critical_access_threshold": run["critical_access_threshold"],
        "is_primary_setting": "1",
        "is_pareto": "0",
        "session_dir": run.get("session_dir", ""),
        "run_dir": run.get("run_dir", ""),
        "sequence": _sequence_key(run["sequence"]),
        "power_sequence": _sequence_key(run.get("power_sequence", [])),
        "road_sequence": _sequence_key(run.get("road_sequence", [])),
    }


def _write_results_discussion_zh(
    *,
    optimized_rows: List[Dict[str, Any]],
    all_rows: List[Dict[str, Any]],
    result_dir: str,
    cfg: TradeoffConfig,
    primary_w_e: float,
    primary_w_a: float,
    primary_cri_threshold: float,
    primary_access_threshold: float,
) -> str:
    report_path = os.path.join(result_dir, "results_discussion_zh.md")
    baseline = next(row for row in optimized_rows if row["experiment_id"] == "baseline_reference")
    baseline_tri = float(baseline["triangle_area"])
    baseline_gini = float(baseline["gini_restore"])
    baseline_p90 = float(baseline["p90_restore"])
    baseline_access = float(baseline["p90_access_restore"])

    sensitivity_rows = [row for row in all_rows if row["row_kind"] == "sensitivity"]
    representative_ids = [
        "baseline_reference",
        "single_triangle",
        "single_gini_restore",
        "single_maximin_time_avg_cri",
        "weighted_gini_restore_l100",
        "guardrail_gini_restore",
    ]

    lines: List[str] = []
    lines.append("# Task B 正式结果表与讨论")
    lines.append("")
    lines.append("## 数据来源")
    lines.append(f"- 主结果文件：`{os.path.basename(os.path.join(result_dir, 'tradeoff_summary.csv'))}`")
    lines.append(f"- 图：`{os.path.basename(os.path.join(result_dir, 'tradeoff_scatter.png'))}`")
    lines.append(f"- 主设定：`w_e={primary_w_e}, w_a={primary_w_a}, cri_threshold={primary_cri_threshold}, critical_access_threshold={primary_access_threshold}`")
    lines.append(f"- 队伍配置：`{cfg.crew_mode}`，power crews=`{cfg.power_crews}`，road crews=`{cfg.road_crews}`")
    lines.append(f"- 调度映射：`{cfg.bus_dispatch_mode}`，使用 `{cfg.bus_to_link_source}`")
    if cfg.enabled_experiment_ids:
        lines.append(f"- 启用实验：`{cfg.enabled_experiment_ids}`")
    lines.append("")
    lines.append("## 表 1 主实验详细结果")
    lines.append("")
    lines.append("| 实验 | 类型 | 描述 | Triangle | ΔTri % | Gini | ΔGini % | Min time-avg CRI | P90 restore | P90 access | Time-avg share access | Pareto |")
    lines.append("|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|")
    for row in optimized_rows:
        lines.append(
            "| {experiment_id} | {rule_type} | {description} | {triangle:.3f} | {dtri:+.1f} | {gini:.3f} | {dgini:+.1f} | {min_final_cri:.3f} | {p90:.3f} | {p90_access:.3f} | {time_avg_share_access:.3f} | {pareto} |".format(
                experiment_id=row["experiment_id"],
                rule_type=row["rule_type"],
                description=row["description"],
                triangle=float(row["triangle_area"]),
                dtri=_pct_change(float(row["triangle_area"]), baseline_tri),
                gini=float(row["gini_restore"]),
                dgini=_pct_change(float(row["gini_restore"]), baseline_gini),
                min_final_cri=float(row["min_time_avg_cri"]),
                p90=float(row["p90_restore"]),
                p90_access=float(row["p90_access_restore"]),
                time_avg_share_access=float(row["time_avg_share_access"]),
                pareto="是" if row.get("is_pareto") == "1" else "否",
            )
        )
    lines.append("")
    lines.append("## 表 2 代表性方案的敏感性范围")
    lines.append("")
    lines.append("| 代表方案 | Gini restore 范围 | Min time-avg CRI 范围 | P90 restore 范围 | P90 access 范围 | Time-avg share access 范围 |")
    lines.append("|---|---:|---:|---:|---:|---:|")
    for exp_id in representative_ids:
        matches = [row for row in sensitivity_rows if row["source_experiment_id"] == exp_id]
        if not matches:
            continue
        gini_vals = [float(row["gini_restore"]) for row in matches]
        min_final_vals = [float(row["min_time_avg_cri"]) for row in matches]
        p90_vals = [float(row["p90_restore"]) for row in matches]
        p90_access_vals = [float(row["p90_access_restore"]) for row in matches]
        time_avg_vals = [float(row["time_avg_share_access"]) for row in matches]
        lines.append(
            f"| {exp_id} | {min(gini_vals):.3f} - {max(gini_vals):.3f} | {min(min_final_vals):.3f} - {max(min_final_vals):.3f} | {min(p90_vals):.3f} - {max(p90_vals):.3f} | {min(p90_access_vals):.3f} - {max(p90_access_vals):.3f} | {min(time_avg_vals):.3f} - {max(time_avg_vals):.3f} |"
        )
    lines.append("")
    lines.append("## 结果讨论")
    lines.append("")
    best_triangle = min(optimized_rows, key=lambda row: float(row["triangle_area"]))
    best_gini = min(optimized_rows, key=lambda row: float(row["gini_restore"]))
    best_maximin = max(optimized_rows, key=lambda row: float(row["min_time_avg_cri"]))
    lines.append(
        f"- 效率最优方案是 `{best_triangle['experiment_id']}`，Triangle 为 `{float(best_triangle['triangle_area']):.3f}`，较基线改善 `{_pct_change(float(best_triangle['triangle_area']), baseline_tri):.1f}%`。"
    )
    lines.append(
        f"- 公平性最优方案是 `{best_gini['experiment_id']}`，Gini restore 为 `{float(best_gini['gini_restore']):.3f}`，较基线改善 `{_pct_change(float(best_gini['gini_restore']), baseline_gini):.1f}%`。"
    )
    lines.append(
        f"- 新增的 maximin 视角下，`{best_maximin['experiment_id']}` 取得最高的最小时间平均 CRI，`min_time_avg_cri={float(best_maximin['min_time_avg_cri']):.3f}`。"
    )
    pareto_ids = [row["experiment_id"] for row in optimized_rows if row.get("is_pareto") == "1"]
    lines.append(f"- `Triangle-Gini` 平面上的 Pareto 解为：`{pareto_ids}`。")
    lines.append(
        f"- `share_access_final` 在当前假设下最终都接近完全恢复，更能区分方案的是 `p90_access_restore` 与 `time_avg_share_access`。基线的 `time_avg_share_access={float(baseline['time_avg_share_access']):.3f}`。"
    )
    lines.append(
        "- 如果后续需要进一步拉开 access 类目标，建议优先调整 critical locations 的布点或扩大受损资产集合；如果要加强搜索，则可以继续增加 SA 迭代次数。"
    )
    lines.append("")
    lines.append("## 建议优先写进正文的结果")
    lines.append("")
    lines.append(f"1. `{best_triangle['experiment_id']}`：代表效率优先。")
    lines.append(f"2. `{best_gini['experiment_id']}`：代表公平性优先。")
    lines.append(f"3. `{best_maximin['experiment_id']}`：代表新增的 maximin equity 规则。")

    with open(report_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    return report_path


def run_tradeoff_study(cfg: TradeoffConfig) -> Dict[str, str]:
    paths = generate_taskb_inputs(destinations=cfg.critical_locations)
    result_dir = os.path.join(cfg.result_root, f"tradeoff_{_ts()}")
    os.makedirs(result_dir, exist_ok=True)

    primary_w_e, primary_w_a = cfg.cri_weight_pairs[0]
    primary_cri_threshold = cfg.cri_thresholds[1] if len(cfg.cri_thresholds) > 1 else cfg.cri_thresholds[0]
    primary_access_threshold = (
        cfg.critical_access_thresholds[1] if len(cfg.critical_access_thresholds) > 1 else cfg.critical_access_thresholds[0]
    )

    reference_run = run_model_multi(
        cfg.base_sequence,
        result_root=result_dir,
        message="Trade-off reference baseline",
        Scenario="tradeoff_reference",
        run_dir=os.path.join(result_dir, "reference"),
        strict=cfg.strict,
        debug=False,
        save_artifacts=cfg.save_reference_artifacts,
        crew_mode=cfg.crew_mode,
        power_crews=cfg.power_crews,
        road_crews=cfg.road_crews,
        multifunction_crews=cfg.multifunction_crews,
        bus_dispatch_mode=cfg.bus_dispatch_mode,
        bus_location_source=cfg.bus_location_source,
        bus_to_link_source=cfg.bus_to_link_source,
        broken_link_factors=cfg.broken_link_factors,
        dest_path=paths["critical_location_path"],
        cri_w_e=primary_w_e,
        cri_w_a=primary_w_a,
        cri_threshold=primary_cri_threshold,
        critical_access_threshold=primary_access_threshold,
        objective="triangle",
    )

    optimized_rows: List[Dict[str, Any]] = [
        _base_row(reference_run, experiment_id="baseline_reference", rule_type="baseline", label="Baseline sequence")
    ]
    sequence_bank: Dict[str, Dict[str, Any]] = {
        _sequence_key(reference_run["sequence"]): {
            "source_experiment_id": "baseline_reference",
            "sequence": list(reference_run["sequence"]),
        }
    }

    experiments = _build_experiments(reference_run["metric_catalog"], cfg)

    for spec in experiments:
        print(f"[tradeoff] running {spec['experiment_id']} ...")
        res = optimize_sequence_sa(
            base_sequence=cfg.base_sequence,
            result_root=result_dir,
            message=spec["label"],
            Scenario=spec["experiment_id"],
            objective=spec["objective"],
            sa=cfg.sa,
            strict=cfg.strict,
            save_baseline=cfg.save_baseline,
            save_best_artifacts=cfg.save_best_artifacts,
            save_best_debug=cfg.save_best_debug,
            crew_mode=cfg.crew_mode,
            power_crews=cfg.power_crews,
            road_crews=cfg.road_crews,
            multifunction_crews=cfg.multifunction_crews,
            bus_dispatch_mode=cfg.bus_dispatch_mode,
            bus_location_source=cfg.bus_location_source,
            bus_to_link_source=cfg.bus_to_link_source,
            broken_link_factors=cfg.broken_link_factors,
            dest_path=paths["critical_location_path"],
            cri_w_e=primary_w_e,
            cri_w_a=primary_w_a,
            cri_threshold=primary_cri_threshold,
            critical_access_threshold=primary_access_threshold,
            objective_weights=spec.get("objective_weights"),
            objective_reference_values=spec.get("objective_reference_values"),
            guardrail_primary_metric=spec.get("guardrail_primary_metric", "triangle"),
            guardrail_metric=spec.get("guardrail_metric", ""),
            guardrail_limit=spec.get("guardrail_limit"),
        )
        best_run = dict(res["best_run"])
        best_run["session_dir"] = res["session_dir"]
        optimized_rows.append(
            _base_row(
                best_run,
                experiment_id=spec["experiment_id"],
                rule_type=spec["rule_type"],
                label=spec["label"],
            )
        )
        sequence_bank.setdefault(
            _sequence_key(best_run["sequence"]),
            {"source_experiment_id": spec["experiment_id"], "sequence": list(best_run["sequence"])},
        )

    _mark_pareto(optimized_rows, x_key="triangle_area", y_key="gini_restore")

    scatter_png = os.path.join(result_dir, "tradeoff_scatter.png")
    scatter_path = _plot_tradeoff_scatter(optimized_rows, scatter_png, x_key="triangle_area", y_key="gini_restore")

    all_rows: List[Dict[str, Any]] = list(optimized_rows)
    if cfg.run_sensitivity:
        sensitivity_tmp = os.path.join(result_dir, "_sensitivity_tmp")
        os.makedirs(sensitivity_tmp, exist_ok=True)

        for seq_key, seq_info in sequence_bank.items():
            seq = list(seq_info["sequence"])
            for (w_e, w_a) in cfg.cri_weight_pairs:
                for cri_threshold in cfg.cri_thresholds:
                    for access_threshold in cfg.critical_access_thresholds:
                        run = run_model_multi(
                            seq,
                            result_root=result_dir,
                            message="Sensitivity re-evaluation",
                            Scenario="sensitivity",
                            run_dir=sensitivity_tmp,
                            strict=cfg.strict,
                            debug=False,
                            save_artifacts=False,
                            crew_mode=cfg.crew_mode,
                            power_crews=cfg.power_crews,
                            road_crews=cfg.road_crews,
                            multifunction_crews=cfg.multifunction_crews,
                            bus_dispatch_mode=cfg.bus_dispatch_mode,
                            bus_location_source=cfg.bus_location_source,
                            bus_to_link_source=cfg.bus_to_link_source,
                            broken_link_factors=cfg.broken_link_factors,
                            dest_path=paths["critical_location_path"],
                            cri_w_e=w_e,
                            cri_w_a=w_a,
                            cri_threshold=cri_threshold,
                            critical_access_threshold=access_threshold,
                            objective="triangle",
                        )
                        all_rows.append(
                            {
                                "row_kind": "sensitivity",
                                "experiment_id": (
                                    f"sens::{seq_info['source_experiment_id']}::"
                                    f"we{w_e:.2f}_wa{w_a:.2f}_ct{cri_threshold:.2f}_at{access_threshold:.2f}"
                                ),
                                "source_experiment_id": seq_info["source_experiment_id"],
                                "rule_type": "sensitivity",
                                "rule_label": "Sensitivity re-evaluation",
                                "description": "Sensitivity re-evaluation of an optimized sequence under alternate weights and thresholds.",
                                "objective": run["objective"],
                                "objective_value": run["objective_value"],
                                "triangle_area": run["triangle_area"],
                                "var_restore": _metric_or_nan(run, "equity:var_restore"),
                                "gini_restore": _metric_or_nan(run, "equity:gini_restore"),
                                "p90_restore": _metric_or_nan(run, "equity:p90_restore"),
                                "min_time_avg_cri": _metric_or_nan(run, "equity:min_time_avg_cri"),
                                "p90_access_restore": _metric_or_nan(run, "critical_access:p90_access_restore"),
                                "share_access_initial": _metric_or_nan(run, "critical_access:share_access_initial"),
                                "share_access_final": _metric_or_nan(run, "critical_access:share_access_final"),
                                "time_avg_share_access": _metric_or_nan(run, "critical_access:time_avg_share_access"),
                                "cri_w_e": w_e,
                                "cri_w_a": w_a,
                                "cri_threshold": cri_threshold,
                                "critical_access_threshold": access_threshold,
                                "is_primary_setting": (
                                    "1"
                                    if (w_e, w_a) == (primary_w_e, primary_w_a)
                                    and cri_threshold == primary_cri_threshold
                                    and access_threshold == primary_access_threshold
                                    else "0"
                                ),
                                "is_pareto": "",
                                "session_dir": "",
                                "run_dir": run["run_dir"],
                                "sequence": seq_key,
                                "power_sequence": _sequence_key(run.get("power_sequence", [])),
                                "road_sequence": _sequence_key(run.get("road_sequence", [])),
                            }
                        )

    csv_path = os.path.join(result_dir, "tradeoff_summary.csv")
    fieldnames = [
        "row_kind",
        "experiment_id",
        "source_experiment_id",
        "rule_type",
        "rule_label",
        "description",
        "objective",
        "objective_value",
        "triangle_area",
        "var_restore",
        "gini_restore",
        "p90_restore",
        "min_time_avg_cri",
        "p90_access_restore",
        "share_access_initial",
        "share_access_final",
        "time_avg_share_access",
        "cri_w_e",
        "cri_w_a",
        "cri_threshold",
        "critical_access_threshold",
        "is_primary_setting",
        "is_pareto",
        "session_dir",
        "run_dir",
        "sequence",
        "power_sequence",
        "road_sequence",
    ]
    with open(csv_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_rows)

    summary_md = os.path.join(result_dir, "tradeoff_summary.md")
    pareto_ids = [row["experiment_id"] for row in optimized_rows if row.get("is_pareto") == "1"]
    with open(summary_md, "w", encoding="utf-8") as f:
        f.write("# Trade-off Summary\n\n")
        f.write(f"- Result directory: `{result_dir}`\n")
        f.write(f"- Critical locations: `{cfg.critical_locations}`\n")
        f.write(f"- Primary CRI weights: `(w_e={primary_w_e}, w_a={primary_w_a})`\n")
        f.write(f"- Primary CRI threshold: `{primary_cri_threshold}`\n")
        f.write(f"- Primary critical-access threshold: `{primary_access_threshold}`\n")
        f.write(f"- Crew mode: `{cfg.crew_mode}` (power=`{cfg.power_crews}`, road=`{cfg.road_crews}`)\n")
        f.write(f"- Bus dispatch mode: `{cfg.bus_dispatch_mode}`\n")
        f.write(f"- Bus-to-link source: `{cfg.bus_to_link_source}`\n")
        if cfg.enabled_experiment_ids:
            f.write(f"- Enabled experiments: `{cfg.enabled_experiment_ids}`\n")
        f.write(f"- Pareto experiments (triangle vs gini_restore): `{pareto_ids}`\n")
        f.write(f"- Combined CSV: `{os.path.basename(csv_path)}`\n")
        f.write(f"- Scatter plot: `{os.path.basename(scatter_path)}`\n")
        f.write("- The combined CSV includes both `power_sequence` and `road_sequence` for specialized-crew runs.\n")
        f.write("\n## Optimized Results\n\n")
        f.write("| Experiment | Type | Description | Triangle | Gini | Min time-avg CRI | P90 access | Time-avg share access | Pareto |\n")
        f.write("|---|---|---|---:|---:|---:|---:|---:|---|\n")
        for row in optimized_rows:
            f.write(
                f"| {row['experiment_id']} | {row['rule_type']} | {row['description']} | "
                f"{float(row['triangle_area']):.3f} | {float(row['gini_restore']):.3f} | "
                f"{float(row['min_time_avg_cri']):.3f} | {float(row['p90_access_restore']):.3f} | "
                f"{float(row['time_avg_share_access']):.3f} | {'Yes' if row.get('is_pareto') == '1' else 'No'} |\n"
            )

    discussion_md = _write_results_discussion_zh(
        optimized_rows=optimized_rows,
        all_rows=all_rows,
        result_dir=result_dir,
        cfg=cfg,
        primary_w_e=primary_w_e,
        primary_w_a=primary_w_a,
        primary_cri_threshold=primary_cri_threshold,
        primary_access_threshold=primary_access_threshold,
    )

    return {
        "result_dir": result_dir,
        "csv_path": csv_path,
        "scatter_png": scatter_path,
        "summary_md": summary_md,
        "discussion_md": discussion_md,
        "critical_location_path": paths["critical_location_path"],
        "destinations": repr(paths["destinations"]),
    }


def main() -> None:
    info = run_tradeoff_study(TradeoffConfig())
    print("Trade-off study complete.")
    print("Result dir:", info["result_dir"])
    print("CSV:", info["csv_path"])
    print("Scatter:", info["scatter_png"])
    print("Summary:", info["summary_md"])
    print("Discussion:", info["discussion_md"])
    print("Critical locations:", info["critical_location_path"], info["destinations"])


if __name__ == "__main__":
    main()
