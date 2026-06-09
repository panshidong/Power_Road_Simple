from __future__ import annotations

import ast
import csv
import json
import math
import os
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple

from resilience_measurement import optimize_sequence_sa, run_model_multi
from simulated_annealing import SAConfig
from taskB_setup import generate_taskb_inputs


SELECTED_RUN_DIR = Path("/home/workenv/results/batch_tradeoff_20260326_151436/scenario_runs/tradeoff_20260326_180346")
PRIMARY_WEIGHTS = (0.133, 0.867)
WEIGHT_SENSITIVITY = [(0.3, 0.7), (0.7, 0.3)]
THRESHOLD_SENSITIVITY = [0.85, 0.95]
WEIGHTED_LAMBDAS = [0.25, 0.75]
WEIGHTED_METRICS = [
    {
        "short": "gini_restore",
        "metric": "equity:gini_restore",
        "label": "restoration-time Gini",
        "interpretation": "lower inequality in zone-level restoration times",
    },
    {
        "short": "p90_restore",
        "metric": "equity:p90_restore",
        "label": "P90 CRI restoration time",
        "interpretation": "faster recovery for at least 90% of zones",
    },
    {
        "short": "maximin_time_avg_cri",
        "metric": "equity:maximin_time_avg_cri_loss",
        "label": "maximin time-averaged CRI loss",
        "interpretation": "higher performance for the worst-served zone",
    },
    {
        "short": "p90_access_restore",
        "metric": "critical_access:p90_access_restore",
        "label": "P90 critical-access restoration time",
        "interpretation": "faster restoration of critical access for at least 90% of zones",
    },
]


def _ts() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def _read_csv(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _write_csv(path: Path, rows: List[Dict[str, Any]], fieldnames: Iterable[str]) -> None:
    field_list = list(fieldnames)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=field_list)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in field_list})


def _metric_or_nan(run: Dict[str, Any], key: str) -> float:
    return float(run["metric_catalog"].get(key, float("nan")))


def _row_from_run(
    run: Dict[str, Any],
    *,
    row_kind: str,
    experiment_id: str,
    source_experiment_id: str,
    case_type: str,
    case_label: str,
    lambda_value: float | str = "",
    weighted_metric: str = "",
) -> Dict[str, Any]:
    return {
        "row_kind": row_kind,
        "experiment_id": experiment_id,
        "source_experiment_id": source_experiment_id,
        "case_type": case_type,
        "case_label": case_label,
        "lambda": lambda_value,
        "weighted_metric": weighted_metric,
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
        "run_dir": run["run_dir"],
        "sequence": repr(list(run["sequence"])),
        "power_sequence": repr(list(run.get("power_sequence", []))),
        "road_sequence": repr(list(run.get("road_sequence", []))),
    }


def _find_selected_scenario(batch_dir: Path, selected_run_dir: Path) -> Dict[str, Any]:
    manifest_path = batch_dir / "scenario_manifest.csv"
    disasters_path = batch_dir / "random_disasters.json"

    selected_rel = os.path.relpath(selected_run_dir, Path.cwd())
    selected_name = selected_run_dir.name
    selected_manifest = None
    for row in _read_csv(manifest_path):
        if row["study_result_dir"] == selected_rel or row["study_result_dir"].endswith(selected_name):
            selected_manifest = row
            break
    if selected_manifest is None:
        raise RuntimeError(f"Could not find {selected_run_dir} in {manifest_path}")

    raw = json.loads(disasters_path.read_text(encoding="utf-8"))
    scenarios = raw["scenarios"] if isinstance(raw, dict) else raw
    for scenario in scenarios:
        if scenario["scenario_id"] == selected_manifest["scenario_id"]:
            return {"manifest": selected_manifest, "scenario": scenario}
    raise RuntimeError(f"Could not find scenario {selected_manifest['scenario_id']} in {disasters_path}")


def _link_factors_from_scenario(scenario: Dict[str, Any]) -> Dict[Tuple[int, int], float]:
    out: Dict[Tuple[int, int], float] = {}
    for item in scenario["broken_links"]:
        out[(int(item["u"]), int(item["v"]))] = float(item["remaining_capacity_factor"])
    return out


def _float(row: Dict[str, Any], key: str) -> float:
    try:
        return float(row[key])
    except Exception:
        return float("nan")


def _range_rows(rows: List[Dict[str, Any]], *, case_type: str) -> List[Dict[str, Any]]:
    metrics = [
        "triangle_area",
        "gini_restore",
        "p90_restore",
        "min_time_avg_cri",
        "p90_access_restore",
        "time_avg_share_access",
    ]
    grouped: Dict[str, List[Dict[str, Any]]] = {}
    for row in rows:
        if row["case_type"] == case_type:
            grouped.setdefault(row["source_experiment_id"], []).append(row)

    out: List[Dict[str, Any]] = []
    for source_id, items in sorted(grouped.items()):
        row: Dict[str, Any] = {
            "source_experiment_id": source_id,
            "case_type": case_type,
            "n_cases": len(items),
        }
        for metric in metrics:
            vals = [_float(item, metric) for item in items]
            vals = [v for v in vals if not math.isnan(v)]
            if vals:
                row[f"{metric}_min"] = min(vals)
                row[f"{metric}_max"] = max(vals)
                row[f"{metric}_range"] = max(vals) - min(vals)
            else:
                row[f"{metric}_min"] = float("nan")
                row[f"{metric}_max"] = float("nan")
                row[f"{metric}_range"] = float("nan")
        out.append(row)
    return out


def _best_weighted_row(rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    # Favor low triangle and low gini after normalizing by the first candidate row.
    if not rows:
        return {}
    ref = rows[0]
    tri_ref = max(_float(ref, "triangle_area"), 1e-9)
    gini_ref = max(_float(ref, "gini_restore"), 1e-9)
    return min(rows, key=lambda row: _float(row, "triangle_area") / tri_ref + _float(row, "gini_restore") / gini_ref)


def _best_weighted_by_metric(rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    grouped: Dict[str, List[Dict[str, Any]]] = {}
    for row in rows:
        grouped.setdefault(str(row.get("weighted_metric", "")), []).append(row)
    for metric, items in sorted(grouped.items()):
        if not metric or not items:
            continue
        out.append(_best_weighted_row(items))
    return out


def main() -> None:
    os.chdir(Path(__file__).resolve().parent)
    selected_run_dir = SELECTED_RUN_DIR
    batch_dir = selected_run_dir.parents[1]
    selected = _find_selected_scenario(batch_dir, selected_run_dir)
    scenario = selected["scenario"]

    paths = generate_taskb_inputs(destinations=[1, 24])
    broken_link_factors = _link_factors_from_scenario(scenario)
    source_rows = [row for row in _read_csv(selected_run_dir / "tradeoff_summary.csv") if row["row_kind"] == "optimized"]

    out_dir = selected_run_dir / f"extra_sensitivity_{_ts()}"
    eval_dir = out_dir / "evaluations"
    weighted_dir = out_dir / "weighted_runs"
    eval_dir.mkdir(parents=True, exist_ok=True)
    weighted_dir.mkdir(parents=True, exist_ok=True)

    cases: List[Dict[str, Any]] = [
        {
            "case_type": "nominal_recheck",
            "case_label": "Current-code recheck at w=(0.133,0.867), threshold=0.9",
            "w_e": PRIMARY_WEIGHTS[0],
            "w_a": PRIMARY_WEIGHTS[1],
            "cri_threshold": 0.9,
            "access_threshold": 0.9,
        }
    ]
    cases.extend(
        {
            "case_type": "weight_sensitivity",
            "case_label": f"Weight sensitivity w=({w_e},{w_a})",
            "w_e": w_e,
            "w_a": w_a,
            "cri_threshold": 0.9,
            "access_threshold": 0.9,
        }
        for w_e, w_a in WEIGHT_SENSITIVITY
    )
    cases.extend(
        {
            "case_type": "threshold_sensitivity",
            "case_label": f"Threshold sensitivity threshold={threshold}",
            "w_e": PRIMARY_WEIGHTS[0],
            "w_a": PRIMARY_WEIGHTS[1],
            "cri_threshold": threshold,
            "access_threshold": threshold,
        }
        for threshold in THRESHOLD_SENSITIVITY
    )

    sensitivity_rows: List[Dict[str, Any]] = []
    nominal_runs: Dict[str, Dict[str, Any]] = {}
    for source_row in source_rows:
        seq = ast.literal_eval(source_row["sequence"])
        source_id = source_row["experiment_id"]
        for case in cases:
            case_id = f"{source_id}__{case['case_type']}__we{case['w_e']:.3f}_wa{case['w_a']:.3f}_ct{case['cri_threshold']:.2f}_at{case['access_threshold']:.2f}"
            print(f"[extra-sensitivity] evaluating {case_id}")
            run = run_model_multi(
                seq,
                result_root=str(out_dir),
                message=f"Extra sensitivity for {source_id}: {case['case_label']}",
                Scenario=case_id,
                run_dir=str(eval_dir / case_id),
                strict=True,
                debug=False,
                save_artifacts=False,
                crew_mode="specialized",
                power_crews=1,
                road_crews=1,
                multifunction_crews=1,
                bus_dispatch_mode="link_only",
                bus_to_link_source="new_bus_to_link.json",
                dest_path=str(paths["critical_location_path"]),
                broken_link_factors=broken_link_factors,
                cri_w_e=case["w_e"],
                cri_w_a=case["w_a"],
                cri_threshold=case["cri_threshold"],
                critical_access_threshold=case["access_threshold"],
                objective="triangle",
            )
            if case["case_type"] == "nominal_recheck":
                nominal_runs[source_id] = run
            sensitivity_rows.append(
                _row_from_run(
                    run,
                    row_kind="extra_sensitivity",
                    experiment_id=case_id,
                    source_experiment_id=source_id,
                    case_type=case["case_type"],
                    case_label=case["case_label"],
                )
            )

    baseline_ref = nominal_runs.get("baseline_reference") or next(iter(nominal_runs.values()))
    reference_values = {
        "triangle": float(baseline_ref["metric_catalog"]["triangle"]),
    }
    for metric_info in WEIGHTED_METRICS:
        metric_key = metric_info["metric"]
        reference_values[metric_key] = float(baseline_ref["metric_catalog"][metric_key])
    base_sequence = ast.literal_eval(next(row for row in source_rows if row["experiment_id"] == "baseline_reference")["sequence"])

    weighted_rows: List[Dict[str, Any]] = []
    for metric_info in WEIGHTED_METRICS:
        metric_key = metric_info["metric"]
        short = metric_info["short"]
        for lam in WEIGHTED_LAMBDAS:
            exp_id = f"weighted_{short}_l{int(round(lam * 100)):03d}_extra"
            print(f"[extra-sensitivity] optimizing {exp_id}")
            res = optimize_sequence_sa(
                base_sequence=base_sequence,
                result_root=str(weighted_dir),
                message=f"Weighted lambda sensitivity: lambda={lam}, metric={metric_key}",
                Scenario=exp_id,
                objective="weighted_sum",
                sa=SAConfig(seed=0, max_iter=80, T0=1.4, alpha=0.97, neighbor="swap"),
                strict=True,
                save_baseline=False,
                save_best_artifacts=False,
                save_best_debug=False,
                crew_mode="specialized",
                power_crews=1,
                road_crews=1,
                multifunction_crews=1,
                bus_dispatch_mode="link_only",
                bus_to_link_source="new_bus_to_link.json",
                dest_path=str(paths["critical_location_path"]),
                broken_link_factors=broken_link_factors,
                cri_w_e=PRIMARY_WEIGHTS[0],
                cri_w_a=PRIMARY_WEIGHTS[1],
                cri_threshold=0.9,
                critical_access_threshold=0.9,
                objective_weights={"triangle": 1.0, metric_key: float(lam)},
                objective_reference_values={
                    "triangle": reference_values["triangle"],
                    metric_key: reference_values[metric_key],
                },
            )
            best_run = dict(res["best_run"])
            best_run["session_dir"] = res["session_dir"]
            weighted_rows.append(
                _row_from_run(
                    best_run,
                    row_kind="weighted_lambda",
                    experiment_id=exp_id,
                    source_experiment_id="weighted_lambda_sensitivity",
                    case_type="weighted_lambda",
                    case_label=f"Weighted objective with lambda={lam} on {metric_key}",
                    lambda_value=lam,
                    weighted_metric=metric_key,
                )
            )

    all_fieldnames = [
        "row_kind",
        "experiment_id",
        "source_experiment_id",
        "case_type",
        "case_label",
        "lambda",
        "weighted_metric",
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
        "run_dir",
        "sequence",
        "power_sequence",
        "road_sequence",
    ]
    _write_csv(out_dir / "extra_sensitivity_rows.csv", sensitivity_rows, all_fieldnames)
    _write_csv(out_dir / "weighted_lambda_rows.csv", weighted_rows, all_fieldnames)

    range_rows = _range_rows(sensitivity_rows, case_type="weight_sensitivity")
    range_rows += _range_rows(sensitivity_rows, case_type="threshold_sensitivity")
    range_fieldnames = sorted({key for row in range_rows for key in row.keys()})
    _write_csv(out_dir / "extra_sensitivity_ranges.csv", range_rows, range_fieldnames)

    chosen_weighted = _best_weighted_row(weighted_rows)
    chosen_weighted_by_metric = _best_weighted_by_metric(weighted_rows)
    summary_path = out_dir / "extra_sensitivity_summary.md"
    with summary_path.open("w", encoding="utf-8") as f:
        f.write("# Extra Sensitivity Summary\n\n")
        f.write(f"- Selected run: `{selected_run_dir}`\n")
        f.write(f"- Scenario: `{scenario['scenario_id']}` seed `{scenario['seed']}`\n")
        f.write("- Important note: these rows re-evaluate the selected sequences with the current codebase.\n")
        f.write("- Weight sensitivity: `(w_e,w_a)=(0.3,0.7)` and `(0.7,0.3)`, thresholds fixed at `0.9`.\n")
        f.write("- Threshold sensitivity: `cri_threshold` and `critical_access_threshold` jointly set to `0.85` and `0.95`, weights fixed at `(0.133,0.867)`.\n")
        f.write("- Weighted lambda sensitivity tests four metrics: restoration-time Gini, P90 CRI restoration time, maximin time-averaged CRI loss, and P90 critical-access restoration time.\n")
        f.write("- Output files: `extra_sensitivity_rows.csv`, `extra_sensitivity_ranges.csv`, `weighted_lambda_rows.csv`.\n\n")

        f.write("## Weighted Lambda Results\n\n")
        f.write("| Experiment | Weighted metric | Lambda | Triangle | Gini | P90 restore | Min time-avg CRI | P90 access | Time-avg share access |\n")
        f.write("|---|---|---:|---:|---:|---:|---:|---:|---:|\n")
        for row in weighted_rows:
            f.write(
                f"| {row['experiment_id']} | {row['weighted_metric']} | {float(row['lambda']):.2f} | {float(row['triangle_area']):.3f} | "
                f"{float(row['gini_restore']):.3f} | {float(row['p90_restore']):.3f} | "
                f"{float(row['min_time_avg_cri']):.3f} | {float(row['p90_access_restore']):.3f} | "
                f"{float(row['time_avg_share_access']):.3f} |\n"
            )
        if chosen_weighted:
            f.write(
                f"\nChosen balanced weighted case by normalized triangle+gini score: "
                f"`{chosen_weighted['experiment_id']}`.\n"
            )
        if chosen_weighted_by_metric:
            f.write("\nBest balanced case within each weighted metric family:\n\n")
            for row in chosen_weighted_by_metric:
                f.write(f"- `{row['weighted_metric']}`: `{row['experiment_id']}`\n")

        f.write("\n## Largest Sensitivity Ranges\n\n")
        for case_type in ["weight_sensitivity", "threshold_sensitivity"]:
            rows = [row for row in range_rows if row["case_type"] == case_type]
            top_gini = sorted(rows, key=lambda row: float(row["gini_restore_range"]), reverse=True)[:3]
            top_min_cri = sorted(rows, key=lambda row: float(row["min_time_avg_cri_range"]), reverse=True)[:3]
            f.write(f"### {case_type}\n\n")
            f.write("Top Gini ranges:\n\n")
            for row in top_gini:
                f.write(
                    f"- `{row['source_experiment_id']}`: "
                    f"{float(row['gini_restore_min']):.3f} - {float(row['gini_restore_max']):.3f}\n"
                )
            f.write("\nTop min time-avg CRI ranges:\n\n")
            for row in top_min_cri:
                f.write(
                    f"- `{row['source_experiment_id']}`: "
                    f"{float(row['min_time_avg_cri_min']):.3f} - {float(row['min_time_avg_cri_max']):.3f}\n"
                )
            f.write("\n")

    print("Extra sensitivity complete.")
    print("Output dir:", out_dir)
    print("Rows CSV:", out_dir / "extra_sensitivity_rows.csv")
    print("Ranges CSV:", out_dir / "extra_sensitivity_ranges.csv")
    print("Weighted CSV:", out_dir / "weighted_lambda_rows.csv")
    print("Summary:", summary_path)


if __name__ == "__main__":
    main()
