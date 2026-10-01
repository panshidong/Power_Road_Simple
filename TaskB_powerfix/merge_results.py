from __future__ import annotations

"""Merge per-scenario rows from all workers into the same four artifacts the
original batch_tradeoff_runner produced (scenario_manifest.csv,
scenario_tradeoff_rows.csv, aggregate_tradeoff_summary.csv,
aggregate_tradeoff_errorbars.png, batch_summary.md) under results/merged/."""

import csv
import glob
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "template"))

from batch_tradeoff_runner import _aggregate_rows, _plot_errorbar_tradeoff  # noqa: E402
from runner_fast100 import build_fast100_config  # noqa: E402


def main() -> None:
    cfg = build_fast100_config()
    files = sorted(glob.glob(os.path.join(HERE, "results", "worker_*", "rows_scenario_*.csv")))
    rows = []
    for fp in files:
        with open(fp, newline="", encoding="utf-8") as f:
            rows.extend(csv.DictReader(f))
    rows.sort(key=lambda r: (r["scenario_id"], r["experiment_id"]))
    scenario_ids = sorted({r["scenario_id"] for r in rows})
    out_dir = os.path.join(HERE, "results", "merged")
    os.makedirs(out_dir, exist_ok=True)

    with open(os.path.join(out_dir, "scenario_tradeoff_rows.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    seen = {}
    for r in rows:
        seen.setdefault(r["scenario_id"], {
            "scenario_id": r["scenario_id"], "scenario_seed": r["scenario_seed"],
            "broken_bus_count": r["broken_bus_count"], "broken_link_count": r["broken_link_count"],
            "study_result_dir": r["study_result_dir"]})
    with open(os.path.join(out_dir, "scenario_manifest.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(next(iter(seen.values())).keys()))
        w.writeheader()
        w.writerows(seen[s] for s in scenario_ids)

    agg = _aggregate_rows(rows, x_metric=cfg.aggregate_x_metric, y_metric=cfg.aggregate_y_metric)
    with open(os.path.join(out_dir, "aggregate_tradeoff_summary.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(agg[0].keys()))
        w.writeheader()
        w.writerows(agg)
    png = _plot_errorbar_tradeoff(agg, out_png=os.path.join(out_dir, "aggregate_tradeoff_errorbars.png"),
                                  x_metric=cfg.aggregate_x_metric, y_metric=cfg.aggregate_y_metric)

    with open(os.path.join(out_dir, "batch_summary.md"), "w", encoding="utf-8") as f:
        f.write("# Batch Trade-off Summary (power-performance rule FIXED; parallel rerun)\n\n")
        f.write(f"- Scenario count: `{len(scenario_ids)}` of `{cfg.n_scenarios}` planned\n")
        f.write(f"- Scenario seed start: `{cfg.seed0}`\n")
        f.write(f"- Experiments: `{cfg.study_cfg.enabled_experiment_ids}` + baseline_reference\n")
        f.write(f"- SA: `{cfg.study_cfg.sa}`\n")
        f.write(f"- Crew mode: `{cfg.study_cfg.crew_mode}` (power=`{cfg.study_cfg.power_crews}`, road=`{cfg.study_cfg.road_crews}`)\n")
        f.write(f"- CRI weights: `{cfg.study_cfg.cri_weight_pairs}`\n\n")
        f.write("## Aggregate Means\n\n| Experiment | Type | Triangle mean ± CI | Gini mean ± CI | Min time-avg CRI ± CI | P90 access mean ± CI |\n|---|---|---:|---:|---:|---:|\n")
        for r in agg:
            f.write(f"| {r['experiment_id']} | {r['rule_type']} | {float(r['triangle_area_mean']):.3f} ± {float(r['triangle_area_ci95']):.3f} | "
                    f"{float(r['gini_restore_mean']):.3f} ± {float(r['gini_restore_ci95']):.3f} | "
                    f"{float(r['min_time_avg_cri_mean']):.3f} ± {float(r['min_time_avg_cri_ci95']):.3f} | "
                    f"{float(r['p90_access_restore_mean']):.3f} ± {float(r['p90_access_restore_ci95']):.3f} |\n")
    print(f"merged {len(rows)} rows from {len(scenario_ids)} scenarios -> {out_dir}")


if __name__ == "__main__":
    main()
