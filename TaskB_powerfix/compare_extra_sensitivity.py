from __future__ import annotations

"""Compare the post-fix representative-case sensitivity (scenario_084) with the
paper's pre-fix run (Tables 5-7 of the equity paper). Writes
results/merged/CHANGE_VS_EXTRA_SENSITIVITY.md."""

import csv
import glob
import os

OLD = "/home/workenv/results/batch_tradeoff_20260326_151436/scenario_runs/tradeoff_20260326_180346/extra_sensitivity_20260604_123941"
HERE = os.path.dirname(os.path.abspath(__file__))


def load(p):
    with open(p, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def main():
    new_dirs = sorted(glob.glob(os.path.join(HERE, "results/worker_11/scenario_runs/tradeoff_20260903_153305/extra_sensitivity_*")))
    NEW = new_dirs[-1]
    L = []
    w = L.append
    w("# Representative case (scenario_084, seed 20260408): pre-fix vs post-fix sensitivity\n")
    w(f"- pre-fix: `{OLD}`\n- post-fix: `{NEW}`\n")
    w("Note: the post-fix sequences are the *re-optimized* sequences from the post-fix batch, not the paper's old sequences re-evaluated; both the sequences and the metric changed.\n")

    for label, fname, keyf in [("Table 5/6-style ranges (weight & threshold sensitivity)", "extra_sensitivity_ranges.csv",
                                lambda r: (r["case_type"], r["source_experiment_id"]))]:
        o = {keyf(r): r for r in load(os.path.join(OLD, fname))}
        n = {keyf(r): r for r in load(os.path.join(NEW, fname))}
        for case_type in ["weight_sensitivity", "threshold_sensitivity"]:
            w(f"## {case_type}\n")
            w("| strategy | Gini range old → new | Min avg CRI old → new | P90 restore old → new |")
            w("|---|---|---|---|")
            for k in sorted(set(o) | set(n)):
                if k[0] != case_type:
                    continue
                ro, rn = o.get(k), n.get(k)
                def rng(r, m):
                    return "—" if r is None else f"{float(r[m+'_min']):.3f}–{float(r[m+'_max']):.3f}"
                def rng1(r, m):
                    return "—" if r is None else f"{float(r[m+'_min']):.1f}–{float(r[m+'_max']):.1f}"
                w(f"| {k[1]} | {rng(ro,'gini_restore')} → {rng(rn,'gini_restore')} | {rng(ro,'min_time_avg_cri')} → {rng(rn,'min_time_avg_cri')} | {rng1(ro,'p90_restore')} → {rng1(rn,'p90_restore')} |")
            w("")

    w("## Table 7-style weighted-lambda results\n")
    o = {r["experiment_id"]: r for r in load(os.path.join(OLD, "weighted_lambda_rows.csv"))}
    n = {r["experiment_id"]: r for r in load(os.path.join(NEW, "weighted_lambda_rows.csv"))}
    w("| experiment | Triangle old → new | Gini old → new | Min avg CRI old → new | P90 access old → new |")
    w("|---|---|---|---|---|")
    for k in sorted(set(o) | set(n)):
        ro, rn = o.get(k), n.get(k)
        f = lambda r, m, nd: "—" if r is None else f"{float(r[m]):.{nd}f}"
        w(f"| {k} | {f(ro,'triangle_area',1)} → {f(rn,'triangle_area',1)} | {f(ro,'gini_restore',3)} → {f(rn,'gini_restore',3)} | "
          f"{f(ro,'min_time_avg_cri',3)} → {f(rn,'min_time_avg_cri',3)} | {f(ro,'p90_access_restore',1)} → {f(rn,'p90_access_restore',1)} |")

    w("\n## Nominal recheck of each strategy's sequence (w=(0.133,0.867), threshold 0.9)\n")
    o = {r["source_experiment_id"]: r for r in load(os.path.join(OLD, "extra_sensitivity_rows.csv")) if r["case_type"] == "nominal_recheck"}
    n = {r["source_experiment_id"]: r for r in load(os.path.join(NEW, "extra_sensitivity_rows.csv")) if r["case_type"] == "nominal_recheck"}
    w("| strategy | Triangle old → new | Gini old → new | Min avg CRI old → new |")
    w("|---|---|---|---|")
    for k in sorted(set(o) | set(n)):
        ro, rn = o.get(k), n.get(k)
        f = lambda r, m, nd: "—" if r is None else f"{float(r[m]):.{nd}f}"
        w(f"| {k} | {f(ro,'triangle_area',1)} → {f(rn,'triangle_area',1)} | {f(ro,'gini_restore',3)} → {f(rn,'gini_restore',3)} | {f(ro,'min_time_avg_cri',3)} → {f(rn,'min_time_avg_cri',3)} |")

    out = os.path.join(HERE, "results/merged/CHANGE_VS_EXTRA_SENSITIVITY.md")
    open(out, "w", encoding="utf-8").write("\n".join(L))
    print("\n".join(L))


if __name__ == "__main__":
    main()
