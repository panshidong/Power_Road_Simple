from __future__ import annotations

"""
Old-vs-new comparison after the power-performance fix.

Reads the pre-fix results from backup_prefix_powerbug_20260903/results/ and the
post-fix results from results/, recomputes every number the paper draft reports
(marginal means, paired differences, significance flags, coverage, sensitivity),
and writes a markdown change list to POSTFIX_CHANGE_LIST.md. Read-only.
"""

import csv
import math
from collections import defaultdict
from typing import Dict, List

OLD = "backup_prefix_powerbug_20260903/results"
NEW = "results"
ORDER = ["S0_centrality", "S1_separate_shapley", "S2_integrated_shapley", "S3_od_representation"]
SHORT = {"S0_centrality": "S0", "S1_separate_shapley": "S1", "S2_integrated_shapley": "S2", "S3_od_representation": "S3"}
BASELINES = ORDER[:3]
METRICS = ["triangle_area", "weighted_triangle_area", "gini_restore", "min_time_avg_cri", "p90_access_restore"]
LOWER_BETTER = {"triangle_area": True, "weighted_triangle_area": True, "gini_restore": True,
                "min_time_avg_cri": False, "p90_access_restore": True}


def load_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def mean(xs):
    return sum(xs) / len(xs) if xs else float("nan")


def ci95(xs):
    if len(xs) <= 1:
        return 0.0
    mu = mean(xs)
    sd = math.sqrt(sum((x - mu) ** 2 for x in xs) / (len(xs) - 1))
    return 1.96 * sd / math.sqrt(len(xs))


class Study:
    def __init__(self, root: str):
        rows = load_csv(f"{root}/task_c_main_100/scenario_strategy_rows.csv")
        self.by: Dict[str, Dict[str, Dict[str, float]]] = defaultdict(dict)
        self.seq: Dict[str, Dict[str, str]] = defaultdict(dict)
        for r in rows:
            self.by[r["scenario_id"]][r["strategy_id"]] = {m: float(r[m]) for m in METRICS}
            self.seq[r["scenario_id"]][r["strategy_id"]] = (r["power_sequence"], r["road_sequence"])
        self.scen = sorted(s for s, d in self.by.items() if all(k in d for k in ORDER))
        rank = load_csv(f"{root}/task_c_main_100/criticality_rankings.csv")
        s3_road = [r for r in rank if r["strategy_id"] == "S3_od_representation" and r["asset_type"] == "road"]
        self.links_total = len(s3_road)
        self.links_nonzero = sum(1 for r in s3_road if float(r["score"]) > 0)
        self.identical = sum(1 for s in self.scen if self.seq[s]["S0_centrality"][1] == self.seq[s]["S3_od_representation"][1])
        try:
            self.sens = sorted(load_csv(f"{root}/task_c_sensitivity/taskC_od_sensitivity_aggregate.csv"),
                               key=lambda r: int(r["sensitivity_tag"][1:]))
        except FileNotFoundError:
            self.sens = []

    def marg(self, sid, m):
        xs = [self.by[s][sid][m] for s in self.scen]
        return mean(xs), ci95(xs)

    def paired(self, m, b):
        d = [self.by[s]["S3_od_representation"][m] - self.by[s][b][m] for s in self.scen]
        lb = LOWER_BETTER[m]
        better = sum(1 for x in d if x != 0 and (x < 0) == lb)
        worse = sum(1 for x in d if x != 0 and (x > 0) == lb)
        tied = sum(1 for x in d if x == 0)
        mu, c = mean(d), ci95(d)
        base = mean([self.by[s][b][m] for s in self.scen])
        return dict(mean=mu, ci=c, sig=abs(mu) > c, better=better, worse=worse, tied=tied, pct=mu / base * 100)


def fmt(m, c, nd=1):
    return f"{m:.{nd}f} ± {c:.{nd}f}"


def main():
    old, new = Study(OLD), Study(NEW)
    out: List[str] = []
    w = out.append
    w("# Post-fix change list (pre-fix → post-fix)\n")
    w(f"Pre-fix scenarios: {len(old.scen)}; post-fix scenarios: {len(new.scen)} (same seed0 → identical damage sets).\n")

    # sequences changed?
    changed = {sid: sum(1 for s in new.scen if old.seq[s][sid] != new.seq[s][sid]) for sid in ORDER}
    w("## 0. Did the repair sequences themselves change?\n")
    w("| strategy | scenarios where (power seq, road seq) differs old vs new |")
    w("|---|---|")
    for sid in ORDER:
        w(f"| {SHORT[sid]} | {changed[sid]} / {len(new.scen)} |")
    w("\nExpected: S0 and S3 unchanged (their scores do not use the fixed function); S1/S2 may change "
      "(their Shapley value functions use power performance).\n")

    # marginal means
    w("## 1. Marginal means (mean ± 95% CI)\n")
    for m in METRICS:
        nd = 3 if m in ("gini_restore", "min_time_avg_cri") else 1
        w(f"### {m}\n")
        w("| strategy | pre-fix | post-fix |")
        w("|---|---|---|")
        for sid in ORDER:
            w(f"| {SHORT[sid]} | {fmt(*old.marg(sid, m), nd)} | {fmt(*new.marg(sid, m), nd)} |")
        w("")
        rank_old = sorted(ORDER, key=lambda s: old.marg(s, m)[0], reverse=not LOWER_BETTER[m])
        rank_new = sorted(ORDER, key=lambda s: new.marg(s, m)[0], reverse=not LOWER_BETTER[m])
        w(f"Best-to-worst order: pre-fix {' < '.join(SHORT[s] for s in rank_old)} → post-fix "
          f"{' < '.join(SHORT[s] for s in rank_new)}\n")

    # paired
    w("## 2. Paired S3 − baseline (mean ± 95% CI, % change, significant?, S3 better/worse/tied)\n")
    for m in METRICS:
        nd = 3 if m in ("gini_restore", "min_time_avg_cri") else 1
        w(f"### {m}\n")
        w("| vs | pre-fix | post-fix | significance change |")
        w("|---|---|---|---|")
        for b in BASELINES:
            po, pn = old.paired(m, b), new.paired(m, b)
            flag = "" if po["sig"] == pn["sig"] else (f"**{'n.s.→sig' if pn['sig'] else 'sig→n.s.'}**")
            w(f"| {SHORT[b]} | {pn and fmt(po['mean'], po['ci'], nd)} ({po['pct']:+.1f}%, "
              f"{'sig' if po['sig'] else 'n.s.'}, {po['better']}/{po['worse']}/{po['tied']}) | "
              f"{fmt(pn['mean'], pn['ci'], nd)} ({pn['pct']:+.1f}%, {'sig' if pn['sig'] else 'n.s.'}, "
              f"{pn['better']}/{pn['worse']}/{pn['tied']}) | {flag} |")
        w("")

    # coverage
    w("## 3. Static-score coverage and degeneracy (should be unchanged)\n")
    w(f"- damaged road links with nonzero S3 score: pre {old.links_nonzero}/{old.links_total} "
      f"({old.links_nonzero/old.links_total*100:.1f}%) → post {new.links_nonzero}/{new.links_total} "
      f"({new.links_nonzero/new.links_total*100:.1f}%)")
    w(f"- S3 road order identical to S0: pre {old.identical}/{len(old.scen)} → post {new.identical}/{len(new.scen)}\n")

    # sensitivity
    w("## 4. k-sensitivity (15-scenario batch, S3 only)\n")
    if old.sens and new.sens:
        w("| k | pre-fix triangle | post-fix triangle | pre-fix weighted | post-fix weighted |")
        w("|---|---|---|---|---|")
        for ro, rn in zip(old.sens, new.sens):
            w(f"| {rn['sensitivity_tag'][1:]} | {float(ro['triangle_area_mean']):.1f} | {float(rn['triangle_area_mean']):.1f} | "
              f"{float(ro['weighted_triangle_area_mean']):.1f} | {float(rn['weighted_triangle_area_mean']):.1f} |")
        tn = [float(r["triangle_area_mean"]) for r in new.sens]
        w(f"\npost-fix spread across k: {(max(tn)-min(tn))/mean(tn)*100:.1f}% of mean\n")
    else:
        w("(post-fix sensitivity not yet available)\n")

    text = "\n".join(out)
    with open("POSTFIX_CHANGE_LIST.md", "w", encoding="utf-8") as f:
        f.write(text)
    print(text)


if __name__ == "__main__":
    main()
