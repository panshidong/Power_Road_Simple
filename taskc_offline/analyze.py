from __future__ import annotations

"""Merge per-scenario rows of the Task C offline-table experiment and compute the
paired statistics used by the OD paper (bootstrap 95% intervals, trimmed means,
medians, wins/losses/ties, exact sign tests)."""

import csv
import glob
import json
import math
import os
import random
import statistics as st
from collections import defaultdict
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
B = 20000
RNG_SEED = 20260907


def load(stage):
    rows = []
    for fp in sorted(glob.glob(os.path.join(HERE, "results", stage, "worker_*", "rows_*.csv"))):
        rows.extend(csv.DictReader(open(fp, newline="", encoding="utf-8")))
    return rows


def by_scenario(rows):
    by = defaultdict(dict)
    for r in rows:
        by[r["scenario_id"]][r["strategy_id"]] = r
    return by


def mean(x):
    return sum(x) / len(x)


def trimmed(x, frac_each=0.025):
    x = sorted(x); k = int(len(x) * frac_each)
    return x[k:len(x) - k] if k else x


def boot_ci(values, stat, rng):
    n = len(values); out = []
    for _ in range(B):
        s = [values[rng.randrange(n)] for _ in range(n)]
        out.append(stat(s))
    out.sort()
    return out[int(0.025 * B)], out[int(0.975 * B) - 1]


def sign_p(d):
    d = [x for x in d if abs(x) > 1e-12]; n = len(d); k = sum(1 for x in d if x < 0)
    if n == 0:
        return 1.0
    return min(1.0, 2 * sum(comb(n, i) for i in range(min(k, n - k) + 1)) / 2 ** n)


def descriptive(by, strategies, metric="triangle_area"):
    rng = random.Random(RNG_SEED); out = []
    for s in strategies:
        v = [float(by[sc][s][metric]) for sc in by]
        lo, hi = boot_ci(v, mean, rng)
        out.append(dict(strategy=s, n=len(v), mean=mean(v), ci_lo=lo, ci_hi=hi, median=st.median(v), sd=st.pstdev(v) * math.sqrt(len(v) / (len(v) - 1)),
                        p90=sorted(v)[int(0.9 * len(v))], p95=sorted(v)[int(0.95 * len(v))], max=max(v)))
    return out


def paired(by, cand, ref, metric="triangle_area", lower_better=True):
    rng = random.Random(RNG_SEED)
    scen = sorted(sc for sc in by if cand in by[sc] and ref in by[sc])
    a = [float(by[sc][ref][metric]) for sc in scen]; b = [float(by[sc][cand][metric]) for sc in scen]
    d = [x - y for x, y in zip(a, b)] if lower_better else [y - x for x, y in zip(a, b)]  # positive favours candidate
    idx = list(range(len(d)))
    def pct(sample_idx):
        ra = mean([a[i] for i in sample_idx]); rb = mean([b[i] for i in sample_idx])
        return (ra - rb) / ra * 100 if lower_better else (rb - ra) / ra * 100
    lo, hi = boot_ci(d, mean, rng)
    plo, phi = boot_ci(idx, pct, rng)
    tlo, thi = boot_ci(d, lambda s: mean(trimmed(s)), rng)
    return dict(candidate=cand, reference=ref, metric=metric, n=len(d), reference_mean=mean(a), candidate_mean=mean(b),
                mean_reduction=mean(d), ci_lo=lo, ci_hi=hi, pct=(mean(d) / mean(a) * 100), pct_lo=plo, pct_hi=phi,
                trimmed=mean(trimmed(d)), trimmed_lo=tlo, trimmed_hi=thi, median=st.median(d),
                wins=sum(1 for x in d if x > 1e-12), losses=sum(1 for x in d if x < -1e-12), ties=sum(1 for x in d if abs(x) <= 1e-12),
                sign_p=sign_p(d))


def fmt_pair(p):
    f = ",.1f" if abs(p["reference_mean"]) > 10 else ".4f"
    return (f"| {p['candidate']} vs {p['reference']} | {p['mean_reduction']:{f}} ({p['ci_lo']:{f}}, {p['ci_hi']:{f}}) | "
            f"{p['pct']:.1f}% ({p['pct_lo']:.1f}, {p['pct_hi']:.1f}) | {p['trimmed']:{f}} ({p['trimmed_lo']:{f}}, {p['trimmed_hi']:{f}}) | "
            f"{p['median']:{f}} | {p['wins']}/{p['losses']}/{p['ties']} | {p['sign_p']:.2g} |")


def main():
    import argparse
    ap = argparse.ArgumentParser(); ap.add_argument("--variant", default="", help="'' (depot-only O-D) or 'allpairs'")
    args = ap.parse_args(); sfx = f"_{args.variant}" if args.variant else ""
    main_rows = load("main_300"); by = by_scenario(main_rows)
    S = ["CEN", "JSH", "IJSH", "OD_CEN", "OD_JSH", "OD_IJSH"]
    if args.variant:
        # baselines from the depot-only run (unchanged), O-D rows from the variant run
        od = by_scenario(load("main_300" + sfx))
        for sc in by:
            for k in list(by[sc]):
                if k.startswith("OD_"):
                    del by[sc][k]
            by[sc].update(od.get(sc, {}))
        S = ["CEN", "JSH", "IJSH", "OD_CEN"] + (["OD_IJSH"] if all("OD_IJSH" in by[sc] for sc in by if "OD_CEN" in by[sc]) else [])
        main_rows = [r for sc in by for r in by[sc].values()]
    complete = {sc for sc in by if all(s in by[sc] for s in S)}
    by = {sc: by[sc] for sc in complete}
    decomp = by_scenario(load("decomp_300"))
    for sc in by:
        if sc in decomp:
            by[sc].update(decomp[sc])
    os.makedirs(os.path.join(HERE, "results", "merged"), exist_ok=True)
    allrows = [r for sc in sorted(by) for r in by[sc].values()]
    with open(os.path.join(HERE, "results", "merged", f"rows{sfx}.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(allrows[0])); w.writeheader(); w.writerows(allrows)

    # reproduction check against Codex's evaluation rows
    ref = by_scenario([r for fp in glob.glob("/home/workenv/TaskA_v4/results/evaluation/worker_*/rows_*.csv") for r in csv.DictReader(open(fp))])
    maxdev = max(abs(float(by[sc][s]["triangle_area"]) - float(ref[sc][s]["triangle_area"])) for sc in by for s in ("CEN", "JSH", "IJSH") if sc in ref)
    fallbacks = sum(int(r["tapb_fallbacks"]) for r in allrows)

    out = {"n": len(by), "max_abs_dev_vs_taskA_rows": maxdev, "tapb_fallbacks": fallbacks}
    out["descriptive"] = descriptive(by, S)
    comps = [("OD_CEN", "CEN"), ("OD_JSH", "JSH"), ("OD_IJSH", "IJSH"), ("OD_CEN", "JSH"), ("OD_CEN", "IJSH"), ("OD_IJSH", "CEN"), ("OD_IJSH", "JSH"), ("IJSH", "CEN"), ("JSH", "CEN"), ("IJSH", "JSH")]
    comps = [(c, r) for c, r in comps if all(c in by[sc] and r in by[sc] for sc in by)]
    out["paired"] = [paired(by, c, r) for c, r in comps]
    supp = [(c, r) for c, r in [("OD_CEN", "CEN"), ("OD_CEN", "JSH"), ("OD_CEN", "IJSH"), ("OD_JSH", "JSH"), ("OD_IJSH", "IJSH")] if all(c in by[sc] for sc in by)]
    out["paired_weighted"] = [paired(by, c, r, "weighted_triangle_area") for c, r in supp]
    out["paired_gini"] = [paired(by, c, r, "gini_restore") for c, r in supp]
    out["paired_mincri"] = [paired(by, c, r, "min_time_avg_cri", lower_better=False) for c, r in supp]
    out["descriptive_weighted"] = descriptive(by, S, "weighted_triangle_area")
    out["descriptive_gini"] = descriptive(by, S, "gini_restore")
    out["descriptive_mincri"] = descriptive(by, S, "min_time_avg_cri")
    out["descriptive_p90"] = descriptive(by, S, "p90_access_restore")
    # road-order agreement and static-score sparsity
    same_as = {s: sum(1 for sc in by if by[sc]["OD_CEN"]["road_sequence"] == by[sc][s]["road_sequence"]) for s in ("CEN", "JSH", "IJSH")}
    out["od_road_order_identical_to"] = same_as
    nz = sum(int(by[sc]["OD_CEN"]["n_road_nonzero_od"]) for sc in by); tot = sum(int(by[sc]["OD_CEN"]["n_road"]) for sc in by)
    out["damaged_links_nonzero_share"] = nz / tot; out["scenarios_with_nonzero"] = sum(1 for sc in by if int(by[sc]["OD_CEN"]["n_road_nonzero_od"]) > 0)
    out["mean_damage"] = dict(power=mean([int(by[sc]["CEN"]["n_power"]) for sc in by]), road=mean([int(by[sc]["CEN"]["n_road"]) for sc in by]))
    if decomp:
        dsc = [sc for sc in by if "MIX_JSHroad_IJSHpower" in by[sc] and "MIX_IJSHroad_JSHpower" in by[sc]]
        sub = {sc: by[sc] for sc in dsc}
        out["decomp_n"] = len(dsc)
        out["decomp"] = [paired(sub, "MIX_IJSHroad_JSHpower", "JSH"), paired(sub, "MIX_JSHroad_IJSHpower", "JSH"),
                         paired(sub, "IJSH", "MIX_IJSHroad_JSHpower"), paired(sub, "IJSH", "MIX_JSHroad_IJSHpower")]
    # sensitivity
    sens = {}
    for k in (1, 2, 3, 5):
        rows = load(f"sens_k{k}{sfx}") if k != 3 else main_rows
        if not rows:
            continue
        b = by_scenario(rows)
        first100 = sorted(b)[:100]
        for s in ("OD_CEN", "OD_IJSH"):
            v = [float(b[sc][s]["triangle_area"]) for sc in first100 if s in b[sc]]
            if v:
                sens[f"k{k}_{s}"] = dict(n=len(v), mean=mean(v), median=st.median(v))
    out["sensitivity_first100"] = sens
    json.dump(out, open(os.path.join(HERE, "results", "merged", f"stats{sfx}.json"), "w"), indent=1)

    L = [f"# Task C (O-D road score{' , all-pairs O-D set' if args.variant else ', depot-origin O-D set'}) vs pre-event Task A tables — {out['n']} evaluation cases (seeds from 20264001)\n",
         f"Reproduction check: max |Δ triangle| vs TaskA_v4 evaluation rows for CEN/JSH/IJSH = {maxdev:.3g}; TAP-B fallbacks = {fallbacks}.\n",
         "## Descriptive (triangle loss)\n", "| strategy | mean (95% CI) | median | SD | P90 | P95 | max |", "|---|---|---|---|---|---|---|"]
    for d in out["descriptive"]:
        L.append(f"| {d['strategy']} | {d['mean']:,.1f} ({d['ci_lo']:,.1f}, {d['ci_hi']:,.1f}) | {d['median']:,.1f} | {d['sd']:,.1f} | {d['p90']:,.1f} | {d['p95']:,.1f} | {d['max']:,.1f} |")
    hdr = ["| comparison | mean reduction (95% CI) | % of reference (95% CI) | trimmed (95% CI) | median | wins/losses/ties | sign p |", "|---|---|---|---|---|---|---|"]
    L += ["\n## Paired: triangle loss (positive favours the candidate)\n"] + hdr + [fmt_pair(p) for p in out["paired"]]
    L += ["\n## Paired: critical-weighted triangle\n"] + hdr + [fmt_pair(p) for p in out["paired_weighted"]]
    L += ["\n## Paired: Gini of restoration times (positive = candidate lower/more even)\n"] + hdr + [fmt_pair(p) for p in out["paired_gini"]]
    L += ["\n## Paired: min time-averaged CRI (positive = candidate higher/better)\n"] + hdr + [fmt_pair(p) for p in out["paired_mincri"]]
    L += ["\n## Supplemental descriptives\n", "| strategy | weighted triangle mean | Gini mean | min CRI mean | P90 access-restore mean |", "|---|---|---|---|---|"]
    for a, g, c, p9 in zip(out["descriptive_weighted"], out["descriptive_gini"], out["descriptive_mincri"], out["descriptive_p90"]):
        L.append(f"| {a['strategy']} | {a['mean']:,.1f} | {g['mean']:.3f} | {c['mean']:.3f} | {p9['mean']:,.1f} |")
    L += ["\n## Structure\n", f"- O-D road order identical to CEN / JSH / IJSH road order in {same_as} of {out['n']} cases",
          f"- damaged links with nonzero O-D score: {100*out['damaged_links_nonzero_share']:.1f}%; scenarios with ≥1 such link: {out['scenarios_with_nonzero']}/{out['n']}",
          f"- mean damage: {out['mean_damage']['power']:.1f} buses, {out['mean_damage']['road']:.1f} roads"]
    if decomp:
        L += [f"\n## Decomposition of IJSH vs JSH (n={out['decomp_n']}): which block carries the access benefit?\n"] + hdr + [fmt_pair(p) for p in out["decomp"]]
    L += ["\n## Sensitivity to k (first 100 evaluation cases)\n", "| variant | n | mean | median |", "|---|---|---|---|"]
    for k, v in sens.items():
        L.append(f"| {k} | {v['n']} | {v['mean']:,.1f} | {v['median']:,.1f} |")
    open(os.path.join(HERE, "results", "merged", f"SUMMARY{sfx}.md"), "w", encoding="utf-8").write("\n".join(L))
    print("\n".join(L))


if __name__ == "__main__":
    main()
