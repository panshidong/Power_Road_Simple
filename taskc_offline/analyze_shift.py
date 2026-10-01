from __future__ import annotations

"""Shifted-damage ensembles: O-D representation vs the three pre-event tables (same seeds and
definitions as TaskA_v4/robustness_worker.py). Reproduction check against Codex's robustness rows."""

import csv, glob, json, os
from analyze import by_scenario, descriptive, fmt_pair, paired, HERE

VARIANTS = ["small", "large", "light", "clustered"]


def load_variant(v, sfx=""):
    rows = []
    for fp in sorted(glob.glob(os.path.join(HERE, "results", f"shift_{v}", "worker_*", "rows_*.csv"))):
        rows.extend(csv.DictReader(open(fp, newline="", encoding="utf-8")))
    if sfx:
        rows = [r for r in rows if not r["strategy_id"].startswith("OD_")]
        for fp in sorted(glob.glob(os.path.join(HERE, "results", f"shift_{v}{sfx}", "worker_*", "rows_*.csv"))):
            rows.extend(csv.DictReader(open(fp, newline="", encoding="utf-8")))
    failed = glob.glob(os.path.join(HERE, "results", f"shift_{v}", "worker_*", "FAILED_*.txt"))
    return rows, failed


def main():
    import argparse
    ap = argparse.ArgumentParser(); ap.add_argument("--variant", default="")
    args = ap.parse_args(); sfx = f"_{args.variant}" if args.variant else ""
    ref_rows = [r for fp in glob.glob("/home/workenv/TaskA_v4/results/robustness/worker_*/rows_*.csv") for r in csv.DictReader(open(fp))]
    ref = by_scenario(ref_rows)
    out = {}
    L = ["# O-D representation under shifted damage distributions (tables built under the main distribution)\n"]
    hdr = ["| comparison | mean reduction (95% CI) | % of reference (95% CI) | trimmed (95% CI) | median | wins/losses/ties | sign p |", "|---|---|---|---|---|---|---|"]
    for v in VARIANTS:
        rows, failed = load_variant(v, sfx)
        by = by_scenario(rows)
        S = ["CEN", "JSH", "IJSH", "OD_CEN"] + (["OD_IJSH"] if any("OD_IJSH" in d for d in by.values()) else [])
        by = {sc: d for sc, d in by.items() if all(s in d for s in S)}
        if not by:
            continue
        dev = max((abs(float(by[sc][s]["triangle_area"]) - float(ref[sc][s]["triangle_area"])) for sc in by for s in ("CEN", "JSH", "IJSH") if sc in ref and s in ref[sc]), default=float("nan"))
        desc = descriptive(by, S)
        pairs = [paired(by, "OD_CEN", "CEN"), paired(by, "OD_CEN", "JSH"), paired(by, "OD_CEN", "IJSH")] + ([paired(by, "OD_IJSH", "IJSH")] if "OD_IJSH" in S else []) + [paired(by, "IJSH", "JSH"), paired(by, "JSH", "CEN")]
        nb = sum(int(by[sc]["CEN"]["n_power"]) for sc in by) / len(by); nr = sum(int(by[sc]["CEN"]["n_road"]) for sc in by) / len(by)
        out[v] = {"n": len(by), "failed": len(failed), "max_abs_dev_vs_taskA": dev, "mean_damage": {"power": nb, "road": nr}, "descriptive": desc, "paired": pairs,
                  "paired_weighted": [paired(by, "OD_CEN", r, "weighted_triangle_area") for r in ("CEN", "JSH", "IJSH")]}
        L += [f"\n## {v} (n = {len(by)}, failed/skipped = {len(failed)}, mean damage {nb:.1f} buses / {nr:.1f} roads; max |Δ| vs TaskA rows = {dev:.3g})\n",
              "| strategy | mean (95% CI) | median |", "|---|---|---|"]
        for d in desc:
            L.append(f"| {d['strategy']} | {d['mean']:,.1f} ({d['ci_lo']:,.1f}, {d['ci_hi']:,.1f}) | {d['median']:,.1f} |")
        L += [""] + hdr + [fmt_pair(p) for p in pairs]
    json.dump(out, open(os.path.join(HERE, "results", "merged", f"shift_stats{sfx}.json"), "w"), indent=1)
    open(os.path.join(HERE, "results", "merged", f"SHIFT_SUMMARY{sfx}.md"), "w", encoding="utf-8").write("\n".join(L))
    print("\n".join(L))


if __name__ == "__main__":
    main()
