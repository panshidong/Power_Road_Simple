from __future__ import annotations

"""Aggregate scenario Shapley scores into pre-event expected criticality tables."""

import argparse
import glob
import json
import os
import statistics
from collections import defaultdict


def aggregate(records, strategy: str):
    values = defaultdict(list)
    for record in records:
        for asset, value in record["scores"][strategy].items():
            values[asset].append(float(value))
    return (
        {asset: statistics.fmean(items) for asset, items in values.items()},
        {asset: len(items) for asset, items in values.items()},
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--construction-root", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--expected", type=int, default=400)
    args = ap.parse_args()

    paths = sorted(glob.glob(os.path.join(args.construction_root, "worker_*", "scores_scenario_*.json")))
    records = [json.load(open(path, encoding="utf-8")) for path in paths]
    records.sort(key=lambda item: int(item["scenario_seed"]))
    ids = [item["scenario_id"] for item in records]
    if len(records) != args.expected or len(set(ids)) != len(ids):
        raise RuntimeError(f"Expected {args.expected} unique construction records; found {len(records)} files / {len(set(ids))} IDs")

    checkpoints = [n for n in (50, 100, 200, 400, len(records)) if n <= len(records)]
    checkpoints = sorted(set(checkpoints))
    payload = {
        "method": "Monte Carlo estimate of E[scenario Shapley | asset damaged]",
        "construction_seed_min": records[0]["scenario_seed"],
        "construction_seed_max": records[-1]["scenario_seed"],
        "n_construction": len(records),
        "physical_road_components": True,
        "alpha": 1.0,
        "checkpoints": {},
    }
    for n in checkpoints:
        tables = {}
        counts = {}
        for strategy in ("CEN", "JSH", "IJSH"):
            tables[strategy], counts[strategy] = aggregate(records[:n], strategy)
        payload["checkpoints"][str(n)] = {"tables": tables, "counts": counts}

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=1)

    final = payload["checkpoints"][str(checkpoints[-1])]
    for strategy in ("CEN", "JSH", "IJSH"):
        counts = final["counts"][strategy]
        power = [value for asset, value in counts.items() if asset.startswith("power:")]
        roads = [value for asset, value in counts.items() if asset.startswith("road:")]
        print(
            f"{strategy}: power {len(power)} components, min/median coverage {min(power)}/{statistics.median(power):.0f}; "
            f"roads {len(roads)} components, min/median coverage {min(roads)}/{statistics.median(roads):.0f}"
        )


if __name__ == "__main__":
    main()
