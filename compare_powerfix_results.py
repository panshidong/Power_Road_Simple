from __future__ import annotations

import ast
import csv
import os
import statistics
from collections import defaultdict
from typing import Any, Dict, List, Mapping, Sequence


OLD_DIR = "results/task_a_criticality_final_100"
NEW_DIR = "results/task_a_criticality_powerfix_100"
OUTPUT_DIR = os.path.join(NEW_DIR, "powerfix_comparison")


def read_csv(path: str) -> List[Dict[str, str]]:
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: str, rows: Sequence[Mapping[str, Any]]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def keyed(rows: Sequence[Mapping[str, str]]) -> Dict[tuple[str, str], Mapping[str, str]]:
    return {(row["scenario_id"], row["strategy_id"]): row for row in rows}


def main() -> None:
    old_rows = read_csv(os.path.join(OLD_DIR, "scenario_strategy_rows.csv"))
    new_rows = read_csv(os.path.join(NEW_DIR, "scenario_strategy_rows.csv"))
    old = keyed(old_rows)
    new = keyed(new_rows)
    if len(old) != 400 or set(old) != set(new):
        raise ValueError("Both result sets must contain the same 400 scenario-strategy keys")

    old_disasters = open(os.path.join(OLD_DIR, "random_disasters.json"), encoding="utf-8").read()
    new_disasters = open(os.path.join(NEW_DIR, "random_disasters.json"), encoding="utf-8").read()
    if old_disasters != new_disasters:
        raise ValueError("Old and corrected result sets do not use identical disaster inputs")

    detail: List[Dict[str, Any]] = []
    grouped: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    for key in sorted(old):
        old_row, new_row = old[key], new[key]
        strategy = key[1]
        item = {
            "scenario_id": key[0],
            "strategy_id": strategy,
            "old_triangle_area": float(old_row["triangle_area"]),
            "new_triangle_area": float(new_row["triangle_area"]),
            "triangle_change_new_minus_old": float(new_row["triangle_area"]) - float(old_row["triangle_area"]),
            "full_sequence_changed": int(ast.literal_eval(old_row["sequence"]) != ast.literal_eval(new_row["sequence"])),
            "power_sequence_changed": int(ast.literal_eval(old_row["power_sequence"]) != ast.literal_eval(new_row["power_sequence"])),
            "road_sequence_changed": int(ast.literal_eval(old_row["road_sequence"]) != ast.literal_eval(new_row["road_sequence"])),
        }
        detail.append(item)
        grouped[strategy].append(item)

    summary: List[Dict[str, Any]] = []
    for strategy, items in sorted(grouped.items()):
        old_values = [float(row["old_triangle_area"]) for row in items]
        new_values = [float(row["new_triangle_area"]) for row in items]
        summary.append(
            {
                "strategy_id": strategy,
                "n": len(items),
                "old_mean_triangle": statistics.fmean(old_values),
                "new_mean_triangle": statistics.fmean(new_values),
                "mean_change_new_minus_old": statistics.fmean(n - o for n, o in zip(new_values, old_values)),
                "full_sequence_changed_scenarios": sum(int(row["full_sequence_changed"]) for row in items),
                "power_sequence_changed_scenarios": sum(int(row["power_sequence_changed"]) for row in items),
                "road_sequence_changed_scenarios": sum(int(row["road_sequence_changed"]) for row in items),
            }
        )

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    write_csv(os.path.join(OUTPUT_DIR, "scenario_comparison.csv"), detail)
    write_csv(os.path.join(OUTPUT_DIR, "summary.csv"), summary)
    print("Power-fix comparison complete")
    for row in summary:
        print(row)


if __name__ == "__main__":
    main()
