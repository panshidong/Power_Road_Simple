from __future__ import annotations

import csv
import math
import os
from typing import Any, Dict, Iterable, List, Mapping, Sequence

from disaster import DisasterScenario, generate_scenarios
from task_a_criticality import (
    Asset,
    FullFunctionalityValue,
    TaskACriticalityConfig,
    asset_key,
    asset_type,
    sampled_shapley_checkpoints,
)


SAMPLE_SIZES = (30, 60, 120, 240)
REFERENCE_SIZE = max(SAMPLE_SIZES)
RESULT_DIR = "results/task_a_shapley_stability_10"


def _write_csv(path: str, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        return
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def _read_csv(path: str) -> List[Dict[str, str]]:
    if not os.path.exists(path):
        return []
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _sort_assets(scores: Mapping[Asset, float], assets: Iterable[Asset]) -> List[Asset]:
    def tie_key(asset: Asset) -> tuple[int, int, int]:
        if isinstance(asset, tuple):
            return (1, int(asset[0]), int(asset[1]))
        return (0, int(asset), -1)

    return sorted(assets, key=lambda asset: (-float(scores[asset]), tie_key(asset)))


def _spearman(order: Sequence[Asset], reference: Sequence[Asset]) -> float:
    if len(order) != len(reference):
        raise ValueError("Rankings must contain the same assets")
    n = len(order)
    if n <= 1:
        return 1.0
    ranks = {asset: idx + 1 for idx, asset in enumerate(order)}
    ref_ranks = {asset: idx + 1 for idx, asset in enumerate(reference)}
    sum_d2 = sum((ranks[asset] - ref_ranks[asset]) ** 2 for asset in order)
    return 1.0 - 6.0 * float(sum_d2) / float(n * (n * n - 1))


def _mean(values: Sequence[float]) -> float:
    return sum(values) / len(values) if values else float("nan")


def _sample_std(values: Sequence[float]) -> float:
    if len(values) <= 1:
        return 0.0
    mean = _mean(values)
    return math.sqrt(sum((value - mean) ** 2 for value in values) / (len(values) - 1))


def _ci95(values: Sequence[float]) -> float:
    if len(values) <= 1:
        return 0.0
    # Student-t 0.975 critical value for n=10; the study fixes ten scenarios.
    return 2.262 * _sample_std(values) / math.sqrt(len(values))


def _scenario_scores(
    scenario: DisasterScenario,
    cfg: TaskACriticalityConfig,
) -> Dict[str, Dict[int, Dict[Asset, float]]]:
    power_assets = [int(asset) for asset in scenario.broken_buses]
    road_assets = [tuple(map(int, link)) for link in scenario.broken_links]
    link_factors = {
        tuple(map(int, link)): float(factor)
        for link, factor in scenario.link_capacity_factors.items()
    }

    power_game = FullFunctionalityValue(
        power_assets=power_assets,
        road_assets=[],
        link_factors={},
        cfg=cfg,
    )
    road_game = FullFunctionalityValue(
        power_assets=[],
        road_assets=road_assets,
        link_factors=link_factors,
        cfg=cfg,
    )
    integrated_game = FullFunctionalityValue(
        power_assets=power_assets,
        road_assets=road_assets,
        link_factors=link_factors,
        cfg=cfg,
    )
    return {
        "S1_separate_power": sampled_shapley_checkpoints(
            power_assets,
            power_game.power_value,
            sample_sizes=SAMPLE_SIZES,
            seed=cfg.shapley_seed + scenario.seed + 101,
        ),
        "S1_separate_road": sampled_shapley_checkpoints(
            road_assets,
            road_game.road_value,
            sample_sizes=SAMPLE_SIZES,
            seed=cfg.shapley_seed + scenario.seed + 202,
        ),
        "S2_integrated": sampled_shapley_checkpoints(
            list(power_assets) + list(road_assets),
            integrated_game.integrated_value,
            sample_sizes=SAMPLE_SIZES,
            seed=cfg.shapley_seed + scenario.seed + 303,
        ),
    }


def _detail_rows(
    scenario: DisasterScenario,
    score_sets: Mapping[str, Mapping[int, Mapping[Asset, float]]],
) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for method, checkpoints in score_sets.items():
        sectors = ("power", "road") if method == "S2_integrated" else (method.rsplit("_", 1)[-1],)
        for sector in sectors:
            reference_scores = checkpoints[REFERENCE_SIZE]
            assets = [asset for asset in reference_scores if asset_type(asset) == sector]
            reference = _sort_assets(reference_scores, assets)
            top_count = min(3, len(assets))
            for sample_size in SAMPLE_SIZES:
                order = _sort_assets(checkpoints[sample_size], assets)
                overlap = len(set(order[:top_count]) & set(reference[:top_count])) / max(top_count, 1)
                rows.append(
                    {
                        "scenario_id": scenario.scenario_id,
                        "scenario_seed": scenario.seed,
                        "method": method,
                        "sector": sector,
                        "sample_size": sample_size,
                        "reference_size": REFERENCE_SIZE,
                        "asset_count": len(assets),
                        "spearman_rho": _spearman(order, reference),
                        "top1_agreement": float(order[:1] == reference[:1]),
                        "top3_overlap": overlap,
                        "ranking": " | ".join(asset_key(asset) for asset in order),
                        "reference_ranking": " | ".join(asset_key(asset) for asset in reference),
                    }
                )
    return rows


def _score_rows(
    scenario: DisasterScenario,
    score_sets: Mapping[str, Mapping[int, Mapping[Asset, float]]],
) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for method, checkpoints in score_sets.items():
        for sample_size, scores in checkpoints.items():
            for asset, score in scores.items():
                rows.append(
                    {
                        "scenario_id": scenario.scenario_id,
                        "scenario_seed": scenario.seed,
                        "method": method,
                        "sample_size": sample_size,
                        "asset": asset_key(asset),
                        "asset_type": asset_type(asset),
                        "score": float(score),
                    }
                )
    return rows


def _aggregate(rows: Sequence[Mapping[str, Any]]) -> List[Dict[str, Any]]:
    groups: Dict[tuple[str, str, int], List[Mapping[str, Any]]] = {}
    for row in rows:
        key = (str(row["method"]), str(row["sector"]), int(row["sample_size"]))
        groups.setdefault(key, []).append(row)
    output: List[Dict[str, Any]] = []
    for (method, sector, sample_size), items in sorted(groups.items()):
        rho = [float(item["spearman_rho"]) for item in items]
        top1 = [float(item["top1_agreement"]) for item in items]
        top3 = [float(item["top3_overlap"]) for item in items]
        output.append(
            {
                "method": method,
                "sector": sector,
                "sample_size": sample_size,
                "reference_size": REFERENCE_SIZE,
                "n_scenarios": len(items),
                "spearman_rho_mean": _mean(rho),
                "spearman_rho_ci95": _ci95(rho),
                "top1_agreement_mean": _mean(top1),
                "top3_overlap_mean": _mean(top3),
            }
        )
    return output


def main() -> None:
    os.makedirs(RESULT_DIR, exist_ok=True)
    details_path = os.path.join(RESULT_DIR, "stability_detail.csv")
    scores_path = os.path.join(RESULT_DIR, "stability_scores.csv")
    aggregate_path = os.path.join(RESULT_DIR, "stability_aggregate.csv")
    details: List[Dict[str, Any]] = list(_read_csv(details_path))
    scores: List[Dict[str, Any]] = list(_read_csv(scores_path))
    detail_scenarios = {str(row["scenario_id"]) for row in details}
    score_scenarios = {str(row["scenario_id"]) for row in scores}
    completed = detail_scenarios & score_scenarios

    scenarios = generate_scenarios(
        n_scenarios=10,
        seed0=20260402,
        bus_count_range=(8, 15),
        link_count_range=(3, 13),
        link_drop_range=(0.5, 1.0),
        out_json=os.path.join(RESULT_DIR, "random_disasters.json"),
    )
    cfg = TaskACriticalityConfig(
        shapley_samples=120,
        shapley_seed=20260402,
        integrated_power_weight=0.5,
        use_full_functionality_shapley=True,
    )
    for index, scenario in enumerate(scenarios, start=1):
        if scenario.scenario_id in completed:
            print(f"[stability] skipping completed {scenario.scenario_id} ({index}/10)")
            continue
        print(f"[stability] running {scenario.scenario_id} ({index}/10)")
        details = [row for row in details if row["scenario_id"] != scenario.scenario_id]
        scores = [row for row in scores if row["scenario_id"] != scenario.scenario_id]
        score_sets = _scenario_scores(scenario, cfg)
        details.extend(_detail_rows(scenario, score_sets))
        scores.extend(_score_rows(scenario, score_sets))
        _write_csv(details_path, details)
        _write_csv(scores_path, scores)
        _write_csv(aggregate_path, _aggregate(details))

    print("Task A Shapley stability study complete.")
    print("Aggregate:", aggregate_path)


if __name__ == "__main__":
    main()
