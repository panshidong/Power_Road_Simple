from __future__ import annotations

import os
from typing import Any, Dict, List, Set

from disaster import generate_scenarios
from resilience_measurement import run_model_multi
from task_a_criticality import (
    FullFunctionalityValue,
    TaskACriticalityConfig,
    build_task_a_strategies,
    exact_shapley_by_permutation,
    sampled_shapley,
    shapley_efficiency_gap,
)


def _assert_close(actual: float, expected: float, *, tol: float, label: str) -> None:
    if abs(float(actual) - float(expected)) > float(tol):
        raise AssertionError(f"{label}: expected {expected}, got {actual}")


def check_known_shapley_game() -> Dict[str, Any]:
    assets = [1, 2, 3]

    def value_fn(coalition: Set[int]) -> float:
        base = 0.0
        base += 2.0 if 1 in coalition else 0.0
        base += 3.0 if 2 in coalition else 0.0
        base += 5.0 if 3 in coalition else 0.0
        base += 4.0 if 1 in coalition and 2 in coalition else 0.0
        return base

    exact = exact_shapley_by_permutation(assets, value_fn)
    expected = {1: 4.0, 2: 5.0, 3: 5.0}
    for asset, expected_value in expected.items():
        _assert_close(exact[asset], expected_value, tol=1e-12, label=f"exact Shapley asset {asset}")

    sampled = sampled_shapley(assets, value_fn, samples=200, seed=19)
    gap = shapley_efficiency_gap(assets, sampled, value_fn)
    _assert_close(gap, 0.0, tol=1e-12, label="sampled Shapley efficiency gap")
    return {"exact": exact, "sampled": sampled, "efficiency_gap": gap}


def check_strategy_generation() -> Dict[str, Any]:
    scenario = generate_scenarios(
        n_scenarios=1,
        seed0=20260405,
        bus_count_range=(4, 4),
        link_count_range=(3, 3),
        out_json="/tmp/task_a_check_disasters.json",
    )[0]
    strategies = build_task_a_strategies(
        scenario,
        cfg=TaskACriticalityConfig(shapley_samples=24, shapley_seed=20260405),
    )
    ids = [strategy.strategy_id for strategy in strategies]
    expected_ids = ["S0_centrality", "S1_separate_shapley", "S2_integrated_shapley"]
    if ids != expected_ids:
        raise AssertionError(f"Unexpected strategy ids: {ids}")

    s0 = strategies[0]
    if "power radial" not in s0.score_method or "road weighted edge" not in s0.score_method:
        raise AssertionError(f"S0 score method should be network-specific; got {s0.score_method!r}")

    s2 = strategies[2]
    if not s2.preserve_sequence_order:
        raise AssertionError("S2 must preserve the integrated full-network ranking")
    if len(s2.sequence) != len(scenario.broken_buses) + len(scenario.broken_links):
        raise AssertionError("S2 sequence does not include every damaged asset exactly once")

    return {
        "scenario_id": scenario.scenario_id,
        "strategy_ids": ids,
        "s2_first_assets": list(s2.sequence[:5]),
    }


def check_full_functionality_shapley_game() -> Dict[str, Any]:
    try:
        TaskACriticalityConfig(use_full_functionality_shapley=False)
    except ValueError:
        pass
    else:
        raise AssertionError("Proxy Shapley mode must fail instead of silently changing the value function")

    scenario = generate_scenarios(
        n_scenarios=1,
        seed0=20260407,
        bus_count_range=(1, 1),
        link_count_range=(1, 1),
        out_json="/tmp/task_a_check_full_functionality_disasters.json",
    )[0]
    cfg = TaskACriticalityConfig(shapley_samples=2, shapley_seed=20260407)
    power_assets = [int(asset) for asset in scenario.broken_buses]
    road_assets = [tuple(map(int, link)) for link in scenario.broken_links]
    link_factors = {tuple(map(int, link)): float(factor) for link, factor in scenario.link_capacity_factors.items()}
    assets = list(power_assets) + list(road_assets)

    full_game = FullFunctionalityValue(
        power_assets=power_assets,
        road_assets=road_assets,
        link_factors=link_factors,
        cfg=cfg,
    )
    exact = exact_shapley_by_permutation(assets, full_game.integrated_value)
    gap = shapley_efficiency_gap(assets, exact, full_game.integrated_value)
    _assert_close(gap, 0.0, tol=1e-12, label="full-functionality exact Shapley efficiency gap")

    empty_value = full_game.integrated_value(set())
    full_value = full_game.integrated_value(set(assets))
    if full_value < empty_value:
        raise AssertionError(f"Expected repaired-state functionality to improve; empty={empty_value}, full={full_value}")

    power_game = FullFunctionalityValue(
        power_assets=power_assets,
        road_assets=[],
        link_factors={},
        cfg=cfg,
    )
    _ = power_game.power_value(set())
    if power_game.road_cache:
        raise AssertionError("Separate power Shapley game should not evaluate road/TAP-B states")

    road_game = FullFunctionalityValue(
        power_assets=[],
        road_assets=road_assets,
        link_factors=link_factors,
        cfg=cfg,
    )
    _ = road_game.road_value(set())
    if any(power_key for power_key, _ in road_game.road_cache):
        raise AssertionError("Separate road Shapley game should not include broken power buses")

    return {
        "scenario_id": scenario.scenario_id,
        "assets": assets,
        "empty_value": empty_value,
        "full_value": full_value,
        "efficiency_gap": gap,
    }


def check_specialized_crews_and_dynamic_tapb() -> Dict[str, Any]:
    scenario = generate_scenarios(
        n_scenarios=1,
        seed0=20260406,
        bus_count_range=(2, 2),
        link_count_range=(2, 2),
        out_json="/tmp/task_a_check_dynamic_disasters.json",
    )[0]
    s2 = build_task_a_strategies(
        scenario,
        cfg=TaskACriticalityConfig(shapley_samples=12, shapley_seed=20260406),
    )[2]
    run = run_model_multi(
        list(s2.sequence),
        result_root="results",
        run_dir="results/task_a_checks_dynamic",
        message="Task A check: specialized crews and dynamic TAP-B states",
        Scenario="task_a_check_s2",
        strict=True,
        save_artifacts=False,
        crew_mode="specialized",
        power_crews=1,
        road_crews=1,
        preserve_sequence_order=True,
        broken_link_factors=dict(scenario.link_capacity_factors),
        objective="triangle",
    )

    if run["crew_mode"] != "specialized":
        raise AssertionError(f"Expected specialized crew mode, got {run['crew_mode']!r}")
    if not run["preserve_sequence_order"]:
        raise AssertionError("Expected S2 run to preserve the integrated ranking")
    if int(run["power_crews"]) < 1 or int(run["road_crews"]) < 1:
        raise AssertionError("Expected separate power and road repair crews")
    if list(run["sequence"]) != list(s2.sequence):
        raise AssertionError("Simulator changed the S2 integrated priority list")

    snapshot_paths: List[str] = []
    for event in run["event_log"]:
        path = event["state"].get("s_txt_path")
        if path:
            snapshot_paths.append(path)
    unique_snapshots = sorted(set(snapshot_paths))
    if len(unique_snapshots) < 2:
        raise AssertionError("Expected multiple TAP-B state snapshots across repair events")
    missing = [path for path in unique_snapshots if not os.path.exists(path)]
    if missing:
        raise AssertionError(f"Missing TAP-B state snapshot files: {missing}")

    return {
        "scenario_id": scenario.scenario_id,
        "event_count": len(run["event_log"]),
        "tapb_state_snapshots": len(unique_snapshots),
        "sequence": list(run["sequence"]),
    }


def main() -> None:
    checks = {
        "known_shapley_game": check_known_shapley_game(),
        "strategy_generation": check_strategy_generation(),
        "full_functionality_shapley_game": check_full_functionality_shapley_game(),
        "specialized_crews_dynamic_tapb": check_specialized_crews_and_dynamic_tapb(),
    }
    for name, payload in checks.items():
        print(f"PASS {name}: {payload}")


if __name__ == "__main__":
    main()
