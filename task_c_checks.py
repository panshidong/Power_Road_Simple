from __future__ import annotations

from typing import Any, Dict

from disaster import generate_scenarios
from resilience_measurement import run_model_multi
from task_a_criticality import TaskACriticalityConfig, build_task_a_strategies, read_road_network
from task_c_od import (
    TaskCODConfig,
    build_critical_od_pairs,
    build_static_road_importance,
    build_task_c_strategies,
    k_shortest_paths,
    load_bus_location,
    route_overlap_structure,
)
from task_c_weighted_triangle import CriticalWeightConfig, compute_weighted_triangle


def check_k_shortest_paths_basic() -> Dict[str, Any]:
    road_net = read_road_network("tap-b/net/SiouxFalls_net.txt")
    paths = k_shortest_paths(road_net.costs, 1, 20, 3)
    if not paths:
        raise AssertionError("Expected at least one baseline path between node 1 and node 20")
    costs = [cost for cost, _ in paths]
    if costs != sorted(costs):
        raise AssertionError(f"k-shortest paths must be non-decreasing in cost, got {costs}")
    node_paths = [tuple(path) for _, path in paths]
    if len(set(node_paths)) != len(node_paths):
        raise AssertionError("k-shortest paths returned duplicate routes")
    for _, path in paths:
        if path[0] != 1 or path[-1] != 20:
            raise AssertionError(f"Path does not connect the requested O-D pair: {path}")
    return {"n_paths": len(paths), "costs": costs}


def check_bus_location_used_not_bus_to_link() -> Dict[str, Any]:
    """Regression guard for the destination-mapping bug in the first S3 design:
    critical-bus destinations must come from bus_location.json (single node per
    bus, matching resilience_measurement's own crew-dispatch mapping), not from
    bus_to_link.json (the unrelated power->road signal-dependency mapping)."""
    cfg = TaskCODConfig()
    bus_location = load_bus_location(cfg.bus_location_path)
    od_pairs = build_critical_od_pairs(cfg, bus_location=bus_location)

    for bus in cfg.critical_buses:
        matching = [od for od in od_pairs if od.destination_kind == "critical_bus" and od.destination_ref == bus]
        if not matching:
            raise AssertionError(f"No O-D pair generated for critical bus {bus}")
        expected_node = bus_location[bus]
        for od in matching:
            if od.destination != expected_node:
                raise AssertionError(
                    f"Critical bus {bus} destination {od.destination} does not match "
                    f"bus_location.json ({expected_node})"
                )
    return {"critical_buses": list(cfg.critical_buses), "n_od_pairs": len(od_pairs)}


def check_od_set_is_scenario_independent() -> Dict[str, Any]:
    """The O-D set (and therefore the resulting road-link importance table)
    must not depend on which assets a particular disruption scenario happens
    to damage."""
    cfg = TaskCODConfig()
    bus_location = load_bus_location(cfg.bus_location_path)
    pairs_a = build_critical_od_pairs(cfg, bus_location=bus_location)
    pairs_b = build_critical_od_pairs(cfg, bus_location=bus_location)
    ids_a = sorted(od.pair_id for od in pairs_a)
    ids_b = sorted(od.pair_id for od in pairs_b)
    if ids_a != ids_b:
        raise AssertionError("O-D pair set is not deterministic/repeatable")

    task_a_cfg = TaskACriticalityConfig()
    importance_1 = build_static_road_importance(cfg, task_a_cfg)
    importance_2 = build_static_road_importance(cfg, task_a_cfg)
    if importance_1 != importance_2:
        raise AssertionError("Static road importance table is not deterministic across calls")

    return {"n_od_pairs": len(pairs_a), "n_scored_links": len(importance_1)}


def check_route_overlap_basic() -> Dict[str, Any]:
    cfg = TaskCODConfig(k_paths=2)
    road_net = read_road_network("tap-b/net/SiouxFalls_net.txt")
    bus_location = load_bus_location(cfg.bus_location_path)
    od_pairs = build_critical_od_pairs(cfg, bus_location=bus_location)
    od_paths, link_to_services = route_overlap_structure(road_net, od_pairs, cfg)

    unreachable = [pair_id for pair_id, paths in od_paths.items() if not paths]
    if unreachable:
        raise AssertionError(f"O-D pairs with no feasible route: {unreachable}")
    if not link_to_services:
        raise AssertionError("Expected at least one link to carry O-D route overlap")

    return {"n_od_pairs": len(od_pairs), "n_links_with_service_overlap": len(link_to_services)}


def check_s3_strategy_generation() -> Dict[str, Any]:
    scenario = generate_scenarios(
        n_scenarios=1,
        seed0=20260901,
        bus_count_range=(5, 5),
        link_count_range=(4, 4),
        out_json="/tmp/task_c_check_strategy_disasters.json",
    )[0]
    task_a_cfg = TaskACriticalityConfig(shapley_samples=20, shapley_seed=20260901)
    od_cfg = TaskCODConfig(k_paths=3)

    baseline_strategies = build_task_a_strategies(scenario, cfg=task_a_cfg)
    s0 = next(s for s in baseline_strategies if s.strategy_id == "S0_centrality")

    s3_strategies = build_task_c_strategies(scenario, od_cfg=od_cfg, task_a_cfg=task_a_cfg)
    strategy_ids = [s.strategy_id for s in s3_strategies]
    if strategy_ids != ["S3_od_representation"]:
        raise AssertionError(f"Unexpected Task C strategy ids: {strategy_ids}")
    s3 = s3_strategies[0]

    if len(s3.sequence) != len(scenario.broken_buses) + len(scenario.broken_links):
        raise AssertionError("S3 sequence does not include every damaged asset exactly once")
    if set(s3.sequence) != set(s0.sequence):
        raise AssertionError("S3 and S0 should rank exactly the same set of damaged assets")
    if set(s3.power_sequence) != set(s0.power_sequence):
        raise AssertionError("S3's power side should reuse S0's power assets/scores unchanged")

    return {
        "scenario_id": scenario.scenario_id,
        "s0_road_sequence": list(s0.road_sequence),
        "s3_road_sequence": list(s3.road_sequence),
        "road_sequences_identical": list(s0.road_sequence) == list(s3.road_sequence),
    }


def check_weighted_triangle_runs_and_matches_functional_logic() -> Dict[str, Any]:
    scenario = generate_scenarios(
        n_scenarios=1,
        seed0=20260902,
        bus_count_range=(2, 2),
        link_count_range=(2, 2),
        out_json="/tmp/task_c_check_weighted_triangle_disasters.json",
    )[0]
    task_a_cfg = TaskACriticalityConfig(shapley_samples=10, shapley_seed=20260902)
    od_cfg = TaskCODConfig(k_paths=2)
    s3 = build_task_c_strategies(scenario, od_cfg=od_cfg, task_a_cfg=task_a_cfg)[0]

    run = run_model_multi(
        s3.sequence,
        result_root="results",
        run_dir="results/task_c_checks_weighted_triangle",
        message="Task C check: weighted triangle",
        Scenario="task_c_check_weighted",
        strict=True,
        save_artifacts=False,
        crew_mode="specialized",
        power_crews=1,
        road_crews=1,
        preserve_sequence_order=True,
        broken_link_factors=dict(scenario.link_capacity_factors),
        objective="triangle",
    )

    w_cfg = CriticalWeightConfig()
    result = compute_weighted_triangle(run, cfg=w_cfg)

    if result["weighted_triangle_area"] < 0:
        raise AssertionError(f"weighted_triangle_area should be non-negative, got {result['weighted_triangle_area']}")
    if not (0.0 <= result["shelter_access_initial"] <= 1.0 + 1e-9):
        raise AssertionError(f"shelter_access_initial out of [0,1]: {result['shelter_access_initial']}")
    if result["shelter_access_final"] < result["shelter_access_initial"] - 1e-9:
        raise AssertionError("Shelter access should not get worse than initial after full repair")
    if result["weighted_power_final"] < result["weighted_power_initial"] - 1e-9:
        raise AssertionError("Weighted power performance should not get worse than initial after full repair")

    return {"scenario_id": scenario.scenario_id, **result}


def main() -> None:
    checks = {
        "k_shortest_paths_basic": check_k_shortest_paths_basic(),
        "bus_location_used_not_bus_to_link": check_bus_location_used_not_bus_to_link(),
        "od_set_is_scenario_independent": check_od_set_is_scenario_independent(),
        "route_overlap_basic": check_route_overlap_basic(),
        "s3_strategy_generation": check_s3_strategy_generation(),
        "weighted_triangle_runs": check_weighted_triangle_runs_and_matches_functional_logic(),
    }
    for name, payload in checks.items():
        print(f"PASS {name}: {payload}")


if __name__ == "__main__":
    main()
