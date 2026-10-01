from __future__ import annotations

import heapq
import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Set, Tuple

from disaster import DisasterScenario
from task_a_criticality import (
    RoadNetwork,
    TaskACriticalityConfig,
    TaskAStrategy,
    _normalize_link,
    _sorted_by_score,
    _split_assets,
    power_radial_service_centrality,
    read_road_network,
)

Asset = Any
Link = Tuple[int, int]

"""
Task C: O-D based interdependence representation (proposal Section 4.4, RQ3/H3).

REDESIGN (superseding the first version): the O-D set used to score road-link
importance is a FIXED, damage-independent set of "critical" destinations (a
small number of designated power buses and a shelter node), not the set of
assets that happen to be damaged in a given random scenario. This mirrors how
S0's topological centrality is a static property of the network, computed once
and then looked up for whichever assets are damaged in a scenario -- so S3 is
now genuinely comparable to S0 as "a different way to score the same road
network," rather than behaving like a third scenario-specific Shapley variant
(which is what the damage-driven version accidentally became). See
TASKC_ASSUMPTIONS.md for why this changed and what was wrong with the first
version (crew-origin/destination construction bug, and the damage-specific
O-D set).

This module only adds new code. It imports from task_a_criticality.py (Task A)
but does not modify it, so Task A's own branch/behavior is unaffected.
"""


@dataclass(frozen=True)
class TaskCODConfig:
    road_net_path: str = "tap-b/net/SiouxFalls_net.txt"
    bus_location_path: str = "bus_location.json"
    depot_nodes: Tuple[int, ...] = (1,)
    # Fixed, damage-independent "critical" destination set (see TASKC_ASSUMPTIONS.md
    # for how these were chosen).
    critical_buses: Tuple[int, ...] = (8, 17)
    shelter_nodes: Tuple[int, ...] = (24,)
    k_paths: int = 3
    path_rank_weights: Tuple[float, ...] = (1.0, 0.6, 0.35)
    allocation_method: str = "weighted_overlap"

    def __post_init__(self) -> None:
        if int(self.k_paths) < 1:
            raise ValueError("k_paths must be >= 1")
        weights = [float(w) for w in self.path_rank_weights]
        if not weights:
            raise ValueError("path_rank_weights must be non-empty")
        if any(w <= 0.0 for w in weights):
            raise ValueError("path_rank_weights must be strictly positive")
        if weights != sorted(weights, reverse=True):
            raise ValueError("path_rank_weights must be non-increasing by path rank")
        if self.allocation_method not in ("weighted_overlap",):
            raise ValueError("allocation_method must be 'weighted_overlap'")
        if not self.depot_nodes:
            raise ValueError("depot_nodes must be non-empty")
        if not self.critical_buses and not self.shelter_nodes:
            raise ValueError("at least one critical bus or shelter node is required")


def load_taskc_od_config(path: str = "taskC_od_config.json") -> TaskCODConfig:
    if not os.path.exists(path):
        return TaskCODConfig()
    raw = json.loads(Path(path).read_text(encoding="utf-8"))
    return TaskCODConfig(
        road_net_path=str(raw.get("road_net_path", "tap-b/net/SiouxFalls_net.txt")),
        bus_location_path=str(raw.get("bus_location_path", "bus_location.json")),
        depot_nodes=tuple(int(n) for n in raw.get("depot_nodes", [1])),
        critical_buses=tuple(int(b) for b in raw.get("critical_buses", [8, 17])),
        shelter_nodes=tuple(int(n) for n in raw.get("shelter_nodes", [24])),
        k_paths=int(raw.get("k_paths", 3)),
        path_rank_weights=tuple(float(w) for w in raw.get("path_rank_weights", [1.0, 0.6, 0.35])),
        allocation_method=str(raw.get("allocation_method", "weighted_overlap")),
    )


def load_bus_location(path: str) -> Dict[int, int]:
    """bus -> single road node where the bus physically sits (NOT bus_to_link.json,
    which is the unrelated power->road signal-dependency mapping). This is the
    same mapping the base simulator itself uses to send a crew to repair a bus
    (resilience_measurement.run_model_multi's bus_location_path)."""
    raw = json.loads(Path(path).read_text(encoding="utf-8"))
    return {int(k): int(v) for k, v in raw.items()}


@dataclass(frozen=True)
class ODPair:
    destination_kind: str  # "critical_bus" or "shelter"
    destination_ref: int  # bus id or shelter node id (for labeling)
    origin: int
    destination: int
    pair_id: str


def build_critical_od_pairs(
    cfg: TaskCODConfig,
    *,
    bus_location: Mapping[int, int],
) -> List[ODPair]:
    """Build the fixed, damage-independent O-D set: depot(s) -> each designated
    critical bus's physical location, and depot(s) -> each designated shelter
    node. This does not depend on any disruption scenario."""
    pairs: List[ODPair] = []
    for origin in cfg.depot_nodes:
        for bus in cfg.critical_buses:
            destination = bus_location.get(int(bus))
            if destination is None:
                raise KeyError(f"No bus_location entry for critical bus {bus}")
            if int(origin) == int(destination):
                continue
            pairs.append(
                ODPair(
                    destination_kind="critical_bus",
                    destination_ref=int(bus),
                    origin=int(origin),
                    destination=int(destination),
                    pair_id=f"critical_bus:{int(bus)}|o{int(origin)}|d{int(destination)}",
                )
            )
        for shelter in cfg.shelter_nodes:
            if int(origin) == int(shelter):
                continue
            pairs.append(
                ODPair(
                    destination_kind="shelter",
                    destination_ref=int(shelter),
                    origin=int(origin),
                    destination=int(shelter),
                    pair_id=f"shelter:{int(shelter)}|o{int(origin)}|d{int(shelter)}",
                )
            )
    return pairs


def _adjacency(costs: Mapping[Link, float]) -> Dict[int, List[Tuple[int, float]]]:
    graph: Dict[int, List[Tuple[int, float]]] = {}
    for (u, v), cost in costs.items():
        graph.setdefault(int(u), []).append((int(v), max(float(cost), 1e-9)))
        graph.setdefault(int(v), [])
    return graph


def _path_cost(costs: Mapping[Link, float], path: Sequence[int]) -> float:
    total = 0.0
    for u, v in zip(path[:-1], path[1:]):
        total += float(costs.get((int(u), int(v)), costs.get((int(v), int(u)), 0.0)))
    return total


def dijkstra_path(
    costs: Mapping[Link, float],
    source: int,
    target: int,
    *,
    excluded_nodes: Optional[Set[int]] = None,
    excluded_links: Optional[Set[Link]] = None,
) -> Optional[Tuple[float, List[int]]]:
    excluded_nodes = excluded_nodes or set()
    excluded_links = excluded_links or set()
    graph = _adjacency(costs)

    dist: Dict[int, float] = {int(source): 0.0}
    prev: Dict[int, int] = {}
    visited: Set[int] = set()
    pq: List[Tuple[float, int]] = [(0.0, int(source))]

    while pq:
        d, u = heapq.heappop(pq)
        if u in visited:
            continue
        visited.add(u)
        if u == int(target):
            break
        for v, w in graph.get(u, []):
            if v in excluded_nodes or (u, v) in excluded_links or (v, u) in excluded_links:
                continue
            nd = d + w
            if nd < dist.get(v, float("inf")):
                dist[v] = nd
                prev[v] = u
                heapq.heappush(pq, (nd, v))

    if int(target) not in dist:
        return None
    path = [int(target)]
    while path[-1] != int(source):
        path.append(prev[path[-1]])
    path.reverse()
    return dist[int(target)], path


def k_shortest_paths(
    costs: Mapping[Link, float],
    source: int,
    target: int,
    k: int,
) -> List[Tuple[float, List[int]]]:
    """Yen's algorithm for k loopless shortest paths, built on dijkstra_path above."""
    source, target = int(source), int(target)
    if source == target:
        return [(0.0, [source])]

    first = dijkstra_path(costs, source, target)
    if first is None:
        return []

    a_paths: List[Tuple[float, List[int]]] = [first]
    b_candidates: List[Tuple[float, List[int]]] = []
    seen_paths: Set[Tuple[int, ...]] = {tuple(first[1])}

    for _ in range(1, int(k)):
        prev_cost, prev_path = a_paths[-1]
        for i in range(len(prev_path) - 1):
            spur_node = prev_path[i]
            root_path = prev_path[: i + 1]

            excluded_links: Set[Link] = set()
            for _, path in a_paths:
                if len(path) > i and path[: i + 1] == root_path:
                    excluded_links.add((path[i], path[i + 1]))

            excluded_nodes = set(root_path[:-1])
            spur = dijkstra_path(
                costs,
                spur_node,
                target,
                excluded_nodes=excluded_nodes,
                excluded_links=excluded_links,
            )
            if spur is None:
                continue
            spur_cost, spur_path = spur
            total_path = root_path[:-1] + spur_path
            key = tuple(total_path)
            if key in seen_paths:
                continue
            total_cost = _path_cost(costs, root_path) + spur_cost
            b_candidates.append((total_cost, total_path))
            seen_paths.add(key)

        if not b_candidates:
            break
        b_candidates.sort(key=lambda item: item[0])
        best = b_candidates.pop(0)
        a_paths.append(best)

    return a_paths


def links_on_path(path: Sequence[int]) -> List[Link]:
    return [(int(path[i]), int(path[i + 1])) for i in range(len(path) - 1)]


def route_overlap_structure(
    road_net: RoadNetwork,
    od_pairs: Sequence[ODPair],
    cfg: TaskCODConfig,
) -> Tuple[Dict[str, List[Tuple[float, List[int]]]], Dict[Link, List[Tuple[str, int, float]]]]:
    """Build the O-D route-overlap ("hyper-edge") structure over the FIXED
    critical O-D set (independent of any disruption scenario).

    Returns:
      od_paths: pair_id -> up to k candidate (cost, node_path) routes, cheapest first
      link_to_services: directed road link -> list of (pair_id, path_rank, weight)
        for every candidate route that traverses it
    """
    weights = list(cfg.path_rank_weights)
    od_paths: Dict[str, List[Tuple[float, List[int]]]] = {}
    link_to_services: Dict[Link, List[Tuple[str, int, float]]] = {}

    for od in od_pairs:
        paths = k_shortest_paths(road_net.costs, od.origin, od.destination, cfg.k_paths)
        od_paths[od.pair_id] = paths
        for rank, (_, path) in enumerate(paths):
            weight = weights[rank] if rank < len(weights) else weights[-1]
            for link in links_on_path(path):
                link_to_services.setdefault(link, []).append((od.pair_id, rank, weight))

    return od_paths, link_to_services


def allocate_link_importance(
    link_to_services: Mapping[Link, List[Tuple[str, int, float]]],
) -> Dict[Link, float]:
    """Explicit allocation rule: importance(link) = sum of path-rank weights of every
    O-D service and candidate route that traverses that (directed) link."""
    return {link: sum(weight for (_, _, weight) in entries) for link, entries in link_to_services.items()}


def undirected_link_importance(
    directed_importance: Mapping[Link, float],
    links: Sequence[Link],
) -> Dict[Link, float]:
    """Sum a directed importance map across both travel directions for the given
    (undirected) links; road damage is treated as physical two-way blockage."""
    out: Dict[Link, float] = {}
    for link in links:
        u, v = int(link[0]), int(link[1])
        out[(u, v)] = float(directed_importance.get((u, v), 0.0)) + float(directed_importance.get((v, u), 0.0))
    return out


def build_static_road_importance(
    od_cfg: TaskCODConfig,
    task_a_cfg: TaskACriticalityConfig,
) -> Dict[Link, float]:
    """Compute the static, damage-independent road-link importance table once
    for the whole network (all links, not just any particular scenario's
    damaged links) -- analogous to S0's weighted_edge_betweenness. Callers
    (e.g. task_c_runner.py) should compute this ONCE and reuse it across all
    scenarios, since it does not depend on scenario damage."""
    road_net = read_road_network(task_a_cfg.road_net_path)
    bus_location = load_bus_location(od_cfg.bus_location_path)
    od_pairs = build_critical_od_pairs(od_cfg, bus_location=bus_location)
    _, link_to_services = route_overlap_structure(road_net, od_pairs, od_cfg)
    directed_importance = allocate_link_importance(link_to_services)
    return undirected_link_importance(directed_importance, list(road_net.costs.keys()))


def build_task_c_strategies(
    scenario: DisasterScenario,
    *,
    od_cfg: Optional[TaskCODConfig] = None,
    task_a_cfg: Optional[TaskACriticalityConfig] = None,
    static_road_importance: Optional[Dict[Link, float]] = None,
) -> List[TaskAStrategy]:
    """Build the S3 (O-D representation) road-repair strategy for a scenario.

    The road-side score comes from a FIXED critical-O-D route-overlap table
    (build_static_road_importance), looked up for whichever road links this
    scenario happens to have damaged -- exactly like S0 looks up its
    precomputed betweenness scores. The power-side score is unchanged from S0
    (power_radial_service_centrality), to isolate the road-side method change.

    Pass a precomputed `static_road_importance` (from
    build_static_road_importance) when scoring many scenarios in a loop, to
    avoid recomputing the same damage-independent table every time.

    Returns a single-element list so callers can treat this exactly like
    task_a_criticality.build_task_a_strategies' output.
    """
    od_cfg = od_cfg or load_taskc_od_config()
    task_a_cfg = task_a_cfg or TaskACriticalityConfig()

    if static_road_importance is None:
        static_road_importance = build_static_road_importance(od_cfg, task_a_cfg)

    power_assets = [int(b) for b in scenario.broken_buses]
    road_assets = [_normalize_link(link) for link in scenario.broken_links]

    road_scores: Dict[Asset, float] = {
        link: float(static_road_importance.get(link, static_road_importance.get((link[1], link[0]), 0.0)))
        for link in road_assets
    }
    power_scores = power_radial_service_centrality()
    s3_scores: Dict[Asset, float] = {
        **{bus: float(power_scores.get(bus, 0.0)) for bus in power_assets},
        **road_scores,
    }
    s3_sequence = _sorted_by_score(power_assets, s3_scores) + _sorted_by_score(road_assets, s3_scores)

    strategy = TaskAStrategy(
        strategy_id="S3_od_representation",
        strategy_label="S3 O-D based interdependence representation",
        score_method=(
            f"static critical-O-D route-overlap ({od_cfg.allocation_method}, k={od_cfg.k_paths}, "
            f"critical_buses={list(od_cfg.critical_buses)}, shelter_nodes={list(od_cfg.shelter_nodes)}) "
            "road importance + power radial downstream service impact"
        ),
        sequence=s3_sequence,
        power_sequence=_split_assets(s3_sequence)[0],
        road_sequence=_split_assets(s3_sequence)[1],
        scores=s3_scores,
    )
    return [strategy]
