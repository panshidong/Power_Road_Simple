from __future__ import annotations

import heapq
import itertools
import random
import shutil
from collections import deque
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, Iterable, List, Mapping, Sequence, Set, Tuple

from disaster import DisasterScenario
from power_util import delete_buses
from resilience_measurement import eval_power_resilience, eval_road_resilience

Asset = Any
Link = Tuple[int, int]


POWER_CONNECTIONS: Dict[int, List[int]] = {
    1: [2],
    2: [3, 19],
    3: [4, 23],
    4: [5],
    5: [6],
    6: [7, 26],
    7: [8],
    8: [9, 21],
    9: [10],
    10: [11],
    11: [12],
    12: [13, 22],
    13: [14],
    14: [15],
    15: [16],
    16: [17],
    17: [18],
    18: [33],
    19: [20],
    20: [21],
    21: [],
    22: [],
    23: [24],
    24: [25],
    25: [29],
    26: [27],
    27: [28],
    28: [29],
    29: [30],
    30: [31],
    31: [32],
    32: [],
    33: [],
}


@dataclass(frozen=True)
class RoadNetwork:
    nodes: Tuple[int, ...]
    costs: Dict[Link, float]


@dataclass
class TaskACriticalityConfig:
    road_net_path: str = "tap-b/net/SiouxFalls_net.txt"
    bus_to_link_path: str = "new_bus_to_link.json"
    shapley_samples: int = 120
    shapley_seed: int = 20260402
    integrated_power_weight: float = 0.5
    use_full_functionality_shapley: bool = True
    power_road_factor: float = 0.5
    baseline_tstt: float = 7475338.0

    def __post_init__(self) -> None:
        if not self.use_full_functionality_shapley:
            raise ValueError("Task A Shapley must use simulator functionality; proxy Shapley is disabled.")


@dataclass
class TaskAStrategy:
    strategy_id: str
    strategy_label: str
    score_method: str
    sequence: List[Asset]
    power_sequence: List[int]
    road_sequence: List[Link]
    scores: Dict[Asset, float] = field(default_factory=dict)
    preserve_sequence_order: bool = True


def asset_type(asset: Asset) -> str:
    return "road" if isinstance(asset, tuple) else "power"


def asset_key(asset: Asset) -> str:
    if isinstance(asset, tuple):
        return f"road:{int(asset[0])}-{int(asset[1])}"
    return f"power:{int(asset)}"


def _asset_sort_key(asset: Asset) -> Tuple[int, int, int]:
    if isinstance(asset, tuple):
        return (1, int(asset[0]), int(asset[1]))
    return (0, int(asset), -1)


def _asset_from_json_like(value: Any) -> Asset:
    if isinstance(value, tuple):
        return (int(value[0]), int(value[1]))
    if isinstance(value, list) and len(value) == 2:
        return (int(value[0]), int(value[1]))
    return int(value)


def _normalize_link(link: Link) -> Link:
    return (int(link[0]), int(link[1]))


def _pair_in(links: Set[Link], link: Link) -> bool:
    u, v = _normalize_link(link)
    return (u, v) in links or (v, u) in links


def _sorted_by_score(assets: Iterable[Asset], scores: Mapping[Asset, float]) -> List[Asset]:
    return sorted(
        list(assets),
        key=lambda asset: (-float(scores.get(asset, 0.0)), _asset_sort_key(asset)),
    )


def _split_assets(sequence: Sequence[Asset]) -> Tuple[List[int], List[Link]]:
    power = [int(asset) for asset in sequence if not isinstance(asset, tuple)]
    road = [_normalize_link(asset) for asset in sequence if isinstance(asset, tuple)]
    return power, road


def _directed_adjacency(edges: Mapping[Link, float]) -> Dict[int, List[int]]:
    graph: Dict[int, List[int]] = {}
    for (u, v) in edges:
        graph.setdefault(int(u), []).append(int(v))
        graph.setdefault(int(v), graph.get(int(v), []))
    return {node: sorted(neighbors) for node, neighbors in graph.items()}


def _weighted_adjacency(edges: Mapping[Link, float]) -> Dict[int, List[Tuple[int, float]]]:
    graph: Dict[int, List[Tuple[int, float]]] = {}
    for (u, v), cost in edges.items():
        graph.setdefault(int(u), []).append((int(v), max(float(cost), 1e-9)))
        graph.setdefault(int(v), graph.get(int(v), []))
    return {node: sorted(neighbors) for node, neighbors in graph.items()}


def read_road_network(path: str) -> RoadNetwork:
    costs: Dict[Link, float] = {}
    nodes: Set[int] = set()
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("~") or s.startswith(";") or s.startswith("<"):
                continue
            parts = s.split()
            if len(parts) < 5:
                continue
            try:
                u = int(parts[0])
                v = int(parts[1])
                cost = float(parts[4])
            except Exception:
                continue
            costs[(u, v)] = max(cost, 1e-9)
            nodes.add(u)
            nodes.add(v)
    return RoadNetwork(nodes=tuple(sorted(nodes)), costs=costs)


def node_betweenness(graph: Mapping[int, Sequence[int]]) -> Dict[int, float]:
    nodes = sorted(int(n) for n in graph.keys())
    scores = {node: 0.0 for node in nodes}
    for source in nodes:
        stack: List[int] = []
        pred: Dict[int, List[int]] = {node: [] for node in nodes}
        sigma = {node: 0.0 for node in nodes}
        dist = {node: -1 for node in nodes}
        sigma[source] = 1.0
        dist[source] = 0
        queue: deque[int] = deque([source])
        while queue:
            v = queue.popleft()
            stack.append(v)
            for w in graph.get(v, []):
                if dist[w] < 0:
                    queue.append(w)
                    dist[w] = dist[v] + 1
                if dist[w] == dist[v] + 1:
                    sigma[w] += sigma[v]
                    pred[w].append(v)
        delta = {node: 0.0 for node in nodes}
        while stack:
            w = stack.pop()
            for v in pred[w]:
                if sigma[w] > 0:
                    delta[v] += (sigma[v] / sigma[w]) * (1.0 + delta[w])
            if w != source:
                scores[w] += delta[w]
    return scores


def edge_betweenness(graph: Mapping[int, Sequence[int]]) -> Dict[Link, float]:
    nodes = sorted(int(n) for n in graph.keys())
    scores: Dict[Link, float] = {(int(u), int(v)): 0.0 for u, vs in graph.items() for v in vs}
    for source in nodes:
        stack: List[int] = []
        pred: Dict[int, List[int]] = {node: [] for node in nodes}
        sigma = {node: 0.0 for node in nodes}
        dist = {node: -1 for node in nodes}
        sigma[source] = 1.0
        dist[source] = 0
        queue: deque[int] = deque([source])
        while queue:
            v = queue.popleft()
            stack.append(v)
            for w in graph.get(v, []):
                if dist[w] < 0:
                    queue.append(w)
                    dist[w] = dist[v] + 1
                if dist[w] == dist[v] + 1:
                    sigma[w] += sigma[v]
                    pred[w].append(v)
        delta = {node: 0.0 for node in nodes}
        while stack:
            w = stack.pop()
            for v in pred[w]:
                if sigma[w] > 0:
                    contrib = (sigma[v] / sigma[w]) * (1.0 + delta[w])
                    scores[(v, w)] = scores.get((v, w), 0.0) + contrib
                    delta[v] += contrib
    return scores


def weighted_edge_betweenness(edges: Mapping[Link, float]) -> Dict[Link, float]:
    graph = _weighted_adjacency(edges)
    nodes = sorted(int(n) for n in graph.keys())
    scores: Dict[Link, float] = {_normalize_link(link): 0.0 for link in edges}
    tol = 1e-12

    for source in nodes:
        stack: List[int] = []
        pred: Dict[int, List[int]] = {node: [] for node in nodes}
        sigma = {node: 0.0 for node in nodes}
        dist = {node: float("inf") for node in nodes}
        sigma[source] = 1.0
        dist[source] = 0.0
        pq: List[Tuple[float, int]] = [(0.0, source)]

        while pq:
            d, v = heapq.heappop(pq)
            if d > dist[v] + tol:
                continue
            stack.append(v)
            for w, weight in graph.get(v, []):
                vw_dist = dist[v] + float(weight)
                if vw_dist < dist[w] - tol:
                    dist[w] = vw_dist
                    heapq.heappush(pq, (vw_dist, w))
                    sigma[w] = sigma[v]
                    pred[w] = [v]
                elif abs(vw_dist - dist[w]) <= tol:
                    sigma[w] += sigma[v]
                    pred[w].append(v)

        delta = {node: 0.0 for node in nodes}
        while stack:
            w = stack.pop()
            for v in pred[w]:
                if sigma[w] > 0:
                    contrib = (sigma[v] / sigma[w]) * (1.0 + delta[w])
                    scores[(v, w)] = scores.get((v, w), 0.0) + contrib
                    delta[v] += contrib
    return scores


def power_radial_service_centrality() -> Dict[int, float]:
    scores: Dict[int, float] = {}
    for bus in sorted(POWER_CONNECTIONS):
        unavailable = set(map(int, delete_buses([bus])))
        scores[int(bus)] = float(len(unavailable)) / 33.0
    return scores


def centrality_scores(road_net: RoadNetwork) -> Tuple[Dict[int, float], Dict[Link, float]]:
    power_scores = power_radial_service_centrality()
    road_scores_directed = weighted_edge_betweenness(road_net.costs)
    road_scores: Dict[Link, float] = {}
    for link in road_net.costs:
        u, v = link
        road_scores[(u, v)] = max(
            float(road_scores_directed.get((u, v), 0.0)),
            float(road_scores_directed.get((v, u), 0.0)),
        )
    return power_scores, road_scores


def sampled_shapley(
    assets: Sequence[Asset],
    value_fn: Callable[[Set[Asset]], float],
    *,
    samples: int,
    seed: int,
) -> Dict[Asset, float]:
    players = [_asset_from_json_like(asset) for asset in assets]
    if not players:
        return {}
    if int(samples) <= 0:
        raise ValueError("samples must be positive for sampled Shapley")
    rng = random.Random(int(seed))
    scores = {asset: 0.0 for asset in players}
    cache: Dict[frozenset[Asset], float] = {}

    def cached_value(coalition: Set[Asset]) -> float:
        key = frozenset(coalition)
        if key not in cache:
            cache[key] = float(value_fn(set(coalition)))
        return cache[key]

    for _ in range(int(samples)):
        order = list(players)
        rng.shuffle(order)
        coalition: Set[Asset] = set()
        prev = cached_value(coalition)
        for asset in order:
            coalition.add(asset)
            curr = cached_value(coalition)
            scores[asset] += curr - prev
            prev = curr
    denom = max(int(samples), 1)
    return {asset: value / denom for asset, value in scores.items()}


def exact_shapley_by_permutation(
    assets: Sequence[Asset],
    value_fn: Callable[[Set[Asset]], float],
    *,
    max_assets: int = 8,
) -> Dict[Asset, float]:
    players = [_asset_from_json_like(asset) for asset in assets]
    if len(players) > int(max_assets):
        raise ValueError(f"exact Shapley is limited to {max_assets} assets; got {len(players)}")
    if not players:
        return {}
    scores = {asset: 0.0 for asset in players}
    n_perm = 0
    cache: Dict[frozenset[Asset], float] = {}

    def cached_value(coalition: Set[Asset]) -> float:
        key = frozenset(coalition)
        if key not in cache:
            cache[key] = float(value_fn(set(coalition)))
        return cache[key]

    for order in itertools.permutations(players):
        n_perm += 1
        coalition: Set[Asset] = set()
        prev = cached_value(coalition)
        for asset in order:
            coalition.add(asset)
            curr = cached_value(coalition)
            scores[asset] += curr - prev
            prev = curr
    return {asset: value / float(n_perm) for asset, value in scores.items()}


def shapley_efficiency_gap(
    assets: Sequence[Asset],
    scores: Mapping[Asset, float],
    value_fn: Callable[[Set[Asset]], float],
) -> float:
    players = [_asset_from_json_like(asset) for asset in assets]
    full = float(value_fn(set(players)))
    empty = float(value_fn(set()))
    return abs(sum(float(scores.get(asset, 0.0)) for asset in players) - (full - empty))


def _damaged_link_factors(scenario: DisasterScenario) -> Dict[Link, float]:
    out: Dict[Link, float] = {}
    for link in scenario.broken_links:
        norm = _normalize_link(link)
        out[norm] = float(scenario.link_capacity_factors.get(norm, scenario.link_capacity_factors.get((norm[1], norm[0]), 0.0)))
    return out


class FullFunctionalityValue:
    def __init__(
        self,
        *,
        power_assets: Sequence[int],
        road_assets: Sequence[Link],
        link_factors: Mapping[Link, float],
        cfg: TaskACriticalityConfig,
    ) -> None:
        self.power_assets = [int(asset) for asset in power_assets]
        self.road_assets = [_normalize_link(link) for link in road_assets]
        self.link_factors = {_normalize_link(link): float(factor) for link, factor in link_factors.items()}
        self.cfg = cfg
        self.power_cache: Dict[Tuple[int, ...], float] = {}
        self.road_cache: Dict[Tuple[Tuple[int, ...], Tuple[Link, ...]], float] = {}

        if self.cfg.bus_to_link_path != "bus_to_link.json":
            shutil.copy2(self.cfg.bus_to_link_path, "bus_to_link.json")

    def _remaining_state(self, repaired: Set[Asset]) -> Tuple[List[int], List[Link]]:
        repaired_buses = {int(asset) for asset in repaired if not isinstance(asset, tuple)}
        repaired_links = {_normalize_link(asset) for asset in repaired if isinstance(asset, tuple)}
        remaining_buses = [bus for bus in self.power_assets if bus not in repaired_buses]
        remaining_links = [link for link in self.road_assets if not _pair_in(repaired_links, link)]
        return remaining_buses, remaining_links

    def power_functionality(self, repaired: Set[Asset]) -> float:
        remaining_buses, remaining_links = self._remaining_state(repaired)
        del remaining_links
        key = tuple(sorted(int(bus) for bus in remaining_buses))
        if key not in self.power_cache:
            self.power_cache[key] = float(eval_power_resilience(list(remaining_buses)))
        return self.power_cache[key]

    def road_functionality(self, repaired: Set[Asset]) -> float:
        remaining_buses, remaining_links = self._remaining_state(repaired)
        key = (
            tuple(sorted(int(bus) for bus in remaining_buses)),
            tuple(sorted(_normalize_link(link) for link in remaining_links)),
        )
        if key not in self.road_cache:
            road_func = eval_road_resilience(
                list(remaining_buses),
                list(remaining_links),
                broken_link_factors=dict(self.link_factors),
                base_net=self.cfg.road_net_path,
                trips="tap-b/net/SiouxFalls_trips.txt",
                net1="work/SiouxFalls_net1.txt",
                net2="work/SiouxFalls_net2.txt",
                broken_link_factor=0.0,
                power_road_factor=self.cfg.power_road_factor,
                baseline_tstt=self.cfg.baseline_tstt,
                strict=True,
            )
            self.road_cache[key] = float(road_func)
        return self.road_cache[key]

    def state_functionality(self, repaired: Set[Asset]) -> Tuple[float, float]:
        return self.power_functionality(repaired), self.road_functionality(repaired)

    def power_value(self, repaired: Set[Asset]) -> float:
        return float(self.power_functionality(repaired))

    def road_value(self, repaired: Set[Asset]) -> float:
        return float(self.road_functionality(repaired))

    def integrated_value(self, repaired: Set[Asset]) -> float:
        power_func, road_func = self.state_functionality(repaired)
        w_power = max(0.0, min(1.0, float(self.cfg.integrated_power_weight)))
        return w_power * float(power_func) + (1.0 - w_power) * float(road_func)


def build_task_a_strategies(
    scenario: DisasterScenario,
    *,
    cfg: TaskACriticalityConfig | None = None,
) -> List[TaskAStrategy]:
    cfg = cfg or TaskACriticalityConfig()
    road_net = read_road_network(cfg.road_net_path)

    power_assets = [int(asset) for asset in scenario.broken_buses]
    road_assets = [_normalize_link(link) for link in scenario.broken_links]
    all_assets: List[Asset] = list(power_assets) + list(road_assets)
    link_factors = _damaged_link_factors(scenario)
    power_full_value = FullFunctionalityValue(
        power_assets=power_assets,
        road_assets=[],
        link_factors={},
        cfg=cfg,
    )
    road_full_value = FullFunctionalityValue(
        power_assets=[],
        road_assets=road_assets,
        link_factors=link_factors,
        cfg=cfg,
    )
    integrated_full_value = FullFunctionalityValue(
        power_assets=power_assets,
        road_assets=road_assets,
        link_factors=link_factors,
        cfg=cfg,
    )

    power_cent, road_cent = centrality_scores(road_net)
    s0_scores: Dict[Asset, float] = {
        **{bus: float(power_cent.get(bus, 0.0)) for bus in power_assets},
        **{link: float(road_cent.get(link, road_cent.get((link[1], link[0]), 0.0))) for link in road_assets},
    }
    s0_sequence = _sorted_by_score(power_assets, s0_scores) + _sorted_by_score(road_assets, s0_scores)

    def separate_power_value(repaired: Set[Asset]) -> float:
        return power_full_value.power_value(repaired)

    def separate_road_value(repaired: Set[Asset]) -> float:
        return road_full_value.road_value(repaired)

    power_shapley = sampled_shapley(
        power_assets,
        separate_power_value,
        samples=cfg.shapley_samples,
        seed=cfg.shapley_seed + scenario.seed + 101,
    )
    road_shapley = sampled_shapley(
        road_assets,
        separate_road_value,
        samples=cfg.shapley_samples,
        seed=cfg.shapley_seed + scenario.seed + 202,
    )
    s1_scores: Dict[Asset, float] = {**power_shapley, **road_shapley}
    s1_sequence = _sorted_by_score(power_assets, s1_scores) + _sorted_by_score(road_assets, s1_scores)

    def integrated_value(repaired: Set[Asset]) -> float:
        return integrated_full_value.integrated_value(repaired)

    s2_scores = sampled_shapley(
        all_assets,
        integrated_value,
        samples=cfg.shapley_samples,
        seed=cfg.shapley_seed + scenario.seed + 303,
    )
    s2_sequence = _sorted_by_score(all_assets, s2_scores)

    strategies = [
        TaskAStrategy(
            strategy_id="S0_centrality",
            strategy_label="S0 separate network centrality",
            score_method="power radial downstream service impact + road weighted edge betweenness",
            sequence=s0_sequence,
            power_sequence=_split_assets(s0_sequence)[0],
            road_sequence=_split_assets(s0_sequence)[1],
            scores=s0_scores,
        ),
        TaskAStrategy(
            strategy_id="S1_separate_shapley",
            strategy_label="S1 separate network sampled Shapley",
            score_method="separate sampled Shapley on full-state power/road functionality",
            sequence=s1_sequence,
            power_sequence=_split_assets(s1_sequence)[0],
            road_sequence=_split_assets(s1_sequence)[1],
            scores=s1_scores,
        ),
        TaskAStrategy(
            strategy_id="S2_integrated_shapley",
            strategy_label="S2 integrated sampled Shapley",
            score_method="sampled Shapley on integrated full-state power-road functionality",
            sequence=s2_sequence,
            power_sequence=_split_assets(s2_sequence)[0],
            road_sequence=_split_assets(s2_sequence)[1],
            scores=s2_scores,
            preserve_sequence_order=True,
        ),
    ]
    return strategies


def strategy_ranking_rows(scenario: DisasterScenario, strategy: TaskAStrategy) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for rank, asset in enumerate(strategy.sequence, start=1):
        rows.append(
            {
                "scenario_id": scenario.scenario_id,
                "scenario_seed": scenario.seed,
                "strategy_id": strategy.strategy_id,
                "strategy_label": strategy.strategy_label,
                "score_method": strategy.score_method,
                "rank": rank,
                "asset": asset_key(asset),
                "asset_type": asset_type(asset),
                "score": float(strategy.scores.get(asset, 0.0)),
            }
        )
    return rows
