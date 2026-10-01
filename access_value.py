from __future__ import annotations

"""
Option B for Task A hypothesis 2: an ACCESS-AWARE characteristic function for the
Shapley criticality games.

The existing Task A value function (task_a_criticality.FullFunctionalityValue) is
static functionality: v(S) = w*P_power(S) + (1-w)*P_road(S) for the state in which
the coalition S of assets is repaired. It carries only the signal channel
(power -> road capacity) and nothing about crews having to reach damaged buses over
the damaged road network -- the interdependency channel that actually drives the
time-explicit simulator. This module adds that channel as a third term:

    A(S) = mean over all damaged buses b of clip( TT0_b / TT_b(S), 0, 1 )

where TT_b(S) is the shortest-path travel time from the crew depot to bus b's
dispatch target on the equilibrium travel times of road state S (same dispatch rule
as the simulator: bus_dispatch_mode="link_only", i.e. the nearer endpoint of the
bus's bus_to_link link), and TT0_b is the same quantity on the undisrupted network.
A(S) = 1 when every damaged bus is as reachable as it would be with no road damage.

    v_alpha(S) = w*P_power(S) + (1-w)*P_road(S) + alpha * A(S)

Both strategies built here keep Task A's execution unchanged (two specialized
crews, priority lists), so any improvement over S1/S2 comes only from the value
function. All simulator functions are imported unmodified.
"""

from dataclasses import dataclass
from typing import Any, Callable, Dict, List, Mapping, Sequence, Set, Tuple
import json

from disaster import DisasterScenario
from resilience_measurement import eval_power_resilience, eval_road_resilience
from road_util import calculate_shortest_path_cost
from task_a_criticality import (
    TaskACriticalityConfig,
    TaskAStrategy,
    _damaged_link_factors,
    _normalize_link,
    _pair_in,
    _sorted_by_score,
    _split_assets,
    sampled_shapley,
)

Asset = Any
Link = Tuple[int, int]


@dataclass
class AccessConfig:
    alpha: float = 0.5
    depot_node: int = 1
    bus_to_link_path: str = "new_bus_to_link.json"
    trips: str = "tap-b/net/SiouxFalls_trips.txt"
    net1: str = "work/SiouxFalls_net1.txt"
    net2: str = "work/SiouxFalls_net2.txt"


class AccessAwareValue:
    """Static functionality (as in Task A) plus the depot->damaged-bus reachability term."""

    def __init__(self, *, power_assets: Sequence[int], road_assets: Sequence[Link], link_factors: Mapping[Link, float],
                 cfg: TaskACriticalityConfig, acc: AccessConfig) -> None:
        self.power_assets = [int(b) for b in power_assets]
        self.road_assets = [_normalize_link(l) for l in road_assets]
        self.link_factors = {_normalize_link(l): float(f) for l, f in link_factors.items()}
        self.cfg = cfg
        self.acc = acc
        raw = json.load(open(acc.bus_to_link_path, encoding="utf-8"))
        self.bus2link = {int(k): (int(v[0]), int(v[1])) for k, v in raw.items()}
        self.power_cache: Dict[Tuple[int, ...], float] = {}
        self.road_cache: Dict[Tuple[Tuple[int, ...], Tuple[Link, ...]], Tuple[float, float]] = {}
        # undisrupted reference travel times to each damaged bus's dispatch target
        eval_road_resilience([], [], broken_link_factors={}, base_net=cfg.road_net_path, trips=acc.trips,
                             net1=acc.net1, net2=acc.net2, broken_link_factor=0.0,
                             power_road_factor=cfg.power_road_factor, baseline_tstt=cfg.baseline_tstt, strict=True)
        self.tt0 = {b: self._tt_to_bus(b, "s.txt") for b in self.power_assets}

    def _tt_to_bus(self, bus: int, s_txt: str) -> float:
        u, v = self.bus2link[int(bus)]
        return min(calculate_shortest_path_cost(s_txt, self.acc.depot_node, u),
                   calculate_shortest_path_cost(s_txt, self.acc.depot_node, v))

    def _remaining(self, repaired: Set[Asset]) -> Tuple[List[int], List[Link]]:
        rb = {int(a) for a in repaired if not isinstance(a, tuple)}
        rl = {_normalize_link(a) for a in repaired if isinstance(a, tuple)}
        return [b for b in self.power_assets if b not in rb], [l for l in self.road_assets if not _pair_in(rl, l)]

    def power(self, repaired: Set[Asset]) -> float:
        rem_b, _ = self._remaining(repaired)
        key = tuple(sorted(rem_b))
        if key not in self.power_cache:
            self.power_cache[key] = float(eval_power_resilience(list(rem_b)))
        return self.power_cache[key]

    def road_and_access(self, repaired: Set[Asset]) -> Tuple[float, float]:
        rem_b, rem_l = self._remaining(repaired)
        key = (tuple(sorted(rem_b)), tuple(sorted(rem_l)))
        if key not in self.road_cache:
            road_func = eval_road_resilience(list(rem_b), list(rem_l), broken_link_factors=dict(self.link_factors),
                                             base_net=self.cfg.road_net_path, trips=self.acc.trips, net1=self.acc.net1,
                                             net2=self.acc.net2, broken_link_factor=0.0,
                                             power_road_factor=self.cfg.power_road_factor,
                                             baseline_tstt=self.cfg.baseline_tstt, strict=True)
            # s.txt now holds the equilibrium for this state -> reachability of every damaged bus
            ratios = []
            for b in self.power_assets:
                tt = self._tt_to_bus(b, "s.txt")
                ratios.append(0.0 if tt == float("inf") or tt != tt else max(0.0, min(1.0, self.tt0[b] / max(tt, 1e-9))))
            access = sum(ratios) / len(ratios) if ratios else 1.0
            self.road_cache[key] = (float(road_func), float(access))
        return self.road_cache[key]

    # --- characteristic functions
    def joint_value(self, repaired: Set[Asset]) -> float:
        w = max(0.0, min(1.0, float(self.cfg.integrated_power_weight)))
        road, access = self.road_and_access(repaired)
        # Normalization leaves the ranking unchanged for a fixed alpha while
        # keeping the characteristic value on an interpretable unit scale.
        numerator = w * self.power(repaired) + (1.0 - w) * road + self.acc.alpha * access
        return numerator / (1.0 + self.acc.alpha)

    def road_value_with_access(self, repaired: Set[Asset]) -> float:
        road, access = self.road_and_access(repaired)
        return road + self.acc.alpha * access

    def power_value(self, repaired: Set[Asset]) -> float:
        return self.power(repaired)


def build_access_strategies(scenario: DisasterScenario, *, cfg: TaskACriticalityConfig, alphas: Sequence[float],
                            depot_node: int = 1) -> List[TaskAStrategy]:
    power_assets = [int(b) for b in scenario.broken_buses]
    road_assets = [_normalize_link(l) for l in scenario.broken_links]
    all_assets: List[Asset] = list(power_assets) + list(road_assets)
    link_factors = _damaged_link_factors(scenario)
    out: List[TaskAStrategy] = []
    for alpha in alphas:
        acc = AccessConfig(alpha=float(alpha), depot_node=depot_node)
        # separate games: power static (as S1), road static + access
        sep = AccessAwareValue(power_assets=power_assets, road_assets=road_assets, link_factors=link_factors, cfg=cfg, acc=acc)
        p_sh = sampled_shapley(power_assets, sep.power_value, samples=cfg.shapley_samples, seed=cfg.shapley_seed + scenario.seed + 101)
        # road sub-game evaluated with all power damage present (as S1's road game does for its own state)
        r_sh = sampled_shapley(road_assets, sep.road_value_with_access, samples=cfg.shapley_samples, seed=cfg.shapley_seed + scenario.seed + 202)
        s1b_scores = {**p_sh, **r_sh}
        s1b_seq = _sorted_by_score(power_assets, s1b_scores) + _sorted_by_score(road_assets, s1b_scores)
        tag = f"a{alpha:g}".replace(".", "p")
        out.append(TaskAStrategy(strategy_id=f"S1b_separate_access_{tag}", strategy_label=f"S1b separate Shapley + access (alpha={alpha:g})",
                                 score_method="separate sampled Shapley; road game = P_road + alpha*A(depot->bus reachability)",
                                 sequence=s1b_seq, power_sequence=_split_assets(s1b_seq)[0], road_sequence=_split_assets(s1b_seq)[1],
                                 scores=s1b_scores))
        # joint game with access term
        joint = AccessAwareValue(power_assets=power_assets, road_assets=road_assets, link_factors=link_factors, cfg=cfg, acc=acc)
        s2b_scores = sampled_shapley(all_assets, joint.joint_value, samples=cfg.shapley_samples, seed=cfg.shapley_seed + scenario.seed + 303)
        s2b_seq = _sorted_by_score(all_assets, s2b_scores)
        out.append(TaskAStrategy(strategy_id=f"S2b_joint_access_{tag}", strategy_label=f"S2b joint Shapley + access (alpha={alpha:g})",
                                 score_method="joint sampled Shapley on 0.5*P_power + 0.5*P_road + alpha*A(depot->bus reachability)",
                                 sequence=s2b_seq, power_sequence=_split_assets(s2b_seq)[0], road_sequence=_split_assets(s2b_seq)[1],
                                 scores=s2b_scores, preserve_sequence_order=True))
    return out
