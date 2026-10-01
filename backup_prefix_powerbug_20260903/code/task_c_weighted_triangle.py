from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Any, Dict, List, Sequence, Tuple

from power_util import get_functional_nodes
from road_util import calculate_shortest_path_cost

"""
Post-hoc, critical-weighted composite resilience-triangle metric.

This does NOT modify resilience_measurement.py or its triangle_area
computation (the "standard" triangle used as the primary efficiency metric
elsewhere in this dissertation). Instead it re-walks the event_log already
returned by run_model_multi() -- which every strategy's run already produces
-- and recomputes an additional composite metric that upweights a small,
designated set of "critical" power buses and adds a shelter-accessibility
term. See TASKC_ASSUMPTIONS.md for why this is a separate, supplemental
metric rather than a change to the primary one.
"""

BUS_COUNT = 33


@dataclass(frozen=True)
class CriticalWeightConfig:
    critical_buses: Tuple[int, ...] = (8, 17)  # same fixed set used for the S3 O-D destinations
    bus_weight_multiplier: float = 3.0
    shelter_node: int = 24
    depot_node: int = 1
    weight_road: float = 1.0
    weight_power: float = 1.0
    weight_shelter: float = 1.0


def _bus_weight_map(cfg: CriticalWeightConfig) -> Dict[int, float]:
    return {b: (cfg.bus_weight_multiplier if b in cfg.critical_buses else 1.0) for b in range(1, BUS_COUNT + 1)}


def weighted_power_func(broken_buses: Sequence[int], cfg: CriticalWeightConfig) -> float:
    """Same functional/broken-bus logic as resilience_measurement.eval_power_resilience
    (via the same get_functional_nodes helper), but each bus is counted by its
    weight instead of uniformly by 1/33."""
    functional = set(get_functional_nodes(set(int(b) for b in broken_buses)))
    weight_map = _bus_weight_map(cfg)
    total_weight = sum(weight_map.values())
    served_weight = sum(weight_map[b] for b in functional if b in weight_map)
    return served_weight / total_weight if total_weight > 0 else 0.0


def shelter_access_func(s_txt_path: str, *, cfg: CriticalWeightConfig, baseline_tt: float) -> float:
    tt = calculate_shortest_path_cost(s_txt_path, cfg.depot_node, cfg.shelter_node)
    return max(0.0, min(1.0, float(baseline_tt) / max(float(tt), 1e-9)))


def compute_weighted_triangle(
    run: Dict[str, Any],
    *,
    cfg: CriticalWeightConfig,
) -> Dict[str, float]:
    """Recompute a critical-weighted composite triangle from an existing
    run_model_multi() result. Reuses that run's own baseline_s.txt (written
    unconditionally by run_model_multi into run_dir) so the shelter-access
    baseline reflects the same undisrupted network state used elsewhere in
    that run, rather than a separately hardcoded constant."""
    baseline_s_txt = os.path.join(run["run_dir"], "baseline_s.txt")
    baseline_shelter_tt = calculate_shortest_path_cost(baseline_s_txt, cfg.depot_node, cfg.shelter_node)

    event_log: List[Dict[str, Any]] = run["event_log"]
    triangle = 0.0
    prev_time = None
    prev_terms = None

    road_terms: List[float] = []
    power_terms: List[float] = []
    shelter_terms: List[float] = []

    for entry in event_log:
        t = float(entry["time"])
        road_func = float(entry["state"]["road_func"])
        w_power = weighted_power_func(entry["broken_buses"], cfg)
        s_func = shelter_access_func(entry["state"]["s_txt_path"], cfg=cfg, baseline_tt=baseline_shelter_tt)

        road_terms.append(road_func)
        power_terms.append(w_power)
        shelter_terms.append(s_func)

        if prev_time is not None:
            dt = t - prev_time
            triangle += (
                cfg.weight_road * (1.0 - prev_terms[0])
                + cfg.weight_power * (1.0 - prev_terms[1])
                + cfg.weight_shelter * (1.0 - prev_terms[2])
            ) * dt

        prev_time = t
        prev_terms = (road_func, w_power, s_func)

    return {
        "weighted_triangle_area": triangle,
        "weighted_power_initial": power_terms[0] if power_terms else float("nan"),
        "weighted_power_final": power_terms[-1] if power_terms else float("nan"),
        "shelter_access_initial": shelter_terms[0] if shelter_terms else float("nan"),
        "shelter_access_final": shelter_terms[-1] if shelter_terms else float("nan"),
        "baseline_shelter_tt": float(baseline_shelter_tt),
    }
