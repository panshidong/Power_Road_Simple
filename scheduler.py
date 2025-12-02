
from __future__ import annotations
from typing import Any, Callable, Dict, Iterable, List, Optional, Tuple
import json
import os

from crews import CrewPool, Crew
from equity import all_metrics
from road_util import calculate_shortest_path_cost

# Types
Asset = Any  # int (bus) or (u,v) tuple for road link
TimelineEntry = Tuple[float, Asset]  # (finish_time, asset)

def _load_bus_to_link(path: str = "bus_to_link.json") -> Dict[int, Tuple[int, int]]:
    """Load optional bus->(u,v) mapping.
    Returns an empty dict if file is missing or malformed (we'll treat unmapped buses as having no travel target).
    """
    if not os.path.exists(path):
        return {}
    try:
        with open(path, "r", encoding="utf-8") as f:
            raw = json.load(f)
    except Exception:
        return {}
    out: Dict[int, Tuple[int, int]] = {}
    if isinstance(raw, dict):
        for k, v in raw.items():
            try:
                ki = int(k)
            except Exception:
                continue
            if isinstance(v, (list, tuple)) and len(v) == 2:
                try:
                    out[ki] = (int(v[0]), int(v[1]))
                except Exception:
                    continue
    return out

def _candidate_nodes_for_asset(asset: Asset, bus2link: Dict[int, Tuple[int, int]]) -> Optional[Tuple[int, int]]:
    """Return (u,v) candidate nodes for an asset, or None if no association exists.
    - Road link: returns (u,v).
    - Bus with mapping: returns its mapped (u,v).
    - Bus without mapping: returns None (no associated road nodes).
    """
    if isinstance(asset, tuple) and len(asset) == 2:
        u, v = int(asset[0]), int(asset[1])
        return (u, v)
    # asset is a bus id (int)
    try:
        a = int(asset)
    except Exception:
        return None
    return bus2link.get(a, None)

def travel_time_from_stxt(asset: Asset, crew: Crew, s_txt_path: str = "s.txt",
                          bus_to_link_path: str = "bus_to_link.json",
                          default_depot: int = 1) -> float:
    """Dispatch travel time using TAP-B output (s.txt).
    If an asset has no associated road endpoints (e.g., an unmapped bus), returns 0.0 and does not change crew location.
    """
    start = crew.location_node if crew.location_node is not None else default_depot
    bus2link = _load_bus_to_link(bus_to_link_path)
    cand = _candidate_nodes_for_asset(asset, bus2link)
    if not cand:
        # No associated nodes — assume work is performed in place; 0 dispatch time.
        return 0.0
    u, v = cand
    tu = calculate_shortest_path_cost(s_txt_path, start, u)
    tv = calculate_shortest_path_cost(s_txt_path, start, v)
    return min(tu, tv)

def default_service_time(asset: Asset, service_time_power: float, service_time_road: float) -> float:
    return service_time_road if isinstance(asset, tuple) else service_time_power

def evaluate_with_crews(
    sequence: List[Asset],
    crew_pool: CrewPool,
    service_time_power: float = 20.0,
    service_time_road: float = 10.0,
    travel_time_fn: Optional[Callable[[Asset, Crew], float]] = None,
    equity_groups: Optional[Iterable[str]] = None,
    equity_svi: Optional[Iterable[float]] = None,
    default_depot: int = 1,
) -> Dict[str, Any]:
    if travel_time_fn is None:
        def travel_time_fn(asset, crew):
            return travel_time_from_stxt(asset, crew, s_txt_path="s.txt", bus_to_link_path="bus_to_link.json", default_depot=default_depot)

    pool = crew_pool.copy()

    finish_times: Dict[Asset, float] = {}
    timeline: List[TimelineEntry] = []

    bus2link = _load_bus_to_link("bus_to_link.json")

    for asset in sequence:
        crew = pool.next_available(asset)

        # Dispatch time (0 if asset has no associated road mapping)
        t_travel = travel_time_fn(asset, crew) / max(crew.speed, 1e-6)

        # Decide arrival node for the crew's next origin.
        # If the asset has candidate nodes, arrive to the nearer endpoint; otherwise keep current location.
        start = crew.location_node if crew.location_node is not None else default_depot
        cand = _candidate_nodes_for_asset(asset, bus2link)
        arrive_node = None
        if cand:
            u, v = cand
            tu = calculate_shortest_path_cost("s.txt", start, u)
            tv = calculate_shortest_path_cost("s.txt", start, v)
            arrive_node = u if tu <= tv else v

        service = default_service_time(asset, service_time_power, service_time_road)
        begin = max(crew.available_time, 0.0) + t_travel
        done = begin + service

        finish_times[asset] = done
        timeline.append((done, asset))

        crew.available_time = done
        if arrive_node is not None:
            crew.location_node = arrive_node  # only update location if we have a valid node

    timeline.sort(key=lambda x: x[0])

    all_times = [finish_times[a] for a in sequence if a in finish_times]
    equity_result = all_metrics(all_times, equity_svi)

    return {
        "timeline": timeline,
        "finish_times": finish_times,
        "equity": equity_result,
    }
