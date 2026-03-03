from __future__ import annotations

import json
import os
from typing import Any, Callable, Dict, Iterable, List, Optional, Tuple

from crews import Crew, CrewPool
from equity import all_metrics
from road_util import calculate_shortest_path_cost

Asset = Any  # int (bus) or (u,v) tuple for road link
TimelineEntry = Tuple[float, Asset]  # (finish_time, asset)


def _require_file(path: str) -> None:
    if not os.path.exists(path):
        raise FileNotFoundError(f"Required file not found: {path!r}")


def _load_bus_to_link(path: str = "bus_to_link.json", *, strict: bool = True) -> Dict[int, Tuple[int, int]]:
    if not os.path.exists(path):
        if strict:
            raise FileNotFoundError(f"bus_to_link.json not found: {path!r}")
        return {}
    with open(path, "r", encoding="utf-8") as f:
        raw = json.load(f)

    out: Dict[int, Tuple[int, int]] = {}
    if isinstance(raw, dict):
        for k, v in raw.items():
            ki = int(k)
            if isinstance(v, (list, tuple)) and len(v) == 2:
                out[ki] = (int(v[0]), int(v[1]))
    return out


def _candidate_nodes_for_asset(asset: Asset, bus2link: Dict[int, Tuple[int, int]]) -> Optional[Tuple[int, int]]:
    if isinstance(asset, tuple) and len(asset) == 2:
        return (int(asset[0]), int(asset[1]))
    try:
        a = int(asset)
    except Exception:
        return None
    return bus2link.get(a, None)


def travel_time_from_stxt(
    asset: Asset,
    crew: Crew,
    *,
    s_txt_path: str = "s.txt",
    bus_to_link_path: str = "bus_to_link.json",
    default_depot: int = 1,
    strict: bool = True,
) -> float:
    """
    Dispatch travel time using TAP-B output (s.txt).

    strict=True policy:
      - missing s.txt -> raise
      - missing bus_to_link.json -> raise
      - bus without mapping -> raise (prevents silent 0 travel time bias)

    strict=False:
      - unmapped bus -> travel time 0.0
    """
    _require_file(s_txt_path)
    bus2link = _load_bus_to_link(bus_to_link_path, strict=strict)

    start = crew.location_node if crew.location_node is not None else int(default_depot)
    cand = _candidate_nodes_for_asset(asset, bus2link)
    if cand is None:
        if strict:
            raise KeyError(f"Asset {asset!r} has no (u,v) mapping; check {bus_to_link_path}.")
        return 0.0

    u, v = cand
    tu = calculate_shortest_path_cost(s_txt_path, start, u)
    tv = calculate_shortest_path_cost(s_txt_path, start, v)
    return float(min(tu, tv))


def default_service_time(asset: Asset, service_time_power: float, service_time_road: float) -> float:
    return float(service_time_road) if isinstance(asset, tuple) else float(service_time_power)


def evaluate_with_crews(
    *,
    sequence: List[Asset],
    crew_pool: CrewPool,
    service_time_power: float = 20.0,
    service_time_road: float = 10.0,
    travel_time_fn: Optional[Callable[[Asset, Crew], float]] = None,
    equity_groups: Optional[Iterable[str]] = None,  # reserved
    equity_svi: Optional[Iterable[float]] = None,
    default_depot: int = 1,
    s_txt_path: str = "s.txt",
    bus_to_link_path: str = "bus_to_link.json",
    strict: bool = True,
) -> Dict[str, Any]:
    """
    Build completion timeline under multiple crews.

    This uses ONE s.txt snapshot for dispatch costs (IO-based), by design.
    """
    if strict:
        _require_file(s_txt_path)
        _require_file(bus_to_link_path)

    if travel_time_fn is None:

        def travel_time_fn(a: Asset, c: Crew) -> float:
            return travel_time_from_stxt(
                a,
                c,
                s_txt_path=s_txt_path,
                bus_to_link_path=bus_to_link_path,
                default_depot=default_depot,
                strict=strict,
            )

    pool = crew_pool.copy()
    finish_times: Dict[Asset, float] = {}
    timeline: List[TimelineEntry] = []

    bus2link = _load_bus_to_link(bus_to_link_path, strict=(False if not strict else True))

    for asset in sequence:
        crew = pool.next_available(asset)

        # dispatch
        t_travel = float(travel_time_fn(asset, crew)) / max(float(crew.speed), 1e-9)

        # arrival node update (nearer endpoint)
        start = crew.location_node if crew.location_node is not None else int(default_depot)
        cand = _candidate_nodes_for_asset(asset, bus2link)
        if cand is None and strict:
            raise KeyError(f"Asset {asset!r} has no mapping for arrival node; check {bus_to_link_path}.")

        arrive_node: Optional[int] = None
        if cand is not None:
            u, v = cand
            tu = calculate_shortest_path_cost(s_txt_path, start, u)
            tv = calculate_shortest_path_cost(s_txt_path, start, v)
            arrive_node = u if tu <= tv else v

        service = default_service_time(asset, service_time_power, service_time_road)
        begin = max(float(crew.available_time), 0.0) + t_travel
        done = begin + float(service)

        finish_times[asset] = done
        timeline.append((done, asset))

        crew.available_time = done
        if arrive_node is not None:
            crew.location_node = arrive_node

    timeline.sort(key=lambda x: x[0])
    all_times = [finish_times[a] for a in sequence if a in finish_times]
    equity_result = all_metrics(all_times, equity_svi)

    return {"timeline": timeline, "finish_times": finish_times, "equity": equity_result}
