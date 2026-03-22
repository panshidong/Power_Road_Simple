from __future__ import annotations

import json
import os
from typing import Any, Callable, Dict, Iterable, List, Optional, Tuple

from crews import Crew, CrewPool
from equity import all_metrics
from road_util import calculate_shortest_path_cost

Asset = Any  # int (bus) or (u,v) tuple for road link
TimelineEntry = Tuple[float, Asset]  # (finish_time, asset)


def _require_file(path: str, *, strict: bool) -> None:
    if strict and not os.path.exists(path):
        raise FileNotFoundError(f"Required file not found: {path!r}")


def _load_bus_to_link(path: str = "bus_to_link.json", *, strict: bool = False) -> Dict[int, Tuple[int, int]]:
    """
    Optional mapping (fallback only). In the new QA policy, bus dispatch does NOT require this file.
    """
    if not os.path.exists(path):
        if strict:
            raise FileNotFoundError(f"bus_to_link.json not found: {path!r}")
        return {}
    with open(path, "r", encoding="utf-8") as f:
        raw = json.load(f)
    out: Dict[int, Tuple[int, int]] = {}
    if isinstance(raw, dict):
        for k, v in raw.items():
            try:
                b = int(k)
            except Exception:
                continue
            if isinstance(v, (list, tuple)) and len(v) == 2:
                out[b] = (int(v[0]), int(v[1]))
    return out


def _load_bus_location(path: str = "bus_location.json", *, strict: bool = True) -> Dict[int, int]:
    """
    Required for bus dispatch in QA policy.
    Expects dict-like JSON: { "bus_id": node_id, ... }
    """
    if not os.path.exists(path):
        if strict:
            raise FileNotFoundError(f"bus_location.json not found: {path!r}")
        return {}
    with open(path, "r", encoding="utf-8") as f:
        raw = json.load(f)
    out: Dict[int, int] = {}
    if isinstance(raw, dict):
        for k, v in raw.items():
            out[int(k)] = int(v)
    return out


def _candidate_nodes_for_asset(
    asset: Asset,
    *,
    busloc: Dict[int, int],
    bus2link: Dict[int, Tuple[int, int]],
    strict: bool,
    bus_dispatch_mode: str = "location_first",
) -> Optional[Tuple[int, int]]:
    """
    Road link: (u,v)
    Bus dispatch modes:
      - "location_first": prefer bus_location -> (node,node), then bus_to_link -> (u,v)
      - "link_first": prefer bus_to_link -> (u,v), then bus_location -> (node,node)
      - "link_only": use bus_to_link only
    """
    if isinstance(asset, tuple) and len(asset) == 2:
        return (int(asset[0]), int(asset[1]))

    try:
        b = int(asset)
    except Exception:
        return None

    if bus_dispatch_mode not in {"location_first", "link_first", "link_only"}:
        raise ValueError(f"Unknown bus_dispatch_mode={bus_dispatch_mode!r}")

    if bus_dispatch_mode == "location_first":
        if b in busloc:
            n = int(busloc[b])
            return (n, n)
        if b in bus2link:
            return bus2link[b]
    elif bus_dispatch_mode == "link_first":
        if b in bus2link:
            return bus2link[b]
        if b in busloc:
            n = int(busloc[b])
            return (n, n)
    elif bus_dispatch_mode == "link_only":
        if b in bus2link:
            return bus2link[b]

    if strict:
        if bus_dispatch_mode == "link_only":
            raise KeyError(f"Bus {b} missing in bus_to_link. Provide new_bus_to_link.json / bus_to_link.json.")
        raise KeyError(f"Bus {b} missing in bus_location and bus_to_link. Provide the configured dispatch mapping.")
    return None


def travel_time_from_stxt(
    asset: Asset,
    crew: Crew,
    *,
    s_txt_path: str = "s.txt",
    bus_location_path: str = "bus_location.json",
    bus_to_link_path: str = "bus_to_link.json",
    default_depot: int = 1,
    strict: bool = True,
    bus_dispatch_mode: str = "location_first",
) -> float:
    """
    Dispatch travel time using TAP-B output (s.txt).

    Policy:
      - s.txt is required (strict=True)
      - bus_location.json is required only when the dispatch mode uses it
      - bus_to_link.json is required in link_only mode
    """
    _require_file(s_txt_path, strict=strict)

    need_busloc = bus_dispatch_mode in {"location_first", "link_first"}
    need_bus2link = bus_dispatch_mode in {"link_first", "link_only"}

    busloc = _load_bus_location(bus_location_path, strict=(strict and need_busloc))
    bus2link = _load_bus_to_link(bus_to_link_path, strict=(strict and bus_dispatch_mode == "link_only"))

    start = crew.location_node if crew.location_node is not None else int(default_depot)
    cand = _candidate_nodes_for_asset(
        asset,
        busloc=busloc,
        bus2link=bus2link,
        strict=strict,
        bus_dispatch_mode=bus_dispatch_mode,
    )
    if cand is None:
        # For strict=False only; treat as 0 travel
        return 0.0

    u, v = cand
    tu = calculate_shortest_path_cost(s_txt_path, start, u)
    tv = calculate_shortest_path_cost(s_txt_path, start, v)
    m = float(min(tu, tv))
    if strict and (m == float("inf") or m != m):
        raise RuntimeError(f"No finite path in s.txt for asset={asset!r} start={start} cand={(u,v)} tu={tu} tv={tv}")
    return m


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
    bus_location_path: str = "bus_location.json",
    bus_to_link_path: str = "bus_to_link.json",
    strict: bool = True,
    bus_dispatch_mode: str = "location_first",
) -> Dict[str, Any]:
    """
    Build completion timeline under multiple crews.

    Dispatch uses ONE s.txt snapshot for travel costs (IO-based).
    Bus dispatch can use bus_location.json, bus_to_link.json, or both.
    """
    if travel_time_fn is None:

        def travel_time_fn(a: Asset, c: Crew) -> float:
            return travel_time_from_stxt(
                a,
                c,
                s_txt_path=s_txt_path,
                bus_location_path=bus_location_path,
                bus_to_link_path=bus_to_link_path,
                default_depot=default_depot,
                strict=strict,
                bus_dispatch_mode=bus_dispatch_mode,
            )

    if strict:
        _require_file(s_txt_path, strict=True)
        if bus_dispatch_mode in {"location_first", "link_first"}:
            _require_file(bus_location_path, strict=True)
        if bus_dispatch_mode == "link_only":
            _require_file(bus_to_link_path, strict=True)

    pool = crew_pool.copy()
    finish_times: Dict[Asset, float] = {}
    timeline: List[TimelineEntry] = []

    busloc = _load_bus_location(bus_location_path, strict=(strict and bus_dispatch_mode in {"location_first", "link_first"}))
    bus2link = _load_bus_to_link(bus_to_link_path, strict=(strict and bus_dispatch_mode == "link_only"))

    for asset in sequence:
        crew = pool.next_available(asset)

        # dispatch
        t_travel_raw = float(travel_time_fn(asset, crew))
        t_travel = t_travel_raw / max(float(crew.speed), 1e-9)

        # arrival node update
        start = crew.location_node if crew.location_node is not None else int(default_depot)
        cand = _candidate_nodes_for_asset(
            asset,
            busloc=busloc,
            bus2link=bus2link,
            strict=strict,
            bus_dispatch_mode=bus_dispatch_mode,
        )

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
