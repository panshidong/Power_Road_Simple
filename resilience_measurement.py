from __future__ import annotations

import json
import os
import shutil
import textwrap
from datetime import datetime
from typing import Any, Dict, List, Tuple, Optional

from crews import CrewPool, make_crews
from scheduler import evaluate_with_crews
from run_tapb import run_tapb
from road_util import capacity_adjustment, eval_tot_OD_travel_time
from power_util import delete_buses, get_functional_nodes
from interdependency import power_to_road

from simulated_annealing import SAConfig, simulated_annealing
from cri import (
    load_json,
    build_zone_lists,
    compute_accessibility_TT,
    compute_accessibility_ratio_by_zone,
    compute_E_by_zone,
    compute_CRI_by_zone,
    critical_access_summary_from_series,
    equity_summary_from_CRI,
)

BUS_COUNT = 33

DEFAULT_BASE_NET = "tap-b/net/SiouxFalls_net.txt"
DEFAULT_TRIPS = "tap-b/net/SiouxFalls_trips.txt"

DEFAULT_WORK_DIR = "work"
DEFAULT_NET1 = os.path.join(DEFAULT_WORK_DIR, "SiouxFalls_net1.txt")
DEFAULT_NET2 = os.path.join(DEFAULT_WORK_DIR, "SiouxFalls_net2.txt")
DEFAULT_CRITICAL_LOCATION_PATH = "critical_location.json"
DEFAULT_BUS_LOCATION_SOURCE = "original_bus_location.json"
DEFAULT_BUS_TO_LINK_SOURCE = "new_bus_to_link.json"


def _ts() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def _require(path: str) -> None:
    if not os.path.exists(path):
        raise FileNotFoundError(f"Required file not found: {path!r}")


def _ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def _fmt(x: float, digits: int = 4) -> str:
    try:
        x = float(x)
    except Exception:
        return str(x)
    if x == 0:
        return "0"
    return f"{x:.{digits}g}"


def _wrap_title(title: str, width: int = 72) -> str:
    lines: List[str] = []
    for chunk in str(title).split("\n"):
        lines.extend(textwrap.wrap(chunk, width=width) or [""])
    return "\n".join(lines)


def _split_sequence_by_skill(sequence: List[Any]) -> Tuple[List[Any], List[Any]]:
    power_seq = [asset for asset in sequence if not isinstance(asset, tuple)]
    road_seq = [asset for asset in sequence if isinstance(asset, tuple)]
    return power_seq, road_seq


def _canonicalize_sequence_for_crews(sequence: List[Any], *, crew_mode: str) -> Tuple[List[Any], List[Any], List[Any]]:
    seq = list(sequence)
    power_seq, road_seq = _split_sequence_by_skill(seq)
    if crew_mode == "specialized":
        seq = list(power_seq) + list(road_seq)
    return seq, power_seq, road_seq


def _specialized_neighbor_factory(*, n_power: int, neighbor_mode: str):
    def _neighbor(seq: List[Any], rng) -> List[Any]:
        blocks: List[Tuple[int, int]] = []
        if n_power >= 2:
            blocks.append((0, n_power))
        n_road = len(seq) - n_power
        if n_road >= 2:
            blocks.append((n_power, len(seq)))
        if not blocks:
            return list(seq)

        start, end = blocks[0] if len(blocks) == 1 else rng.choice(blocks)
        cand = list(seq)
        if neighbor_mode == "insert":
            i, j = rng.sample(range(start, end), 2)
            x = cand.pop(i)
            cand.insert(j, x)
        else:
            i, j = rng.sample(range(start, end), 2)
            cand[i], cand[j] = cand[j], cand[i]
        return cand

    return _neighbor


def eval_power_resilience(broken_buses: List[int]) -> float:
    functional = set(get_functional_nodes(set(map(int, broken_buses))))
    return float(len(functional)) / float(BUS_COUNT)


def _prepare_state_and_run_tapb(
    *,
    broken_buses: List[int],
    broken_links: List[Tuple[int, int]],
    broken_link_factors: Optional[Dict[Tuple[int, int], float]],
    base_net: str,
    trips: str,
    net1: str,
    net2: str,
    broken_link_factor: float,
    power_road_factor: float,
    strict: bool,
) -> None:
    _require(base_net)
    _require(trips)
    _ensure_dir(os.path.dirname(net1) or ".")

    capacity_adjustment(base_net, net1, broken_links, broken_link_factor, link_factors=broken_link_factors)
    _require(net1)

    unfunctional_nodes = delete_buses(list(map(int, broken_buses)))
    if strict and unfunctional_nodes is None:
        raise RuntimeError("delete_buses returned None; check power_util.delete_buses implementation.")

    power_to_road(unfunctional_nodes, net1, net2, power_road_factor)
    _require(net2)

    if os.path.exists("s.txt"):
        os.remove("s.txt")
    run_tapb(net2, trips, check=True)
    if strict:
        _require("s.txt")


def eval_road_resilience(
    broken_buses: List[int],
    broken_links: List[Tuple[int, int]],
    *,
    broken_link_factors: Optional[Dict[Tuple[int, int], float]],
    base_net: str,
    trips: str,
    net1: str,
    net2: str,
    broken_link_factor: float,
    power_road_factor: float,
    baseline_tstt: float,
    strict: bool,
) -> float:
    _prepare_state_and_run_tapb(
        broken_buses=broken_buses,
        broken_links=broken_links,
        broken_link_factors=broken_link_factors,
        base_net=base_net,
        trips=trips,
        net1=net1,
        net2=net2,
        broken_link_factor=broken_link_factor,
        power_road_factor=power_road_factor,
        strict=strict,
    )
    current_tstt = float(eval_tot_OD_travel_time("s.txt"))
    if strict and current_tstt <= 0:
        raise RuntimeError(f"current_tstt <= 0 from s.txt; got {current_tstt}")
    func = float(baseline_tstt) / max(current_tstt, 1e-9)
    return max(0.0, min(1.0, func))


def _plot_triangle_with_shading(
    *,
    time_series: List[float],
    road_series: List[float],
    power_series: List[float],
    out_png: str,
    title: str,
) -> str:
    try:
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        out_svg = os.path.splitext(out_png)[0] + ".svg"
        _write_triangle_svg(
            time_series=time_series,
            road_series=road_series,
            power_series=power_series,
            out_svg=out_svg,
            title=title,
        )
        return out_svg

    if len(time_series) != len(road_series) or len(time_series) != len(power_series):
        raise ValueError("time_series/road_series/power_series length mismatch")

    fig, ax = plt.subplots(figsize=(11, 7), constrained_layout=True)
    ax.plot(time_series, road_series, label="Road functionality", linewidth=2.0)
    ax.plot(time_series, power_series, label="Power functionality", linewidth=2.0)
    ones = [1.0 for _ in time_series]
    ax.fill_between(time_series, road_series, ones, alpha=0.2, label="Road complement")
    ax.fill_between(time_series, power_series, ones, alpha=0.2, label="Power complement")
    ax.set_ylim(0.0, 1.05)
    ax.set_xlabel("Time")
    ax.set_ylabel("Functionality")
    ax.set_title(_wrap_title(title), pad=14)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=2, frameon=True)
    fig.savefig(out_png, dpi=200, bbox_inches="tight", pad_inches=0.3)
    plt.close(fig)
    return out_png


def _write_triangle_svg(
    *,
    time_series: List[float],
    road_series: List[float],
    power_series: List[float],
    out_svg: str,
    title: str,
) -> None:
    if len(time_series) != len(road_series) or len(time_series) != len(power_series):
        raise ValueError("time_series/road_series/power_series length mismatch")

    width = 960
    height = 620
    margin_left = 70
    margin_right = 30
    margin_top = 90
    margin_bottom = 120

    xmin = min(time_series) if time_series else 0.0
    xmax = max(time_series) if time_series else 1.0
    if xmax <= xmin:
        xmax = xmin + 1.0
    ymin = 0.0
    ymax = 1.05

    def sx(x: float) -> float:
        return margin_left + (float(x) - xmin) / (xmax - xmin) * (width - margin_left - margin_right)

    def sy(y: float) -> float:
        return height - margin_bottom - (float(y) - ymin) / (ymax - ymin) * (height - margin_top - margin_bottom)

    def polyline(xs: List[float], ys: List[float]) -> str:
        return " ".join(f"{sx(x):.2f},{sy(y):.2f}" for x, y in zip(xs, ys))

    def fill_polygon(xs: List[float], ys: List[float], top_y: float) -> str:
        top = [(sx(x), sy(top_y)) for x in reversed(xs)]
        bottom = [(sx(x), sy(y)) for x, y in zip(xs, ys)]
        pts = bottom + top
        return " ".join(f"{px:.2f},{py:.2f}" for px, py in pts)

    road_fill = fill_polygon(time_series, road_series, 1.0)
    power_fill = fill_polygon(time_series, power_series, 1.0)
    road_line = polyline(time_series, road_series)
    power_line = polyline(time_series, power_series)

    title_lines = _wrap_title(title, width=64).splitlines() or [title]
    title_svg = "".join(
        f'<tspan x="{width/2:.0f}" dy="{0 if idx == 0 else 20}">{line}</tspan>'
        for idx, line in enumerate(title_lines)
    )

    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
  <rect width="100%" height="100%" fill="white"/>
  <text x="{width/2:.0f}" y="28" text-anchor="middle" font-family="Arial" font-size="18">{title_svg}</text>
  <line x1="{margin_left}" y1="{height-margin_bottom}" x2="{width-margin_right}" y2="{height-margin_bottom}" stroke="#222" stroke-width="2"/>
  <line x1="{margin_left}" y1="{margin_top}" x2="{margin_left}" y2="{height-margin_bottom}" stroke="#222" stroke-width="2"/>
  <text x="{width/2:.0f}" y="{height-18}" text-anchor="middle" font-family="Arial" font-size="14">Time</text>
  <text x="20" y="{height/2:.0f}" text-anchor="middle" font-family="Arial" font-size="14" transform="rotate(-90 20,{height/2:.0f})">Functionality</text>
  <polygon points="{road_fill}" fill="#1f77b4" opacity="0.18"/>
  <polygon points="{power_fill}" fill="#ff7f0e" opacity="0.18"/>
  <polyline points="{road_line}" fill="none" stroke="#1f77b4" stroke-width="3"/>
  <polyline points="{power_line}" fill="none" stroke="#ff7f0e" stroke-width="3"/>
  <rect x="{width-300}" y="{height-88}" width="250" height="48" fill="white" stroke="#ccc"/>
  <line x1="{width-285}" y1="{height-68}" x2="{width-255}" y2="{height-68}" stroke="#1f77b4" stroke-width="3"/>
  <text x="{width-245}" y="{height-63}" font-family="Arial" font-size="13">Road functionality</text>
  <line x1="{width-285}" y1="{height-46}" x2="{width-255}" y2="{height-46}" stroke="#ff7f0e" stroke-width="3"/>
  <text x="{width-245}" y="{height-41}" font-family="Arial" font-size="13">Power functionality</text>
</svg>
"""
    with open(out_svg, "w", encoding="utf-8") as f:
        f.write(svg)


def _load_destinations(dest_path: str) -> List[int]:
    if not os.path.exists(dest_path) and dest_path == DEFAULT_CRITICAL_LOCATION_PATH and os.path.exists("taskB_essential_destinations.json"):
        dest_path = "taskB_essential_destinations.json"
    dest_obj = load_json(dest_path)
    if not isinstance(dest_obj, dict) or "destinations" not in dest_obj:
        raise KeyError(f"Destination file must contain 'destinations': {dest_path!r}")
    return [int(x) for x in dest_obj["destinations"]]


def _metric_catalog(
    *,
    triangle_area: float,
    equity_summary: Dict[str, float],
    critical_access_summary: Dict[str, float],
) -> Dict[str, float]:
    out = {"triangle": float(triangle_area)}
    for k, v in equity_summary.items():
        out[f"equity:{k}"] = float(v)
    for k, v in critical_access_summary.items():
        out[f"critical_access:{k}"] = float(v)
    return out


def _metric_value(metric_key: str, *, metric_catalog: Dict[str, float]) -> float:
    if metric_key not in metric_catalog:
        raise KeyError(f"Unknown metric {metric_key!r}. Available: {sorted(metric_catalog.keys())}")
    return float(metric_catalog[metric_key])


def _objective_value(
    objective: str,
    *,
    metric_catalog: Dict[str, float],
    objective_weights: Optional[Dict[str, float]] = None,
    objective_reference_values: Optional[Dict[str, float]] = None,
    guardrail_primary_metric: str = "triangle",
    guardrail_metric: str = "",
    guardrail_limit: Optional[float] = None,
    guardrail_penalty: float = 1e6,
) -> Tuple[float, Dict[str, Any]]:
    """
    Minimization objective.
      - single metric: "triangle", "equity:*", "critical_access:*"
      - "weighted_sum" using objective_weights and optional references
      - "guardrail" minimizing guardrail_primary_metric with penalty if guardrail_metric > guardrail_limit
    """
    if objective in metric_catalog:
        value = _metric_value(objective, metric_catalog=metric_catalog)
        return value, {"mode": "single", "metric": objective}

    if objective == "weighted_sum":
        if not objective_weights:
            raise ValueError("objective_weights is required for objective='weighted_sum'")
        refs = objective_reference_values or {}
        total = 0.0
        parts: Dict[str, float] = {}
        for metric_key, weight in objective_weights.items():
            raw = _metric_value(metric_key, metric_catalog=metric_catalog)
            ref = float(refs.get(metric_key, 1.0))
            norm = raw / ref if ref not in (0.0, -0.0) else raw
            contrib = float(weight) * norm
            total += contrib
            parts[metric_key] = contrib
        return total, {"mode": "weighted_sum", "weights": dict(objective_weights), "references": dict(refs), "parts": parts}

    if objective == "guardrail":
        if not guardrail_metric:
            raise ValueError("guardrail_metric is required for objective='guardrail'")
        if guardrail_limit is None:
            raise ValueError("guardrail_limit is required for objective='guardrail'")
        primary = _metric_value(guardrail_primary_metric, metric_catalog=metric_catalog)
        guard = _metric_value(guardrail_metric, metric_catalog=metric_catalog)
        violation = max(0.0, guard - float(guardrail_limit))
        total = primary + float(guardrail_penalty) * violation
        return total, {
            "mode": "guardrail",
            "primary_metric": guardrail_primary_metric,
            "guardrail_metric": guardrail_metric,
            "guardrail_limit": float(guardrail_limit),
            "guardrail_penalty": float(guardrail_penalty),
            "guardrail_value": guard,
            "guardrail_violation": violation,
        }

    raise ValueError(f"Unknown objective {objective!r}")


def run_model_multi(
    sequence: List[Any],
    *,
    result_root: str,
    message: str,
    Scenario: str,
    run_dir: Optional[str] = None,
    strict: bool = True,
    debug: bool = False,
    save_artifacts: bool = True,
    # crew
    crew_mode: str = "specialized",
    power_crews: int = 1,
    road_crews: int = 1,
    multifunction_crews: int = 1,
    depot_node: int = 1,
    crew_speed: float = 1.0,
    service_time_power: float = 20.0,
    service_time_road: float = 10.0,
    bus_dispatch_mode: str = "link_only",
    bus_location_source: str = DEFAULT_BUS_LOCATION_SOURCE,
    bus_to_link_source: str = DEFAULT_BUS_TO_LINK_SOURCE,
    # network
    base_net: str = DEFAULT_BASE_NET,
    trips: str = DEFAULT_TRIPS,
    net1: str = DEFAULT_NET1,
    net2: str = DEFAULT_NET2,
    broken_link_factor: float = 0.0,
    broken_link_factors: Optional[Dict[Tuple[int, int], float]] = None,
    power_road_factor: float = 0.5,
    baseline_tstt: float = 7475338.0,
    # Task B / CRI
    node_to_zone_path: str = "taskB_node_to_zone.json",
    bus_to_zone_path: str = "taskB_bus_to_zone.json",
    dest_path: str = DEFAULT_CRITICAL_LOCATION_PATH,
    cri_w_e: float = 0.133,
    cri_w_a: float = 0.867,
    cri_threshold: float = 0.9,
    critical_access_threshold: float = 0.9,
    # objective
    objective: str = "triangle",
    objective_weights: Optional[Dict[str, float]] = None,
    objective_reference_values: Optional[Dict[str, float]] = None,
    guardrail_primary_metric: str = "triangle",
    guardrail_metric: str = "",
    guardrail_limit: Optional[float] = None,
    guardrail_penalty: float = 1e6,
) -> Dict[str, Any]:
    timestamp = _ts()
    if run_dir is None:
        run_dir = os.path.join(result_root, f"{Scenario}_{timestamp}")
    _ensure_dir(run_dir)

    # Ensure mapping files are in place (your existing behavior)
    if os.path.exists("bus_location.json"):
        os.remove("bus_location.json")
    if os.path.exists("bus_to_link.json"):
        os.remove("bus_to_link.json")
    shutil.copy2(bus_to_link_source, "bus_to_link.json")
    if os.path.exists(bus_location_source):
        shutil.copy2(bus_location_source, "bus_location.json")

    # Load Task B mappings
    node_to_zone = load_json(node_to_zone_path)
    bus_to_zone = load_json(bus_to_zone_path)
    destinations = _load_destinations(dest_path)

    zones = build_zone_lists(node_to_zone)

    # Baseline accessibility TT0: run TAP-B once on no-disruption state
    _prepare_state_and_run_tapb(
        broken_buses=[],
        broken_links=[],
        broken_link_factors=None,
        base_net=base_net,
        trips=trips,
        net1=net1,
        net2=net2,
        broken_link_factor=broken_link_factor,
        power_road_factor=power_road_factor,
        strict=strict,
    )
    baseline_s = os.path.join(run_dir, "baseline_s.txt")
    shutil.copy2("s.txt", baseline_s)
    TT0 = compute_accessibility_TT(s_txt_path=baseline_s, zones=zones, destinations=destinations)

    seq, power_sequence, road_sequence = _canonicalize_sequence_for_crews(list(sequence), crew_mode=crew_mode)
    broken_buses = {a for a in seq if not isinstance(a, tuple)}
    broken_links = {a for a in seq if isinstance(a, tuple)}

    # Dispatch snapshot (initial damaged state) for scheduling
    _prepare_state_and_run_tapb(
        broken_buses=list(broken_buses),
        broken_links=list(broken_links),
        broken_link_factors=broken_link_factors,
        base_net=base_net,
        trips=trips,
        net1=net1,
        net2=net2,
        broken_link_factor=broken_link_factor,
        power_road_factor=power_road_factor,
        strict=strict,
    )
    dispatch_s = os.path.join(run_dir, "dispatch_s.txt")
    shutil.copy2("s.txt", dispatch_s)

    pool: CrewPool = make_crews(
        mode=crew_mode,
        power_crews=power_crews,
        road_crews=road_crews,
        multifunction_crews=multifunction_crews,
        depot_node=depot_node,
        speed=crew_speed,
    )

    sim = evaluate_with_crews(
        sequence=seq,
        crew_pool=pool,
        service_time_power=service_time_power,
        service_time_road=service_time_road,
        default_depot=depot_node,
        s_txt_path=dispatch_s,
        bus_to_link_path="bus_to_link.json",
        strict=strict,
        bus_dispatch_mode=bus_dispatch_mode,
    )
    events = sorted(sim["timeline"], key=lambda x: x[0])

    # initial functionality + CRI
    road_func = eval_road_resilience(
        list(broken_buses),
        list(broken_links),
        broken_link_factors=broken_link_factors,
        base_net=base_net,
        trips=trips,
        net1=net1,
        net2=net2,
        broken_link_factor=broken_link_factor,
        power_road_factor=power_road_factor,
        baseline_tstt=baseline_tstt,
        strict=strict,
    )
    power_func = eval_power_resilience(list(broken_buses))

    # CRI series (aligned with time_series)
    time_series = [0.0]
    road_series = [road_func]
    power_series = [power_func]

    # build CRI_z(t=0)
    func_buses0 = set(get_functional_nodes(set(map(int, broken_buses))))
    E0 = compute_E_by_zone(zones=zones, bus_to_zone=bus_to_zone, functional_buses=func_buses0, empty_zone_policy="error")
    TT = compute_accessibility_TT(s_txt_path=dispatch_s, zones=zones, destinations=destinations)
    A0 = compute_accessibility_ratio_by_zone(zones=zones, TT0=TT0, TT=TT)
    CRI0 = compute_CRI_by_zone(zones=zones, E=E0, TT0=TT0, TT=TT, w_e=cri_w_e, w_a=cri_w_a)
    cri_series: List[Dict[int, float]] = [CRI0]
    access_series: List[Dict[int, float]] = [A0]

    triangle_area = 0.0
    t_prev = 0.0

    debug_lines: List[str] = []
    if debug:
        debug_lines.append(f"[INIT] t=0 road={road_func} power={power_func}")

    for (t, asset) in events:
        dt = float(t - t_prev)
        triangle_area += ((1 - road_func) + (1 - power_func)) * dt
        t_prev = float(t)

        if isinstance(asset, tuple):
            broken_links.discard(asset)
        else:
            broken_buses.discard(asset)

        # Re-evaluate system-level functions (updates s.txt)
        road_func = eval_road_resilience(
            list(broken_buses),
            list(broken_links),
            broken_link_factors=broken_link_factors,
            base_net=base_net,
            trips=trips,
            net1=net1,
            net2=net2,
            broken_link_factor=broken_link_factor,
            power_road_factor=power_road_factor,
            baseline_tstt=baseline_tstt,
            strict=strict,
        )
        power_func = eval_power_resilience(list(broken_buses))

        time_series.append(t_prev)
        road_series.append(road_func)
        power_series.append(power_func)

        # CRI at this event time, using current s.txt
        func_buses = set(get_functional_nodes(set(map(int, broken_buses))))
        E = compute_E_by_zone(zones=zones, bus_to_zone=bus_to_zone, functional_buses=func_buses, empty_zone_policy="error")
        current_s = os.path.abspath("s.txt")
        TT = compute_accessibility_TT(s_txt_path=current_s, zones=zones, destinations=destinations)
        A = compute_accessibility_ratio_by_zone(zones=zones, TT0=TT0, TT=TT)
        CRI = compute_CRI_by_zone(zones=zones, E=E, TT0=TT0, TT=TT, w_e=cri_w_e, w_a=cri_w_a)
        cri_series.append(CRI)
        access_series.append(A)

        if debug:
            debug_lines.append(f"[EVENT] t={t_prev} asset={asset} road={road_func} power={power_func}")

    equity_summary = equity_summary_from_CRI(zones=zones, time_series=time_series, cri_series=cri_series, threshold=cri_threshold)
    critical_access_summary = critical_access_summary_from_series(
        zones=zones,
        time_series=time_series,
        access_series=access_series,
        threshold=critical_access_threshold,
    )
    metric_catalog = _metric_catalog(
        triangle_area=triangle_area,
        equity_summary=equity_summary,
        critical_access_summary=critical_access_summary,
    )
    obj_value, objective_detail = _objective_value(
        objective,
        metric_catalog=metric_catalog,
        objective_weights=objective_weights,
        objective_reference_values=objective_reference_values,
        guardrail_primary_metric=guardrail_primary_metric,
        guardrail_metric=guardrail_metric,
        guardrail_limit=guardrail_limit,
        guardrail_penalty=guardrail_penalty,
    )

    tri_png = os.path.join(run_dir, f"{Scenario}_triangle.png")
    tri_plot_path = ""

    if save_artifacts:
        # write summary files
        raw_path = os.path.join(run_dir, f"{Scenario}_raw.txt")
        en_path = os.path.join(run_dir, f"{Scenario}_report_en.txt")
        zh_path = os.path.join(run_dir, f"{Scenario}_report_zh.txt")
        dbg_path = os.path.join(run_dir, f"{Scenario}_debug.txt") if debug_lines else ""

        with open(raw_path, "w", encoding="utf-8") as f:
            print(f"timestamp: {timestamp}", file=f)
            print(message, file=f)
            print("sequence:", seq, file=f)
            print("power_sequence:", power_sequence, file=f)
            print("road_sequence:", road_sequence, file=f)
            print("timeline:", events, file=f)
            print("triangle_area:", triangle_area, file=f)
            print("equity_summary:", equity_summary, file=f)
            print("critical_access_summary:", critical_access_summary, file=f)
            print("metric_catalog:", metric_catalog, file=f)
            print("objective:", objective, file=f)
            print("objective_detail:", objective_detail, file=f)
            print("objective_value:", obj_value, file=f)
            print("bus_dispatch_mode:", bus_dispatch_mode, file=f)
            print("bus_location_source:", bus_location_source, file=f)
            print("bus_to_link_source:", bus_to_link_source, file=f)
            print("crew_mode:", crew_mode, file=f)
            print("power_crews:", power_crews, file=f)
            print("road_crews:", road_crews, file=f)
            print("multifunction_crews:", multifunction_crews, file=f)
            print("broken_link_factors:", broken_link_factors, file=f)

        with open(en_path, "w", encoding="utf-8") as f:
            print("Task B run report (CRI + equity)", file=f)
            print(f"Timestamp: {timestamp}", file=f)
            print(f"Scenario: {Scenario}", file=f)
            print(f"Objective: {objective}  value={_fmt(obj_value)}", file=f)
            print(f"CRI weights: w_e={_fmt(cri_w_e)}  w_a={_fmt(cri_w_a)}", file=f)
            print(f"CRI threshold: {_fmt(cri_threshold)}", file=f)
            print(f"Critical access threshold: {_fmt(critical_access_threshold)}", file=f)
            print(f"Bus dispatch mode: {bus_dispatch_mode}", file=f)
            print(f"Bus-to-link source: {bus_to_link_source}", file=f)
            print(f"Crew mode: {crew_mode}", file=f)
            print(f"Power crews: {power_crews}  Road crews: {road_crews}  Multifunction crews: {multifunction_crews}", file=f)
            print(f"Power sequence: {power_sequence}", file=f)
            print(f"Road sequence: {road_sequence}", file=f)
            print(f"Triangle area: {_fmt(triangle_area)}", file=f)
            print("", file=f)
            print("Equity summary (CRI-based):", file=f)
            for k, v in equity_summary.items():
                print(f"- {k}: {_fmt(v)}", file=f)
            print("", file=f)
            print("Critical-access summary (A_z-based):", file=f)
            for k, v in critical_access_summary.items():
                print(f"- {k}: {_fmt(v)}", file=f)

        with open(zh_path, "w", encoding="utf-8") as f:
            print("Task B 运行报告（CRI + 公平性）", file=f)
            print(f"时间戳：{timestamp}", file=f)
            print(f"场景：{Scenario}", file=f)
            print(f"目标函数：{objective}  值={_fmt(obj_value)}", file=f)
            print(f"CRI 权重：w_e={_fmt(cri_w_e)}  w_a={_fmt(cri_w_a)}", file=f)
            print(f"CRI 阈值：{_fmt(cri_threshold)}", file=f)
            print(f"关键可达性阈值：{_fmt(critical_access_threshold)}", file=f)
            print(f"调度映射模式：{bus_dispatch_mode}", file=f)
            print(f"Bus-to-link 来源：{bus_to_link_source}", file=f)
            print(f"队伍模式：{crew_mode}", file=f)
            print(f"电力队伍：{power_crews}  路网队伍：{road_crews}  多功能队伍：{multifunction_crews}", file=f)
            print(f"电力修复序列：{power_sequence}", file=f)
            print(f"路网修复序列：{road_sequence}", file=f)
            print(f"Triangle 补面积：{_fmt(triangle_area)}", file=f)
            print("", file=f)
            print("公平性汇总指标（基于 CRI_z 的分布）:", file=f)
            for k, v in equity_summary.items():
                print(f"- {k}：{_fmt(v)}", file=f)
            print("", file=f)
            print("关键可达性汇总指标（基于 A_z 的分布）:", file=f)
            for k, v in critical_access_summary.items():
                print(f"- {k}：{_fmt(v)}", file=f)

        if debug_lines:
            with open(dbg_path, "w", encoding="utf-8") as f:
                for line in debug_lines:
                    f.write(line + "\n")

        tri_plot_path = _plot_triangle_with_shading(
            time_series=time_series,
            road_series=road_series,
            power_series=power_series,
            out_png=tri_png,
            title=f"{Scenario} | obj={objective}={_fmt(obj_value)} | area={_fmt(triangle_area)} | {timestamp}",
        )

        # save CRI series as json for later analysis
        cri_path = os.path.join(run_dir, "cri_series.json")
        with open(cri_path, "w", encoding="utf-8") as f:
            f.write(
                __import__("json").dumps(
                    {
                        "time_series": time_series,
                        "sequence": seq,
                        "power_sequence": power_sequence,
                        "road_sequence": road_sequence,
                        "zones": zones,
                        "destinations": destinations,
                        "cri_weights": {"w_e": cri_w_e, "w_a": cri_w_a},
                        "cri_threshold": cri_threshold,
                        "critical_access_threshold": critical_access_threshold,
                        "broken_link_factors": {
                            f"{u}-{v}": factor for (u, v), factor in (broken_link_factors or {}).items()
                        },
                        "cri_series": cri_series,
                        "access_series": access_series,
                    },
                    indent=2,
                )
            )

    return {
        "run_dir": run_dir,
        "timestamp": timestamp,
        "sequence": seq,
        "power_sequence": list(power_sequence),
        "road_sequence": list(road_sequence),
        "triangle_area": triangle_area,
        "equity_summary": equity_summary,
        "critical_access_summary": critical_access_summary,
        "metric_catalog": metric_catalog,
        "objective": objective,
        "objective_detail": objective_detail,
        "objective_value": obj_value,
        "cri_weights": {"w_e": float(cri_w_e), "w_a": float(cri_w_a)},
        "cri_threshold": float(cri_threshold),
        "critical_access_threshold": float(critical_access_threshold),
        "bus_dispatch_mode": bus_dispatch_mode,
        "bus_location_source": bus_location_source,
        "bus_to_link_source": bus_to_link_source,
        "crew_mode": crew_mode,
        "power_crews": int(power_crews),
        "road_crews": int(road_crews),
        "multifunction_crews": int(multifunction_crews),
        "broken_link_factors": dict(broken_link_factors or {}),
        "triangle_png": (tri_plot_path if save_artifacts else ""),
    }


def optimize_sequence_sa(
    *,
    base_sequence: List[Any],
    result_root: str,
    message: str,
    Scenario: str,
    objective: str = "triangle",
    sa: Optional[SAConfig] = None,
    strict: bool = True,
    save_baseline: bool = True,
    save_best_artifacts: bool = True,
    save_best_debug: bool = True,
    # pass-through
    crew_mode: str = "specialized",
    power_crews: int = 1,
    road_crews: int = 1,
    multifunction_crews: int = 1,
    bus_dispatch_mode: str = "link_only",
    bus_location_source: str = DEFAULT_BUS_LOCATION_SOURCE,
    bus_to_link_source: str = DEFAULT_BUS_TO_LINK_SOURCE,
    dest_path: str = DEFAULT_CRITICAL_LOCATION_PATH,
    cri_w_e: float = 0.133,
    cri_w_a: float = 0.867,
    cri_threshold: float = 0.9,
    critical_access_threshold: float = 0.9,
    objective_weights: Optional[Dict[str, float]] = None,
    objective_reference_values: Optional[Dict[str, float]] = None,
    guardrail_primary_metric: str = "triangle",
    guardrail_metric: str = "",
    guardrail_limit: Optional[float] = None,
    guardrail_penalty: float = 1e6,
    broken_link_factors: Optional[Dict[Tuple[int, int], float]] = None,
) -> Dict[str, Any]:
    if sa is None:
        sa = SAConfig()

    session_ts = _ts()
    session_dir = os.path.join(result_root, f"{Scenario}_{session_ts}")
    _ensure_dir(session_dir)

    log_path = os.path.join(session_dir, "sa_log.csv")
    with open(log_path, "w", encoding="utf-8") as f:
        f.write(
            "iter,tag,objective_value,triangle_area,var_restore,gini_restore,p90_restore,"
            "p90_access_restore,sequence\n"
        )

    canonical_base_sequence, power_base_sequence, road_base_sequence = _canonicalize_sequence_for_crews(
        list(base_sequence),
        crew_mode=crew_mode,
    )
    neighbor_fn = None
    if crew_mode == "specialized":
        neighbor_fn = _specialized_neighbor_factory(
            n_power=len(power_base_sequence),
            neighbor_mode=sa.neighbor,
        )

    def _append_log(it: int, tag: str, run: Dict[str, Any]) -> None:
        es = run["equity_summary"]
        cas = run["critical_access_summary"]
        with open(log_path, "a", encoding="utf-8") as f:
            f.write(
                f"{it},{tag},{run['objective_value']},{run['triangle_area']},"
                f"{es.get('var_restore')},{es.get('gini_restore')},{es.get('p90_restore')},"
                f"{cas.get('p90_access_restore')},"
                f"{repr(run['sequence'])}\n"
            )

    baseline_run = None
    if save_baseline:
        baseline_dir = os.path.join(session_dir, "baseline")
        baseline_run = run_model_multi(
            canonical_base_sequence,
            result_root=result_root,
            message=message + " | baseline",
            Scenario=Scenario,
            run_dir=baseline_dir,
            strict=strict,
            debug=False,
            save_artifacts=True,
            crew_mode=crew_mode,
            power_crews=power_crews,
            road_crews=road_crews,
            multifunction_crews=multifunction_crews,
            bus_dispatch_mode=bus_dispatch_mode,
            bus_location_source=bus_location_source,
            bus_to_link_source=bus_to_link_source,
            dest_path=dest_path,
            cri_w_e=cri_w_e,
            cri_w_a=cri_w_a,
            objective=objective,
            cri_threshold=cri_threshold,
            critical_access_threshold=critical_access_threshold,
            objective_weights=objective_weights,
            objective_reference_values=objective_reference_values,
            guardrail_primary_metric=guardrail_primary_metric,
            guardrail_metric=guardrail_metric,
            guardrail_limit=guardrail_limit,
            guardrail_penalty=guardrail_penalty,
            broken_link_factors=broken_link_factors,
        )
        _append_log(-1, "baseline", baseline_run)

    it_counter = {"i": 0}

    def _eval(seq: List[Any], run_tag: str) -> Dict[str, Any]:
        i = it_counter["i"]
        it_counter["i"] += 1
        tmp_dir = os.path.join(session_dir, "_tmp")
        _ensure_dir(tmp_dir)

        run = run_model_multi(
            seq,
            result_root=result_root,
            message=message + f" | {run_tag}",
            Scenario=Scenario,
            run_dir=tmp_dir,
            strict=strict,
            debug=False,
            save_artifacts=False,  # no per-iteration files
            crew_mode=crew_mode,
            power_crews=power_crews,
            road_crews=road_crews,
            multifunction_crews=multifunction_crews,
            bus_dispatch_mode=bus_dispatch_mode,
            bus_location_source=bus_location_source,
            bus_to_link_source=bus_to_link_source,
            dest_path=dest_path,
            cri_w_e=cri_w_e,
            cri_w_a=cri_w_a,
            objective=objective,
            cri_threshold=cri_threshold,
            critical_access_threshold=critical_access_threshold,
            objective_weights=objective_weights,
            objective_reference_values=objective_reference_values,
            guardrail_primary_metric=guardrail_primary_metric,
            guardrail_metric=guardrail_metric,
            guardrail_limit=guardrail_limit,
            guardrail_penalty=guardrail_penalty,
            broken_link_factors=broken_link_factors,
        )
        _append_log(i, run_tag, run)
        return run

    sa_res = simulated_annealing(
        initial=list(canonical_base_sequence),
        evaluate=_eval,
        neighbor_fn=neighbor_fn,
        config=sa,
        scenario_prefix=f"{Scenario}_sa",
    )

    best_dir = os.path.join(session_dir, "best")
    best_run = run_model_multi(
        sa_res["best_sequence"],
        result_root=result_root,
        message=message + " | BEST",
        Scenario=Scenario,
        run_dir=best_dir,
        strict=strict,
        debug=False,
        save_artifacts=save_best_artifacts,
        crew_mode=crew_mode,
        power_crews=power_crews,
        road_crews=road_crews,
        multifunction_crews=multifunction_crews,
        bus_dispatch_mode=bus_dispatch_mode,
        bus_location_source=bus_location_source,
        bus_to_link_source=bus_to_link_source,
        dest_path=dest_path,
        cri_w_e=cri_w_e,
        cri_w_a=cri_w_a,
        objective=objective,
        cri_threshold=cri_threshold,
        critical_access_threshold=critical_access_threshold,
        objective_weights=objective_weights,
        objective_reference_values=objective_reference_values,
        guardrail_primary_metric=guardrail_primary_metric,
        guardrail_metric=guardrail_metric,
        guardrail_limit=guardrail_limit,
        guardrail_penalty=guardrail_penalty,
        broken_link_factors=broken_link_factors,
    )

    best_debug_run = None
    if save_best_debug:
        best_debug_dir = os.path.join(session_dir, "best_debug")
        best_debug_run = run_model_multi(
            sa_res["best_sequence"],
            result_root=result_root,
            message=message + " | BEST_DEBUG",
            Scenario=Scenario,
            run_dir=best_debug_dir,
            strict=strict,
            debug=True,
            save_artifacts=True,
            crew_mode=crew_mode,
            power_crews=power_crews,
            road_crews=road_crews,
            multifunction_crews=multifunction_crews,
            bus_dispatch_mode=bus_dispatch_mode,
            bus_location_source=bus_location_source,
            bus_to_link_source=bus_to_link_source,
            dest_path=dest_path,
            cri_w_e=cri_w_e,
            cri_w_a=cri_w_a,
            objective=objective,
            cri_threshold=cri_threshold,
            critical_access_threshold=critical_access_threshold,
            objective_weights=objective_weights,
            objective_reference_values=objective_reference_values,
            guardrail_primary_metric=guardrail_primary_metric,
            guardrail_metric=guardrail_metric,
            guardrail_limit=guardrail_limit,
            guardrail_penalty=guardrail_penalty,
            broken_link_factors=broken_link_factors,
        )

    return {
        "session_dir": session_dir,
        "sa_log": log_path,
        "baseline_run": baseline_run,
        "best_sequence": sa_res["best_sequence"],
        "best_objective_value": sa_res["best_objective_value"],
        "best_run": best_run,
        "best_debug_run": best_debug_run,
        "sa_config": sa,
        "objective": objective,
        "objective_weights": objective_weights,
        "objective_reference_values": objective_reference_values,
        "guardrail_primary_metric": guardrail_primary_metric,
        "guardrail_metric": guardrail_metric,
        "guardrail_limit": guardrail_limit,
        "bus_dispatch_mode": bus_dispatch_mode,
        "bus_location_source": bus_location_source,
        "bus_to_link_source": bus_to_link_source,
        "crew_mode": crew_mode,
        "power_crews": int(power_crews),
        "road_crews": int(road_crews),
        "multifunction_crews": int(multifunction_crews),
        "broken_link_factors": dict(broken_link_factors or {}),
    }
