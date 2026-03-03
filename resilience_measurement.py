"""
Resilience evaluation (research / strict-by-default).

Enhancements:
  1) Timestamped outputs + formatted reports (EN + ZH) + optional debug log
  2) Resilience triangle plotting (PNG)
  3) Optional optimization: objective can be triangle or equity metric
"""

from __future__ import annotations

import os
import shutil
import math
import random
from datetime import datetime
from typing import Any, Dict, List, Tuple, Optional

from crews import CrewPool, make_crews
from scheduler import evaluate_with_crews
from run_tapb import run_tapb
from road_util import capacity_adjustment, eval_tot_OD_travel_time
from power_util import delete_buses, get_functional_nodes
from interdependency import power_to_road

BUS_COUNT = 33

# Default TAP-B inputs (adjust if your paths differ)
DEFAULT_BASE_NET = "tap-b/net/SiouxFalls_net.txt"
DEFAULT_TRIPS = "tap-b/net/SiouxFalls_trips.txt"

# Working IO files
DEFAULT_WORK_DIR = "work"
DEFAULT_NET1 = os.path.join(DEFAULT_WORK_DIR, "SiouxFalls_net1.txt")
DEFAULT_NET2 = os.path.join(DEFAULT_WORK_DIR, "SiouxFalls_net2.txt")


def _ts() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def _require(path: str) -> None:
    if not os.path.exists(path):
        raise FileNotFoundError(f"Required file not found: {path!r}")


def _ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def _fmt(x: float, digits: int = 4) -> str:
    # significant digits formatting
    if x is None:
        return "NA"
    if isinstance(x, (int,)):
        return str(x)
    try:
        x = float(x)
    except Exception:
        return str(x)
    if x == 0:
        return "0"
    # use general format with significant digits
    return f"{x:.{digits}g}"


def eval_power_resilience(broken_buses: List[int]) -> float:
    functional = set(get_functional_nodes(set(map(int, broken_buses))))
    return float(len(functional)) / float(BUS_COUNT)


def _prepare_state_and_run_tapb(
    *,
    broken_buses: List[int],
    broken_links: List[Tuple[int, int]],
    base_net: str,
    trips: str,
    net1: str,
    net2: str,
    broken_link_factor: float,
    power_road_factor: float,
    strict: bool,
) -> None:
    """
    IO workflow:
      base_net --(capacity_adjustment on broken_links)--> net1
      broken_buses --(delete_buses)--> unfunctional_nodes
      net1 + unfunctional_nodes --(power_to_road)--> net2
      run_tapb(net2, trips) -> produces s.txt
    """
    _require(base_net)
    _require(trips)
    _ensure_dir(os.path.dirname(net1) or ".")

    capacity_adjustment(base_net, net1, broken_links, broken_link_factor)
    _require(net1)

    unfunctional_nodes = delete_buses(list(map(int, broken_buses)))
    if strict and unfunctional_nodes is None:
        raise RuntimeError("delete_buses returned None; check power_util.delete_buses implementation.")

    power_to_road(unfunctional_nodes, net1, net2, power_road_factor)
    _require(net2)

    # fresh s.txt
    if os.path.exists("s.txt"):
        os.remove("s.txt")
    run_tapb(net2, trips, check=True)

    if strict:
        _require("s.txt")


def eval_road_resilience(
    broken_buses: List[int],
    broken_links: List[Tuple[int, int]],
    *,
    base_net: str = DEFAULT_BASE_NET,
    trips: str = DEFAULT_TRIPS,
    net1: str = DEFAULT_NET1,
    net2: str = DEFAULT_NET2,
    broken_link_factor: float = 0.0,
    power_road_factor: float = 0.5,
    baseline_tstt: float = 7475338.0,
    strict: bool = True,
) -> float:
    _prepare_state_and_run_tapb(
        broken_buses=broken_buses,
        broken_links=broken_links,
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


def _plot_triangle(
    *,
    time_series: List[float],
    road_series: List[float],
    power_series: List[float],
    out_png: str,
    title: str,
) -> None:
    import matplotlib.pyplot as plt

    if len(time_series) != len(road_series) or len(time_series) != len(power_series):
        raise ValueError("time_series/road_series/power_series length mismatch")

    plt.figure()
    plt.plot(time_series, road_series, label="Road functionality", linewidth=2.0)
    plt.plot(time_series, power_series, label="Power functionality", linewidth=2.0)
    plt.ylim(0.0, 1.05)
    plt.xlabel("Time")
    plt.ylabel("Functionality")
    plt.title(title)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


def _write_reports(
    *,
    result_folder: str,
    scenario: str,
    timestamp: str,
    message: str,
    params: Dict[str, Any],
    sequence: List[Any],
    timeline: List[Tuple[float, Any]],
    equity: Dict[str, Any],
    time_series: List[float],
    road_series: List[float],
    power_series: List[float],
    triangle_area: float,
    objective: str,
    objective_value: float,
    debug_lines: Optional[List[str]] = None,
) -> Dict[str, str]:
    _ensure_dir(result_folder)

    # filenames (no overwrite)
    raw_path = os.path.join(result_folder, f"{scenario}_multi_summary_{timestamp}.txt")
    en_path = os.path.join(result_folder, f"{scenario}_report_en_{timestamp}.txt")
    zh_path = os.path.join(result_folder, f"{scenario}_report_zh_{timestamp}.txt")
    dbg_path = os.path.join(result_folder, f"{scenario}_debug_{timestamp}.txt") if debug_lines else ""

    # raw (still readable, but keep simple)
    with open(raw_path, "w", encoding="utf-8") as f:
        print(f"timestamp: {timestamp}", file=f)
        print(message, file=f)
        print("sequence:", list(sequence), file=f)
        print("params:", params, file=f)
        print("timeline:", timeline, file=f)
        print("equity:", equity, file=f)
        print("time:", time_series, file=f)
        print("road functionality:", road_series, file=f)
        print("power functionality:", power_series, file=f)
        print("triangle_area:", triangle_area, file=f)
        print("objective:", objective, file=f)
        print("objective_value:", objective_value, file=f)

    # English report
    with open(en_path, "w", encoding="utf-8") as f:
        print("Run report (Power–Road interdependency)", file=f)
        print(f"Timestamp: {timestamp}", file=f)
        print("", file=f)
        print("Summary", file=f)
        print(f"- Scenario: {scenario}", file=f)
        print(f"- Objective: {objective}", file=f)
        print(f"- Objective value: {_fmt(objective_value)}", file=f)
        print(f"- Triangle area (road+power complement integral): {_fmt(triangle_area)}", file=f)
        print("", file=f)

        print("Inputs", file=f)
        for k, v in params.items():
            print(f"- {k}: {v}", file=f)
        print("", file=f)

        print("Sequence", file=f)
        print(list(sequence), file=f)
        print("", file=f)

        print("Crew completion timeline (finish_time, asset)", file=f)
        for t, a in timeline:
            print(f"- t={_fmt(t)}  asset={a}", file=f)
        print("", file=f)

        print("Equity metrics (on asset completion times)", file=f)
        for k, v in (equity or {}).items():
            print(f"- {k}: {_fmt(v)}", file=f)
        print("", file=f)

        print("Resilience series", file=f)
        print(f"- time points: {len(time_series)}", file=f)
        print("  (time, road, power)", file=f)
        for i in range(len(time_series)):
            print(f"  {_fmt(time_series[i])}\t{_fmt(road_series[i])}\t{_fmt(power_series[i])}", file=f)

    # Chinese report
    with open(zh_path, "w", encoding="utf-8") as f:
        print("运行报告（电力–道路互依恢复仿真）", file=f)
        print(f"时间戳：{timestamp}", file=f)
        print("", file=f)
        print("摘要", file=f)
        print(f"- 场景：{scenario}", file=f)
        print(f"- 目标函数：{objective}", file=f)
        print(f"- 目标函数值：{_fmt(objective_value)}", file=f)
        print(f"- Resilience triangle（路+电补面积积分）：{_fmt(triangle_area)}", file=f)
        print("", file=f)

        print("输入参数", file=f)
        for k, v in params.items():
            print(f"- {k}：{v}", file=f)
        print("", file=f)

        print("修复序列（int=母线，tuple=道路link）", file=f)
        print(list(sequence), file=f)
        print("", file=f)

        print("维修完成时间线（finish_time, asset）", file=f)
        for t, a in timeline:
            print(f"- t={_fmt(t)}  资产={a}", file=f)
        print("", file=f)

        print("公平性指标（基于资产完成时间分布）", file=f)
        for k, v in (equity or {}).items():
            print(f"- {k}：{_fmt(v)}", file=f)
        print("", file=f)

        print("功能性曲线数据（用于画 triangle）", file=f)
        print("（time, road, power）", file=f)
        for i in range(len(time_series)):
            print(f"{_fmt(time_series[i])}\t{_fmt(road_series[i])}\t{_fmt(power_series[i])}", file=f)

    # debug log
    if debug_lines:
        with open(dbg_path, "w", encoding="utf-8") as f:
            print(f"timestamp: {timestamp}", file=f)
            for line in debug_lines:
                f.write(line.rstrip("\n") + "\n")

    return {"raw": raw_path, "en": en_path, "zh": zh_path, "debug": dbg_path}


def _objective_from_run(
    *,
    objective: str,
    triangle_area: float,
    equity: Dict[str, Any],
) -> float:
    """
    objective:
      - "triangle"
      - "equity:<metric>" where metric in {"gini","theil","atkinson","jain", ...}
    Minimization objective.
    """
    if objective == "triangle":
        return float(triangle_area)

    if objective.startswith("equity:"):
        metric = objective.split(":", 1)[1].strip().lower()
        if not equity or metric not in equity:
            raise KeyError(f"Equity metric {metric!r} not found in equity dict keys={list((equity or {}).keys())}")

        val = float(equity[metric])

        # Jain is usually higher is better (maximization). Convert to minimization.
        if metric == "jain":
            return 1.0 - val
        return val

    raise ValueError(f"Unknown objective: {objective!r}")


def run_model_multi(
    sequence: List[Any],
    result_folder: str,
    message: str,
    Scenario: str,
    plot_control: bool,
    focus: bool,
    *,
    strict: bool = True,
    debug: bool = False,
    crew_mode: str = "specialized",  # "specialized" or "multifunction"
    power_crews: int = 1,
    road_crews: int = 1,
    multifunction_crews: int = 1,
    depot_node: int = 1,
    crew_speed: float = 1.0,
    service_time_power: float = 20.0,
    service_time_road: float = 10.0,
    # network params
    base_net: str = DEFAULT_BASE_NET,
    trips: str = DEFAULT_TRIPS,
    net1: str = DEFAULT_NET1,
    net2: str = DEFAULT_NET2,
    broken_link_factor: float = 0.0,
    power_road_factor: float = 0.5,
    baseline_tstt: float = 7475338.0,
    # objective
    objective: str = "triangle",
) -> Dict[str, Any]:
    """
    Multi-crew evaluator with event-driven integration.
    Dispatch costs come from ONE s.txt snapshot built from the initial damaged state.
    Returns a dict with run artifacts and paths.
    """
    timestamp = _ts()

    # scenario mapping files (keep same as before)
    if os.path.exists("bus_location.json"):
        os.remove("bus_location.json")
    if os.path.exists("bus_to_link.json"):
        os.remove("bus_to_link.json")
    shutil.copy2("original_bus_to_link.json", "bus_to_link.json")
    shutil.copy2("original_bus_location.json", "bus_location.json")

    seq = list(sequence)
    broken_buses = {a for a in seq if not isinstance(a, tuple)}
    broken_links = {a for a in seq if isinstance(a, tuple)}

    # build dispatch snapshot
    _prepare_state_and_run_tapb(
        broken_buses=list(broken_buses),
        broken_links=list(broken_links),
        base_net=base_net,
        trips=trips,
        net1=net1,
        net2=net2,
        broken_link_factor=broken_link_factor,
        power_road_factor=power_road_factor,
        strict=strict,
    )

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
        s_txt_path="s.txt",
        bus_to_link_path="bus_to_link.json",
        strict=strict,
    )
    events = sorted(sim["timeline"], key=lambda x: x[0])

    # initial functionality
    road_func = eval_road_resilience(
        list(broken_buses),
        list(broken_links),
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

    # event-driven integration
    t_prev = 0.0
    time_series = [0.0]
    road_series = [road_func]
    power_series = [power_func]
    triangle_area = 0.0

    debug_lines: List[str] = []
    if debug:
        debug_lines.append(f"[INIT] t=0 road={road_func} power={power_func} broken_buses={sorted(list(broken_buses))} broken_links={sorted(list(broken_links))}")

    for (t, asset) in events:
        dt = float(t - t_prev)
        triangle_area += ((1 - road_func) + (1 - power_func)) * dt
        if debug:
            debug_lines.append(f"[INTERVAL] ({t_prev} -> {t}) dt={dt} road={road_func} power={power_func} add_area={((1-road_func)+(1-power_func))*dt}")

        t_prev = float(t)

        # apply repair
        if isinstance(asset, tuple):
            broken_links.discard(asset)
        else:
            broken_buses.discard(asset)

        # re-evaluate functionality after repair
        road_func = eval_road_resilience(
            list(broken_buses),
            list(broken_links),
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

        if debug:
            debug_lines.append(f"[EVENT] t={t_prev} asset={asset} -> road={road_func} power={power_func} broken_buses={sorted(list(broken_buses))} broken_links={sorted(list(broken_links))}")

    equity = sim.get("equity", {}) or {}
    obj_value = _objective_from_run(objective=objective, triangle_area=triangle_area, equity=equity)

    params = {
        "strict": strict,
        "crew_mode": crew_mode,
        "power_crews": power_crews,
        "road_crews": road_crews,
        "multifunction_crews": multifunction_crews,
        "depot_node": depot_node,
        "crew_speed": crew_speed,
        "service_time_power": service_time_power,
        "service_time_road": service_time_road,
        "broken_link_factor": broken_link_factor,
        "power_road_factor": power_road_factor,
        "baseline_tstt": baseline_tstt,
        "base_net": base_net,
        "trips": trips,
        "net1": net1,
        "net2": net2,
    }

    paths = _write_reports(
        result_folder=result_folder,
        scenario=Scenario,
        timestamp=timestamp,
        message=message,
        params=params,
        sequence=seq,
        timeline=events,
        equity=equity,
        time_series=time_series,
        road_series=road_series,
        power_series=power_series,
        triangle_area=triangle_area,
        objective=objective,
        objective_value=obj_value,
        debug_lines=(debug_lines if debug else None),
    )

    # plot triangle
    tri_png = os.path.join(result_folder, f"{Scenario}_triangle_{timestamp}.png")
    _plot_triangle(
        time_series=time_series,
        road_series=road_series,
        power_series=power_series,
        out_png=tri_png,
        title=f"{Scenario} | objective={objective} | obj={_fmt(obj_value)} | area={_fmt(triangle_area)} | {timestamp}",
    )

    return {
        "timestamp": timestamp,
        "sequence": seq,
        "timeline": events,
        "equity": equity,
        "time_series": time_series,
        "road_series": road_series,
        "power_series": power_series,
        "triangle_area": triangle_area,
        "objective": objective,
        "objective_value": obj_value,
        "paths": paths,
        "triangle_png": tri_png,
    }


def optimize_sequence(
    *,
    base_sequence: List[Any],
    result_folder: str,
    message: str,
    Scenario: str,
    objective: str = "triangle",
    method: str = "bruteforce",  # "bruteforce" or "random"
    max_iter: int = 200,
    seed: int = 0,
    strict: bool = True,
    debug: bool = False,
    # pass-through run params
    crew_mode: str = "multifunction",
    multifunction_crews: int = 1,
    power_crews: int = 1,
    road_crews: int = 1,
    depot_node: int = 1,
    crew_speed: float = 1.0,
    service_time_power: float = 20.0,
    service_time_road: float = 10.0,
    base_net: str = DEFAULT_BASE_NET,
    trips: str = DEFAULT_TRIPS,
    net1: str = DEFAULT_NET1,
    net2: str = DEFAULT_NET2,
    broken_link_factor: float = 0.0,
    power_road_factor: float = 0.5,
    baseline_tstt: float = 7475338.0,
) -> Dict[str, Any]:
    """
    Optimization wrapper.

    bruteforce: exact best for small n (recommended for QA / small experiments).
    random: simple local search for larger n (correct but not guaranteed optimal).
    """
    random.seed(seed)
    seq0 = list(base_sequence)
    n = len(seq0)

    best_seq = seq0
    best_run = run_model_multi(
        sequence=best_seq,
        result_folder=result_folder,
        message=message + " | baseline",
        Scenario=Scenario + "_baseline",
        plot_control=False,
        focus=False,
        strict=strict,
        debug=debug,
        crew_mode=crew_mode,
        multifunction_crews=multifunction_crews,
        power_crews=power_crews,
        road_crews=road_crews,
        depot_node=depot_node,
        crew_speed=crew_speed,
        service_time_power=service_time_power,
        service_time_road=service_time_road,
        base_net=base_net,
        trips=trips,
        net1=net1,
        net2=net2,
        broken_link_factor=broken_link_factor,
        power_road_factor=power_road_factor,
        baseline_tstt=baseline_tstt,
        objective=objective,
    )
    best_val = float(best_run["objective_value"])

    if method == "bruteforce":
        if n > 8:
            raise ValueError("bruteforce too large; use method='random' or reduce n (suggest n<=8).")
        from itertools import permutations
        for perm in permutations(seq0, n):
            cand = list(perm)
            run = run_model_multi(
                sequence=cand,
                result_folder=result_folder,
                message=message + " | candidate",
                Scenario=Scenario + "_cand",
                plot_control=False,
                focus=False,
                strict=strict,
                debug=False,
                crew_mode=crew_mode,
                multifunction_crews=multifunction_crews,
                power_crews=power_crews,
                road_crews=road_crews,
                depot_node=depot_node,
                crew_speed=crew_speed,
                service_time_power=service_time_power,
                service_time_road=service_time_road,
                base_net=base_net,
                trips=trips,
                net1=net1,
                net2=net2,
                broken_link_factor=broken_link_factor,
                power_road_factor=power_road_factor,
                baseline_tstt=baseline_tstt,
                objective=objective,
            )
            val = float(run["objective_value"])
            if val < best_val:
                best_val = val
                best_seq = cand
                best_run = run

    elif method == "random":
        # simple swap-based local search
        cand = seq0[:]
        for it in range(int(max_iter)):
            i, j = random.sample(range(n), 2)
            new = cand[:]
            new[i], new[j] = new[j], new[i]
            run = run_model_multi(
                sequence=new,
                result_folder=result_folder,
                message=message + f" | iter={it}",
                Scenario=Scenario + "_iter",
                plot_control=False,
                focus=False,
                strict=strict,
                debug=False,
                crew_mode=crew_mode,
                multifunction_crews=multifunction_crews,
                power_crews=power_crews,
                road_crews=road_crews,
                depot_node=depot_node,
                crew_speed=crew_speed,
                service_time_power=service_time_power,
                service_time_road=service_time_road,
                base_net=base_net,
                trips=trips,
                net1=net1,
                net2=net2,
                broken_link_factor=broken_link_factor,
                power_road_factor=power_road_factor,
                baseline_tstt=baseline_tstt,
                objective=objective,
            )
            val = float(run["objective_value"])
            if val < best_val:
                best_val = val
                best_seq = new
                best_run = run
                cand = new

    else:
        raise ValueError(f"Unknown method: {method!r}")

    # final record
    best_run["best_sequence"] = best_seq
    best_run["best_objective_value"] = best_val
    return best_run
