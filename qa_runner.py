from __future__ import annotations

import json
import os
import random
import shutil
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Tuple, Set, Optional

from road_util import eval_tot_OD_travel_time
from power_util import get_functional_nodes
from resilience_measurement import (
    _prepare_state_and_run_tapb,
    eval_road_resilience,
    eval_power_resilience,
    _plot_triangle_with_shading,
    DEFAULT_BASE_NET,
    DEFAULT_TRIPS,
    DEFAULT_NET1,
    DEFAULT_NET2,
)

# ----------------------------
# Config
# ----------------------------

@dataclass
class QAConfig:
    seed: int = 0
    n_assets: int = 6           # number of assets in random sequence
    n_buses: int = 3            # how many buses among selected assets
    n_links: int = 3            # how many road links among selected assets
    depot_node: int = 1

    strict: bool = True

    # model params
    broken_link_factor: float = 0.0
    power_road_factor: float = 0.5
    baseline_tstt: float = 7475338.0

    # file paths
    base_net: str = DEFAULT_BASE_NET
    trips: str = DEFAULT_TRIPS
    net1: str = DEFAULT_NET1
    net2: str = DEFAULT_NET2

    # TaskB mapping files
    node_to_zone_path: str = "taskB_node_to_zone.json"
    bus_to_zone_path: str = "taskB_bus_to_zone.json"
    dest_path: str = "taskB_essential_destinations.json"
    bus_location_path: str = "bus_location.json"


def ts() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def load_json(path: str) -> Any:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Missing file: {path}")
    return json.loads(p.read_text(encoding="utf-8"))


def parse_nodes_from_net(netfile: str) -> Set[int]:
    nodes: Set[int] = set()
    with open(netfile, "r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("~") or s.startswith(";"):
                continue
            parts = s.split()
            if len(parts) < 2:
                continue
            try:
                u = int(parts[0]); v = int(parts[1])
            except Exception:
                continue
            nodes.add(u); nodes.add(v)
    if not nodes:
        raise RuntimeError(f"No nodes parsed from {netfile}")
    return nodes


def parse_links_from_net(netfile: str) -> List[Tuple[int,int]]:
    links: Set[Tuple[int,int]] = set()
    with open(netfile, "r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("~") or s.startswith(";"):
                continue
            parts = s.split()
            if len(parts) < 2:
                continue
            try:
                u = int(parts[0]); v = int(parts[1])
            except Exception:
                continue
            links.add((u,v))
    if not links:
        raise RuntimeError(f"No links parsed from {netfile}")
    return sorted(list(links))


def validate_taskB_inputs(cfg: QAConfig) -> Dict[str, Any]:
    """
    Fail-fast validation:
      - node_to_zone keys are valid nodes
      - bus_location nodes are valid nodes
      - bus_to_zone agrees with bus_location (optional but strong check)
      - destinations are valid nodes
      - every zone has >=1 bus (per your current assumption)
    """
    nodes = parse_nodes_from_net(cfg.base_net)
    node_to_zone = load_json(cfg.node_to_zone_path)
    bus_to_zone = load_json(cfg.bus_to_zone_path)
    dest_obj = load_json(cfg.dest_path)
    bus_loc = load_json(cfg.bus_location_path)

    # 1) node_to_zone validity
    bad_nodes = [int(k) for k in node_to_zone.keys() if int(k) not in nodes]
    if bad_nodes:
        raise RuntimeError(f"node_to_zone contains nodes not in network: {bad_nodes[:20]} ...")

    zones = sorted(set(int(z) for z in node_to_zone.values()))

    # 2) bus_location validity
    bad_busloc = []
    for b, n in bus_loc.items():
        if int(n) not in nodes:
            bad_busloc.append((int(b), int(n)))
    if bad_busloc:
        raise RuntimeError(f"bus_location has invalid node(s): {bad_busloc[:20]} ...")

    # 3) destinations validity
    dests = list(dest_obj["destinations"])
    bad_dest = [int(d) for d in dests if int(d) not in nodes]
    if bad_dest:
        raise RuntimeError(f"essential destinations invalid: {bad_dest}")

    # 4) bus_to_zone vs bus_location (strong consistency check)
    # since your current rule is bus_to_zone generated from bus_location, enforce it
    mismatch = []
    for b, n in bus_loc.items():
        z_expected = int(n)  # zone == node
        z = bus_to_zone.get(str(int(b)))
        if z is None:
            mismatch.append((int(b), "missing_in_bus_to_zone"))
        elif int(z) != z_expected:
            mismatch.append((int(b), int(z), z_expected))
    if mismatch:
        raise RuntimeError(f"bus_to_zone inconsistent with bus_location (first 20): {mismatch[:20]}")

    # 5) each zone has >=1 bus
    bus_counts = {z: 0 for z in zones}
    for b, z in bus_to_zone.items():
        z = int(z)
        if z in bus_counts:
            bus_counts[z] += 1
    empty = [z for z, c in bus_counts.items() if c == 0]
    if empty:
        raise RuntimeError(
            f"Found zone(s) with zero buses: {empty[:20]} ... "
            f"Per current assumption this must error. Fix bus_location/bus_to_zone."
        )

    return {
        "nodes_count": len(nodes),
        "zones_count": len(zones),
        "destinations": dests,
        "bus_count": len(bus_loc),
    }


def random_sequence(cfg: QAConfig) -> List[Any]:
    rng = random.Random(cfg.seed)
    # buses come from bus_location keys
    bus_loc = load_json(cfg.bus_location_path)
    all_buses = sorted([int(k) for k in bus_loc.keys()])

    # links come from network file
    links = parse_links_from_net(cfg.base_net)

    if cfg.n_buses + cfg.n_links != cfg.n_assets:
        raise ValueError("n_assets must equal n_buses + n_links for this QA generator.")

    if cfg.n_buses > len(all_buses):
        raise ValueError("n_buses larger than bus pool.")
    if cfg.n_links > len(links):
        raise ValueError("n_links larger than link pool.")

    buses = rng.sample(all_buses, cfg.n_buses)
    lks = rng.sample(links, cfg.n_links)
    seq: List[Any] = buses + lks
    rng.shuffle(seq)
    return seq


def qa_evaluate_sequence(cfg: QAConfig, sequence: List[Any], out_dir: str) -> str:
    """
    Run one full evaluation (event-driven), log each step:
      - repaired asset
      - TSTT after re-run
      - # functional buses
      - road/power functionality
      - time
    Also saves shaded resilience triangle plot.
    """
    os.makedirs(out_dir, exist_ok=True)

    # Ensure mapping files used by your main code are present (as in your pipeline)
    if os.path.exists("bus_location.json"):
        os.remove("bus_location.json")
    if os.path.exists("bus_to_link.json"):
        os.remove("bus_to_link.json")
    shutil.copy2("original_bus_to_link.json", "bus_to_link.json")
    shutil.copy2("original_bus_location.json", "bus_location.json")

    # initial broken sets derived from sequence
    seq = list(sequence)
    broken_buses = {a for a in seq if not isinstance(a, tuple)}
    broken_links = {a for a in seq if isinstance(a, tuple)}

    # We need a dispatch snapshot to build a consistent schedule (same as main run)
    _prepare_state_and_run_tapb(
        broken_buses=list(broken_buses),
        broken_links=list(broken_links),
        base_net=cfg.base_net,
        trips=cfg.trips,
        net1=cfg.net1,
        net2=cfg.net2,
        broken_link_factor=cfg.broken_link_factor,
        power_road_factor=cfg.power_road_factor,
        strict=cfg.strict,
    )
    dispatch_s = os.path.join(out_dir, "dispatch_s.txt")
    shutil.copy2("s.txt", dispatch_s)

    # schedule (uses dispatch_s)
    from scheduler import evaluate_with_crews
    from crews import make_crews

    pool = make_crews(mode="multifunction", multifunction_crews=1, depot_node=cfg.depot_node, speed=1.0)

    sim = evaluate_with_crews(
        sequence=seq,
        crew_pool=pool,
        service_time_power=20.0,
        service_time_road=10.0,
        default_depot=cfg.depot_node,
        s_txt_path=dispatch_s,
        bus_location_path="bus_location.json",
        bus_to_link_path="bus_to_link.json",  # optional fallback
        strict=cfg.strict,
    )
    events = sorted(sim["timeline"], key=lambda x: x[0])

    # step-by-step log
    log_csv = os.path.join(out_dir, "qa_step_log.csv")
    with open(log_csv, "w", encoding="utf-8") as f:
        f.write("step,time,repair_asset,tstt,functional_bus_count,road_func,power_func,broken_bus_count,broken_link_count\n")

    # initial evaluation
    road_func = eval_road_resilience(
        list(broken_buses),
        list(broken_links),
        base_net=cfg.base_net,
        trips=cfg.trips,
        net1=cfg.net1,
        net2=cfg.net2,
        broken_link_factor=cfg.broken_link_factor,
        power_road_factor=cfg.power_road_factor,
        baseline_tstt=cfg.baseline_tstt,
        strict=cfg.strict,
    )
    power_func = eval_power_resilience(list(broken_buses))
    tstt = float(eval_tot_OD_travel_time("s.txt"))
    func_buses = set(get_functional_nodes(set(map(int, broken_buses))))
    func_bus_count = len(func_buses)

    time_series = [0.0]
    road_series = [road_func]
    power_series = [power_func]

    with open(log_csv, "a", encoding="utf-8") as f:
        f.write(f"0,0.0,INIT,{tstt},{func_bus_count},{road_func},{power_func},{len(broken_buses)},{len(broken_links)}\n")

    # iterate events
    for idx, (t, asset) in enumerate(events, start=1):
        # apply repair
        if isinstance(asset, tuple):
            broken_links.discard(asset)
        else:
            broken_buses.discard(asset)

        # re-evaluate after repair (this updates s.txt)
        road_func = eval_road_resilience(
            list(broken_buses),
            list(broken_links),
            base_net=cfg.base_net,
            trips=cfg.trips,
            net1=cfg.net1,
            net2=cfg.net2,
            broken_link_factor=cfg.broken_link_factor,
            power_road_factor=cfg.power_road_factor,
            baseline_tstt=cfg.baseline_tstt,
            strict=cfg.strict,
        )
        power_func = eval_power_resilience(list(broken_buses))
        tstt = float(eval_tot_OD_travel_time("s.txt"))
        func_buses = set(get_functional_nodes(set(map(int, broken_buses))))
        func_bus_count = len(func_buses)

        time_series.append(float(t))
        road_series.append(road_func)
        power_series.append(power_func)

        with open(log_csv, "a", encoding="utf-8") as f:
            f.write(f"{idx},{float(t)},{repr(asset)},{tstt},{func_bus_count},{road_func},{power_func},{len(broken_buses)},{len(broken_links)}\n")

    # plot shaded triangle
    tri_png = os.path.join(out_dir, "qa_triangle.png")
    _plot_triangle_with_shading(
        time_series=time_series,
        road_series=road_series,
        power_series=power_series,
        out_png=tri_png,
        title=f"QA Triangle | steps={len(events)} | seed={cfg.seed} | {ts()}",
    )

    # simple QA summary text
    summary_txt = os.path.join(out_dir, "qa_summary.txt")
    with open(summary_txt, "w", encoding="utf-8") as f:
        f.write("QA Summary\n")
        f.write(f"- timestamp: {ts()}\n")
        f.write(f"- seed: {cfg.seed}\n")
        f.write(f"- sequence: {seq}\n")
        f.write(f"- steps: {len(events)}\n")
        f.write(f"- outputs:\n  - {log_csv}\n  - {tri_png}\n")

    return out_dir


def main():
    cfg = QAConfig(seed=0, n_assets=6, n_buses=3, n_links=3)

    # 1) Input QA
    info = validate_taskB_inputs(cfg)

    # 2) Random QA sequence
    seq = random_sequence(cfg)

    # 3) Run QA evaluation
    out_dir = os.path.join("results", f"QA_{ts()}")
    qa_evaluate_sequence(cfg, seq, out_dir)

    print("QA done.")
    print("Input check:", info)
    print("Sequence:", seq)
    print("Output dir:", out_dir)
    print("See qa_step_log.csv and qa_triangle.png")


if __name__ == "__main__":
    main()
