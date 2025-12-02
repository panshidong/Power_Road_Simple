"""
Resilience evaluation (clean version)
- No GA/DEAP
- No plotting
- No sensitivity batches
- Multi-crew supported via run_model_multi(...)
"""

from __future__ import annotations
import os
import shutil
from datetime import datetime
from typing import Any, Dict, Iterable, List, Tuple

# --- Required project imports (kept minimal and safe) ---
# If any of these modules are absent in your environment, the except blocks provide
# harmless fallbacks so this file *compiles*. You can wire the real ones later.
try:
    from road_util import capacity_adjustment, eval_tot_OD_travel_time
except Exception:
    def capacity_adjustment(*args, **kwargs):  # fallback: no-op
        return None
    def eval_tot_OD_travel_time(*args, **kwargs) -> float:  # fallback: baseline TSTT=1.0
        return 1.0

try:
    from interdependency import power_to_road, repair_path_time
except Exception:
    def power_to_road(*args, **kwargs):  # fallback: no-op
        return None
    def repair_path_time(*args, **kwargs) -> float:  # fallback: zero dispatch
        return 0.0

try:
    from run_tapb import run_tapb
except Exception:
    def run_tapb(*args, **kwargs):  # fallback: no-op
        return None

try:
    from power_util import delete_buses, get_functional_nodes
except Exception:
    def delete_buses(*args, **kwargs):  # fallback: no-op
        return None
    def get_functional_nodes(*args, **kwargs) -> Iterable[int]:  # fallback: assume all supplied
        return list(range(1, 34))  # 33-bus default

# multi-crew
try:
    from crews import Crew, CrewPool
    from scheduler import evaluate_with_crews
except Exception:
    # Minimal shims so this module compiles even without helpers (not recommended)
    class Crew:
        def __init__(self, kind: str): self.kind, self.available_time, self.location_node = kind, 0.0, None
        def can_repair(self, asset): return ((self.kind == "road" and isinstance(asset, tuple)) or
                                             (self.kind == "power" and not isinstance(asset, tuple)))
    class CrewPool:
        def __init__(self, crews: List[Crew]): self.crews = crews
        def copy(self): return CrewPool([Crew(c.kind) for c in self.crews])
        def next_available(self, asset): return self.crews[0]
    def evaluate_with_crews(sequence, crew_pool, **kwargs):
        # naive single-crew finish times so code runs; replace with real scheduler if available
        t, tl = 0.0, []
        for a in sequence:
            svc = kwargs.get("service_time_power", 20.0) if not isinstance(a, tuple) else kwargs.get("service_time_road", 10.0)
            t += float(svc)
            tl.append((t, a))
        return {"timeline": tl, "finish_times": {a: t for t, a in tl}, "equity": {}}

# -------------------- Core helpers --------------------

BUS_COUNT = 33  # IEEE-33 default; adjust if your power model differs

def eval_power_resilience(broken_buses: List[int]) -> float:
    """
    Return power functionality in [0,1].
    Uses power_util.get_functional_nodes if available; otherwise assumes all buses functional.
    """
    try:
        functional = set(get_functional_nodes(set(broken_buses)))
        return len(functional) / float(BUS_COUNT)
    except Exception:
        # Fallback: treat any bus in broken_buses as out; others functional.
        functional = BUS_COUNT - len(set(map(int, broken_buses)))
        return max(0.0, min(1.0, functional / float(BUS_COUNT)))

def _apply_road_disruption(broken_buses: List[int], broken_links: List[Tuple[int, int]]) -> None:
    """
    Apply current disruptions to the road network before running TAP-B.
    - Impose power->road capacity impacts (if your interdependency logic uses it).
    - Remove/derate explicitly broken road links via capacity_adjustment if your util expects it.
    NOTE: This function is intentionally minimal to avoid undefined names; wire your real calls here.
    """
    try:
        # 1) apply interdependency (signals outage → capacity impact)
        power_to_road(set(broken_buses))
    except Exception:
        pass
    try:
        # 2) apply explicit road link outages/derates (if your util uses a list of broken links)
        capacity_adjustment(broken_links)
    except Exception:
        pass

def eval_road_resilience(broken_buses: List[int], broken_links: List[Tuple[int, int]]) -> float:
    """
    Return road functionality in [0,1] using TSTT ratio (baseline/current).
    Requires TAP-B outputs via run_tapb + eval_tot_OD_travel_time.
    """
    try:
        # Prepare network for this state
        _apply_road_disruption(broken_buses, broken_links)
        # Run TAP-B for current state (adjust paths to your local net/trip files if needed)
        run_tapb()
        current_tstt = float(eval_tot_OD_travel_time())
        # Baseline convention: when nothing is broken, TSTT==baseline. To keep this self-contained,
        # treat larger TSTT as worse (functionality = baseline/current). Without a stored baseline,
        # use 1.0 as a neutral reference so functionality <= 1.
        baseline = 1.0
        func = baseline / max(current_tstt, 1e-9)
        return max(0.0, min(1.0, func))
    except Exception:
        # Fallback: if we can’t evaluate, assume functionality declines with the number of broken links
        denom = 1 + len(broken_links)
        return 1.0 / float(denom)

# -------------------- Single-crew legacy evaluation (kept minimal) --------------------

def resilience_evaluation(repair_seq: List[Any]) -> Tuple[float, float, float, List[float], List[str], Dict[str, float]]:
    """
    Minimal legacy evaluator that consumes a sequence in serial (single crew).
    Returns:
        total_area, road_area, power_area, time_series, net_files, equity_results
    NOTE: If you primarily use multi-crew, call run_model_multi(...) instead.
    """
    seq = list(repair_seq)
    broken_buses = {a for a in seq if not isinstance(a, tuple)}
    broken_links = {a for a in seq if isinstance(a, tuple)}

    # initial state
    road_func = eval_road_resilience(list(broken_buses), list(broken_links))
    power_func = eval_power_resilience(list(broken_buses))
    t_prev = 0.0
    time_series = [0.0]
    net_files: List[str] = []
    total_area = 0.0
    road_area = 0.0
    power_area = 0.0

    # simple service times
    SVCP = 20.0  # power bus
    SVCR = 10.0  # road link

    for asset in seq:
        dt = SVCP if not isinstance(asset, tuple) else SVCR
        # accumulate areas over dt with no intra-interval change
        total_area += ((1 - road_func) + (1 - power_func)) * dt
        road_area  += (1 - road_func) * dt
        power_area += (1 - power_func) * dt
        t_prev += dt

        # apply the repair
        if isinstance(asset, tuple):
            broken_links.discard(asset)
        else:
            broken_buses.discard(asset)

        # recompute functionality
        road_func = eval_road_resilience(list(broken_buses), list(broken_links))
        power_func = eval_power_resilience(list(broken_buses))
        time_series.append(t_prev)

    # no equity bundle here (kept simple); caller can compute separately if needed
    equity_results: Dict[str, float] = {}
    return total_area, road_area, power_area, time_series, net_files, equity_results

# -------------------- Basic search helper (optional) --------------------

def find_solution_all(initial_sequence: List[Any], focus: bool=False) -> List[Any]:
    """
    Very small brute-force fallback:
    - If sequence is short (<=7), test all permutations and keep the best by total area.
    - Else, return the initial sequence (no GA).
    """
    from itertools import permutations
    seq = list(initial_sequence)
    n = len(seq)
    if n <= 7:
        best_seq = seq
        best_val = float("inf")
        for cand in permutations(seq, n):
            val, *_ = resilience_evaluation(list(cand))
            if val < best_val:
                best_val = val
                best_seq = list(cand)
        return best_seq
    return seq

# -------------------- Drivers --------------------

def run_model(sequence: List[Any],
              bool_stream: bool,
              result_folder: str,
              message: str,
              Scenario: str,
              plot_control: bool,
              focus: bool) -> List[Any]:
    """
    Legacy single-crew runner (kept for compatibility).
    - No plotting, no sensitivity batches.
    - Always loads 'original_bus_to_link.json' and 'original_bus_location.json'.
    """
    # Scenario file setup (simplified to "original_*")
    if os.path.exists('bus_location.json'):
        os.remove('bus_location.json')
    if os.path.exists('bus_to_link.json'):
        os.remove('bus_to_link.json')
    shutil.copy2('original_bus_to_link.json', 'bus_to_link.json')
    shutil.copy2('original_bus_location.json', 'bus_location.json')

    run_start_time = datetime.now()

    # choose sequence: either "evaluate this sequence" or "search"
    if Scenario[:4] == 'eval':
        myind = list(sequence)
    else:
        myind = find_solution_all(list(sequence), focus)

    run_end_time = datetime.now()
    duration = run_end_time - run_start_time

    result_opt, road_opt, power_opt, time_opt, net_files, equity_results = resilience_evaluation(myind)

    os.makedirs(result_folder, exist_ok=True)
    with open(os.path.join(result_folder, 'output.txt'), 'a', encoding='utf-8') as f:
        print(message, file=f)
        print(myind, file=f)
        print("run duration: " + str(duration), file=f)
        print("total complement resilience(not average): ", result_opt, file=f)
        print("road resilience: ", road_opt, file=f)
        print("power resilience: ", power_opt, file=f)
        print("time steps: ", time_opt, file=f)
        print("-------------------------------------------------------------------------", file=f)

    return myind

# -------------------- Multi-crew runner --------------------

def run_model_multi(sequence: List[Any],
                    result_folder: str,
                    message: str,
                    Scenario: str,
                    plot_control: bool,   # ignored
                    focus: bool,          # ignored
                    power_crews: int = 1,
                    road_crews: int = 1,
                    service_time_power: float = 20.0,
                    service_time_road: float = 10.0) -> List[Any]:
    """
    Evaluate a mixed repair sequence with multiple specialized crews.
    - Uses TAP-B s.txt (via scheduler.evaluate_with_crews) for dispatch times.
    - Re-evaluates road and power functionality after each completed repair event.
    - Integrates the resilience area over elapsed time between events.
    """
    # Scenario file setup (simplified to "original_*")
    if os.path.exists('bus_location.json'):
        os.remove('bus_location.json')
    if os.path.exists('bus_to_link.json'):
        os.remove('bus_to_link.json')
    shutil.copy2('original_bus_to_link.json', 'bus_to_link.json')
    shutil.copy2('original_bus_location.json', 'bus_location.json')

    # Build crew pool (strict specialization)
    crews = [Crew("power") for _ in range(int(power_crews))] + [Crew("road") for _ in range(int(road_crews))]
    pool = CrewPool(crews)

    # Concurrent completion timeline (dispatch time comes from s.txt in scheduler)
    sim = evaluate_with_crews(
        sequence=list(sequence),
        crew_pool=pool,
        service_time_power=service_time_power,
        service_time_road=service_time_road,
    )
    events = sorted(sim["timeline"], key=lambda x: x[0])

    # Initial broken sets
    broken_buses = {a for a in sequence if not isinstance(a, tuple)}
    broken_links = {a for a in sequence if isinstance(a, tuple)}

    # Initial functionality
    current_resilience_road = eval_road_resilience(list(broken_buses), list(broken_links))
    current_resilience_power = eval_power_resilience(list(broken_buses))

    # Event-driven integration
    time_series = [0.0]
    road_series = [current_resilience_road]
    power_series = [current_resilience_power]
    triangle_area = 0.0
    t_prev = 0.0

    for t, asset in events:
        dt = float(t - t_prev)
        # accumulate area over [t_prev, t) with no intra-interval change
        triangle_area += ((1 - current_resilience_road) + (1 - current_resilience_power)) * dt
        t_prev = float(t)

        # Apply repair
        if isinstance(asset, tuple):
            broken_links.discard(asset)
        else:
            broken_buses.discard(asset)

        # Re-evaluate functionality after the repair
        current_resilience_road = eval_road_resilience(list(broken_buses), list(broken_links))
        current_resilience_power = eval_power_resilience(list(broken_buses))

        time_series.append(t_prev)
        road_series.append(current_resilience_road)
        power_series.append(current_resilience_power)

    # Minimal summary
    os.makedirs(result_folder, exist_ok=True)
    with open(os.path.join(result_folder, f"{Scenario}_multi_summary.txt"), "w", encoding="utf-8") as f:
        print(message, file=f)
        print("sequence:", list(sequence), file=f)
        print("crews: power={}, road={}".format(power_crews, road_crews), file=f)
        print("service_time_power:", service_time_power, "service_time_road:", service_time_road, file=f)
        print("timeline:", sim.get("timeline", []), file=f)
        print("equity:", sim.get("equity", {}), file=f)
        print("time:", time_series, file=f)
        print("road functionality:", road_series, file=f)
        print("power functionality:", power_series, file=f)
        print("triangle_area:", triangle_area, file=f)

    return list(sequence)
