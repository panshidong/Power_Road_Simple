from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Dict, List, Tuple, Iterable, Any

from road_util import calculate_shortest_path_cost


def load_json(path: str | Path) -> Any:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Missing file: {p}")
    return json.loads(p.read_text(encoding="utf-8"))


def gini(values: List[float]) -> float:
    xs = [float(x) for x in values if x is not None and not math.isnan(float(x))]
    n = len(xs)
    if n == 0:
        return float("nan")
    xs.sort()
    s = sum(xs)
    if s == 0:
        return 0.0
    # Gini = (2*sum(i*x_i)/(n*sum x)) - (n+1)/n
    cum = 0.0
    for i, x in enumerate(xs, start=1):
        cum += i * x
    return (2.0 * cum) / (n * s) - (n + 1.0) / n


def percentile(values: List[float], p: float) -> float:
    xs = sorted([float(x) for x in values if x is not None and not math.isnan(float(x))])
    if not xs:
        return float("nan")
    if p <= 0:
        return xs[0]
    if p >= 100:
        return xs[-1]
    k = (len(xs) - 1) * (p / 100.0)
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return xs[int(k)]
    return xs[f] * (c - k) + xs[c] * (k - f)


def build_zone_lists(node_to_zone: Dict[str, int]) -> List[int]:
    # zones are ints; node_to_zone maps node(str)->zone(int)
    zones = sorted(set(int(z) for z in node_to_zone.values()))
    return zones


def compute_accessibility_TT(
    *,
    s_txt_path: str,
    zones: List[int],
    destinations: List[int],
) -> Dict[int, float]:
    """
    TT_z = min_{d in destinations} dist(z, d) using s.txt-derived directed graph.
    Here zone==node.
    """
    out: Dict[int, float] = {}
    for z in zones:
        best = float("inf")
        for d in destinations:
            cost = calculate_shortest_path_cost(s_txt_path, int(z), int(d))
            if cost < best:
                best = cost
        out[int(z)] = float(best)
    return out


def compute_accessibility_ratio_by_zone(
    *,
    zones: List[int],
    TT0: Dict[int, float],
    TT: Dict[int, float],
) -> Dict[int, float]:
    """
    A_z(t) = min(1, TT0_z / TT_z) with safety:
      - if TT is inf or <=0: A=0
      - if TT0 is inf: A=0
    """
    out: Dict[int, float] = {}
    for z in zones:
        tt0 = float(TT0.get(z, float("inf")))
        tt = float(TT.get(z, float("inf")))

        if tt0 == 0 and tt == 0:
            a = 1.0
        elif math.isinf(tt0) or math.isinf(tt) or tt <= 0:
            a = 0.0
        else:
            a = min(1.0, tt0 / tt)
        out[int(z)] = float(a)
    return out


def compute_E_by_zone(
    *,
    zones: List[int],
    bus_to_zone: Dict[str, int],
    functional_buses: Iterable[int],
    empty_zone_policy: str = "exclude",  # "exclude" or "zero"
) -> Dict[int, float]:
    """
    E_z(t) = (# functional buses in zone z) / (# buses in zone z).
    If a zone has 0 buses:
      - "exclude": not included in E output
      - "zero": E_z = 0
    """
    func_set = set(int(b) for b in functional_buses)
    total: Dict[int, int] = {z: 0 for z in zones}
    good: Dict[int, int] = {z: 0 for z in zones}

    for b_str, z in bus_to_zone.items():
        z = int(z)
        if z not in total:
            continue
        total[z] += 1
        if int(b_str) in func_set:
            good[z] += 1

    out: Dict[int, float] = {}
    for z in zones:
        if total[z] == 0:
            if empty_zone_policy == "zero":
                out[z] = 0.0
            else:
                continue
        else:
            out[z] = good[z] / float(total[z])
    return out


def compute_CRI_by_zone(
    *,
    zones: List[int],
    E: Dict[int, float],
    TT0: Dict[int, float],
    TT: Dict[int, float],
    w_e: float = 0.5,
    w_a: float = 0.5,
) -> Dict[int, float]:
    """
    CRI = w_e*E + w_a*A. Missing E defaults to 0.
    """
    A = compute_accessibility_ratio_by_zone(zones=zones, TT0=TT0, TT=TT)
    out: Dict[int, float] = {}
    for z in zones:
        e = float(E.get(z, 0.0))
        a = float(A.get(z, 0.0))

        cri = w_e * e + w_a * a
        # clamp
        cri = max(0.0, min(1.0, cri))
        out[z] = cri
    return out


def first_hit_times_from_series(
    *,
    zones: List[int],
    time_series: List[float],
    value_series: List[Dict[int, float]],  # aligned with time_series
    threshold: float = 0.9,
) -> Dict[int, float]:
    """
    First time when value_z(t) >= threshold.
    If never reaches, returns final time.
    """
    out: Dict[int, float] = {}
    Tfinal = float(time_series[-1]) if time_series else 0.0
    for z in zones:
        t_hit = None
        for t, val_map in zip(time_series, value_series):
            if float(val_map.get(z, 0.0)) >= float(threshold):
                t_hit = float(t)
                break
        out[z] = t_hit if t_hit is not None else Tfinal
    return out


def restoration_times_from_series(
    *,
    zones: List[int],
    time_series: List[float],
    cri_series: List[Dict[int, float]],
    threshold: float = 0.9,
) -> Dict[int, float]:
    return first_hit_times_from_series(
        zones=zones,
        time_series=time_series,
        value_series=cri_series,
        threshold=threshold,
    )


def time_average_by_zone(
    *,
    zones: List[int],
    time_series: List[float],
    value_series: List[Dict[int, float]],
) -> Dict[int, float]:
    """
    Piecewise-constant time average for each zone over the recovery horizon.
    """
    if not zones:
        return {}
    if not time_series or not value_series:
        return {z: float("nan") for z in zones}
    if len(time_series) == 1:
        return {z: float(value_series[0].get(z, 0.0)) for z in zones}

    total_time = float(time_series[-1] - time_series[0])
    if total_time <= 0:
        return {z: float(value_series[-1].get(z, 0.0)) for z in zones}

    out: Dict[int, float] = {}
    for z in zones:
        acc = 0.0
        for i in range(len(time_series) - 1):
            dt = float(time_series[i + 1] - time_series[i])
            acc += float(value_series[i].get(z, 0.0)) * dt
        out[z] = acc / total_time
    return out


def equity_summary_from_CRI(
    *,
    zones: List[int],
    time_series: List[float],
    cri_series: List[Dict[int, float]],
    threshold: float = 0.9,
) -> Dict[str, float]:
    """
    Returns equity metrics based on:
      - restoration times T_z(threshold)
      - final CRI distribution
    """
    Tz = restoration_times_from_series(zones=zones, time_series=time_series, cri_series=cri_series, threshold=threshold)
    tvals = [Tz[z] for z in zones]

    final_cri = cri_series[-1] if cri_series else {z: 0.0 for z in zones}
    cvals = [float(final_cri.get(z, 0.0)) for z in zones]
    avg_cri = time_average_by_zone(zones=zones, time_series=time_series, value_series=cri_series)
    avg_vals = [float(avg_cri.get(z, 0.0)) for z in zones]

    out = {
        # restoration-time dispersion
        "var_restore": float(_variance(tvals)),
        "gini_restore": float(gini(tvals)),
        "p90_restore": float(percentile(tvals, 90)),
        "p95_restore": float(percentile(tvals, 95)),
        # maximin-style CRI summary
        "min_time_avg_cri": float(min(avg_vals)) if avg_vals else float("nan"),
        "maximin_time_avg_cri_loss": float(1.0 - min(avg_vals)) if avg_vals else float("nan"),
        # final CRI dispersion (optional but useful)
        "min_final_cri": float(min(cvals)) if cvals else float("nan"),
        "maximin_final_cri_loss": float(1.0 - min(cvals)) if cvals else float("nan"),
        "var_final_cri": float(_variance(cvals)),
        "gini_final_cri": float(gini(cvals)),
        "p10_final_cri": float(percentile(cvals, 10)),
        "p90_final_cri": float(percentile(cvals, 90)),
    }
    return out


def critical_access_summary_from_series(
    *,
    zones: List[int],
    time_series: List[float],
    access_series: List[Dict[int, float]],
    threshold: float = 0.9,
) -> Dict[str, float]:
    """
    Returns critical-access metrics based on A_z(t) and an adequacy threshold.
    """
    Tz = first_hit_times_from_series(
        zones=zones,
        time_series=time_series,
        value_series=access_series,
        threshold=threshold,
    )
    tvals = [Tz[z] for z in zones]

    final_access = access_series[-1] if access_series else {z: 0.0 for z in zones}
    avals = [float(final_access.get(z, 0.0)) for z in zones]
    reached = [1.0 if float(final_access.get(z, 0.0)) >= float(threshold) else 0.0 for z in zones]
    shares = [
        sum(1.0 if float(access_map.get(z, 0.0)) >= float(threshold) else 0.0 for z in zones) / float(len(zones))
        for access_map in access_series
    ] if zones else []

    if len(time_series) >= 2 and shares:
        total_time = float(time_series[-1] - time_series[0])
        if total_time > 0:
            time_avg_share = sum(
                float(shares[i]) * float(time_series[i + 1] - time_series[i])
                for i in range(len(time_series) - 1)
            ) / total_time
        else:
            time_avg_share = float(shares[-1])
    elif shares:
        time_avg_share = float(shares[-1])
    else:
        time_avg_share = float("nan")

    out = {
        "var_access_restore": float(_variance(tvals)),
        "gini_access_restore": float(gini(tvals)),
        "p90_access_restore": float(percentile(tvals, 90)),
        "p95_access_restore": float(percentile(tvals, 95)),
        "share_access_initial": float(shares[0]) if shares else float("nan"),
        "mean_access_final": float(sum(avals) / len(avals)) if avals else float("nan"),
        "share_access_final": float(sum(reached) / len(reached)) if reached else float("nan"),
        "time_avg_share_access": float(time_avg_share),
    }
    return out


def _variance(xs: List[float]) -> float:
    ys = [float(x) for x in xs if x is not None and not math.isnan(float(x))]
    n = len(ys)
    if n == 0:
        return float("nan")
    mu = sum(ys) / n
    return sum((y - mu) ** 2 for y in ys) / n
