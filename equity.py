
from __future__ import annotations
import math
from typing import Iterable, List, Optional, Dict

def _as_list(values: Iterable[float]) -> List[float]:
    xs = [float(x) for x in values]
    return xs

def gini(values: Iterable[float]) -> Optional[float]:
    xs = _as_list(values)
    n = len(xs)
    if n == 0:
        return None
    if all(v == 0 for v in xs):
        return 0.0
    xs = sorted(xs)
    cum = 0.0
    for i, x in enumerate(xs, start=1):
        cum += i * x
    s = sum(xs)
    return (2 * cum) / (n * s) - (n + 1) / n

def theil_T(values: Iterable[float]) -> Optional[float]:
    xs = _as_list(values)
    n = len(xs)
    if n == 0:
        return None
    mean = sum(xs) / n if n else 0.0
    if mean <= 0:
        return 0.0
    T = 0.0
    for x in xs:
        if x > 0:
            r = x / mean
            T += r * math.log(r)
    return T / n

def atkinson(values: Iterable[float], epsilon: float = 0.5) -> Optional[float]:
    xs = _as_list(values)
    n = len(xs)
    if n == 0:
        return None
    if epsilon == 1.0:
        # limit case
        geo = 1.0
        cnt = 0
        for x in xs:
            if x > 0:
                geo *= x
                cnt += 1
        if cnt == 0:
            return 0.0
        geo = geo ** (1.0 / cnt)
        mean = sum(xs) / n
        return 1.0 - geo / mean if mean > 0 else 0.0
    mean = sum(xs) / n
    if mean <= 0:
        return 0.0
    s = sum((x ** (1.0 - epsilon) for x in xs if x >= 0))
    Ae = (s / n) ** (1.0 / (1.0 - epsilon))
    return 1.0 - Ae / mean

def jain(values: Iterable[float]) -> Optional[float]:
    xs = _as_list(values)
    n = len(xs)
    if n == 0:
        return None
    num = (sum(xs)) ** 2
    den = n * sum(x * x for x in xs)
    return num / den if den > 0 else 1.0

def p90_minus_p10(values: Iterable[float]) -> Optional[float]:
    xs = sorted(_as_list(values))
    n = len(xs)
    if n == 0:
        return None
    def q(p):
        k = (n - 1) * p
        f, c = math.floor(k), math.ceil(k)
        if f == c:
            return xs[int(k)]
        return xs[f] * (c - k) + xs[c] * (k - f)
    return q(0.9) - q(0.1)

def svi_weighted_mean(values: Iterable[float], svi_weights: Iterable[float]) -> Optional[float]:
    xs = _as_list(values)
    ws = _as_list(svi_weights)
    if len(xs) == 0 or len(xs) != len(ws):
        return None
    sw = sum(ws)
    if sw == 0:
        return sum(xs) / len(xs)
    return sum(x * w for x, w in zip(xs, ws)) / sw

def all_metrics(values: Iterable[float], svi_weights: Optional[Iterable[float]] = None) -> Dict[str, Optional[float]]:
    return {
        "gini": gini(values),
        "theil_T": theil_T(values),
        "atkinson_e0.5": atkinson(values, 0.5),
        "jain": jain(values),
        "p90_p10_gap": p90_minus_p10(values),
        "svi_weighted_mean": None if svi_weights is None else svi_weighted_mean(values, svi_weights),
    }
