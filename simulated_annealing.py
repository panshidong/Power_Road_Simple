from __future__ import annotations

import math
import random
from dataclasses import dataclass
from typing import Any, Callable, Dict, List, Tuple


@dataclass
class SAConfig:
    seed: int = 0
    max_iter: int = 200
    T0: float = 1.0
    alpha: float = 0.98          # temperature *= alpha
    neighbor: str = "swap"       # "swap" or "insert"
    keep_best: bool = True
    accept_equal: bool = True    # if True, allow accept when delta==0


def _neighbor_swap(seq: List[Any], rng: random.Random) -> List[Any]:
    n = len(seq)
    i, j = rng.sample(range(n), 2)
    cand = seq[:]
    cand[i], cand[j] = cand[j], cand[i]
    return cand


def _neighbor_insert(seq: List[Any], rng: random.Random) -> List[Any]:
    n = len(seq)
    i, j = rng.sample(range(n), 2)
    cand = seq[:]
    x = cand.pop(i)
    cand.insert(j, x)
    return cand


def simulated_annealing(
    *,
    initial: List[Any],
    evaluate: Callable[[List[Any], str], Dict[str, Any]],
    # evaluate(seq, run_tag) -> dict containing at least:
    #   dict["objective_value"] (float)
    config: SAConfig,
    scenario_prefix: str = "sa",
) -> Dict[str, Any]:
    """
    Generic SA over permutations.
    - evaluate is user-supplied and can decide whether to write files (recommended).
    - scenario_prefix is passed to evaluate to build per-iteration run tags.
    """
    rng = random.Random(config.seed)
    curr = list(initial)

    curr_run = evaluate(curr, f"{scenario_prefix}_init")
    curr_val = float(curr_run["objective_value"])

    best = curr[:]
    best_run = curr_run
    best_val = curr_val

    T = float(config.T0)

    for it in range(int(config.max_iter)):
        if config.neighbor == "swap":
            cand = _neighbor_swap(curr, rng)
        elif config.neighbor == "insert":
            cand = _neighbor_insert(curr, rng)
        else:
            raise ValueError(f"Unknown neighbor={config.neighbor!r}")

        cand_run = evaluate(cand, f"{scenario_prefix}_it{it:04d}")
        cand_val = float(cand_run["objective_value"])

        delta = cand_val - curr_val

        accept = False
        if delta < 0:
            accept = True
        elif delta == 0 and config.accept_equal:
            accept = True
        else:
            # minimize objective
            prob = math.exp(-delta / max(T, 1e-12))
            if rng.random() < prob:
                accept = True

        if accept:
            curr, curr_run, curr_val = cand, cand_run, cand_val

        if config.keep_best and cand_val < best_val:
            best, best_run, best_val = cand, cand_run, cand_val

        T *= float(config.alpha)
        if T < 1e-12:
            T = 1e-12

    return {
        "best_sequence": best,
        "best_objective_value": best_val,
        "best_run": best_run,
        "final_sequence": curr,
        "final_objective_value": curr_val,
    }
