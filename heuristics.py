
from __future__ import annotations
import math
import random
from typing import Callable, List, Tuple, Any

Sequence = List[Any]

def swap_or_relocate(seq: Sequence) -> Sequence:
    """Neighbor generator: with 50% prob do a swap, else relocate one item."""
    n = len(seq)
    if n < 2:
        return seq[:]
    s = seq[:]
    if random.random() < 0.5:
        i, j = random.sample(range(n), 2)
        s[i], s[j] = s[j], s[i]
    else:
        i, j = random.sample(range(n), 2)
        item = s.pop(i)
        j = j if j < len(s) else len(s)
        s.insert(j, item)
    return s

def simulated_annealing(initial: Sequence,
                        eval_fn: Callable[[Sequence], float],
                        neighbor_fn: Callable[[Sequence], Sequence] = swap_or_relocate,
                        T0: float = 1.0,
                        alpha: float = 0.95,
                        iters: int = 2000,
                        seed: int = 42) -> Tuple[Sequence, float]:
    """Simple, dependency-free SA. Returns best sequence and its objective (lower is better)."""
    rng = random.Random(seed)
    curr = initial[:]
    curr_cost = eval_fn(curr)
    best, best_cost = curr[:], curr_cost
    T = T0 if T0 > 1e-9 else 1.0

    for k in range(1, iters + 1):
        cand = neighbor_fn(curr)
        cand_cost = eval_fn(cand)
        d = cand_cost - curr_cost
        if d <= 0 or rng.random() < math.exp(-d / max(T, 1e-12)):
            curr, curr_cost = cand, cand_cost
            if curr_cost < best_cost:
                best, best_cost = curr[:], curr_cost
        T *= alpha
    return best, best_cost
