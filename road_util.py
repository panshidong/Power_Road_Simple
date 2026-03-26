from __future__ import annotations

import heapq
import os
from typing import Dict, Iterable, List, Mapping, Optional, Tuple

Link = Tuple[int, int]


def _normalized_link_factors(link_factors: Optional[Mapping[Link, float]]) -> Dict[Link, float]:
    out: Dict[Link, float] = {}
    if not link_factors:
        return out
    for (u, v), factor in link_factors.items():
        out[(int(u), int(v))] = float(factor)
    return out


def capacity_adjustment(
    input_file: str,
    output_file: str,
    links: Iterable[Link],
    adj_factor: float,
    link_factors: Optional[Mapping[Link, float]] = None,
) -> None:
    """
    Edit a TAP-B network file by derating (or effectively disabling) specified links.

    - adj_factor < 0.1: set a very large cost-like field (index 4) to approximate removal
    - else: scale capacity field (index 2) by adj_factor

    This is intentionally file-IO based.
    """
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"input_file not found: {input_file!r}")

    links_set = {(int(u), int(v)) for (u, v) in links}
    factor_map = _normalized_link_factors(link_factors)

    with open(input_file, "r", encoding="utf-8") as f:
        lines = f.readlines()

    # keep metadata (assumed first 8 lines; same convention as your earlier code)
    output_lines = lines[:8]
    links_start = 8

    for line in lines[links_start:]:
        s = line.strip()
        if not s or s.startswith("~") or s.startswith(";"):
            continue
        parts = s.split()
        if len(parts) < 10:
            continue

        u = int(parts[0])
        v = int(parts[1])
        key = (u, v)
        key_rev = (v, u)

        if key in links_set or key_rev in links_set:
            factor = float(factor_map.get(key, factor_map.get(key_rev, adj_factor)))
            if factor < 0.1:
                parts[4] = "9999"
            else:
                cap = float(parts[2]) * factor
                parts[2] = f"{cap:.8f}"

        output_lines.append("\t".join(parts) + " \n")

    os.makedirs(os.path.dirname(output_file) or ".", exist_ok=True)
    with open(output_file, "w", encoding="utf-8") as f:
        f.writelines(output_lines)


def eval_tot_OD_travel_time(s_txt_path: str = "s.txt") -> float:
    """
    Read s.txt and compute total system travel time: sum(flow * cost)
    Expected line format:
      (i,j)  flow  cost
    """
    if not os.path.exists(s_txt_path):
        raise FileNotFoundError(f"s.txt not found: {s_txt_path!r}")

    tstt = 0.0
    with open(s_txt_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if not line.startswith("("):
                continue
            parts = line.split()
            if len(parts) != 3:
                continue
            tstt += float(parts[1]) * float(parts[2])

    return float(tstt)


def _read_cost_graph_from_s_txt(file_path: str) -> Dict[Tuple[int, int], float]:
    graph: Dict[Tuple[int, int], float] = {}
    with open(file_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) != 3:
                continue
            if not parts[0].startswith("("):
                continue
            a, b = map(int, parts[0][1:-1].split(","))
            cost = float(parts[2])
            graph[(a, b)] = cost
    return graph


def calculate_shortest_path_cost(file_path: str, start: int, end: int) -> float:
    """
    Dijkstra on directed graph extracted from s.txt.
    Edge weight = cost column.
    """
    if start == end:
        return 0.0
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"file not found: {file_path!r}")

    graph = _read_cost_graph_from_s_txt(file_path)

    pq: List[Tuple[float, int]] = [(0.0, int(start))]
    visited = set()
    dist: Dict[int, float] = {int(start): 0.0}

    while pq:
        d, u = heapq.heappop(pq)
        if u in visited:
            continue
        visited.add(u)
        if u == int(end):
            return float(d)

        for (a, b), w in graph.items():
            if a != u:
                continue
            if b in visited:
                continue
            nd = d + float(w)
            if b not in dist or nd < dist[b]:
                dist[b] = nd
                heapq.heappush(pq, (nd, b))

    return float("inf")
