from __future__ import annotations

import json
import random
from pathlib import Path
from typing import List, Tuple, Dict, Any


def load_candidates(
    *,
    bus_location_path: str = "bus_location.json",
    base_net_path: str = "tap-b/net/SiouxFalls_net.txt",
) -> Dict[str, Any]:
    # buses from bus_location keys
    bus_loc = json.loads(Path(bus_location_path).read_text(encoding="utf-8"))
    buses = sorted([int(k) for k in bus_loc.keys()])

    # links from TAP net file (u,v)
    links = set()
    with open(base_net_path, "r", encoding="utf-8") as f:
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
            links.add((u, v))
    links = sorted(list(links))

    return {"buses": buses, "links": links}


def sample_disaster(
    *,
    n_buses: int,
    n_links: int,
    seed: int,
    candidates: Dict[str, Any],
) -> Dict[str, Any]:
    rng = random.Random(seed)
    buses = candidates["buses"]
    links = candidates["links"]
    if n_buses > len(buses) or n_links > len(links):
        raise ValueError("Requested failures exceed candidate pool.")
    return {
        "seed": seed,
        "broken_buses": rng.sample(buses, n_buses),
        "broken_links": rng.sample(links, n_links),
    }


def generate_scenarios(
    *,
    n_scenarios: int,
    n_buses: int,
    n_links: int,
    seed0: int = 0,
    out_json: str = "random_disasters.json",
) -> str:
    cand = load_candidates()
    scenarios = []
    for i in range(n_scenarios):
        scenarios.append(sample_disaster(n_buses=n_buses, n_links=n_links, seed=seed0 + i, candidates=cand))
    Path(out_json).write_text(json.dumps({"n": n_scenarios, "n_buses": n_buses, "n_links": n_links, "scenarios": scenarios}, indent=2), encoding="utf-8")
    return out_json
