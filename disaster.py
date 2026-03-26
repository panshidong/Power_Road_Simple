from __future__ import annotations

import json
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Tuple

Link = Tuple[int, int]


@dataclass(frozen=True)
class DisasterScenario:
    scenario_id: str
    seed: int
    broken_buses: List[int]
    broken_links: List[Link]
    link_capacity_factors: Dict[Link, float]

    def to_json_dict(self) -> Dict[str, Any]:
        return {
            "scenario_id": self.scenario_id,
            "seed": self.seed,
            "broken_buses": list(self.broken_buses),
            "broken_links": [
                {
                    "u": int(u),
                    "v": int(v),
                    "remaining_capacity_factor": float(self.link_capacity_factors[(u, v)]),
                    "capacity_drop_fraction": float(1.0 - self.link_capacity_factors[(u, v)]),
                }
                for (u, v) in self.broken_links
            ],
        }


def load_candidates(
    *,
    bus_candidate_path: str = "new_bus_to_link.json",
    base_net_path: str = "tap-b/net/SiouxFalls_net.txt",
) -> Dict[str, Any]:
    bus_source = json.loads(Path(bus_candidate_path).read_text(encoding="utf-8"))
    buses = sorted(int(k) for k in bus_source.keys())

    links = set()
    with open(base_net_path, "r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("~") or s.startswith(";"):
                continue
            parts = s.split()
            if len(parts) < 10:
                continue
            try:
                u = int(parts[0])
                v = int(parts[1])
            except Exception:
                continue
            links.add((u, v))

    return {"buses": buses, "links": sorted(links)}


def sample_disaster(
    *,
    scenario_id: str,
    seed: int,
    bus_count_range: Tuple[int, int],
    link_count_range: Tuple[int, int],
    link_drop_range: Tuple[float, float],
    candidates: Dict[str, Any],
) -> DisasterScenario:
    rng = random.Random(seed)
    buses = list(candidates["buses"])
    links = list(candidates["links"])

    n_buses = rng.randint(int(bus_count_range[0]), int(bus_count_range[1]))
    n_links = rng.randint(int(link_count_range[0]), int(link_count_range[1]))

    if n_buses > len(buses) or n_links > len(links):
        raise ValueError("Requested failures exceed candidate pool.")

    broken_buses = rng.sample(buses, n_buses)
    broken_links = rng.sample(links, n_links)

    low_drop = float(link_drop_range[0])
    high_drop = float(link_drop_range[1])
    if not (0.0 <= low_drop <= high_drop <= 1.0):
        raise ValueError("link_drop_range must be within [0, 1] and low<=high")

    link_capacity_factors: Dict[Link, float] = {}
    for link in broken_links:
        drop = rng.uniform(low_drop, high_drop)
        link_capacity_factors[link] = max(0.0, 1.0 - float(drop))

    return DisasterScenario(
        scenario_id=scenario_id,
        seed=seed,
        broken_buses=broken_buses,
        broken_links=broken_links,
        link_capacity_factors=link_capacity_factors,
    )


def generate_scenarios(
    *,
    n_scenarios: int,
    seed0: int = 0,
    bus_count_range: Tuple[int, int] = (8, 15),
    link_count_range: Tuple[int, int] = (3, 13),
    link_drop_range: Tuple[float, float] = (0.5, 1.0),
    bus_candidate_path: str = "new_bus_to_link.json",
    base_net_path: str = "tap-b/net/SiouxFalls_net.txt",
    out_json: str = "random_disasters.json",
) -> List[DisasterScenario]:
    candidates = load_candidates(bus_candidate_path=bus_candidate_path, base_net_path=base_net_path)
    scenarios: List[DisasterScenario] = []
    for idx in range(int(n_scenarios)):
        scenario_id = f"scenario_{idx + 1:03d}"
        scenarios.append(
            sample_disaster(
                scenario_id=scenario_id,
                seed=int(seed0 + idx),
                bus_count_range=bus_count_range,
                link_count_range=link_count_range,
                link_drop_range=link_drop_range,
                candidates=candidates,
            )
        )

    payload = {
        "n_scenarios": int(n_scenarios),
        "seed0": int(seed0),
        "bus_count_range": list(bus_count_range),
        "link_count_range": list(link_count_range),
        "link_drop_range": list(link_drop_range),
        "scenarios": [scenario.to_json_dict() for scenario in scenarios],
    }
    Path(out_json).write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return scenarios


def scenario_sequence(scenario: DisasterScenario) -> List[Any]:
    seq: List[Any] = list(scenario.broken_buses)
    seq.extend(list(scenario.broken_links))
    return seq
