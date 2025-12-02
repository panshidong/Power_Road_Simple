
from __future__ import annotations
from dataclasses import dataclass, field
from typing import List, Optional

@dataclass
class Crew:
    """Specialized repair crew.
    kind: "power" or "road" (strict — a power crew cannot repair road, and vice versa)
    available_time: when this crew becomes free
    location_node: current road node id used for dispatch; None means "use depot"
    speed: scales travel time from the network (1.0 = real time)
    """
    kind: str  # "power" or "road"
    available_time: float = 0.0
    location_node: Optional[int] = None
    speed: float = 1.0

    def can_repair(self, asset) -> bool:
        # Buses are ints; road links are 2-tuples like (u, v)
        is_road = isinstance(asset, tuple) and len(asset) == 2
        return (self.kind == "road" and is_road) or (self.kind == "power" and not is_road)

@dataclass
class CrewPool:
    crews: List[Crew] = field(default_factory=list)

    def copy(self) -> "CrewPool":
        return CrewPool([Crew(c.kind, c.available_time, c.location_node, c.speed) for c in self.crews])

    def next_available(self, asset) -> Crew:
        eligible = [c for c in self.crews if c.can_repair(asset)]
        if not eligible:
            raise ValueError(f"No eligible crew for asset {asset!r}. Ensure specialized crews exist for both types.")
        # Choose earliest available; tie-break by order
        return min(eligible, key=lambda c: c.available_time)

    def reset(self, depot_node: Optional[int] = None):
        for c in self.crews:
            c.available_time = 0.0
            c.location_node = depot_node
