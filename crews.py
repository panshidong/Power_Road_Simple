from __future__ import annotations

from dataclasses import dataclass, field
from typing import FrozenSet, List, Optional


def asset_required_skill(asset) -> str:
    """
    Asset types:
      - bus: int-like
      - road link: (u, v) tuple
    """
    is_road = isinstance(asset, tuple) and len(asset) == 2
    return "road" if is_road else "power"


@dataclass
class Crew:
    """
    Repair crew (supports both specialized and multifunction via skills).

    skills:
        {"power"}, {"road"}, or {"power","road"}
    available_time:
        time when crew becomes free
    location_node:
        current road node id for dispatch origin (None -> depot)
    speed:
        scales travel time derived from network (1.0 = real time)
    """

    skills: FrozenSet[str] = field(default_factory=lambda: frozenset({"power"}))
    available_time: float = 0.0
    location_node: Optional[int] = None
    speed: float = 1.0

    @classmethod
    def specialized(cls, kind: str, **kwargs) -> "Crew":
        if kind not in {"power", "road"}:
            raise ValueError(f"Invalid crew kind: {kind!r}")
        return cls(skills=frozenset({kind}), **kwargs)

    @classmethod
    def multifunction(cls, **kwargs) -> "Crew":
        return cls(skills=frozenset({"power", "road"}), **kwargs)

    def can_repair(self, asset) -> bool:
        return asset_required_skill(asset) in self.skills


@dataclass
class CrewPool:
    crews: List[Crew] = field(default_factory=list)

    def copy(self) -> "CrewPool":
        return CrewPool(
            [
                Crew(
                    skills=c.skills,
                    available_time=float(c.available_time),
                    location_node=c.location_node,
                    speed=float(c.speed),
                )
                for c in self.crews
            ]
        )

    def next_available(self, asset) -> Crew:
        eligible = [c for c in self.crews if c.can_repair(asset)]
        if not eligible:
            req = asset_required_skill(asset)
            raise ValueError(
                f"No eligible crew for asset {asset!r} (requires skill={req!r}). "
                f"Existing crews skills={[sorted(list(c.skills)) for c in self.crews]!r}"
            )
        return min(eligible, key=lambda c: c.available_time)

    def reset(self, depot_node: Optional[int] = None) -> None:
        for c in self.crews:
            c.available_time = 0.0
            c.location_node = depot_node


def make_crews(
    mode: str = "specialized",
    *,
    power_crews: int = 1,
    road_crews: int = 1,
    multifunction_crews: int = 1,
    depot_node: Optional[int] = None,
    speed: float = 1.0,
) -> CrewPool:
    """
    mode:
      - "specialized": power_crews of {"power"} + road_crews of {"road"}
      - "multifunction": multifunction_crews of {"power","road"}
    """
    if speed <= 0:
        raise ValueError("speed must be > 0")

    crews: List[Crew] = []
    if mode == "specialized":
        if power_crews <= 0 or road_crews <= 0:
            raise ValueError("specialized mode requires power_crews>=1 and road_crews>=1")
        crews += [Crew.specialized("power", location_node=depot_node, speed=speed) for _ in range(int(power_crews))]
        crews += [Crew.specialized("road", location_node=depot_node, speed=speed) for _ in range(int(road_crews))]
    elif mode == "multifunction":
        if multifunction_crews <= 0:
            raise ValueError("multifunction mode requires multifunction_crews>=1")
        crews += [Crew.multifunction(location_node=depot_node, speed=speed) for _ in range(int(multifunction_crews))]
    else:
        raise ValueError(f"Unknown crew mode: {mode!r}")

    return CrewPool(crews=crews)
