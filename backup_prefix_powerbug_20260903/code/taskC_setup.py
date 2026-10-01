from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Optional

DEFAULT_OD_CONFIG: Dict[str, Any] = {
    "depot_nodes": [1],
    "k_paths": 3,
    "path_rank_weights": [1.0, 0.6, 0.35],
    "destination_endpoint_mode": "both",
    "allocation_method": "weighted_overlap",
    "bus_to_link_path": "new_bus_to_link.json",
    "shapley_samples": 60,
    "shapley_seed": 20260701,
}


def write_json(path: Path, obj: object) -> None:
    path.write_text(json.dumps(obj, indent=2), encoding="utf-8")


def generate_taskc_inputs(overrides: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """Write taskC_od_config.json with the placeholder O-D assumptions.

    See TASKC_ASSUMPTIONS.md for what each field means and why these particular
    placeholder values were chosen.
    """
    cfg = dict(DEFAULT_OD_CONFIG)
    if overrides:
        cfg.update(overrides)
    write_json(Path("taskC_od_config.json"), cfg)
    return cfg


def main() -> None:
    cfg = generate_taskc_inputs()
    print("Wrote taskC_od_config.json:")
    print(json.dumps(cfg, indent=2))


if __name__ == "__main__":
    main()
