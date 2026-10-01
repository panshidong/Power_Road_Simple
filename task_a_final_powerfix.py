from __future__ import annotations

import time

from simulated_annealing import SAConfig
from task_a_criticality import TaskACriticalityConfig
from task_a_runner import TaskAConfig, run_task_a_study


RESULT_DIR = "results/task_a_criticality_powerfix_100"


def main() -> None:
    """Re-run the production Task A design after the downstream-power fix."""
    cfg = TaskAConfig(
        n_scenarios=100,
        seed0=20260402,
        result_dir=RESULT_DIR,
        resume=True,
        bus_count_range=(8, 15),
        link_count_range=(3, 13),
        link_drop_range=(0.5, 1.0),
        include_heuristic=True,
        strict=True,
        crew_mode="specialized",
        power_crews=1,
        road_crews=1,
        criticality=TaskACriticalityConfig(
            shapley_samples=120,
            shapley_seed=20260402,
            integrated_power_weight=0.5,
            use_full_functionality_shapley=True,
        ),
        heuristic_sa=SAConfig(
            seed=20260402,
            max_iter=60,
            T0=2.0,
            alpha=0.98,
            neighbor="swap",
        ),
    )
    info = None
    for attempt in range(1, 21):
        try:
            info = run_task_a_study(cfg)
            break
        except FileNotFoundError as exc:
            if "s.txt" not in str(exc) or attempt == 20:
                raise
            print(
                f"[task-a-powerfix] transient s.txt failure; "
                f"resuming from checkpoint (retry {attempt}/20)",
                flush=True,
            )
            time.sleep(1.0)
    if info is None:
        raise RuntimeError("Task A power-fix study did not complete")
    print("Task A power-fix study complete.")
    for key, value in info.items():
        print(f"{key}: {value}")


if __name__ == "__main__":
    main()
