from __future__ import annotations

import task_a_stability


RESULT_DIR = "results/task_a_shapley_stability_powerfix_10"


def main() -> None:
    """Run the nested Shapley stability check under the fixed power metric."""
    task_a_stability.RESULT_DIR = RESULT_DIR
    task_a_stability.main()


if __name__ == "__main__":
    main()
