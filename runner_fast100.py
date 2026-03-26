from batch_tradeoff_runner import BatchTradeoffConfig, run_batch_tradeoff_study
from simulated_annealing import SAConfig
from tradeoff_runner import TradeoffConfig


def build_fast100_config() -> BatchTradeoffConfig:
    """
    Fast 100-disaster preset.

    Design goal:
    - keep 100 random disasters
    - keep specialized crews and the current CRI weights
    - reduce wall-clock time by shrinking the experiment set and SA budget
    """
    study_cfg = TradeoffConfig(
        crew_mode="specialized",
        power_crews=1,
        road_crews=1,
        multifunction_crews=1,
        cri_weight_pairs=[(0.133, 0.867)],
        weighted_lambdas=[0.5],
        enabled_experiment_ids=[
            "single_triangle",
            "single_gini_restore",
            "single_maximin_time_avg_cri",
            "single_p90_access_restore",
            "weighted_gini_restore_l050",
            "weighted_maximin_time_avg_cri_l050",
            "guardrail_gini_restore",
        ],
        run_sensitivity=False,
        save_reference_artifacts=False,
        save_baseline=False,
        save_best_artifacts=False,
        save_best_debug=False,
        sa=SAConfig(seed=0, max_iter=80, T0=1.4, alpha=0.97, neighbor="swap"),
    )
    return BatchTradeoffConfig(
        n_scenarios=100,
        seed0=20260325,
        study_cfg=study_cfg,
    )


def main() -> None:
    info = run_batch_tradeoff_study(build_fast100_config())
    print("Fast 100-disaster batch run complete.")
    print("Result dir:", info["result_dir"])
    print("Scenario manifest:", info["scenario_manifest_csv"])
    print("Scenario rows:", info["scenario_rows_csv"])
    print("Aggregate CSV:", info["aggregate_csv"])
    print("Error-bar plot:", info["errorbar_png"])
    print("Summary:", info["summary_md"])


if __name__ == "__main__":
    main()
