from resilience_measurement import run_model_multi, optimize_sequence

def main():
    seq = [11, (1, 3), 15]

    # 1) 只评估（带中文/英文解释报告 + triangle 图 + debug log）
    out = run_model_multi(
        sequence=seq,
        result_folder="results",
        message="test multifunction crew strict",
        Scenario="eval_multifunc",
        plot_control=False,
        focus=False,
        crew_mode="multifunction",
        multifunction_crews=1,
        strict=True,
        debug=True,
        objective="triangle",   # 或 "equity:gini"
    )
    print("Done evaluation:", out["paths"], out["triangle_png"])

    # 2) 小规模优化（3 个资产建议用 bruteforce）
    best = optimize_sequence(
        base_sequence=seq,
        result_folder="results",
        message="optimize demo",
        Scenario="opt_demo",
        objective="triangle",   # 或 "equity:gini"
        method="bruteforce",
        strict=True,
        crew_mode="multifunction",
        multifunction_crews=1,
    )
    print("Best seq:", best["best_sequence"])
    print("Best objective:", best["best_objective_value"])
    print("Best reports:", best["paths"], best["triangle_png"])

if __name__ == "__main__":
    main()