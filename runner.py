from batch_tradeoff_runner import BatchTradeoffConfig, run_batch_tradeoff_study


def main() -> None:
    info = run_batch_tradeoff_study(BatchTradeoffConfig())
    print("Batch trade-off run complete.")
    print("Result dir:", info["result_dir"])
    print("Scenario manifest:", info["scenario_manifest_csv"])
    print("Scenario rows:", info["scenario_rows_csv"])
    print("Aggregate CSV:", info["aggregate_csv"])
    print("Error-bar plot:", info["errorbar_png"])
    print("Summary:", info["summary_md"])


if __name__ == "__main__":
    main()
