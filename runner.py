from tradeoff_runner import TradeoffConfig, run_tradeoff_study


def main() -> None:
    info = run_tradeoff_study(TradeoffConfig())
    print("Task B trade-off run complete.")
    print("Result dir:", info["result_dir"])
    print("CSV:", info["csv_path"])
    print("Scatter:", info["scatter_png"])
    print("Summary:", info["summary_md"])
    print("Discussion:", info["discussion_md"])


if __name__ == "__main__":
    main()
