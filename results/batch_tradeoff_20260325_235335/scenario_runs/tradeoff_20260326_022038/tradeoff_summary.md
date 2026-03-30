# Trade-off Summary

- Result directory: `results/batch_tradeoff_20260325_235335/scenario_runs/tradeoff_20260326_022038`
- Critical locations: `[1, 24]`
- Primary CRI weights: `(w_e=0.133, w_a=0.867)`
- Primary CRI threshold: `0.9`
- Primary critical-access threshold: `0.9`
- Crew mode: `specialized` (power=`1`, road=`1`)
- Bus dispatch mode: `link_only`
- Bus-to-link source: `new_bus_to_link.json`
- Enabled experiments: `['single_triangle', 'single_gini_restore', 'single_maximin_time_avg_cri', 'single_p90_access_restore', 'weighted_gini_restore_l050', 'weighted_maximin_time_avg_cri_l050', 'guardrail_gini_restore']`
- Pareto experiments (triangle vs gini_restore): `['single_triangle', 'single_p90_access_restore', 'guardrail_gini_restore']`
- Combined CSV: `tradeoff_summary.csv`
- Scatter plot: `tradeoff_scatter.png`
- The combined CSV includes both `power_sequence` and `road_sequence` for specialized-crew runs.

## Optimized Results

| Experiment | Type | Description | Triangle | Gini | Min time-avg CRI | P90 access | Time-avg share access | Pareto |
|---|---|---|---:|---:|---:|---:|---:|---|
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 5814.199 | 0.114 | 0.485 | 5346.846 | 0.427 | No |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 2526.492 | 0.265 | 0.523 | 3352.034 | 0.466 | Yes |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 5814.199 | 0.114 | 0.485 | 5346.846 | 0.427 | No |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 5170.920 | 0.240 | 0.555 | 6638.414 | 0.509 | No |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 3065.334 | 0.157 | 0.305 | 2728.835 | 0.290 | Yes |
| weighted_gini_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:gini_restore. | 5088.274 | 0.109 | 0.300 | 4943.567 | 0.259 | No |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 3914.542 | 0.076 | 0.278 | 3483.042 | 0.262 | Yes |
| weighted_maximin_time_avg_cri_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:maximin_time_avg_cri_loss. | 5195.996 | 0.163 | 0.543 | 5253.946 | 0.473 | No |
