# Trade-off Summary

- Result directory: `results/batch_tradeoff_20260325_213339/scenario_runs/tradeoff_20260325_214405`
- Critical locations: `[1, 24]`
- Primary CRI weights: `(w_e=0.133, w_a=0.867)`
- Primary CRI threshold: `0.9`
- Primary critical-access threshold: `0.9`
- Crew mode: `specialized` (power=`1`, road=`1`)
- Bus dispatch mode: `link_only`
- Bus-to-link source: `new_bus_to_link.json`
- Enabled experiments: `['single_triangle', 'single_gini_restore', 'single_maximin_time_avg_cri', 'single_p90_access_restore', 'weighted_gini_restore_l050', 'weighted_maximin_time_avg_cri_l050', 'guardrail_gini_restore']`
- Pareto experiments (triangle vs gini_restore): `['single_triangle', 'single_p90_access_restore', 'weighted_gini_restore_l050', 'guardrail_gini_restore', 'weighted_maximin_time_avg_cri_l050']`
- Combined CSV: `tradeoff_summary.csv`
- Scatter plot: `tradeoff_scatter.png`
- The combined CSV includes both `power_sequence` and `road_sequence` for specialized-crew runs.

## Optimized Results

| Experiment | Type | Description | Triangle | Gini | Min time-avg CRI | P90 access | Time-avg share access | Pareto |
|---|---|---|---:|---:|---:|---:|---:|---|
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 948.829 | 0.296 | 0.519 | 1412.523 | 0.353 | No |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 701.816 | 0.283 | 0.443 | 1037.306 | 0.348 | Yes |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 982.661 | 0.282 | 0.258 | 1075.001 | 0.285 | No |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 1327.948 | 0.331 | 0.541 | 1815.893 | 0.422 | No |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 697.655 | 0.284 | 0.318 | 898.368 | 0.335 | Yes |
| weighted_gini_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:gini_restore. | 944.834 | 0.259 | 0.279 | 1075.001 | 0.286 | Yes |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 974.145 | 0.254 | 0.361 | 1083.093 | 0.285 | Yes |
| weighted_maximin_time_avg_cri_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:maximin_time_avg_cri_loss. | 666.391 | 0.302 | 0.423 | 1032.002 | 0.363 | Yes |
