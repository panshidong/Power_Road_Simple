# Trade-off Summary

- Result directory: `results/batch_tradeoff_20260325_213339/scenario_runs/tradeoff_20260325_213339`
- Critical locations: `[1, 24]`
- Primary CRI weights: `(w_e=0.133, w_a=0.867)`
- Primary CRI threshold: `0.9`
- Primary critical-access threshold: `0.9`
- Crew mode: `specialized` (power=`1`, road=`1`)
- Bus dispatch mode: `link_only`
- Bus-to-link source: `new_bus_to_link.json`
- Enabled experiments: `['single_triangle', 'single_gini_restore', 'single_maximin_time_avg_cri', 'single_p90_access_restore', 'weighted_gini_restore_l050', 'weighted_maximin_time_avg_cri_l050', 'guardrail_gini_restore']`
- Pareto experiments (triangle vs gini_restore): `['single_gini_restore', 'single_p90_access_restore', 'weighted_maximin_time_avg_cri_l050']`
- Combined CSV: `tradeoff_summary.csv`
- Scatter plot: `tradeoff_scatter.png`
- The combined CSV includes both `power_sequence` and `road_sequence` for specialized-crew runs.

## Optimized Results

| Experiment | Type | Description | Triangle | Gini | Min time-avg CRI | P90 access | Time-avg share access | Pareto |
|---|---|---|---:|---:|---:|---:|---:|---|
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 4949.903 | 0.195 | 0.302 | 5338.734 | 0.313 | No |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 4747.240 | 0.210 | 0.314 | 4336.497 | 0.384 | No |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 4951.666 | 0.148 | 0.358 | 4804.604 | 0.360 | Yes |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 5224.817 | 0.245 | 0.377 | 5648.027 | 0.362 | No |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 4837.246 | 0.152 | 0.337 | 4621.908 | 0.356 | Yes |
| weighted_gini_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:gini_restore. | 5019.855 | 0.160 | 0.307 | 4990.501 | 0.315 | No |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 4935.458 | 0.162 | 0.310 | 4941.414 | 0.315 | No |
| weighted_maximin_time_avg_cri_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:maximin_time_avg_cri_loss. | 3856.620 | 0.187 | 0.274 | 3641.838 | 0.317 | Yes |
