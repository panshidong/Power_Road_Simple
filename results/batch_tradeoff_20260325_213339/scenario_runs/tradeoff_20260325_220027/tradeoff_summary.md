# Trade-off Summary

- Result directory: `results/batch_tradeoff_20260325_213339/scenario_runs/tradeoff_20260325_220027`
- Critical locations: `[1, 24]`
- Primary CRI weights: `(w_e=0.133, w_a=0.867)`
- Primary CRI threshold: `0.9`
- Primary critical-access threshold: `0.9`
- Crew mode: `specialized` (power=`1`, road=`1`)
- Bus dispatch mode: `link_only`
- Bus-to-link source: `new_bus_to_link.json`
- Enabled experiments: `['single_triangle', 'single_gini_restore', 'single_maximin_time_avg_cri', 'single_p90_access_restore', 'weighted_gini_restore_l050', 'weighted_maximin_time_avg_cri_l050', 'guardrail_gini_restore']`
- Pareto experiments (triangle vs gini_restore): `['single_gini_restore', 'guardrail_gini_restore']`
- Combined CSV: `tradeoff_summary.csv`
- Scatter plot: `tradeoff_scatter.png`
- The combined CSV includes both `power_sequence` and `road_sequence` for specialized-crew runs.

## Optimized Results

| Experiment | Type | Description | Triangle | Gini | Min time-avg CRI | P90 access | Time-avg share access | Pareto |
|---|---|---|---:|---:|---:|---:|---:|---|
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 1773.020 | 0.207 | 0.266 | 1820.625 | 0.293 | No |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 1525.589 | 0.203 | 0.267 | 1564.898 | 0.289 | No |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 1599.998 | 0.170 | 0.239 | 1562.454 | 0.293 | Yes |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 1436.318 | 0.234 | 0.319 | 1552.480 | 0.311 | No |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 1329.024 | 0.206 | 0.252 | 1314.402 | 0.298 | No |
| weighted_gini_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:gini_restore. | 1421.668 | 0.183 | 0.271 | 1450.238 | 0.292 | No |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 1298.084 | 0.177 | 0.276 | 1314.402 | 0.293 | Yes |
| weighted_maximin_time_avg_cri_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:maximin_time_avg_cri_loss. | 1480.311 | 0.217 | 0.291 | 1551.307 | 0.293 | No |
