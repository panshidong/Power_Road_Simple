# Trade-off Summary

- Result directory: `results/batch_tradeoff_20260325_213339/scenario_runs/tradeoff_20260325_230611`
- Critical locations: `[1, 24]`
- Primary CRI weights: `(w_e=0.133, w_a=0.867)`
- Primary CRI threshold: `0.9`
- Primary critical-access threshold: `0.9`
- Crew mode: `specialized` (power=`1`, road=`1`)
- Bus dispatch mode: `link_only`
- Bus-to-link source: `new_bus_to_link.json`
- Enabled experiments: `['single_triangle', 'single_gini_restore', 'single_maximin_time_avg_cri', 'single_p90_access_restore', 'weighted_gini_restore_l050', 'weighted_maximin_time_avg_cri_l050', 'guardrail_gini_restore']`
- Pareto experiments (triangle vs gini_restore): `['guardrail_gini_restore', 'weighted_maximin_time_avg_cri_l050']`
- Combined CSV: `tradeoff_summary.csv`
- Scatter plot: `tradeoff_scatter.png`
- The combined CSV includes both `power_sequence` and `road_sequence` for specialized-crew runs.

## Optimized Results

| Experiment | Type | Description | Triangle | Gini | Min time-avg CRI | P90 access | Time-avg share access | Pareto |
|---|---|---|---:|---:|---:|---:|---:|---|
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 29373.055 | 0.202 | 0.476 | 33685.718 | 0.360 | No |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 26852.222 | 0.168 | 0.491 | 23028.731 | 0.454 | No |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 31948.962 | 0.172 | 0.360 | 33619.240 | 0.264 | No |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 26880.386 | 0.182 | 0.492 | 22950.981 | 0.454 | No |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 27956.272 | 0.173 | 0.484 | 22702.233 | 0.453 | No |
| weighted_gini_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:gini_restore. | 27949.982 | 0.190 | 0.483 | 22702.233 | 0.453 | No |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 26738.429 | 0.164 | 0.492 | 22756.900 | 0.453 | Yes |
| weighted_maximin_time_avg_cri_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:maximin_time_avg_cri_loss. | 25327.035 | 0.207 | 0.490 | 22974.243 | 0.469 | Yes |
