# Trade-off Summary

- Result directory: `results/batch_tradeoff_20260325_213339/scenario_runs/tradeoff_20260325_220820`
- Critical locations: `[1, 24]`
- Primary CRI weights: `(w_e=0.133, w_a=0.867)`
- Primary CRI threshold: `0.9`
- Primary critical-access threshold: `0.9`
- Crew mode: `specialized` (power=`1`, road=`1`)
- Bus dispatch mode: `link_only`
- Bus-to-link source: `new_bus_to_link.json`
- Enabled experiments: `['single_triangle', 'single_gini_restore', 'single_maximin_time_avg_cri', 'single_p90_access_restore', 'weighted_gini_restore_l050', 'weighted_maximin_time_avg_cri_l050', 'guardrail_gini_restore']`
- Pareto experiments (triangle vs gini_restore): `['baseline_reference', 'single_triangle', 'single_gini_restore', 'single_p90_access_restore', 'weighted_gini_restore_l050', 'guardrail_gini_restore']`
- Combined CSV: `tradeoff_summary.csv`
- Scatter plot: `tradeoff_scatter.png`
- The combined CSV includes both `power_sequence` and `road_sequence` for specialized-crew runs.

## Optimized Results

| Experiment | Type | Description | Triangle | Gini | Min time-avg CRI | P90 access | Time-avg share access | Pareto |
|---|---|---|---:|---:|---:|---:|---:|---|
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 1407.414 | 0.128 | 0.243 | 1386.923 | 0.219 | Yes |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 852.241 | 0.184 | 0.190 | 997.982 | 0.255 | Yes |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 1407.414 | 0.128 | 0.243 | 1386.923 | 0.219 | Yes |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 2084.476 | 0.139 | 0.247 | 2092.293 | 0.248 | No |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 899.863 | 0.137 | 0.290 | 805.040 | 0.341 | Yes |
| weighted_gini_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:gini_restore. | 870.471 | 0.138 | 0.424 | 805.040 | 0.454 | Yes |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 1407.414 | 0.128 | 0.243 | 1386.923 | 0.219 | Yes |
| weighted_maximin_time_avg_cri_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:maximin_time_avg_cri_loss. | 892.692 | 0.199 | 0.454 | 990.467 | 0.475 | No |
