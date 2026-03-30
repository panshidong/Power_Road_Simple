# Trade-off Summary

- Result directory: `results/batch_tradeoff_20260325_235335/scenario_runs/tradeoff_20260326_030530`
- Critical locations: `[1, 24]`
- Primary CRI weights: `(w_e=0.133, w_a=0.867)`
- Primary CRI threshold: `0.9`
- Primary critical-access threshold: `0.9`
- Crew mode: `specialized` (power=`1`, road=`1`)
- Bus dispatch mode: `link_only`
- Bus-to-link source: `new_bus_to_link.json`
- Enabled experiments: `['single_triangle', 'single_gini_restore', 'single_maximin_time_avg_cri', 'single_p90_access_restore', 'weighted_gini_restore_l050', 'weighted_maximin_time_avg_cri_l050', 'guardrail_gini_restore']`
- Pareto experiments (triangle vs gini_restore): `['single_triangle', 'weighted_gini_restore_l050', 'guardrail_gini_restore', 'weighted_maximin_time_avg_cri_l050']`
- Combined CSV: `tradeoff_summary.csv`
- Scatter plot: `tradeoff_scatter.png`
- The combined CSV includes both `power_sequence` and `road_sequence` for specialized-crew runs.

## Optimized Results

| Experiment | Type | Description | Triangle | Gini | Min time-avg CRI | P90 access | Time-avg share access | Pareto |
|---|---|---|---:|---:|---:|---:|---:|---|
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 2044.636 | 0.119 | 0.178 | 1906.293 | 0.188 | No |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 819.288 | 0.128 | 0.512 | 784.871 | 0.491 | Yes |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 1813.834 | 0.074 | 0.374 | 1839.995 | 0.295 | No |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 1415.025 | 0.221 | 0.486 | 1854.486 | 0.289 | No |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 946.387 | 0.151 | 0.266 | 841.667 | 0.246 | No |
| weighted_gini_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:gini_restore. | 891.654 | 0.081 | 0.496 | 894.951 | 0.397 | Yes |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 856.286 | 0.091 | 0.543 | 842.132 | 0.492 | Yes |
| weighted_maximin_time_avg_cri_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:maximin_time_avg_cri_loss. | 1525.300 | 0.068 | 0.385 | 1536.764 | 0.368 | Yes |
