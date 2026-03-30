# Trade-off Summary

- Result directory: `results/batch_tradeoff_20260325_235335/scenario_runs/tradeoff_20260326_032559`
- Critical locations: `[1, 24]`
- Primary CRI weights: `(w_e=0.133, w_a=0.867)`
- Primary CRI threshold: `0.9`
- Primary critical-access threshold: `0.9`
- Crew mode: `specialized` (power=`1`, road=`1`)
- Bus dispatch mode: `link_only`
- Bus-to-link source: `new_bus_to_link.json`
- Enabled experiments: `['single_triangle', 'single_gini_restore', 'single_maximin_time_avg_cri', 'single_p90_access_restore', 'weighted_gini_restore_l050', 'weighted_maximin_time_avg_cri_l050', 'guardrail_gini_restore']`
- Pareto experiments (triangle vs gini_restore): `['single_triangle', 'single_p90_access_restore']`
- Combined CSV: `tradeoff_summary.csv`
- Scatter plot: `tradeoff_scatter.png`
- The combined CSV includes both `power_sequence` and `road_sequence` for specialized-crew runs.

## Optimized Results

| Experiment | Type | Description | Triangle | Gini | Min time-avg CRI | P90 access | Time-avg share access | Pareto |
|---|---|---|---:|---:|---:|---:|---:|---|
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 1876.068 | 0.131 | 0.346 | 1429.912 | 0.326 | No |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 869.400 | 0.267 | 0.433 | 1028.376 | 0.355 | Yes |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 2097.666 | 0.124 | 0.260 | 1852.253 | 0.212 | No |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 1532.436 | 0.283 | 0.560 | 1862.133 | 0.460 | No |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 1295.678 | 0.107 | 0.428 | 1063.890 | 0.376 | Yes |
| weighted_gini_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:gini_restore. | 1409.735 | 0.142 | 0.257 | 1264.441 | 0.245 | No |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 1583.838 | 0.123 | 0.438 | 1279.762 | 0.347 | No |
| weighted_maximin_time_avg_cri_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:maximin_time_avg_cri_loss. | 950.241 | 0.307 | 0.580 | 1015.745 | 0.492 | No |
