# Trade-off Summary

- Result directory: `results/batch_tradeoff_20260325_235335/scenario_runs/tradeoff_20260326_050918`
- Critical locations: `[1, 24]`
- Primary CRI weights: `(w_e=0.133, w_a=0.867)`
- Primary CRI threshold: `0.9`
- Primary critical-access threshold: `0.9`
- Crew mode: `specialized` (power=`1`, road=`1`)
- Bus dispatch mode: `link_only`
- Bus-to-link source: `new_bus_to_link.json`
- Enabled experiments: `['single_triangle', 'single_gini_restore', 'single_maximin_time_avg_cri', 'single_p90_access_restore', 'weighted_gini_restore_l050', 'weighted_maximin_time_avg_cri_l050', 'guardrail_gini_restore']`
- Pareto experiments (triangle vs gini_restore): `['single_triangle', 'single_gini_restore', 'weighted_gini_restore_l050', 'guardrail_gini_restore']`
- Combined CSV: `tradeoff_summary.csv`
- Scatter plot: `tradeoff_scatter.png`
- The combined CSV includes both `power_sequence` and `road_sequence` for specialized-crew runs.

## Optimized Results

| Experiment | Type | Description | Triangle | Gini | Min time-avg CRI | P90 access | Time-avg share access | Pareto |
|---|---|---|---:|---:|---:|---:|---:|---|
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 2581.115 | 0.148 | 0.271 | 2454.541 | 0.291 | No |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 1370.658 | 0.159 | 0.294 | 1386.868 | 0.310 | Yes |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 2157.096 | 0.076 | 0.249 | 1906.154 | 0.270 | Yes |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 2385.738 | 0.153 | 0.366 | 2518.243 | 0.325 | No |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 1527.724 | 0.142 | 0.228 | 1387.458 | 0.300 | No |
| weighted_gini_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:gini_restore. | 1594.918 | 0.087 | 0.254 | 1444.375 | 0.278 | Yes |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 1431.930 | 0.094 | 0.285 | 1359.741 | 0.309 | Yes |
| weighted_maximin_time_avg_cri_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:maximin_time_avg_cri_loss. | 1730.161 | 0.197 | 0.382 | 2103.396 | 0.314 | No |
