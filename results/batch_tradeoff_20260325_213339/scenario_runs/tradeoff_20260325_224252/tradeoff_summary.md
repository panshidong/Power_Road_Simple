# Trade-off Summary

- Result directory: `results/batch_tradeoff_20260325_213339/scenario_runs/tradeoff_20260325_224252`
- Critical locations: `[1, 24]`
- Primary CRI weights: `(w_e=0.133, w_a=0.867)`
- Primary CRI threshold: `0.9`
- Primary critical-access threshold: `0.9`
- Crew mode: `specialized` (power=`1`, road=`1`)
- Bus dispatch mode: `link_only`
- Bus-to-link source: `new_bus_to_link.json`
- Enabled experiments: `['single_triangle', 'single_gini_restore', 'single_maximin_time_avg_cri', 'single_p90_access_restore', 'weighted_gini_restore_l050', 'weighted_maximin_time_avg_cri_l050', 'guardrail_gini_restore']`
- Pareto experiments (triangle vs gini_restore): `['single_gini_restore', 'weighted_gini_restore_l050', 'guardrail_gini_restore']`
- Combined CSV: `tradeoff_summary.csv`
- Scatter plot: `tradeoff_scatter.png`
- The combined CSV includes both `power_sequence` and `road_sequence` for specialized-crew runs.

## Optimized Results

| Experiment | Type | Description | Triangle | Gini | Min time-avg CRI | P90 access | Time-avg share access | Pareto |
|---|---|---|---:|---:|---:|---:|---:|---|
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 1024.258 | 0.293 | 0.303 | 1193.540 | 0.295 | No |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 739.817 | 0.292 | 0.424 | 1058.742 | 0.312 | No |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 1082.167 | 0.216 | 0.306 | 1294.833 | 0.283 | Yes |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 609.215 | 0.336 | 0.595 | 1111.550 | 0.426 | No |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 976.647 | 0.284 | 0.312 | 1152.990 | 0.287 | No |
| weighted_gini_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:gini_restore. | 559.341 | 0.276 | 0.530 | 835.685 | 0.421 | Yes |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 771.576 | 0.238 | 0.461 | 1146.790 | 0.343 | Yes |
| weighted_maximin_time_avg_cri_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:maximin_time_avg_cri_loss. | 566.142 | 0.276 | 0.555 | 822.927 | 0.405 | No |
