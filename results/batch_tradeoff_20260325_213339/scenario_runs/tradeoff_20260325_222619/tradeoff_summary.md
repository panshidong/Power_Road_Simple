# Trade-off Summary

- Result directory: `results/batch_tradeoff_20260325_213339/scenario_runs/tradeoff_20260325_222619`
- Critical locations: `[1, 24]`
- Primary CRI weights: `(w_e=0.133, w_a=0.867)`
- Primary CRI threshold: `0.9`
- Primary critical-access threshold: `0.9`
- Crew mode: `specialized` (power=`1`, road=`1`)
- Bus dispatch mode: `link_only`
- Bus-to-link source: `new_bus_to_link.json`
- Enabled experiments: `['single_triangle', 'single_gini_restore', 'single_maximin_time_avg_cri', 'single_p90_access_restore', 'weighted_gini_restore_l050', 'weighted_maximin_time_avg_cri_l050', 'guardrail_gini_restore']`
- Pareto experiments (triangle vs gini_restore): `['weighted_gini_restore_l050', 'weighted_maximin_time_avg_cri_l050']`
- Combined CSV: `tradeoff_summary.csv`
- Scatter plot: `tradeoff_scatter.png`
- The combined CSV includes both `power_sequence` and `road_sequence` for specialized-crew runs.

## Optimized Results

| Experiment | Type | Description | Triangle | Gini | Min time-avg CRI | P90 access | Time-avg share access | Pareto |
|---|---|---|---:|---:|---:|---:|---:|---|
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 2855.976 | 0.232 | 0.288 | 3107.042 | 0.278 | No |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 2329.310 | 0.269 | 0.216 | 2193.789 | 0.279 | No |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 2677.779 | 0.196 | 0.240 | 2698.340 | 0.309 | No |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 2800.987 | 0.279 | 0.324 | 3107.042 | 0.308 | No |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 2352.561 | 0.192 | 0.208 | 2193.789 | 0.251 | No |
| weighted_gini_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:gini_restore. | 1818.460 | 0.187 | 0.289 | 2007.534 | 0.245 | Yes |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 2350.397 | 0.186 | 0.189 | 2193.789 | 0.217 | No |
| weighted_maximin_time_avg_cri_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:maximin_time_avg_cri_loss. | 1848.806 | 0.182 | 0.279 | 2007.534 | 0.243 | Yes |
