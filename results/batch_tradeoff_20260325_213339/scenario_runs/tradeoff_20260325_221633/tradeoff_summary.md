# Trade-off Summary

- Result directory: `results/batch_tradeoff_20260325_213339/scenario_runs/tradeoff_20260325_221633`
- Critical locations: `[1, 24]`
- Primary CRI weights: `(w_e=0.133, w_a=0.867)`
- Primary CRI threshold: `0.9`
- Primary critical-access threshold: `0.9`
- Crew mode: `specialized` (power=`1`, road=`1`)
- Bus dispatch mode: `link_only`
- Bus-to-link source: `new_bus_to_link.json`
- Enabled experiments: `['single_triangle', 'single_gini_restore', 'single_maximin_time_avg_cri', 'single_p90_access_restore', 'weighted_gini_restore_l050', 'weighted_maximin_time_avg_cri_l050', 'guardrail_gini_restore']`
- Pareto experiments (triangle vs gini_restore): `['single_gini_restore', 'single_maximin_time_avg_cri', 'guardrail_gini_restore']`
- Combined CSV: `tradeoff_summary.csv`
- Scatter plot: `tradeoff_scatter.png`
- The combined CSV includes both `power_sequence` and `road_sequence` for specialized-crew runs.

## Optimized Results

| Experiment | Type | Description | Triangle | Gini | Min time-avg CRI | P90 access | Time-avg share access | Pareto |
|---|---|---|---:|---:|---:|---:|---:|---|
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 3871.971 | 0.135 | 0.228 | 3481.358 | 0.305 | No |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 3398.229 | 0.138 | 0.245 | 3096.330 | 0.303 | No |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 3990.237 | 0.106 | 0.223 | 3286.004 | 0.292 | Yes |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 3357.038 | 0.137 | 0.280 | 3227.528 | 0.310 | Yes |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 3856.381 | 0.138 | 0.234 | 3463.823 | 0.303 | No |
| weighted_gini_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:gini_restore. | 3828.931 | 0.129 | 0.229 | 3481.358 | 0.307 | No |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 3470.234 | 0.128 | 0.223 | 3096.330 | 0.294 | Yes |
| weighted_maximin_time_avg_cri_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:maximin_time_avg_cri_loss. | 3722.196 | 0.138 | 0.212 | 3235.111 | 0.292 | No |
