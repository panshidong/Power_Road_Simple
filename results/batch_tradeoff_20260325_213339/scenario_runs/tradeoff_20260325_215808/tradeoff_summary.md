# Trade-off Summary

- Result directory: `results/batch_tradeoff_20260325_213339/scenario_runs/tradeoff_20260325_215808`
- Critical locations: `[1, 24]`
- Primary CRI weights: `(w_e=0.133, w_a=0.867)`
- Primary CRI threshold: `0.9`
- Primary critical-access threshold: `0.9`
- Crew mode: `specialized` (power=`1`, road=`1`)
- Bus dispatch mode: `link_only`
- Bus-to-link source: `new_bus_to_link.json`
- Enabled experiments: `['single_triangle', 'single_gini_restore', 'single_maximin_time_avg_cri', 'single_p90_access_restore', 'weighted_gini_restore_l050', 'weighted_maximin_time_avg_cri_l050', 'guardrail_gini_restore']`
- Pareto experiments (triangle vs gini_restore): `['baseline_reference', 'single_triangle', 'single_p90_access_restore', 'weighted_gini_restore_l050', 'guardrail_gini_restore', 'weighted_maximin_time_avg_cri_l050']`
- Combined CSV: `tradeoff_summary.csv`
- Scatter plot: `tradeoff_scatter.png`
- The combined CSV includes both `power_sequence` and `road_sequence` for specialized-crew runs.

## Optimized Results

| Experiment | Type | Description | Triangle | Gini | Min time-avg CRI | P90 access | Time-avg share access | Pareto |
|---|---|---|---:|---:|---:|---:|---:|---|
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 3380.129 | 0.154 | 0.191 | 2991.247 | 0.237 | Yes |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 2985.544 | 0.231 | 0.235 | 2936.967 | 0.337 | Yes |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 3470.974 | 0.144 | 0.195 | 2957.204 | 0.221 | No |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 4101.080 | 0.244 | 0.284 | 4641.008 | 0.321 | No |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 3299.625 | 0.194 | 0.195 | 2957.204 | 0.290 | Yes |
| weighted_gini_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:gini_restore. | 3380.129 | 0.154 | 0.191 | 2991.247 | 0.237 | Yes |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 3445.194 | 0.143 | 0.184 | 2936.967 | 0.230 | Yes |
| weighted_maximin_time_avg_cri_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:maximin_time_avg_cri_loss. | 3118.740 | 0.209 | 0.237 | 2957.289 | 0.336 | Yes |
