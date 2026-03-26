# Trade-off Summary

- Result directory: `results/batch_tradeoff_20260325_213339/scenario_runs/tradeoff_20260325_231400`
- Critical locations: `[1, 24]`
- Primary CRI weights: `(w_e=0.133, w_a=0.867)`
- Primary CRI threshold: `0.9`
- Primary critical-access threshold: `0.9`
- Crew mode: `specialized` (power=`1`, road=`1`)
- Bus dispatch mode: `link_only`
- Bus-to-link source: `new_bus_to_link.json`
- Enabled experiments: `['single_triangle', 'single_gini_restore', 'single_maximin_time_avg_cri', 'single_p90_access_restore', 'weighted_gini_restore_l050', 'weighted_maximin_time_avg_cri_l050', 'guardrail_gini_restore']`
- Pareto experiments (triangle vs gini_restore): `['single_triangle', 'single_gini_restore', 'single_maximin_time_avg_cri', 'weighted_gini_restore_l050', 'guardrail_gini_restore']`
- Combined CSV: `tradeoff_summary.csv`
- Scatter plot: `tradeoff_scatter.png`
- The combined CSV includes both `power_sequence` and `road_sequence` for specialized-crew runs.

## Optimized Results

| Experiment | Type | Description | Triangle | Gini | Min time-avg CRI | P90 access | Time-avg share access | Pareto |
|---|---|---|---:|---:|---:|---:|---:|---|
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 1068.645 | 0.200 | 0.292 | 1212.157 | 0.298 | No |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 562.347 | 0.198 | 0.429 | 804.295 | 0.231 | Yes |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 761.234 | 0.124 | 0.253 | 806.673 | 0.205 | Yes |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 608.013 | 0.194 | 0.368 | 788.874 | 0.219 | Yes |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 764.351 | 0.172 | 0.312 | 814.791 | 0.321 | No |
| weighted_gini_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:gini_restore. | 621.674 | 0.149 | 0.266 | 721.528 | 0.194 | Yes |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 648.881 | 0.148 | 0.267 | 748.170 | 0.198 | Yes |
| weighted_maximin_time_avg_cri_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:maximin_time_avg_cri_loss. | 608.475 | 0.195 | 0.367 | 788.874 | 0.220 | No |
