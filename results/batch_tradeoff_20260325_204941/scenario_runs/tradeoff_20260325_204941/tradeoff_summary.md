# Trade-off Summary

- Result directory: `results/batch_tradeoff_20260325_204941/scenario_runs/tradeoff_20260325_204941`
- Critical locations: `[1, 24]`
- Primary CRI weights: `(w_e=0.133, w_a=0.867)`
- Primary CRI threshold: `0.9`
- Primary critical-access threshold: `0.9`
- Crew mode: `specialized` (power=`1`, road=`1`)
- Bus dispatch mode: `link_only`
- Bus-to-link source: `new_bus_to_link.json`
- Pareto experiments (triangle vs gini_restore): `['baseline_reference', 'single_triangle', 'single_var_restore', 'single_gini_restore', 'single_p90_restore', 'single_maximin_time_avg_cri', 'single_p90_access_restore', 'weighted_gini_restore_l025', 'guardrail_gini_restore', 'weighted_p90_restore_l025', 'guardrail_p90_restore', 'weighted_maximin_time_avg_cri_l025', 'guardrail_maximin_time_avg_cri', 'weighted_p90_access_restore_l025', 'guardrail_p90_access_restore']`
- Combined CSV: `tradeoff_summary.csv`
- Scatter plot: `tradeoff_scatter.png`

## Optimized Results

| Experiment | Type | Description | Triangle | Gini | Min time-avg CRI | P90 access | Time-avg share access | Pareto |
|---|---|---|---:|---:|---:|---:|---:|---|
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 2619.690 | 0.136 | 0.281 | 2720.009 | 0.220 | Yes |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 2619.690 | 0.136 | 0.281 | 2720.009 | 0.220 | Yes |
| single_var_restore | single | Equity objective that reduces dispersion in restoration times. | 2619.690 | 0.136 | 0.281 | 2720.009 | 0.220 | Yes |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 2619.690 | 0.136 | 0.281 | 2720.009 | 0.220 | Yes |
| single_p90_restore | single | Equity objective that shortens the time by which 90% of zones recover. | 2619.690 | 0.136 | 0.281 | 2720.009 | 0.220 | Yes |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 2619.690 | 0.136 | 0.281 | 2720.009 | 0.220 | Yes |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 2619.690 | 0.136 | 0.281 | 2720.009 | 0.220 | Yes |
| weighted_gini_restore_l025 | weighted_sum | Weighted trade-off objective. triangle + 0.25 * equity:gini_restore. | 2619.690 | 0.136 | 0.281 | 2720.009 | 0.220 | Yes |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 2619.690 | 0.136 | 0.281 | 2720.009 | 0.220 | Yes |
| weighted_p90_restore_l025 | weighted_sum | Weighted trade-off objective. triangle + 0.25 * equity:p90_restore. | 2619.690 | 0.136 | 0.281 | 2720.009 | 0.220 | Yes |
| guardrail_p90_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:p90_restore. | 2619.690 | 0.136 | 0.281 | 2720.009 | 0.220 | Yes |
| weighted_maximin_time_avg_cri_l025 | weighted_sum | Weighted trade-off objective. triangle + 0.25 * equity:maximin_time_avg_cri_loss. | 2619.690 | 0.136 | 0.281 | 2720.009 | 0.220 | Yes |
| guardrail_maximin_time_avg_cri | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:maximin_time_avg_cri_loss. | 2619.690 | 0.136 | 0.281 | 2720.009 | 0.220 | Yes |
| weighted_p90_access_restore_l025 | weighted_sum | Weighted trade-off objective. triangle + 0.25 * critical_access:p90_access_restore. | 2619.690 | 0.136 | 0.281 | 2720.009 | 0.220 | Yes |
| guardrail_p90_access_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on critical_access:p90_access_restore. | 2619.690 | 0.136 | 0.281 | 2720.009 | 0.220 | Yes |
