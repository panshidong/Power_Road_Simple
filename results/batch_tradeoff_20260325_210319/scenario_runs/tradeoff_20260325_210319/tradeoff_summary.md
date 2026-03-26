# Trade-off Summary

- Result directory: `results/batch_tradeoff_20260325_210319/scenario_runs/tradeoff_20260325_210319`
- Critical locations: `[1, 24]`
- Primary CRI weights: `(w_e=0.133, w_a=0.867)`
- Primary CRI threshold: `0.9`
- Primary critical-access threshold: `0.9`
- Crew mode: `specialized` (power=`1`, road=`1`)
- Bus dispatch mode: `link_only`
- Bus-to-link source: `new_bus_to_link.json`
- Pareto experiments (triangle vs gini_restore): `['single_triangle', 'single_var_restore', 'single_gini_restore', 'weighted_gini_restore_l025', 'guardrail_gini_restore', 'weighted_p90_restore_l025', 'guardrail_p90_restore', 'weighted_maximin_time_avg_cri_l025', 'weighted_p90_access_restore_l025', 'guardrail_p90_access_restore']`
- Combined CSV: `tradeoff_summary.csv`
- Scatter plot: `tradeoff_scatter.png`
- The combined CSV includes both `power_sequence` and `road_sequence` for specialized-crew runs.

## Optimized Results

| Experiment | Type | Description | Triangle | Gini | Min time-avg CRI | P90 access | Time-avg share access | Pareto |
|---|---|---|---:|---:|---:|---:|---:|---|
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 3922.024 | 0.165 | 0.248 | 3970.883 | 0.285 | No |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 3847.844 | 0.163 | 0.240 | 3970.883 | 0.285 | Yes |
| single_var_restore | single | Equity objective that reduces dispersion in restoration times. | 3847.844 | 0.163 | 0.240 | 3970.883 | 0.285 | Yes |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 3847.844 | 0.163 | 0.240 | 3970.883 | 0.285 | Yes |
| single_p90_restore | single | Equity objective that shortens the time by which 90% of zones recover. | 3922.024 | 0.165 | 0.248 | 3970.883 | 0.285 | No |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 3922.024 | 0.165 | 0.248 | 3970.883 | 0.285 | No |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 3922.024 | 0.165 | 0.248 | 3970.883 | 0.285 | No |
| weighted_gini_restore_l025 | weighted_sum | Weighted trade-off objective. triangle + 0.25 * equity:gini_restore. | 3847.844 | 0.163 | 0.240 | 3970.883 | 0.285 | Yes |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 3847.844 | 0.163 | 0.240 | 3970.883 | 0.285 | Yes |
| weighted_p90_restore_l025 | weighted_sum | Weighted trade-off objective. triangle + 0.25 * equity:p90_restore. | 3847.844 | 0.163 | 0.240 | 3970.883 | 0.285 | Yes |
| guardrail_p90_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:p90_restore. | 3847.844 | 0.163 | 0.240 | 3970.883 | 0.285 | Yes |
| weighted_maximin_time_avg_cri_l025 | weighted_sum | Weighted trade-off objective. triangle + 0.25 * equity:maximin_time_avg_cri_loss. | 3847.844 | 0.163 | 0.240 | 3970.883 | 0.285 | Yes |
| guardrail_maximin_time_avg_cri | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:maximin_time_avg_cri_loss. | 3955.780 | 0.165 | 0.252 | 3970.883 | 0.283 | No |
| weighted_p90_access_restore_l025 | weighted_sum | Weighted trade-off objective. triangle + 0.25 * critical_access:p90_access_restore. | 3847.844 | 0.163 | 0.240 | 3970.883 | 0.285 | Yes |
| guardrail_p90_access_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on critical_access:p90_access_restore. | 3847.844 | 0.163 | 0.240 | 3970.883 | 0.285 | Yes |
