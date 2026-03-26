# Trade-off Summary

- Result directory: `results/batch_tradeoff_20260325_210830/scenario_runs/tradeoff_20260325_210830`
- Critical locations: `[1, 24]`
- Primary CRI weights: `(w_e=0.133, w_a=0.867)`
- Primary CRI threshold: `0.9`
- Primary critical-access threshold: `0.9`
- Crew mode: `specialized` (power=`1`, road=`1`)
- Bus dispatch mode: `link_only`
- Bus-to-link source: `new_bus_to_link.json`
- Pareto experiments (triangle vs gini_restore): `['single_gini_restore', 'weighted_gini_restore_l100', 'guardrail_gini_restore']`
- Combined CSV: `tradeoff_summary.csv`
- Scatter plot: `tradeoff_scatter.png`
- The combined CSV includes both `power_sequence` and `road_sequence` for specialized-crew runs.

## Optimized Results

| Experiment | Type | Description | Triangle | Gini | Min time-avg CRI | P90 access | Time-avg share access | Pareto |
|---|---|---|---:|---:|---:|---:|---:|---|
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 4949.903 | 0.195 | 0.302 | 5338.734 | 0.313 | No |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 2945.067 | 0.218 | 0.340 | 3169.093 | 0.358 | No |
| single_var_restore | single | Equity objective that reduces dispersion in restoration times. | 3303.609 | 0.135 | 0.173 | 2882.162 | 0.279 | No |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 4041.660 | 0.089 | 0.191 | 3312.686 | 0.272 | Yes |
| single_p90_restore | single | Equity objective that shortens the time by which 90% of zones recover. | 2855.560 | 0.150 | 0.192 | 2474.804 | 0.269 | No |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 4779.494 | 0.249 | 0.431 | 5380.297 | 0.448 | No |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 2933.237 | 0.149 | 0.189 | 2546.103 | 0.267 | No |
| weighted_gini_restore_l025 | weighted_sum | Weighted trade-off objective. triangle + 0.25 * equity:gini_restore. | 2719.448 | 0.197 | 0.255 | 2560.400 | 0.260 | No |
| weighted_gini_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:gini_restore. | 4350.384 | 0.118 | 0.171 | 3753.112 | 0.265 | No |
| weighted_gini_restore_l100 | weighted_sum | Weighted trade-off objective. triangle + 1.0 * equity:gini_restore. | 3257.802 | 0.122 | 0.192 | 2960.539 | 0.260 | Yes |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 1918.458 | 0.125 | 0.307 | 1572.585 | 0.404 | Yes |
| weighted_p90_restore_l025 | weighted_sum | Weighted trade-off objective. triangle + 0.25 * equity:p90_restore. | 2183.418 | 0.180 | 0.271 | 2071.851 | 0.360 | No |
| weighted_p90_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:p90_restore. | 2501.725 | 0.304 | 0.526 | 2687.583 | 0.545 | No |
| weighted_p90_restore_l100 | weighted_sum | Weighted trade-off objective. triangle + 1.0 * equity:p90_restore. | 2675.153 | 0.201 | 0.269 | 2560.400 | 0.270 | No |
| guardrail_p90_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:p90_restore. | 2045.354 | 0.233 | 0.469 | 2345.740 | 0.452 | No |
| weighted_maximin_time_avg_cri_l025 | weighted_sum | Weighted trade-off objective. triangle + 0.25 * equity:maximin_time_avg_cri_loss. | 3829.087 | 0.187 | 0.277 | 3641.838 | 0.319 | No |
| weighted_maximin_time_avg_cri_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:maximin_time_avg_cri_loss. | 3829.087 | 0.187 | 0.277 | 3641.838 | 0.319 | No |
| weighted_maximin_time_avg_cri_l100 | weighted_sum | Weighted trade-off objective. triangle + 1.0 * equity:maximin_time_avg_cri_loss. | 3605.462 | 0.230 | 0.414 | 3001.010 | 0.446 | No |
| guardrail_maximin_time_avg_cri | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:maximin_time_avg_cri_loss. | 2547.544 | 0.241 | 0.495 | 3007.190 | 0.452 | No |
| weighted_p90_access_restore_l025 | weighted_sum | Weighted trade-off objective. triangle + 0.25 * critical_access:p90_access_restore. | 3203.870 | 0.215 | 0.279 | 3170.029 | 0.307 | No |
| weighted_p90_access_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * critical_access:p90_access_restore. | 3361.085 | 0.262 | 0.438 | 3866.967 | 0.420 | No |
| weighted_p90_access_restore_l100 | weighted_sum | Weighted trade-off objective. triangle + 1.0 * critical_access:p90_access_restore. | 2695.554 | 0.189 | 0.227 | 2447.421 | 0.254 | No |
| guardrail_p90_access_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on critical_access:p90_access_restore. | 2045.354 | 0.233 | 0.469 | 2345.740 | 0.452 | No |
