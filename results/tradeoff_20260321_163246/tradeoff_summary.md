# Trade-off Summary

- Result directory: `results/tradeoff_20260321_163246`
- Critical locations: `[1, 24]`
- Primary CRI weights: `(w_e=0.5, w_a=0.5)`
- Primary CRI threshold: `0.9`
- Primary critical-access threshold: `0.9`
- Bus dispatch mode: `link_only`
- Bus-to-link source: `new_bus_to_link.json`
- Pareto experiments (triangle vs gini_restore): `['single_triangle', 'single_gini_restore', 'weighted_gini_restore_l050', 'weighted_gini_restore_l100', 'weighted_p90_restore_l025', 'weighted_p90_restore_l050', 'weighted_p90_restore_l100', 'weighted_maximin_time_avg_cri_l050', 'weighted_maximin_time_avg_cri_l100', 'guardrail_p90_access_restore']`
- Combined CSV: `tradeoff_summary.csv`
- Scatter plot: `tradeoff_scatter.png`

## Optimized Results

| Experiment | Type | Description | Triangle | Gini | Min time-avg CRI | P90 access | Time-avg share access | Pareto |
|---|---|---|---:|---:|---:|---:|---:|---|
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 310.830 | 0.351 | 0.525 | 629.845 | 0.481 | No |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 137.096 | 0.356 | 0.433 | 296.655 | 0.573 | Yes |
| single_var_restore | single | Equity objective that reduces dispersion in restoration times. | 137.876 | 0.356 | 0.430 | 296.655 | 0.573 | No |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 351.295 | 0.261 | 0.185 | 496.014 | 0.323 | Yes |
| single_p90_restore | single | Equity objective that shortens the time by which 90% of zones recover. | 137.876 | 0.356 | 0.430 | 296.655 | 0.573 | No |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 309.409 | 0.284 | 0.575 | 479.939 | 0.472 | No |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 162.114 | 0.339 | 0.282 | 256.110 | 0.391 | No |
| weighted_gini_restore_l025 | weighted_sum | Weighted trade-off objective. triangle + 0.25 * equity:gini_restore. | 187.331 | 0.310 | 0.485 | 379.325 | 0.488 | No |
| weighted_gini_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:gini_restore. | 201.552 | 0.282 | 0.430 | 303.187 | 0.443 | Yes |
| weighted_gini_restore_l100 | weighted_sum | Weighted trade-off objective. triangle + 1.0 * equity:gini_restore. | 201.552 | 0.282 | 0.430 | 303.187 | 0.443 | Yes |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 199.949 | 0.330 | 0.266 | 308.973 | 0.384 | No |
| weighted_p90_restore_l025 | weighted_sum | Weighted trade-off objective. triangle + 0.25 * equity:p90_restore. | 201.552 | 0.282 | 0.430 | 303.187 | 0.443 | Yes |
| weighted_p90_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:p90_restore. | 201.552 | 0.282 | 0.430 | 303.187 | 0.443 | Yes |
| weighted_p90_restore_l100 | weighted_sum | Weighted trade-off objective. triangle + 1.0 * equity:p90_restore. | 155.021 | 0.325 | 0.488 | 315.031 | 0.482 | Yes |
| guardrail_p90_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:p90_restore. | 160.214 | 0.376 | 0.318 | 300.723 | 0.413 | No |
| weighted_maximin_time_avg_cri_l025 | weighted_sum | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 182.525 | 0.380 | 0.316 | 356.276 | 0.500 | No |
| weighted_maximin_time_avg_cri_l050 | weighted_sum | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 185.462 | 0.295 | 0.460 | 379.315 | 0.468 | Yes |
| weighted_maximin_time_avg_cri_l100 | weighted_sum | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 292.388 | 0.282 | 0.539 | 430.867 | 0.461 | Yes |
| guardrail_maximin_time_avg_cri | guardrail | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 347.328 | 0.313 | 0.543 | 531.105 | 0.507 | No |
| weighted_p90_access_restore_l025 | weighted_sum | Weighted trade-off objective. triangle + 0.25 * critical_access:p90_access_restore. | 199.203 | 0.362 | 0.258 | 367.980 | 0.484 | No |
| weighted_p90_access_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * critical_access:p90_access_restore. | 199.203 | 0.362 | 0.258 | 367.980 | 0.484 | No |
| weighted_p90_access_restore_l100 | weighted_sum | Weighted trade-off objective. triangle + 1.0 * critical_access:p90_access_restore. | 160.567 | 0.368 | 0.500 | 295.031 | 0.547 | No |
| guardrail_p90_access_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on critical_access:p90_access_restore. | 137.096 | 0.356 | 0.433 | 296.655 | 0.573 | Yes |
