# Trade-off Summary

- Result directory: `results/tradeoff_20260321_163115`
- Critical locations: `[1, 24]`
- Primary CRI weights: `(w_e=0.5, w_a=0.5)`
- Primary CRI threshold: `0.9`
- Primary critical-access threshold: `0.9`
- Bus dispatch mode: `link_only`
- Bus-to-link source: `new_bus_to_link.json`
- Pareto experiments (triangle vs gini_restore): `['baseline_reference', 'single_triangle', 'single_var_restore', 'single_gini_restore', 'single_p90_restore', 'single_maximin_time_avg_cri', 'single_p90_access_restore', 'weighted_gini_restore_l025', 'guardrail_gini_restore', 'weighted_p90_restore_l025', 'guardrail_p90_restore', 'weighted_maximin_time_avg_cri_l025', 'guardrail_maximin_time_avg_cri', 'weighted_p90_access_restore_l025', 'guardrail_p90_access_restore']`
- Combined CSV: `tradeoff_summary.csv`
- Scatter plot: `tradeoff_scatter.png`

## Optimized Results

| Experiment | Type | Description | Triangle | Gini | Min time-avg CRI | P90 access | Time-avg share access | Pareto |
|---|---|---|---:|---:|---:|---:|---:|---|
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 310.830 | 0.351 | 0.525 | 629.845 | 0.481 | Yes |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 305.779 | 0.381 | 0.391 | 617.807 | 0.405 | Yes |
| single_var_restore | single | Equity objective that reduces dispersion in restoration times. | 310.830 | 0.351 | 0.525 | 629.845 | 0.481 | Yes |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 310.830 | 0.351 | 0.525 | 629.845 | 0.481 | Yes |
| single_p90_restore | single | Equity objective that shortens the time by which 90% of zones recover. | 310.830 | 0.351 | 0.525 | 629.845 | 0.481 | Yes |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 310.830 | 0.351 | 0.525 | 629.845 | 0.481 | Yes |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 305.779 | 0.381 | 0.391 | 617.807 | 0.405 | Yes |
| weighted_gini_restore_l025 | weighted_sum | Weighted trade-off objective. triangle + 0.25 * equity:gini_restore. | 310.830 | 0.351 | 0.525 | 629.845 | 0.481 | Yes |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 310.830 | 0.351 | 0.525 | 629.845 | 0.481 | Yes |
| weighted_p90_restore_l025 | weighted_sum | Weighted trade-off objective. triangle + 0.25 * equity:p90_restore. | 310.830 | 0.351 | 0.525 | 629.845 | 0.481 | Yes |
| guardrail_p90_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:p90_restore. | 310.830 | 0.351 | 0.525 | 629.845 | 0.481 | Yes |
| weighted_maximin_time_avg_cri_l025 | weighted_sum | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 310.830 | 0.351 | 0.525 | 629.845 | 0.481 | Yes |
| guardrail_maximin_time_avg_cri | guardrail | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 310.830 | 0.351 | 0.525 | 629.845 | 0.481 | Yes |
| weighted_p90_access_restore_l025 | weighted_sum | Weighted trade-off objective. triangle + 0.25 * critical_access:p90_access_restore. | 305.779 | 0.381 | 0.391 | 617.807 | 0.405 | Yes |
| guardrail_p90_access_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on critical_access:p90_access_restore. | 305.779 | 0.381 | 0.391 | 617.807 | 0.405 | Yes |
