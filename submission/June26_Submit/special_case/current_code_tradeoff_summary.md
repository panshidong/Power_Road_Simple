# Trade-off Summary

- Result directory: `/home/workenv/results/batch_tradeoff_20260326_151436/scenario_runs/tradeoff_20260326_180346/current_code_special_case_20260604/tradeoff_20260604_141137`
- Critical locations: `[1, 24]`
- Primary CRI weights: `(w_e=0.133, w_a=0.867)`
- Primary CRI threshold: `0.9`
- Primary critical-access threshold: `0.9`
- Crew mode: `specialized` (power=`1`, road=`1`)
- Bus dispatch mode: `link_only`
- Bus-to-link source: `new_bus_to_link.json`
- Enabled experiments: `['single_triangle', 'single_gini_restore', 'single_maximin_time_avg_cri', 'single_p90_access_restore', 'weighted_gini_restore_l050', 'weighted_maximin_time_avg_cri_l050', 'guardrail_gini_restore']`
- Pareto experiments (triangle vs gini_restore): `['single_triangle', 'single_gini_restore', 'weighted_gini_restore_l050', 'guardrail_gini_restore']`
- Combined CSV: `tradeoff_summary.csv`
- Scatter plot: `tradeoff_scatter.png`
- The combined CSV includes both `power_sequence` and `road_sequence` for specialized-crew runs.

## Optimized Results

| Experiment | Type | Description | Triangle | Gini | Min time-avg CRI | P90 access | Time-avg share access | Pareto |
|---|---|---|---:|---:|---:|---:|---:|---|
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 2105.151 | 0.186 | 0.094 | 1871.258 | 0.237 | No |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 420.560 | 0.243 | 0.319 | 467.914 | 0.347 | Yes |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 1664.082 | 0.134 | 0.129 | 1458.784 | 0.223 | Yes |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 1707.724 | 0.184 | 0.268 | 1574.362 | 0.284 | No |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 488.809 | 0.182 | 0.215 | 484.936 | 0.293 | No |
| weighted_gini_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:gini_restore. | 732.274 | 0.154 | 0.150 | 689.447 | 0.257 | Yes |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 477.108 | 0.177 | 0.251 | 509.054 | 0.297 | Yes |
| weighted_maximin_time_avg_cri_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:maximin_time_avg_cri_loss. | 632.533 | 0.196 | 0.194 | 663.166 | 0.283 | No |
