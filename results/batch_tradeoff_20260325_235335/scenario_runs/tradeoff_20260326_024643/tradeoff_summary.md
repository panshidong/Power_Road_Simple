# Trade-off Summary

- Result directory: `results/batch_tradeoff_20260325_235335/scenario_runs/tradeoff_20260326_024643`
- Critical locations: `[1, 24]`
- Primary CRI weights: `(w_e=0.133, w_a=0.867)`
- Primary CRI threshold: `0.9`
- Primary critical-access threshold: `0.9`
- Crew mode: `specialized` (power=`1`, road=`1`)
- Bus dispatch mode: `link_only`
- Bus-to-link source: `new_bus_to_link.json`
- Enabled experiments: `['single_triangle', 'single_gini_restore', 'single_maximin_time_avg_cri', 'single_p90_access_restore', 'weighted_gini_restore_l050', 'weighted_maximin_time_avg_cri_l050', 'guardrail_gini_restore']`
- Pareto experiments (triangle vs gini_restore): `['single_gini_restore', 'guardrail_gini_restore', 'weighted_maximin_time_avg_cri_l050']`
- Combined CSV: `tradeoff_summary.csv`
- Scatter plot: `tradeoff_scatter.png`
- The combined CSV includes both `power_sequence` and `road_sequence` for specialized-crew runs.

## Optimized Results

| Experiment | Type | Description | Triangle | Gini | Min time-avg CRI | P90 access | Time-avg share access | Pareto |
|---|---|---|---:|---:|---:|---:|---:|---|
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 527.787 | 0.190 | 0.483 | 682.167 | 0.350 | No |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 431.955 | 0.191 | 0.402 | 530.531 | 0.358 | No |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 741.044 | 0.121 | 0.262 | 778.450 | 0.201 | Yes |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 717.873 | 0.218 | 0.506 | 1000.053 | 0.330 | No |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 488.990 | 0.187 | 0.416 | 505.079 | 0.373 | No |
| weighted_gini_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:gini_restore. | 483.875 | 0.167 | 0.443 | 587.730 | 0.306 | No |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 405.612 | 0.163 | 0.402 | 464.744 | 0.325 | Yes |
| weighted_maximin_time_avg_cri_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:maximin_time_avg_cri_loss. | 437.929 | 0.157 | 0.498 | 633.923 | 0.244 | Yes |
