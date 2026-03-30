# Trade-off Summary

- Result directory: `results/dynamic_trial/tradeoff_20260330_142814`
- Critical locations: `[1, 24]`
- Primary CRI weights: `(w_e=0.133, w_a=0.867)`
- Primary CRI threshold: `0.9`
- Primary critical-access threshold: `0.9`
- Crew mode: `specialized` (power=`1`, road=`1`)
- Bus dispatch mode: `link_only`
- Bus-to-link source: `new_bus_to_link.json`
- Enabled experiments: `['single_triangle', 'single_gini_restore', 'single_maximin_time_avg_cri', 'single_p90_access_restore']`
- Pareto experiments (triangle vs gini_restore): `['single_gini_restore', 'single_p90_access_restore']`
- Combined CSV: `tradeoff_summary.csv`
- Scatter plot: `tradeoff_scatter.png`
- The combined CSV includes both `power_sequence` and `road_sequence` for specialized-crew runs.

## Optimized Results

| Experiment | Type | Description | Triangle | Gini | Min time-avg CRI | P90 access | Time-avg share access | Pareto |
|---|---|---|---:|---:|---:|---:|---:|---|
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 170.863 | 0.313 | 0.583 | 352.803 | 0.451 | No |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 170.863 | 0.313 | 0.583 | 352.803 | 0.451 | No |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 167.645 | 0.283 | 0.434 | 272.442 | 0.335 | Yes |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 170.863 | 0.313 | 0.583 | 352.803 | 0.451 | No |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 167.645 | 0.283 | 0.434 | 272.442 | 0.335 | Yes |
