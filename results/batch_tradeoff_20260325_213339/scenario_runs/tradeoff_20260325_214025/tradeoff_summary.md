# Trade-off Summary

- Result directory: `results/batch_tradeoff_20260325_213339/scenario_runs/tradeoff_20260325_214025`
- Critical locations: `[1, 24]`
- Primary CRI weights: `(w_e=0.133, w_a=0.867)`
- Primary CRI threshold: `0.9`
- Primary critical-access threshold: `0.9`
- Crew mode: `specialized` (power=`1`, road=`1`)
- Bus dispatch mode: `link_only`
- Bus-to-link source: `new_bus_to_link.json`
- Enabled experiments: `['single_triangle', 'single_gini_restore', 'single_maximin_time_avg_cri', 'single_p90_access_restore', 'weighted_gini_restore_l050', 'weighted_maximin_time_avg_cri_l050', 'guardrail_gini_restore']`
- Pareto experiments (triangle vs gini_restore): `['single_triangle', 'single_gini_restore', 'single_maximin_time_avg_cri', 'single_p90_access_restore', 'weighted_gini_restore_l050', 'weighted_maximin_time_avg_cri_l050']`
- Combined CSV: `tradeoff_summary.csv`
- Scatter plot: `tradeoff_scatter.png`
- The combined CSV includes both `power_sequence` and `road_sequence` for specialized-crew runs.

## Optimized Results

| Experiment | Type | Description | Triangle | Gini | Min time-avg CRI | P90 access | Time-avg share access | Pareto |
|---|---|---|---:|---:|---:|---:|---:|---|
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 2862.481 | 0.127 | 0.271 | 2721.633 | 0.324 | No |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 2643.604 | 0.100 | 0.354 | 2598.389 | 0.331 | Yes |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 3184.807 | 0.098 | 0.270 | 2952.890 | 0.290 | Yes |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 2641.383 | 0.127 | 0.315 | 2581.490 | 0.332 | Yes |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 2613.284 | 0.127 | 0.299 | 2540.080 | 0.332 | Yes |
| weighted_gini_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:gini_restore. | 2641.383 | 0.127 | 0.315 | 2581.490 | 0.332 | Yes |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 3088.195 | 0.104 | 0.288 | 2948.310 | 0.292 | No |
| weighted_maximin_time_avg_cri_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:maximin_time_avg_cri_loss. | 2641.383 | 0.127 | 0.315 | 2581.490 | 0.332 | Yes |
