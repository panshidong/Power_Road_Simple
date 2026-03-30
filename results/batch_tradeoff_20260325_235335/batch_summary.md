# Batch Trade-off Summary

- Scenario count: `100`
- Scenario seed start: `20260325`
- Power failures per disaster: `(8, 15)`
- Road failures per disaster: `(3, 13)`
- Road capacity-drop fraction range: `(0.5, 1.0)`
- Crew mode: `specialized` (power=`1`, road=`1`)
- CRI weights: `[(0.133, 0.867)]`
- Scenario manifest: `scenario_manifest.csv`
- Per-scenario rows: `scenario_tradeoff_rows.csv`
- Aggregate summary: `aggregate_tradeoff_summary.csv`
- Error-bar trade-off figure: `aggregate_tradeoff_errorbars.png`

## Aggregate Means

| Experiment | Type | Triangle mean ± CI | Gini mean ± CI | P90 access mean ± CI |
|---|---|---:|---:|---:|
| single_triangle | single | 2293.025 ± 588.180 | 0.231 ± 0.013 | 2572.446 ± 622.862 |
| single_p90_access_restore | single | 2935.141 ± 904.549 | 0.204 ± 0.012 | 2704.707 ± 803.412 |
| guardrail_gini_restore | guardrail | 3165.478 ± 979.687 | 0.169 ± 0.010 | 3185.827 ± 944.604 |
| weighted_maximin_time_avg_cri_l050 | weighted_sum | 3173.215 ± 982.977 | 0.233 ± 0.013 | 3565.511 ± 1013.386 |
| weighted_gini_restore_l050 | weighted_sum | 3364.255 ± 954.768 | 0.171 ± 0.011 | 3213.980 ± 836.290 |
| single_maximin_time_avg_cri | single | 4305.178 ± 1448.442 | 0.242 ± 0.011 | 4916.668 ± 1591.224 |
| baseline_reference | baseline | 4988.354 ± 1464.373 | 0.198 ± 0.011 | 4899.384 ± 1402.608 |
| single_gini_restore | single | 5029.536 ± 1539.790 | 0.139 ± 0.009 | 4382.582 ± 1250.266 |
