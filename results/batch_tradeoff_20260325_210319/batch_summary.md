# Batch Trade-off Summary

- Scenario count: `1`
- Scenario seed start: `41`
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
| guardrail_gini_restore | guardrail | 3847.844 ± 0.000 | 0.163 ± 0.000 | 3970.883 ± 0.000 |
| guardrail_p90_access_restore | guardrail | 3847.844 ± 0.000 | 0.163 ± 0.000 | 3970.883 ± 0.000 |
| guardrail_p90_restore | guardrail | 3847.844 ± 0.000 | 0.163 ± 0.000 | 3970.883 ± 0.000 |
| single_gini_restore | single | 3847.844 ± 0.000 | 0.163 ± 0.000 | 3970.883 ± 0.000 |
| single_triangle | single | 3847.844 ± 0.000 | 0.163 ± 0.000 | 3970.883 ± 0.000 |
| single_var_restore | single | 3847.844 ± 0.000 | 0.163 ± 0.000 | 3970.883 ± 0.000 |
| weighted_gini_restore_l025 | weighted_sum | 3847.844 ± 0.000 | 0.163 ± 0.000 | 3970.883 ± 0.000 |
| weighted_maximin_time_avg_cri_l025 | weighted_sum | 3847.844 ± 0.000 | 0.163 ± 0.000 | 3970.883 ± 0.000 |
| weighted_p90_access_restore_l025 | weighted_sum | 3847.844 ± 0.000 | 0.163 ± 0.000 | 3970.883 ± 0.000 |
| weighted_p90_restore_l025 | weighted_sum | 3847.844 ± 0.000 | 0.163 ± 0.000 | 3970.883 ± 0.000 |
| baseline_reference | baseline | 3922.024 ± 0.000 | 0.165 ± 0.000 | 3970.883 ± 0.000 |
| single_maximin_time_avg_cri | single | 3922.024 ± 0.000 | 0.165 ± 0.000 | 3970.883 ± 0.000 |
| single_p90_access_restore | single | 3922.024 ± 0.000 | 0.165 ± 0.000 | 3970.883 ± 0.000 |
| single_p90_restore | single | 3922.024 ± 0.000 | 0.165 ± 0.000 | 3970.883 ± 0.000 |
| guardrail_maximin_time_avg_cri | guardrail | 3955.780 ± 0.000 | 0.165 ± 0.000 | 3970.883 ± 0.000 |
