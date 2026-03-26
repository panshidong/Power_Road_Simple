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
| single_triangle | single | 4011.070 ± 1245.509 | 0.207 ± 0.012 | 4205.712 ± 1213.427 |
| weighted_maximin_time_avg_cri_l050 | weighted_sum | 4232.822 ± 1273.882 | 0.208 ± 0.011 | 4397.740 ± 1263.929 |
| single_p90_access_restore | single | 4414.369 ± 1340.080 | 0.189 ± 0.012 | 4250.376 ± 1278.813 |
| weighted_gini_restore_l050 | weighted_sum | 4529.141 ± 1440.116 | 0.177 ± 0.011 | 4368.587 ± 1290.881 |
| guardrail_gini_restore | guardrail | 4597.564 ± 1394.071 | 0.171 ± 0.010 | 4431.102 ± 1290.004 |
| single_maximin_time_avg_cri | single | 4856.323 ± 1642.090 | 0.216 ± 0.012 | 5027.652 ± 1579.847 |
| single_gini_restore | single | 4942.842 ± 1479.265 | 0.165 ± 0.010 | 4730.717 ± 1390.779 |
| baseline_reference | baseline | 4988.354 ± 1464.373 | 0.198 ± 0.011 | 4899.384 ± 1402.608 |
