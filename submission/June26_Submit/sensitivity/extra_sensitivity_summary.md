# Extra Sensitivity Summary

- Selected run: `/home/workenv/results/batch_tradeoff_20260326_151436/scenario_runs/tradeoff_20260326_180346`
- Scenario: `scenario_084` seed `20260408`
- Important note: these rows re-evaluate the selected sequences with the current codebase.
- Weight sensitivity: `(w_e,w_a)=(0.3,0.7)` and `(0.7,0.3)`, thresholds fixed at `0.9`.
- Threshold sensitivity: `cri_threshold` and `critical_access_threshold` jointly set to `0.85` and `0.95`, weights fixed at `(0.133,0.867)`.
- Weighted lambda sensitivity tests four metrics: restoration-time Gini, P90 CRI restoration time, maximin time-averaged CRI loss, and P90 critical-access restoration time.
- Output files: `extra_sensitivity_rows.csv`, `extra_sensitivity_ranges.csv`, `weighted_lambda_rows.csv`.

## Weighted Lambda Results

| Experiment | Weighted metric | Lambda | Triangle | Gini | P90 restore | Min time-avg CRI | P90 access | Time-avg share access |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| weighted_gini_restore_l025_extra | equity:gini_restore | 0.25 | 504.346 | 0.167 | 509.304 | 0.178 | 509.304 | 0.280 |
| weighted_gini_restore_l075_extra | equity:gini_restore | 0.75 | 801.725 | 0.142 | 746.438 | 0.123 | 746.438 | 0.210 |
| weighted_p90_restore_l025_extra | equity:p90_restore | 0.25 | 641.690 | 0.176 | 633.482 | 0.170 | 633.482 | 0.273 |
| weighted_p90_restore_l075_extra | equity:p90_restore | 0.75 | 641.690 | 0.176 | 633.482 | 0.170 | 633.482 | 0.273 |
| weighted_maximin_time_avg_cri_l025_extra | equity:maximin_time_avg_cri_loss | 0.25 | 641.690 | 0.176 | 633.482 | 0.170 | 633.482 | 0.273 |
| weighted_maximin_time_avg_cri_l075_extra | equity:maximin_time_avg_cri_loss | 0.75 | 543.247 | 0.220 | 554.560 | 0.171 | 554.560 | 0.315 |
| weighted_p90_access_restore_l025_extra | critical_access:p90_access_restore | 0.25 | 641.690 | 0.176 | 633.482 | 0.170 | 633.482 | 0.273 |
| weighted_p90_access_restore_l075_extra | critical_access:p90_access_restore | 0.75 | 641.690 | 0.176 | 633.482 | 0.170 | 633.482 | 0.273 |

Chosen balanced weighted case by normalized triangle+gini score: `weighted_gini_restore_l025_extra`.

Best balanced case within each weighted metric family:

- `critical_access:p90_access_restore`: `weighted_p90_access_restore_l025_extra`
- `equity:gini_restore`: `weighted_gini_restore_l025_extra`
- `equity:maximin_time_avg_cri_loss`: `weighted_maximin_time_avg_cri_l025_extra`
- `equity:p90_restore`: `weighted_p90_restore_l025_extra`

## Largest Sensitivity Ranges

### weight_sensitivity

Top Gini ranges:

- `single_var_restore`: 0.175 - 0.225
- `single_p90_access_restore`: 0.225 - 0.270
- `guardrail_gini_restore`: 0.161 - 0.196

Top min time-avg CRI ranges:

- `guardrail_gini_restore`: 0.121 - 0.249
- `single_triangle`: 0.080 - 0.173
- `single_var_restore`: 0.058 - 0.136

### threshold_sensitivity

Top Gini ranges:

- `guardrail_gini_restore`: 0.153 - 0.309
- `single_var_restore`: 0.170 - 0.294
- `single_p90_access_restore`: 0.212 - 0.335

Top min time-avg CRI ranges:

- `baseline_reference`: 0.094 - 0.094
- `guardrail_gini_restore`: 0.254 - 0.254
- `single_gini_restore`: 0.128 - 0.128

