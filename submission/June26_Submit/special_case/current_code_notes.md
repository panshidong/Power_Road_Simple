# Current-Code Special Case Notes

This folder re-runs the representative special case with the current dynamic evaluator.

## Why the old scatter and sensitivity values differed

The original special-case result folder was generated with older artifacts. For example, the old weighted Gini run contains a single `dispatch_s.txt` file under:

`weighted_gini_restore_l050_20260326_180500/best/dispatch_s.txt`

The current-code re-evaluation folders contain per-state TAP-B snapshots, for example:

`extra_sensitivity_20260604_123941/evaluations/weighted_gini_restore_l050__nominal_recheck__we0.133_wa0.867_ct0.90_at0.90/state_snapshots/*`

This difference matches the later dynamic-simulation correction: travel times are updated after repair events, and crew dispatch is based on the latest network state.

## Current-code optimized special case

The current run uses the same disaster scenario (`scenario_084`, seed `20260408`) and the same CRI weights `(w_e, w_a) = (0.133, 0.867)`.

Variance was intentionally excluded from this figure, following the latest reporting preference.

| Experiment | Triangle | Gini restore | P90 restore | Pareto |
|---|---:|---:|---:|---|
| baseline_reference | 2105.151 | 0.186 | 1871.258 | No |
| single_triangle | 420.560 | 0.243 | 467.914 | Yes |
| single_gini_restore | 1664.082 | 0.134 | 1458.784 | Yes |
| single_maximin_time_avg_cri | 1707.724 | 0.184 | 1574.362 | No |
| single_p90_access_restore | 488.809 | 0.182 | 484.936 | No |
| weighted_gini_restore_l050 | 732.274 | 0.154 | 689.447 | Yes |
| guardrail_gini_restore | 477.108 | 0.177 | 509.054 | Yes |
| weighted_maximin_time_avg_cri_l050 | 632.533 | 0.196 | 663.166 | No |

The weighted Gini result now aligns with the current-code weighted sensitivity scale (`triangle = 732.274`), rather than the original old-result scatter scale (`triangle = 1846.213`).
