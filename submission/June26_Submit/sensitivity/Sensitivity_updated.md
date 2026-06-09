# 5.3 Representative Case Illustration

This section revisits one representative scenario from the 100-scenario batch run and uses it as a focused sensitivity case. The selected run is tradeoff_20260326_180346, corresponding to scenario_084 with random seed 20260408. The disruption contains 9 failed power buses and 10 damaged road links; the road capacity-drop fractions range from 60.5% to 92.1%. All sensitivity results below use the same damaged assets and the same critical locations, nodes 1 and 24.

The purpose of this representative case is not to replace the aggregate 100-scenario results. Instead, it illustrates how the interpretation of a repair strategy changes when the CRI definition, adequacy threshold, and weighted equity objective are perturbed. The analysis was re-run with the current event-driven simulation logic, in which TAP-B travel times are updated after repair events and subsequent crew dispatch uses the latest network state.

## Sensitivity to CRI Weights and Adequacy Thresholds

Two forms of post-optimization sensitivity were applied to the optimized sequences from the selected scenario. First, CRI weight sensitivity re-evaluated each sequence under (w_e, w_a) = (0.3, 0.7) and (0.7, 0.3), while keeping the adequacy thresholds at 0.9. Second, threshold sensitivity jointly changed the CRI adequacy threshold and the critical-access threshold from the nominal value 0.9 to 0.85 and 0.95, while keeping the nominal CRI weights (0.133, 0.867). These runs do not re-optimize the sequence; they ask whether the conclusions remain recognizable when the definition of community recovery is changed.

**Table 5. Sensitivity design for the representative case**

| Dimension | Values tested | What changes in the metrics |
| --- | --- | --- |
| CRI weights | (0.3, 0.7), (0.7, 0.3) | Changes the balance between electric service E_z(t) and accessibility A_z(t) in CRI_z(t). |
| Adequacy thresholds | 0.85, 0.95 | Changes the first-hit restoration times used by Gini, variance, P90 restoration, and access summaries. |
| Weighted objective lambda | 0.25, 0.75 | Changes how strongly a selected equity or access metric is weighted against the resilience triangle around the original lambda = 0.5 setting. |

A useful distinction appears in the results. The two sensitivity dimensions should not be pooled into a single range, because they affect different parts of the calculation. Changing the CRI weights changes the trajectory of CRI_z(t) itself, while changing the adequacy threshold only changes the first-hit time extracted from an already fixed trajectory. The results are therefore separated into weight sensitivity and threshold sensitivity.

**Table 6. Sensitivity ranges under CRI-weight perturbations, thresholds fixed at 0.9**

| Strategy | Gini range | Min avg CRI | P90 rest. | P90 access | Avg share access |
| --- | --- | --- | --- | --- | --- |
| Baseline | 0.186 - 0.192 | 0.094 - 0.193 | 1871.258 - 1871.258 | 1871.258 - 1871.258 | 0.237 - 0.237 |
| Triangle | 0.208 - 0.217 | 0.080 - 0.173 | 474.170 - 507.498 | 507.498 - 507.498 | 0.302 - 0.302 |
| Variance | 0.171 - 0.225 | 0.058 - 0.149 | 699.402 - 735.411 | 699.402 - 699.402 | 0.309 - 0.309 |
| Gini | 0.135 - 0.157 | 0.044 - 0.128 | 1697.868 - 1705.758 | 1679.459 - 1679.459 | 0.200 - 0.200 |
| Maximin CRI | 0.163 - 0.181 | 0.075 - 0.149 | 797.056 - 805.670 | 805.670 - 805.670 | 0.215 - 0.215 |
| P90 Access | 0.224 - 0.270 | 0.161 - 0.237 | 586.755 - 649.204 | 649.204 - 649.204 | 0.301 - 0.301 |
| Weighted Gini | 0.133 - 0.144 | 0.151 - 0.194 | 675.190 - 690.993 | 690.993 - 690.993 | 0.203 - 0.203 |
| Guardrail Gini | 0.161 - 0.196 | 0.121 - 0.254 | 848.889 - 859.099 | 859.099 - 859.099 | 0.302 - 0.302 |
| Weighted Maximin | 0.202 - 0.215 | 0.062 - 0.157 | 1362.951 - 1406.189 | 1362.951 - 1362.951 | 0.311 - 0.311 |

Table 6 isolates the effect of changing w_e and w_a while holding both thresholds at 0.9. The access-only columns remain fixed because A_z(t) and the access threshold are unchanged. The min time-averaged CRI, however, can move noticeably because the whole CRI trajectory is reweighted. This is most visible for the guardrail sequence, where the min time-averaged CRI ranges from 0.121 to 0.254, and for the baseline sequence, where it ranges from 0.094 to 0.193. The Gini-oriented sequence remains one of the strongest restoration-time inequality options under weight changes, with a Gini range of 0.135 to 0.157.

**Table 7. Sensitivity ranges under adequacy-threshold perturbations, CRI weights fixed at (0.133, 0.867)**

| Strategy | Gini range | Min avg CRI | P90 rest. | P90 access | Avg share access |
| --- | --- | --- | --- | --- | --- |
| Baseline | 0.176 - 0.261 | 0.094 - 0.094 | 1871.258 - 1871.258 | 1871.258 - 1880.389 | 0.220 - 0.275 |
| Triangle | 0.185 - 0.299 | 0.173 - 0.173 | 499.093 - 507.498 | 507.498 - 507.498 | 0.241 - 0.352 |
| Variance | 0.170 - 0.294 | 0.149 - 0.149 | 699.402 - 735.411 | 699.402 - 724.609 | 0.283 - 0.332 |
| Gini | 0.135 - 0.218 | 0.128 - 0.128 | 1679.459 - 1705.758 | 1679.459 - 1679.459 | 0.193 - 0.234 |
| Maximin CRI | 0.152 - 0.260 | 0.149 - 0.149 | 805.670 - 805.670 | 805.670 - 805.670 | 0.191 - 0.240 |
| P90 Access | 0.212 - 0.335 | 0.237 - 0.237 | 649.204 - 649.204 | 649.204 - 649.204 | 0.246 - 0.340 |
| Weighted Gini | 0.132 - 0.221 | 0.151 - 0.151 | 690.993 - 690.993 | 690.993 - 690.993 | 0.169 - 0.214 |
| Guardrail Gini | 0.153 - 0.309 | 0.254 - 0.254 | 859.099 - 859.099 | 859.099 - 859.099 | 0.249 - 0.317 |
| Weighted Maximin | 0.175 - 0.284 | 0.157 - 0.157 | 1362.951 - 1406.189 | 1362.951 - 1393.218 | 0.273 - 0.369 |

Table 7 isolates the effect of changing the adequacy thresholds while holding the CRI weights fixed. In this table, min time-averaged CRI is constant for each strategy because it is an area-under-trajectory metric rather than a first-hit threshold metric. By contrast, Gini restoration, P90 restoration, P90 access, and time-averaged share access can change because they depend on when each zone first crosses a selected threshold. This split makes the interpretation cleaner: threshold sensitivity mainly tests how demanding the recovery definition is, while weight sensitivity tests how electric service and accessibility are valued inside CRI.

The equality between P90 restoration and P90 access in several rows is also interpretable rather than simply an output error. In this experiment, the nominal CRI places high weight on accessibility, w_a = 0.867. When access is the dominant component and electric-service recovery is not the limiting condition for the 90th percentile zone, the first time at which CRI_z(t) crosses 0.9 can coincide with the first time at which A_z(t) crosses 0.9.

## Weighted Equity Objective Sensitivity

The original batch run used lambda = 0.5 for the weighted objectives. Therefore, the lambda sensitivity is centered around that setting by testing lambda = 0.25 and lambda = 0.75, rather than comparing against a much larger value. The original weighted test used only restoration-time Gini as the equity index. To avoid over-interpreting a single fairness definition, the weighted objective was expanded to four alternatives: restoration-time Gini, P90 CRI restoration time, maximin time-averaged CRI loss, and P90 critical-access restoration time. Each metric was normalized by the baseline reference value before being combined with the resilience triangle, so the lambda values are comparable across metrics.

**Table 8. Weighted objective sensitivity across equity and access metrics**

| Weighted metric | Lambda | Triangle | Gini | P90 rest. | Min avg CRI | P90 access | Avg share access |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Gini restoration inequality | 0.25 | 504.346 | 0.167 | 509.304 | 0.178 | 509.304 | 0.280 |
| Gini restoration inequality | 0.75 | 801.725 | 0.142 | 746.438 | 0.123 | 746.438 | 0.210 |
| P90 CRI restoration time | 0.25 | 641.690 | 0.176 | 633.482 | 0.170 | 633.482 | 0.273 |
| P90 CRI restoration time | 0.75 | 641.690 | 0.176 | 633.482 | 0.170 | 633.482 | 0.273 |
| Maximin time-avg CRI | 0.25 | 641.690 | 0.176 | 633.482 | 0.170 | 633.482 | 0.273 |
| Maximin time-avg CRI | 0.75 | 543.247 | 0.220 | 554.560 | 0.171 | 554.560 | 0.315 |
| P90 critical access time | 0.25 | 641.690 | 0.176 | 633.482 | 0.170 | 633.482 | 0.273 |
| P90 critical access time | 0.75 | 641.690 | 0.176 | 633.482 | 0.170 | 633.482 | 0.273 |

Table 8 shows that modest changes around the original lambda = 0.5 setting can still change the selected sequence, but the effect depends strongly on which fairness metric is weighted. The lowest triangle among the weighted sensitivity runs is obtained by the Gini-weighted lambda = 0.25 case (triangle = 504.346), while the lowest Gini is obtained by the Gini-weighted lambda = 0.75 case (Gini = 0.142). This pattern is consistent with the intended trade-off: a stronger Gini weight improves restoration-time equality but accepts a larger aggregate triangle loss.

The P90 CRI restoration and P90 critical-access objectives converge to the same sequence at both lambda values in this representative case. This reinforces the earlier interpretation that the P90 restoration result is largely driven by accessibility under the current CRI weights. The maximin-weighted runs show a different behavior: lambda = 0.75 improves time-averaged share access and lowers the triangle relative to the lambda = 0.25 case, but it does not substantially improve the minimum time-averaged CRI. As with the other SA-based comparisons, these results should be interpreted as representative-case search outcomes rather than theoretical dominance relationships.

## Interpretation

Overall, the representative case supports the same qualitative conclusion as the larger trade-off analysis: efficiency-oriented objectives reduce aggregate loss quickly, while equity-oriented objectives change which zones are served earlier and how evenly recovery times are distributed. The sensitivity analysis adds two refinements. First, threshold-based fairness metrics should be reported with their threshold assumptions because their values can move noticeably between 0.85 and 0.95. Second, weighted equity objectives should not rely only on Gini; Gini, P90 restoration, and maximin CRI emphasize different ethical claims.

For reporting, the most defensible weighted comparison is therefore not a single universal weighted solution, but a small set of weighted variants. A Gini-weighted run is appropriate when restoration-time inequality is the focus, a P90-weighted run is appropriate when the planner wants most zones to cross the CRI adequacy threshold quickly, and a maximin run is appropriate when the worst-served zone is the central equity concern.
