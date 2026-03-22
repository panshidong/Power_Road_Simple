# Task B 正式结果表与讨论

## 数据来源
- 主结果文件：`tradeoff_summary.csv`
- 图：`tradeoff_scatter.png`
- 主设定：`w_e=0.5, w_a=0.5, cri_threshold=0.9, critical_access_threshold=0.9`
- 调度映射：`link_only`，使用 `new_bus_to_link.json`

## 表 1 主实验详细结果

| 实验 | 类型 | 描述 | Triangle | ΔTri % | Gini | ΔGini % | Min final CRI | P90 restore | P90 access | Time-avg share access | Pareto |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 310.830 | +0.0 | 0.351 | +0.0 | 1.000 | 521.950 | 629.845 | 0.481 | 是 |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 305.779 | -1.6 | 0.381 | +8.5 | 1.000 | 617.807 | 617.807 | 0.405 | 是 |
| single_var_restore | single | Equity objective that reduces dispersion in restoration times. | 310.830 | +0.0 | 0.351 | +0.0 | 1.000 | 521.950 | 629.845 | 0.481 | 是 |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 310.830 | +0.0 | 0.351 | +0.0 | 1.000 | 521.950 | 629.845 | 0.481 | 是 |
| single_p90_restore | single | Equity objective that shortens the time by which 90% of zones recover. | 310.830 | +0.0 | 0.351 | +0.0 | 1.000 | 521.950 | 629.845 | 0.481 | 是 |
| single_maximin_final_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum final CRI. | 310.830 | +0.0 | 0.351 | +0.0 | 1.000 | 521.950 | 629.845 | 0.481 | 是 |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 305.779 | -1.6 | 0.381 | +8.5 | 1.000 | 617.807 | 617.807 | 0.405 | 是 |
| weighted_gini_restore_l025 | weighted_sum | Weighted trade-off objective. triangle + 0.25 * equity:gini_restore. | 310.830 | +0.0 | 0.351 | +0.0 | 1.000 | 521.950 | 629.845 | 0.481 | 是 |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 310.830 | +0.0 | 0.351 | +0.0 | 1.000 | 521.950 | 629.845 | 0.481 | 是 |
| weighted_p90_restore_l025 | weighted_sum | Weighted trade-off objective. triangle + 0.25 * equity:p90_restore. | 310.830 | +0.0 | 0.351 | +0.0 | 1.000 | 521.950 | 629.845 | 0.481 | 是 |
| guardrail_p90_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:p90_restore. | 310.830 | +0.0 | 0.351 | +0.0 | 1.000 | 521.950 | 629.845 | 0.481 | 是 |
| weighted_maximin_final_cri_l025 | weighted_sum | Equity via maximin: improve the worst-served zone by maximizing the minimum final CRI. | 305.779 | -1.6 | 0.381 | +8.5 | 1.000 | 617.807 | 617.807 | 0.405 | 是 |
| guardrail_maximin_final_cri | guardrail | Equity via maximin: improve the worst-served zone by maximizing the minimum final CRI. | 305.779 | -1.6 | 0.381 | +8.5 | 1.000 | 617.807 | 617.807 | 0.405 | 是 |
| weighted_p90_access_restore_l025 | weighted_sum | Weighted trade-off objective. triangle + 0.25 * critical_access:p90_access_restore. | 305.779 | -1.6 | 0.381 | +8.5 | 1.000 | 617.807 | 617.807 | 0.405 | 是 |
| guardrail_p90_access_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on critical_access:p90_access_restore. | 305.779 | -1.6 | 0.381 | +8.5 | 1.000 | 617.807 | 617.807 | 0.405 | 是 |

## 表 2 代表性方案的敏感性范围

| 代表方案 | Gini restore 范围 | Min final CRI 范围 | P90 restore 范围 | P90 access 范围 | Time-avg share access 范围 |
|---|---:|---:|---:|---:|---:|
| baseline_reference | 0.351 - 0.351 | 1.000 - 1.000 | 521.950 - 521.950 | 629.845 - 629.845 | 0.481 - 0.481 |
| single_triangle | 0.381 - 0.381 | 1.000 - 1.000 | 617.807 - 617.807 | 617.807 - 617.807 | 0.405 - 0.405 |

## 结果讨论

- 效率最优方案是 `single_triangle`，Triangle 为 `305.779`，较基线改善 `-1.6%`。
- 公平性最优方案是 `baseline_reference`，Gini restore 为 `0.351`，较基线改善 `0.0%`。
- 新增的 maximin 视角下，`baseline_reference` 取得最高的最小最终 CRI，`min_final_cri=1.000`。
- `Triangle-Gini` 平面上的 Pareto 解为：`['baseline_reference', 'single_triangle', 'single_var_restore', 'single_gini_restore', 'single_p90_restore', 'single_maximin_final_cri', 'single_p90_access_restore', 'weighted_gini_restore_l025', 'guardrail_gini_restore', 'weighted_p90_restore_l025', 'guardrail_p90_restore', 'weighted_maximin_final_cri_l025', 'guardrail_maximin_final_cri', 'weighted_p90_access_restore_l025', 'guardrail_p90_access_restore']`。
- `share_access_final` 在当前假设下最终都接近完全恢复，更能区分方案的是 `p90_access_restore` 与 `time_avg_share_access`。基线的 `time_avg_share_access=0.481`。
- 如果后续需要进一步拉开 access 类目标，建议优先调整 critical locations 的布点或扩大受损资产集合；如果要加强搜索，则可以继续增加 SA 迭代次数。

## 建议优先写进正文的结果

1. `single_triangle`：代表效率优先。
2. `baseline_reference`：代表公平性优先。
3. `baseline_reference`：代表新增的 maximin equity 规则。
