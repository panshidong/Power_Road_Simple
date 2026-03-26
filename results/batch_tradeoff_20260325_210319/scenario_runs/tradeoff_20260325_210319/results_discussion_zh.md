# Task B 正式结果表与讨论

## 数据来源
- 主结果文件：`tradeoff_summary.csv`
- 图：`tradeoff_scatter.png`
- 主设定：`w_e=0.133, w_a=0.867, cri_threshold=0.9, critical_access_threshold=0.9`
- 队伍配置：`specialized`，power crews=`1`，road crews=`1`
- 调度映射：`link_only`，使用 `new_bus_to_link.json`

## 表 1 主实验详细结果

| 实验 | 类型 | 描述 | Triangle | ΔTri % | Gini | ΔGini % | Min time-avg CRI | P90 restore | P90 access | Time-avg share access | Pareto |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 3922.024 | +0.0 | 0.165 | +0.0 | 0.248 | 3970.883 | 3970.883 | 0.285 | 否 |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 3847.844 | -1.9 | 0.163 | -1.0 | 0.240 | 3970.883 | 3970.883 | 0.285 | 是 |
| single_var_restore | single | Equity objective that reduces dispersion in restoration times. | 3847.844 | -1.9 | 0.163 | -1.0 | 0.240 | 3970.883 | 3970.883 | 0.285 | 是 |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 3847.844 | -1.9 | 0.163 | -1.0 | 0.240 | 3970.883 | 3970.883 | 0.285 | 是 |
| single_p90_restore | single | Equity objective that shortens the time by which 90% of zones recover. | 3922.024 | +0.0 | 0.165 | +0.0 | 0.248 | 3970.883 | 3970.883 | 0.285 | 否 |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 3922.024 | +0.0 | 0.165 | +0.0 | 0.248 | 3970.883 | 3970.883 | 0.285 | 否 |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 3922.024 | +0.0 | 0.165 | +0.0 | 0.248 | 3970.883 | 3970.883 | 0.285 | 否 |
| weighted_gini_restore_l025 | weighted_sum | Weighted trade-off objective. triangle + 0.25 * equity:gini_restore. | 3847.844 | -1.9 | 0.163 | -1.0 | 0.240 | 3970.883 | 3970.883 | 0.285 | 是 |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 3847.844 | -1.9 | 0.163 | -1.0 | 0.240 | 3970.883 | 3970.883 | 0.285 | 是 |
| weighted_p90_restore_l025 | weighted_sum | Weighted trade-off objective. triangle + 0.25 * equity:p90_restore. | 3847.844 | -1.9 | 0.163 | -1.0 | 0.240 | 3970.883 | 3970.883 | 0.285 | 是 |
| guardrail_p90_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:p90_restore. | 3847.844 | -1.9 | 0.163 | -1.0 | 0.240 | 3970.883 | 3970.883 | 0.285 | 是 |
| weighted_maximin_time_avg_cri_l025 | weighted_sum | Weighted trade-off objective. triangle + 0.25 * equity:maximin_time_avg_cri_loss. | 3847.844 | -1.9 | 0.163 | -1.0 | 0.240 | 3970.883 | 3970.883 | 0.285 | 是 |
| guardrail_maximin_time_avg_cri | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:maximin_time_avg_cri_loss. | 3955.780 | +0.9 | 0.165 | +0.1 | 0.252 | 3970.883 | 3970.883 | 0.283 | 否 |
| weighted_p90_access_restore_l025 | weighted_sum | Weighted trade-off objective. triangle + 0.25 * critical_access:p90_access_restore. | 3847.844 | -1.9 | 0.163 | -1.0 | 0.240 | 3970.883 | 3970.883 | 0.285 | 是 |
| guardrail_p90_access_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on critical_access:p90_access_restore. | 3847.844 | -1.9 | 0.163 | -1.0 | 0.240 | 3970.883 | 3970.883 | 0.285 | 是 |

## 表 2 代表性方案的敏感性范围

| 代表方案 | Gini restore 范围 | Min time-avg CRI 范围 | P90 restore 范围 | P90 access 范围 | Time-avg share access 范围 |
|---|---:|---:|---:|---:|---:|

## 结果讨论

- 效率最优方案是 `single_triangle`，Triangle 为 `3847.844`，较基线改善 `-1.9%`。
- 公平性最优方案是 `single_triangle`，Gini restore 为 `0.163`，较基线改善 `-1.0%`。
- 新增的 maximin 视角下，`guardrail_maximin_time_avg_cri` 取得最高的最小时间平均 CRI，`min_time_avg_cri=0.252`。
- `Triangle-Gini` 平面上的 Pareto 解为：`['single_triangle', 'single_var_restore', 'single_gini_restore', 'weighted_gini_restore_l025', 'guardrail_gini_restore', 'weighted_p90_restore_l025', 'guardrail_p90_restore', 'weighted_maximin_time_avg_cri_l025', 'weighted_p90_access_restore_l025', 'guardrail_p90_access_restore']`。
- `share_access_final` 在当前假设下最终都接近完全恢复，更能区分方案的是 `p90_access_restore` 与 `time_avg_share_access`。基线的 `time_avg_share_access=0.285`。
- 如果后续需要进一步拉开 access 类目标，建议优先调整 critical locations 的布点或扩大受损资产集合；如果要加强搜索，则可以继续增加 SA 迭代次数。

## 建议优先写进正文的结果

1. `single_triangle`：代表效率优先。
2. `single_triangle`：代表公平性优先。
3. `guardrail_maximin_time_avg_cri`：代表新增的 maximin equity 规则。
