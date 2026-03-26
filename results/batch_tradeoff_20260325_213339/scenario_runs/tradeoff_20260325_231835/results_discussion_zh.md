# Task B 正式结果表与讨论

## 数据来源
- 主结果文件：`tradeoff_summary.csv`
- 图：`tradeoff_scatter.png`
- 主设定：`w_e=0.133, w_a=0.867, cri_threshold=0.9, critical_access_threshold=0.9`
- 队伍配置：`specialized`，power crews=`1`，road crews=`1`
- 调度映射：`link_only`，使用 `new_bus_to_link.json`
- 启用实验：`['single_triangle', 'single_gini_restore', 'single_maximin_time_avg_cri', 'single_p90_access_restore', 'weighted_gini_restore_l050', 'weighted_maximin_time_avg_cri_l050', 'guardrail_gini_restore']`

## 表 1 主实验详细结果

| 实验 | 类型 | 描述 | Triangle | ΔTri % | Gini | ΔGini % | Min time-avg CRI | P90 restore | P90 access | Time-avg share access | Pareto |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 40278.848 | +0.0 | 0.099 | +0.0 | 0.144 | 36019.379 | 36019.379 | 0.173 | 是 |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 34747.070 | -13.7 | 0.099 | +0.0 | 0.293 | 33919.636 | 28980.655 | 0.324 | 是 |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 43261.832 | +7.4 | 0.058 | -41.1 | 0.070 | 36072.198 | 36072.198 | 0.103 | 是 |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 37391.157 | -7.2 | 0.123 | +24.1 | 0.211 | 37238.779 | 37238.779 | 0.224 | 否 |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 35183.118 | -12.7 | 0.116 | +16.5 | 0.190 | 34342.687 | 34342.687 | 0.194 | 否 |
| weighted_gini_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:gini_restore. | 45348.801 | +12.6 | 0.059 | -40.6 | 0.170 | 36803.661 | 36803.661 | 0.217 | 否 |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 40337.291 | +0.1 | 0.068 | -31.6 | 0.095 | 34587.932 | 34587.932 | 0.117 | 是 |
| weighted_maximin_time_avg_cri_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:maximin_time_avg_cri_loss. | 35327.114 | -12.3 | 0.115 | +15.5 | 0.201 | 34342.687 | 34342.687 | 0.192 | 否 |

## 表 2 代表性方案的敏感性范围

| 代表方案 | Gini restore 范围 | Min time-avg CRI 范围 | P90 restore 范围 | P90 access 范围 | Time-avg share access 范围 |
|---|---:|---:|---:|---:|---:|

## 结果讨论

- 效率最优方案是 `single_triangle`，Triangle 为 `34747.070`，较基线改善 `-13.7%`。
- 公平性最优方案是 `single_gini_restore`，Gini restore 为 `0.058`，较基线改善 `-41.1%`。
- 新增的 maximin 视角下，`single_triangle` 取得最高的最小时间平均 CRI，`min_time_avg_cri=0.293`。
- `Triangle-Gini` 平面上的 Pareto 解为：`['baseline_reference', 'single_triangle', 'single_gini_restore', 'guardrail_gini_restore']`。
- `share_access_final` 在当前假设下最终都接近完全恢复，更能区分方案的是 `p90_access_restore` 与 `time_avg_share_access`。基线的 `time_avg_share_access=0.173`。
- 如果后续需要进一步拉开 access 类目标，建议优先调整 critical locations 的布点或扩大受损资产集合；如果要加强搜索，则可以继续增加 SA 迭代次数。

## 建议优先写进正文的结果

1. `single_triangle`：代表效率优先。
2. `single_gini_restore`：代表公平性优先。
3. `single_triangle`：代表新增的 maximin equity 规则。
