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
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 8139.514 | +0.0 | 0.134 | +0.0 | 0.224 | 6944.065 | 6944.065 | 0.283 | 否 |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 3257.384 | -60.0 | 0.185 | +37.4 | 0.415 | 3591.496 | 3591.496 | 0.335 | 是 |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 6619.940 | -18.7 | 0.093 | -31.0 | 0.261 | 5993.676 | 5993.676 | 0.280 | 是 |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 5272.196 | -35.2 | 0.164 | +22.1 | 0.384 | 5540.559 | 5540.559 | 0.310 | 否 |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 4477.878 | -45.0 | 0.139 | +3.2 | 0.201 | 3690.602 | 3435.684 | 0.340 | 否 |
| weighted_gini_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:gini_restore. | 6201.619 | -23.8 | 0.124 | -7.5 | 0.349 | 6117.453 | 6117.453 | 0.335 | 否 |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 3633.296 | -55.4 | 0.121 | -9.9 | 0.262 | 3282.972 | 3282.972 | 0.256 | 是 |
| weighted_maximin_time_avg_cri_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:maximin_time_avg_cri_loss. | 5021.748 | -38.3 | 0.226 | +68.5 | 0.397 | 5666.670 | 5666.670 | 0.332 | 否 |

## 表 2 代表性方案的敏感性范围

| 代表方案 | Gini restore 范围 | Min time-avg CRI 范围 | P90 restore 范围 | P90 access 范围 | Time-avg share access 范围 |
|---|---:|---:|---:|---:|---:|

## 结果讨论

- 效率最优方案是 `single_triangle`，Triangle 为 `3257.384`，较基线改善 `-60.0%`。
- 公平性最优方案是 `single_gini_restore`，Gini restore 为 `0.093`，较基线改善 `-31.0%`。
- 新增的 maximin 视角下，`single_triangle` 取得最高的最小时间平均 CRI，`min_time_avg_cri=0.415`。
- `Triangle-Gini` 平面上的 Pareto 解为：`['single_triangle', 'single_gini_restore', 'guardrail_gini_restore']`。
- `share_access_final` 在当前假设下最终都接近完全恢复，更能区分方案的是 `p90_access_restore` 与 `time_avg_share_access`。基线的 `time_avg_share_access=0.283`。
- 如果后续需要进一步拉开 access 类目标，建议优先调整 critical locations 的布点或扩大受损资产集合；如果要加强搜索，则可以继续增加 SA 迭代次数。

## 建议优先写进正文的结果

1. `single_triangle`：代表效率优先。
2. `single_gini_restore`：代表公平性优先。
3. `single_triangle`：代表新增的 maximin equity 规则。
