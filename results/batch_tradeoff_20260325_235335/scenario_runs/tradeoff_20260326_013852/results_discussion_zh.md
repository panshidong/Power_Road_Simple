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
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 1629.171 | +0.0 | 0.182 | +0.0 | 0.371 | 1832.644 | 1832.644 | 0.326 | 否 |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 728.627 | -55.3 | 0.171 | -6.1 | 0.354 | 829.916 | 829.916 | 0.353 | 是 |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 2113.159 | +29.7 | 0.118 | -35.3 | 0.320 | 2080.858 | 2080.858 | 0.296 | 是 |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 1201.749 | -26.2 | 0.217 | +18.9 | 0.530 | 1335.604 | 1335.604 | 0.490 | 否 |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 736.947 | -54.8 | 0.149 | -18.4 | 0.335 | 799.189 | 799.189 | 0.302 | 是 |
| weighted_gini_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:gini_restore. | 1390.842 | -14.6 | 0.144 | -20.7 | 0.178 | 1245.665 | 1245.665 | 0.290 | 是 |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 1948.886 | +19.6 | 0.170 | -6.6 | 0.272 | 2072.995 | 2072.995 | 0.298 | 否 |
| weighted_maximin_time_avg_cri_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:maximin_time_avg_cri_loss. | 1093.926 | -32.9 | 0.190 | +4.5 | 0.342 | 1237.340 | 1237.340 | 0.346 | 否 |

## 表 2 代表性方案的敏感性范围

| 代表方案 | Gini restore 范围 | Min time-avg CRI 范围 | P90 restore 范围 | P90 access 范围 | Time-avg share access 范围 |
|---|---:|---:|---:|---:|---:|

## 结果讨论

- 效率最优方案是 `single_triangle`，Triangle 为 `728.627`，较基线改善 `-55.3%`。
- 公平性最优方案是 `single_gini_restore`，Gini restore 为 `0.118`，较基线改善 `-35.3%`。
- 新增的 maximin 视角下，`single_maximin_time_avg_cri` 取得最高的最小时间平均 CRI，`min_time_avg_cri=0.530`。
- `Triangle-Gini` 平面上的 Pareto 解为：`['single_triangle', 'single_gini_restore', 'single_p90_access_restore', 'weighted_gini_restore_l050']`。
- `share_access_final` 在当前假设下最终都接近完全恢复，更能区分方案的是 `p90_access_restore` 与 `time_avg_share_access`。基线的 `time_avg_share_access=0.326`。
- 如果后续需要进一步拉开 access 类目标，建议优先调整 critical locations 的布点或扩大受损资产集合；如果要加强搜索，则可以继续增加 SA 迭代次数。

## 建议优先写进正文的结果

1. `single_triangle`：代表效率优先。
2. `single_gini_restore`：代表公平性优先。
3. `single_maximin_time_avg_cri`：代表新增的 maximin equity 规则。
