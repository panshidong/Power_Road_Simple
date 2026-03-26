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
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 4424.499 | +0.0 | 0.172 | +0.0 | 0.299 | 4232.422 | 4100.951 | 0.316 | 否 |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 2969.852 | -32.9 | 0.255 | +48.3 | 0.375 | 2934.495 | 2934.495 | 0.392 | 否 |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 2699.129 | -39.0 | 0.112 | -34.8 | 0.250 | 2395.133 | 2363.701 | 0.194 | 是 |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 3140.652 | -29.0 | 0.321 | +86.7 | 0.368 | 3313.743 | 3121.181 | 0.393 | 否 |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 3304.829 | -25.3 | 0.141 | -18.3 | 0.358 | 3321.471 | 3190.001 | 0.276 | 否 |
| weighted_gini_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:gini_restore. | 3181.799 | -28.1 | 0.180 | +5.0 | 0.324 | 3066.667 | 3066.667 | 0.270 | 否 |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 3015.413 | -31.8 | 0.157 | -8.9 | 0.385 | 2986.287 | 2986.287 | 0.317 | 否 |
| weighted_maximin_time_avg_cri_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:maximin_time_avg_cri_loss. | 2996.590 | -32.3 | 0.309 | +79.9 | 0.336 | 2948.273 | 2948.273 | 0.373 | 否 |

## 表 2 代表性方案的敏感性范围

| 代表方案 | Gini restore 范围 | Min time-avg CRI 范围 | P90 restore 范围 | P90 access 范围 | Time-avg share access 范围 |
|---|---:|---:|---:|---:|---:|

## 结果讨论

- 效率最优方案是 `single_gini_restore`，Triangle 为 `2699.129`，较基线改善 `-39.0%`。
- 公平性最优方案是 `single_gini_restore`，Gini restore 为 `0.112`，较基线改善 `-34.8%`。
- 新增的 maximin 视角下，`guardrail_gini_restore` 取得最高的最小时间平均 CRI，`min_time_avg_cri=0.385`。
- `Triangle-Gini` 平面上的 Pareto 解为：`['single_gini_restore']`。
- `share_access_final` 在当前假设下最终都接近完全恢复，更能区分方案的是 `p90_access_restore` 与 `time_avg_share_access`。基线的 `time_avg_share_access=0.316`。
- 如果后续需要进一步拉开 access 类目标，建议优先调整 critical locations 的布点或扩大受损资产集合；如果要加强搜索，则可以继续增加 SA 迭代次数。

## 建议优先写进正文的结果

1. `single_gini_restore`：代表效率优先。
2. `single_gini_restore`：代表公平性优先。
3. `guardrail_gini_restore`：代表新增的 maximin equity 规则。
