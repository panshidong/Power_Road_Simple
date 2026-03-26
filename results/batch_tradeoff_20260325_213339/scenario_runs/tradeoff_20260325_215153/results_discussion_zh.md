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
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 1617.261 | -0.7 | 0.182 | +0.2 | 0.376 | 1832.644 | 1832.644 | 0.328 | 否 |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 2275.394 | +39.7 | 0.164 | -9.8 | 0.281 | 2387.905 | 2387.905 | 0.309 | 否 |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 2210.381 | +35.7 | 0.186 | +1.9 | 0.435 | 2599.296 | 2599.296 | 0.347 | 否 |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 912.345 | -44.0 | 0.145 | -20.6 | 0.317 | 933.992 | 933.992 | 0.301 | 是 |
| weighted_gini_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:gini_restore. | 1571.153 | -3.6 | 0.179 | -1.8 | 0.386 | 1744.899 | 1744.899 | 0.327 | 否 |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 1963.989 | +20.6 | 0.171 | -6.0 | 0.271 | 2090.635 | 2090.635 | 0.300 | 否 |
| weighted_maximin_time_avg_cri_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:maximin_time_avg_cri_loss. | 1093.926 | -32.9 | 0.190 | +4.5 | 0.342 | 1237.340 | 1237.340 | 0.346 | 否 |

## 表 2 代表性方案的敏感性范围

| 代表方案 | Gini restore 范围 | Min time-avg CRI 范围 | P90 restore 范围 | P90 access 范围 | Time-avg share access 范围 |
|---|---:|---:|---:|---:|---:|

## 结果讨论

- 效率最优方案是 `single_p90_access_restore`，Triangle 为 `912.345`，较基线改善 `-44.0%`。
- 公平性最优方案是 `single_p90_access_restore`，Gini restore 为 `0.145`，较基线改善 `-20.6%`。
- 新增的 maximin 视角下，`single_maximin_time_avg_cri` 取得最高的最小时间平均 CRI，`min_time_avg_cri=0.435`。
- `Triangle-Gini` 平面上的 Pareto 解为：`['single_p90_access_restore']`。
- `share_access_final` 在当前假设下最终都接近完全恢复，更能区分方案的是 `p90_access_restore` 与 `time_avg_share_access`。基线的 `time_avg_share_access=0.326`。
- 如果后续需要进一步拉开 access 类目标，建议优先调整 critical locations 的布点或扩大受损资产集合；如果要加强搜索，则可以继续增加 SA 迭代次数。

## 建议优先写进正文的结果

1. `single_p90_access_restore`：代表效率优先。
2. `single_p90_access_restore`：代表公平性优先。
3. `single_maximin_time_avg_cri`：代表新增的 maximin equity 规则。
