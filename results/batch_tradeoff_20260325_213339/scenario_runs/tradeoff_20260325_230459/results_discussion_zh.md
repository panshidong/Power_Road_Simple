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
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 2119.795 | +0.0 | 0.249 | +0.0 | 0.435 | 2469.551 | 2469.551 | 0.331 | 是 |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 1316.061 | -37.9 | 0.273 | +9.5 | 0.540 | 2166.924 | 2166.924 | 0.335 | 是 |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 2391.670 | +12.8 | 0.228 | -8.4 | 0.362 | 2481.160 | 2481.160 | 0.328 | 是 |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 1570.230 | -25.9 | 0.269 | +8.0 | 0.556 | 2548.166 | 2548.166 | 0.348 | 否 |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 1543.494 | -27.2 | 0.274 | +10.1 | 0.507 | 2150.025 | 2096.345 | 0.364 | 否 |
| weighted_gini_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:gini_restore. | 1545.180 | -27.1 | 0.261 | +4.9 | 0.503 | 2166.924 | 2166.924 | 0.336 | 是 |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 2209.893 | +4.3 | 0.230 | -7.4 | 0.407 | 2464.009 | 2464.009 | 0.326 | 是 |
| weighted_maximin_time_avg_cri_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:maximin_time_avg_cri_loss. | 1332.835 | -37.1 | 0.285 | +14.6 | 0.536 | 2150.025 | 2096.345 | 0.363 | 否 |

## 表 2 代表性方案的敏感性范围

| 代表方案 | Gini restore 范围 | Min time-avg CRI 范围 | P90 restore 范围 | P90 access 范围 | Time-avg share access 范围 |
|---|---:|---:|---:|---:|---:|

## 结果讨论

- 效率最优方案是 `single_triangle`，Triangle 为 `1316.061`，较基线改善 `-37.9%`。
- 公平性最优方案是 `single_gini_restore`，Gini restore 为 `0.228`，较基线改善 `-8.4%`。
- 新增的 maximin 视角下，`single_maximin_time_avg_cri` 取得最高的最小时间平均 CRI，`min_time_avg_cri=0.556`。
- `Triangle-Gini` 平面上的 Pareto 解为：`['baseline_reference', 'single_triangle', 'single_gini_restore', 'weighted_gini_restore_l050', 'guardrail_gini_restore']`。
- `share_access_final` 在当前假设下最终都接近完全恢复，更能区分方案的是 `p90_access_restore` 与 `time_avg_share_access`。基线的 `time_avg_share_access=0.331`。
- 如果后续需要进一步拉开 access 类目标，建议优先调整 critical locations 的布点或扩大受损资产集合；如果要加强搜索，则可以继续增加 SA 迭代次数。

## 建议优先写进正文的结果

1. `single_triangle`：代表效率优先。
2. `single_gini_restore`：代表公平性优先。
3. `single_maximin_time_avg_cri`：代表新增的 maximin equity 规则。
