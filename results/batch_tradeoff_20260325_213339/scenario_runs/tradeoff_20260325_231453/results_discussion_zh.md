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
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 2983.244 | +0.0 | 0.249 | +0.0 | 0.368 | 3149.939 | 3149.939 | 0.294 | 否 |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 2116.801 | -29.0 | 0.250 | +0.3 | 0.372 | 2158.594 | 2144.594 | 0.321 | 是 |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 2368.060 | -20.6 | 0.177 | -29.1 | 0.228 | 1899.009 | 1899.009 | 0.269 | 是 |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 1345.712 | -54.9 | 0.285 | +14.2 | 0.468 | 1592.465 | 1549.382 | 0.374 | 是 |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 2406.936 | -19.3 | 0.199 | -20.1 | 0.249 | 2052.659 | 2052.659 | 0.269 | 否 |
| weighted_gini_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:gini_restore. | 3038.290 | +1.8 | 0.223 | -10.7 | 0.272 | 2722.138 | 2722.138 | 0.284 | 否 |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 3041.401 | +1.9 | 0.236 | -5.6 | 0.369 | 3157.263 | 3143.263 | 0.303 | 否 |
| weighted_maximin_time_avg_cri_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:maximin_time_avg_cri_loss. | 2165.037 | -27.4 | 0.233 | -6.7 | 0.360 | 2158.594 | 2144.594 | 0.307 | 是 |

## 表 2 代表性方案的敏感性范围

| 代表方案 | Gini restore 范围 | Min time-avg CRI 范围 | P90 restore 范围 | P90 access 范围 | Time-avg share access 范围 |
|---|---:|---:|---:|---:|---:|

## 结果讨论

- 效率最优方案是 `single_maximin_time_avg_cri`，Triangle 为 `1345.712`，较基线改善 `-54.9%`。
- 公平性最优方案是 `single_gini_restore`，Gini restore 为 `0.177`，较基线改善 `-29.1%`。
- 新增的 maximin 视角下，`single_maximin_time_avg_cri` 取得最高的最小时间平均 CRI，`min_time_avg_cri=0.468`。
- `Triangle-Gini` 平面上的 Pareto 解为：`['single_triangle', 'single_gini_restore', 'single_maximin_time_avg_cri', 'weighted_maximin_time_avg_cri_l050']`。
- `share_access_final` 在当前假设下最终都接近完全恢复，更能区分方案的是 `p90_access_restore` 与 `time_avg_share_access`。基线的 `time_avg_share_access=0.294`。
- 如果后续需要进一步拉开 access 类目标，建议优先调整 critical locations 的布点或扩大受损资产集合；如果要加强搜索，则可以继续增加 SA 迭代次数。

## 建议优先写进正文的结果

1. `single_maximin_time_avg_cri`：代表效率优先。
2. `single_gini_restore`：代表公平性优先。
3. `single_maximin_time_avg_cri`：代表新增的 maximin equity 规则。
