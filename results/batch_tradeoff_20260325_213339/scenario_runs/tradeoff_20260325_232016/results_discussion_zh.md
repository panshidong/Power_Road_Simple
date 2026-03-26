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
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 1273.470 | +0.0 | 0.512 | +0.0 | 0.545 | 1868.471 | 2002.814 | 0.632 | 否 |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 1191.203 | -6.5 | 0.471 | -8.1 | 0.568 | 1784.791 | 1919.133 | 0.584 | 否 |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 1267.570 | -0.5 | 0.423 | -17.5 | 0.563 | 1933.143 | 1933.143 | 0.562 | 否 |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 1269.059 | -0.3 | 0.449 | -12.3 | 0.649 | 2197.992 | 2197.992 | 0.597 | 否 |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 1152.537 | -9.5 | 0.510 | -0.5 | 0.564 | 1793.239 | 1830.634 | 0.651 | 否 |
| weighted_gini_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:gini_restore. | 1185.991 | -6.9 | 0.445 | -13.2 | 0.556 | 1782.029 | 1916.371 | 0.574 | 否 |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 1191.203 | -6.5 | 0.471 | -8.1 | 0.568 | 1784.791 | 1919.133 | 0.584 | 否 |
| weighted_maximin_time_avg_cri_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:maximin_time_avg_cri_loss. | 1005.608 | -21.0 | 0.410 | -19.9 | 0.526 | 1084.143 | 1084.143 | 0.584 | 是 |

## 表 2 代表性方案的敏感性范围

| 代表方案 | Gini restore 范围 | Min time-avg CRI 范围 | P90 restore 范围 | P90 access 范围 | Time-avg share access 范围 |
|---|---:|---:|---:|---:|---:|

## 结果讨论

- 效率最优方案是 `weighted_maximin_time_avg_cri_l050`，Triangle 为 `1005.608`，较基线改善 `-21.0%`。
- 公平性最优方案是 `weighted_maximin_time_avg_cri_l050`，Gini restore 为 `0.410`，较基线改善 `-19.9%`。
- 新增的 maximin 视角下，`single_maximin_time_avg_cri` 取得最高的最小时间平均 CRI，`min_time_avg_cri=0.649`。
- `Triangle-Gini` 平面上的 Pareto 解为：`['weighted_maximin_time_avg_cri_l050']`。
- `share_access_final` 在当前假设下最终都接近完全恢复，更能区分方案的是 `p90_access_restore` 与 `time_avg_share_access`。基线的 `time_avg_share_access=0.632`。
- 如果后续需要进一步拉开 access 类目标，建议优先调整 critical locations 的布点或扩大受损资产集合；如果要加强搜索，则可以继续增加 SA 迭代次数。

## 建议优先写进正文的结果

1. `weighted_maximin_time_avg_cri_l050`：代表效率优先。
2. `weighted_maximin_time_avg_cri_l050`：代表公平性优先。
3. `single_maximin_time_avg_cri`：代表新增的 maximin equity 规则。
