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
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 3259.153 | +0.0 | 0.160 | +0.0 | 0.248 | 3558.751 | 3558.751 | 0.240 | 否 |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 1728.357 | -47.0 | 0.190 | +18.4 | 0.383 | 2169.558 | 2169.558 | 0.274 | 是 |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 3834.995 | +17.7 | 0.108 | -32.5 | 0.282 | 3713.685 | 3645.613 | 0.250 | 是 |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 2972.022 | -8.8 | 0.169 | +5.8 | 0.378 | 3284.891 | 3284.891 | 0.367 | 否 |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 2361.834 | -27.5 | 0.262 | +63.7 | 0.490 | 2426.472 | 2426.472 | 0.494 | 否 |
| weighted_gini_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:gini_restore. | 2166.405 | -33.5 | 0.132 | -17.5 | 0.352 | 2751.439 | 2812.890 | 0.298 | 是 |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 2758.298 | -15.4 | 0.143 | -10.7 | 0.234 | 3054.014 | 3054.014 | 0.311 | 否 |
| weighted_maximin_time_avg_cri_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:maximin_time_avg_cri_loss. | 2668.300 | -18.1 | 0.198 | +23.9 | 0.428 | 2981.888 | 2981.888 | 0.369 | 否 |

## 表 2 代表性方案的敏感性范围

| 代表方案 | Gini restore 范围 | Min time-avg CRI 范围 | P90 restore 范围 | P90 access 范围 | Time-avg share access 范围 |
|---|---:|---:|---:|---:|---:|

## 结果讨论

- 效率最优方案是 `single_triangle`，Triangle 为 `1728.357`，较基线改善 `-47.0%`。
- 公平性最优方案是 `single_gini_restore`，Gini restore 为 `0.108`，较基线改善 `-32.5%`。
- 新增的 maximin 视角下，`single_p90_access_restore` 取得最高的最小时间平均 CRI，`min_time_avg_cri=0.490`。
- `Triangle-Gini` 平面上的 Pareto 解为：`['single_triangle', 'single_gini_restore', 'weighted_gini_restore_l050']`。
- `share_access_final` 在当前假设下最终都接近完全恢复，更能区分方案的是 `p90_access_restore` 与 `time_avg_share_access`。基线的 `time_avg_share_access=0.240`。
- 如果后续需要进一步拉开 access 类目标，建议优先调整 critical locations 的布点或扩大受损资产集合；如果要加强搜索，则可以继续增加 SA 迭代次数。

## 建议优先写进正文的结果

1. `single_triangle`：代表效率优先。
2. `single_gini_restore`：代表公平性优先。
3. `single_p90_access_restore`：代表新增的 maximin equity 规则。
