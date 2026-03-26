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
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 2279.451 | +0.0 | 0.183 | +0.0 | 0.265 | 2299.985 | 2256.435 | 0.313 | 否 |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 2088.582 | -8.4 | 0.166 | -9.4 | 0.313 | 2177.078 | 2133.529 | 0.320 | 是 |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 2108.889 | -7.5 | 0.166 | -9.2 | 0.310 | 2144.533 | 2100.984 | 0.320 | 否 |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 1561.099 | -31.5 | 0.201 | +9.4 | 0.402 | 1715.797 | 1497.634 | 0.409 | 是 |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 1780.116 | -21.9 | 0.195 | +6.6 | 0.304 | 1877.892 | 1834.343 | 0.314 | 是 |
| weighted_gini_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:gini_restore. | 1791.659 | -21.4 | 0.198 | +7.9 | 0.350 | 2018.545 | 1938.446 | 0.365 | 否 |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 2088.582 | -8.4 | 0.166 | -9.4 | 0.313 | 2177.078 | 2133.529 | 0.320 | 是 |
| weighted_maximin_time_avg_cri_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:maximin_time_avg_cri_loss. | 2014.391 | -11.6 | 0.181 | -1.3 | 0.337 | 2017.538 | 1973.989 | 0.322 | 是 |

## 表 2 代表性方案的敏感性范围

| 代表方案 | Gini restore 范围 | Min time-avg CRI 范围 | P90 restore 范围 | P90 access 范围 | Time-avg share access 范围 |
|---|---:|---:|---:|---:|---:|

## 结果讨论

- 效率最优方案是 `single_maximin_time_avg_cri`，Triangle 为 `1561.099`，较基线改善 `-31.5%`。
- 公平性最优方案是 `single_triangle`，Gini restore 为 `0.166`，较基线改善 `-9.4%`。
- 新增的 maximin 视角下，`single_maximin_time_avg_cri` 取得最高的最小时间平均 CRI，`min_time_avg_cri=0.402`。
- `Triangle-Gini` 平面上的 Pareto 解为：`['single_triangle', 'single_maximin_time_avg_cri', 'single_p90_access_restore', 'guardrail_gini_restore', 'weighted_maximin_time_avg_cri_l050']`。
- `share_access_final` 在当前假设下最终都接近完全恢复，更能区分方案的是 `p90_access_restore` 与 `time_avg_share_access`。基线的 `time_avg_share_access=0.313`。
- 如果后续需要进一步拉开 access 类目标，建议优先调整 critical locations 的布点或扩大受损资产集合；如果要加强搜索，则可以继续增加 SA 迭代次数。

## 建议优先写进正文的结果

1. `single_maximin_time_avg_cri`：代表效率优先。
2. `single_triangle`：代表公平性优先。
3. `single_maximin_time_avg_cri`：代表新增的 maximin equity 规则。
