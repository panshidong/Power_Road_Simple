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
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 2270.782 | +0.0 | 0.163 | +0.0 | 0.304 | 2149.919 | 2149.919 | 0.285 | 否 |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 1484.289 | -34.6 | 0.205 | +26.4 | 0.401 | 1503.067 | 1503.067 | 0.293 | 否 |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 2325.409 | +2.4 | 0.151 | -7.0 | 0.321 | 2149.713 | 2149.713 | 0.296 | 否 |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 2043.156 | -10.0 | 0.243 | +49.6 | 0.442 | 2189.220 | 2189.220 | 0.352 | 否 |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 1430.925 | -37.0 | 0.153 | -5.8 | 0.320 | 1235.620 | 1235.620 | 0.313 | 否 |
| weighted_gini_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:gini_restore. | 1301.931 | -42.7 | 0.168 | +3.3 | 0.362 | 1414.838 | 1414.838 | 0.326 | 否 |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 1335.143 | -41.2 | 0.151 | -7.1 | 0.324 | 1350.255 | 1350.255 | 0.276 | 是 |
| weighted_maximin_time_avg_cri_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:maximin_time_avg_cri_loss. | 1066.265 | -53.0 | 0.167 | +2.5 | 0.320 | 1150.276 | 1130.026 | 0.316 | 是 |

## 表 2 代表性方案的敏感性范围

| 代表方案 | Gini restore 范围 | Min time-avg CRI 范围 | P90 restore 范围 | P90 access 范围 | Time-avg share access 范围 |
|---|---:|---:|---:|---:|---:|

## 结果讨论

- 效率最优方案是 `weighted_maximin_time_avg_cri_l050`，Triangle 为 `1066.265`，较基线改善 `-53.0%`。
- 公平性最优方案是 `guardrail_gini_restore`，Gini restore 为 `0.151`，较基线改善 `-7.1%`。
- 新增的 maximin 视角下，`single_maximin_time_avg_cri` 取得最高的最小时间平均 CRI，`min_time_avg_cri=0.442`。
- `Triangle-Gini` 平面上的 Pareto 解为：`['guardrail_gini_restore', 'weighted_maximin_time_avg_cri_l050']`。
- `share_access_final` 在当前假设下最终都接近完全恢复，更能区分方案的是 `p90_access_restore` 与 `time_avg_share_access`。基线的 `time_avg_share_access=0.285`。
- 如果后续需要进一步拉开 access 类目标，建议优先调整 critical locations 的布点或扩大受损资产集合；如果要加强搜索，则可以继续增加 SA 迭代次数。

## 建议优先写进正文的结果

1. `weighted_maximin_time_avg_cri_l050`：代表效率优先。
2. `guardrail_gini_restore`：代表公平性优先。
3. `single_maximin_time_avg_cri`：代表新增的 maximin equity 规则。
