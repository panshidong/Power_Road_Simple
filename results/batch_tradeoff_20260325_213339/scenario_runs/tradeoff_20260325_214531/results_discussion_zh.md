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
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 800.852 | +0.0 | 0.205 | +0.0 | 0.186 | 875.309 | 875.309 | 0.273 | 否 |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 763.642 | -4.6 | 0.194 | -5.1 | 0.231 | 867.152 | 867.152 | 0.277 | 是 |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 866.266 | +8.2 | 0.197 | -3.7 | 0.121 | 875.309 | 875.309 | 0.297 | 否 |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 608.883 | -24.0 | 0.247 | +20.9 | 0.322 | 681.269 | 661.269 | 0.350 | 是 |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 825.777 | +3.1 | 0.194 | -5.1 | 0.178 | 867.152 | 867.152 | 0.282 | 否 |
| weighted_gini_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:gini_restore. | 800.852 | +0.0 | 0.205 | +0.0 | 0.186 | 875.309 | 875.309 | 0.273 | 否 |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 1284.708 | +60.4 | 0.189 | -7.5 | 0.208 | 1379.434 | 1379.434 | 0.290 | 是 |
| weighted_maximin_time_avg_cri_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:maximin_time_avg_cri_loss. | 783.027 | -2.2 | 0.217 | +5.8 | 0.205 | 875.309 | 875.309 | 0.286 | 否 |

## 表 2 代表性方案的敏感性范围

| 代表方案 | Gini restore 范围 | Min time-avg CRI 范围 | P90 restore 范围 | P90 access 范围 | Time-avg share access 范围 |
|---|---:|---:|---:|---:|---:|

## 结果讨论

- 效率最优方案是 `single_maximin_time_avg_cri`，Triangle 为 `608.883`，较基线改善 `-24.0%`。
- 公平性最优方案是 `guardrail_gini_restore`，Gini restore 为 `0.189`，较基线改善 `-7.5%`。
- 新增的 maximin 视角下，`single_maximin_time_avg_cri` 取得最高的最小时间平均 CRI，`min_time_avg_cri=0.322`。
- `Triangle-Gini` 平面上的 Pareto 解为：`['single_triangle', 'single_maximin_time_avg_cri', 'guardrail_gini_restore']`。
- `share_access_final` 在当前假设下最终都接近完全恢复，更能区分方案的是 `p90_access_restore` 与 `time_avg_share_access`。基线的 `time_avg_share_access=0.273`。
- 如果后续需要进一步拉开 access 类目标，建议优先调整 critical locations 的布点或扩大受损资产集合；如果要加强搜索，则可以继续增加 SA 迭代次数。

## 建议优先写进正文的结果

1. `single_maximin_time_avg_cri`：代表效率优先。
2. `guardrail_gini_restore`：代表公平性优先。
3. `single_maximin_time_avg_cri`：代表新增的 maximin equity 规则。
