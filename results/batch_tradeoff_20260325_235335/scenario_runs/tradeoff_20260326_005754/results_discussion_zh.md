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
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 812.246 | +0.0 | 0.233 | +0.0 | 0.381 | 1019.143 | 1019.143 | 0.368 | 否 |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 398.111 | -51.0 | 0.271 | +16.4 | 0.478 | 661.146 | 661.146 | 0.374 | 否 |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 1098.606 | +35.3 | 0.155 | -33.4 | 0.270 | 1212.467 | 1212.467 | 0.292 | 是 |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 726.089 | -10.6 | 0.268 | +15.2 | 0.501 | 1139.365 | 1139.365 | 0.391 | 否 |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 498.692 | -38.6 | 0.251 | +7.7 | 0.415 | 655.521 | 489.050 | 0.380 | 否 |
| weighted_gini_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:gini_restore. | 397.702 | -51.0 | 0.232 | -0.5 | 0.617 | 479.129 | 479.129 | 0.539 | 否 |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 373.205 | -54.1 | 0.215 | -7.4 | 0.400 | 497.581 | 458.885 | 0.382 | 是 |
| weighted_maximin_time_avg_cri_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:maximin_time_avg_cri_loss. | 501.509 | -38.3 | 0.283 | +21.4 | 0.608 | 767.214 | 723.465 | 0.469 | 否 |

## 表 2 代表性方案的敏感性范围

| 代表方案 | Gini restore 范围 | Min time-avg CRI 范围 | P90 restore 范围 | P90 access 范围 | Time-avg share access 范围 |
|---|---:|---:|---:|---:|---:|

## 结果讨论

- 效率最优方案是 `guardrail_gini_restore`，Triangle 为 `373.205`，较基线改善 `-54.1%`。
- 公平性最优方案是 `single_gini_restore`，Gini restore 为 `0.155`，较基线改善 `-33.4%`。
- 新增的 maximin 视角下，`weighted_gini_restore_l050` 取得最高的最小时间平均 CRI，`min_time_avg_cri=0.617`。
- `Triangle-Gini` 平面上的 Pareto 解为：`['single_gini_restore', 'guardrail_gini_restore']`。
- `share_access_final` 在当前假设下最终都接近完全恢复，更能区分方案的是 `p90_access_restore` 与 `time_avg_share_access`。基线的 `time_avg_share_access=0.368`。
- 如果后续需要进一步拉开 access 类目标，建议优先调整 critical locations 的布点或扩大受损资产集合；如果要加强搜索，则可以继续增加 SA 迭代次数。

## 建议优先写进正文的结果

1. `guardrail_gini_restore`：代表效率优先。
2. `single_gini_restore`：代表公平性优先。
3. `weighted_gini_restore_l050`：代表新增的 maximin equity 规则。
