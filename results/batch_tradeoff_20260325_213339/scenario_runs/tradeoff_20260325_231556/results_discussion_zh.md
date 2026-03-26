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
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 1314.682 | +0.0 | 0.241 | +0.0 | 0.457 | 1399.773 | 1399.773 | 0.309 | 否 |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 1220.278 | -7.2 | 0.163 | -32.4 | 0.417 | 1326.369 | 1326.369 | 0.283 | 否 |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 1234.921 | -6.1 | 0.155 | -35.8 | 0.390 | 1326.562 | 1326.562 | 0.304 | 是 |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 1314.682 | +0.0 | 0.241 | +0.0 | 0.457 | 1399.773 | 1399.773 | 0.309 | 否 |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 1319.463 | +0.4 | 0.225 | -6.7 | 0.439 | 1381.734 | 1381.734 | 0.323 | 否 |
| weighted_gini_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:gini_restore. | 1209.544 | -8.0 | 0.170 | -29.6 | 0.448 | 1465.683 | 1465.683 | 0.307 | 否 |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 1193.805 | -9.2 | 0.159 | -34.0 | 0.406 | 1326.562 | 1326.562 | 0.291 | 是 |
| weighted_maximin_time_avg_cri_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:maximin_time_avg_cri_loss. | 1176.574 | -10.5 | 0.160 | -33.6 | 0.412 | 1281.885 | 1281.885 | 0.295 | 是 |

## 表 2 代表性方案的敏感性范围

| 代表方案 | Gini restore 范围 | Min time-avg CRI 范围 | P90 restore 范围 | P90 access 范围 | Time-avg share access 范围 |
|---|---:|---:|---:|---:|---:|

## 结果讨论

- 效率最优方案是 `weighted_maximin_time_avg_cri_l050`，Triangle 为 `1176.574`，较基线改善 `-10.5%`。
- 公平性最优方案是 `single_gini_restore`，Gini restore 为 `0.155`，较基线改善 `-35.8%`。
- 新增的 maximin 视角下，`baseline_reference` 取得最高的最小时间平均 CRI，`min_time_avg_cri=0.457`。
- `Triangle-Gini` 平面上的 Pareto 解为：`['single_gini_restore', 'guardrail_gini_restore', 'weighted_maximin_time_avg_cri_l050']`。
- `share_access_final` 在当前假设下最终都接近完全恢复，更能区分方案的是 `p90_access_restore` 与 `time_avg_share_access`。基线的 `time_avg_share_access=0.309`。
- 如果后续需要进一步拉开 access 类目标，建议优先调整 critical locations 的布点或扩大受损资产集合；如果要加强搜索，则可以继续增加 SA 迭代次数。

## 建议优先写进正文的结果

1. `weighted_maximin_time_avg_cri_l050`：代表效率优先。
2. `single_gini_restore`：代表公平性优先。
3. `baseline_reference`：代表新增的 maximin equity 规则。
