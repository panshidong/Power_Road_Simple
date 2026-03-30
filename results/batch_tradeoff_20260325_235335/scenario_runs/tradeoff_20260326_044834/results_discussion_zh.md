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
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 21868.877 | +0.0 | 0.147 | +0.0 | 0.181 | 17335.237 | 17335.237 | 0.237 | 否 |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 4755.988 | -78.3 | 0.235 | +60.0 | 0.528 | 7964.079 | 7964.079 | 0.343 | 是 |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 15679.717 | -28.3 | 0.097 | -34.3 | 0.244 | 12853.359 | 12853.359 | 0.226 | 是 |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 13349.618 | -39.0 | 0.191 | +29.7 | 0.403 | 16216.509 | 16216.509 | 0.309 | 否 |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 9028.629 | -58.7 | 0.140 | -4.7 | 0.418 | 7964.333 | 7964.333 | 0.408 | 否 |
| weighted_gini_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:gini_restore. | 12408.923 | -43.3 | 0.120 | -18.3 | 0.283 | 12786.910 | 12786.910 | 0.275 | 是 |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 8818.931 | -59.7 | 0.128 | -13.2 | 0.178 | 7984.998 | 7984.998 | 0.220 | 是 |
| weighted_maximin_time_avg_cri_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:maximin_time_avg_cri_loss. | 9228.861 | -57.8 | 0.128 | -12.7 | 0.435 | 7929.597 | 7929.597 | 0.488 | 否 |

## 表 2 代表性方案的敏感性范围

| 代表方案 | Gini restore 范围 | Min time-avg CRI 范围 | P90 restore 范围 | P90 access 范围 | Time-avg share access 范围 |
|---|---:|---:|---:|---:|---:|

## 结果讨论

- 效率最优方案是 `single_triangle`，Triangle 为 `4755.988`，较基线改善 `-78.3%`。
- 公平性最优方案是 `single_gini_restore`，Gini restore 为 `0.097`，较基线改善 `-34.3%`。
- 新增的 maximin 视角下，`single_triangle` 取得最高的最小时间平均 CRI，`min_time_avg_cri=0.528`。
- `Triangle-Gini` 平面上的 Pareto 解为：`['single_triangle', 'single_gini_restore', 'weighted_gini_restore_l050', 'guardrail_gini_restore']`。
- `share_access_final` 在当前假设下最终都接近完全恢复，更能区分方案的是 `p90_access_restore` 与 `time_avg_share_access`。基线的 `time_avg_share_access=0.237`。
- 如果后续需要进一步拉开 access 类目标，建议优先调整 critical locations 的布点或扩大受损资产集合；如果要加强搜索，则可以继续增加 SA 迭代次数。

## 建议优先写进正文的结果

1. `single_triangle`：代表效率优先。
2. `single_gini_restore`：代表公平性优先。
3. `single_triangle`：代表新增的 maximin equity 规则。
