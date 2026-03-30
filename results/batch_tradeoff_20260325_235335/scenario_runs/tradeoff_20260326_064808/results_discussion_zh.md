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
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 3860.604 | +0.0 | 0.258 | +0.0 | 0.333 | 4531.222 | 4531.222 | 0.294 | 否 |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 1805.225 | -53.2 | 0.307 | +19.1 | 0.495 | 2322.764 | 1871.536 | 0.473 | 否 |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 3425.732 | -11.3 | 0.162 | -37.4 | 0.168 | 3171.426 | 3171.426 | 0.202 | 是 |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 3574.446 | -7.4 | 0.216 | -16.2 | 0.368 | 3637.388 | 3200.636 | 0.376 | 否 |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 2375.874 | -38.5 | 0.201 | -22.1 | 0.340 | 2189.973 | 1736.152 | 0.424 | 是 |
| weighted_gini_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:gini_restore. | 1214.432 | -68.5 | 0.203 | -21.1 | 0.358 | 1263.214 | 1263.214 | 0.350 | 是 |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 1147.033 | -70.3 | 0.234 | -9.5 | 0.415 | 1762.756 | 1762.756 | 0.228 | 是 |
| weighted_maximin_time_avg_cri_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:maximin_time_avg_cri_loss. | 1338.857 | -65.3 | 0.391 | +51.5 | 0.507 | 2246.160 | 2246.160 | 0.488 | 否 |

## 表 2 代表性方案的敏感性范围

| 代表方案 | Gini restore 范围 | Min time-avg CRI 范围 | P90 restore 范围 | P90 access 范围 | Time-avg share access 范围 |
|---|---:|---:|---:|---:|---:|

## 结果讨论

- 效率最优方案是 `guardrail_gini_restore`，Triangle 为 `1147.033`，较基线改善 `-70.3%`。
- 公平性最优方案是 `single_gini_restore`，Gini restore 为 `0.162`，较基线改善 `-37.4%`。
- 新增的 maximin 视角下，`weighted_maximin_time_avg_cri_l050` 取得最高的最小时间平均 CRI，`min_time_avg_cri=0.507`。
- `Triangle-Gini` 平面上的 Pareto 解为：`['single_gini_restore', 'single_p90_access_restore', 'weighted_gini_restore_l050', 'guardrail_gini_restore']`。
- `share_access_final` 在当前假设下最终都接近完全恢复，更能区分方案的是 `p90_access_restore` 与 `time_avg_share_access`。基线的 `time_avg_share_access=0.294`。
- 如果后续需要进一步拉开 access 类目标，建议优先调整 critical locations 的布点或扩大受损资产集合；如果要加强搜索，则可以继续增加 SA 迭代次数。

## 建议优先写进正文的结果

1. `guardrail_gini_restore`：代表效率优先。
2. `single_gini_restore`：代表公平性优先。
3. `weighted_maximin_time_avg_cri_l050`：代表新增的 maximin equity 规则。
