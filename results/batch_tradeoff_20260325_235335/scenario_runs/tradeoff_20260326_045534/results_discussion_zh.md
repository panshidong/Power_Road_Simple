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
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 2855.976 | +0.0 | 0.232 | +0.0 | 0.288 | 3107.042 | 3107.042 | 0.278 | 否 |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 1464.594 | -48.7 | 0.349 | +50.7 | 0.590 | 2746.981 | 2746.981 | 0.381 | 否 |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 2371.214 | -17.0 | 0.158 | -31.8 | 0.158 | 2208.486 | 2208.486 | 0.216 | 否 |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 1525.521 | -46.6 | 0.256 | +10.7 | 0.543 | 2171.525 | 2171.525 | 0.370 | 否 |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 907.169 | -68.2 | 0.324 | +39.9 | 0.446 | 1251.593 | 1251.593 | 0.418 | 是 |
| weighted_gini_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:gini_restore. | 1822.028 | -36.2 | 0.145 | -37.5 | 0.340 | 2007.534 | 2007.534 | 0.211 | 是 |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 1025.371 | -64.1 | 0.219 | -5.6 | 0.437 | 1485.209 | 1485.209 | 0.329 | 是 |
| weighted_maximin_time_avg_cri_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:maximin_time_avg_cri_loss. | 1528.427 | -46.5 | 0.311 | +34.4 | 0.413 | 2042.471 | 2042.471 | 0.299 | 否 |

## 表 2 代表性方案的敏感性范围

| 代表方案 | Gini restore 范围 | Min time-avg CRI 范围 | P90 restore 范围 | P90 access 范围 | Time-avg share access 范围 |
|---|---:|---:|---:|---:|---:|

## 结果讨论

- 效率最优方案是 `single_p90_access_restore`，Triangle 为 `907.169`，较基线改善 `-68.2%`。
- 公平性最优方案是 `weighted_gini_restore_l050`，Gini restore 为 `0.145`，较基线改善 `-37.5%`。
- 新增的 maximin 视角下，`single_triangle` 取得最高的最小时间平均 CRI，`min_time_avg_cri=0.590`。
- `Triangle-Gini` 平面上的 Pareto 解为：`['single_p90_access_restore', 'weighted_gini_restore_l050', 'guardrail_gini_restore']`。
- `share_access_final` 在当前假设下最终都接近完全恢复，更能区分方案的是 `p90_access_restore` 与 `time_avg_share_access`。基线的 `time_avg_share_access=0.278`。
- 如果后续需要进一步拉开 access 类目标，建议优先调整 critical locations 的布点或扩大受损资产集合；如果要加强搜索，则可以继续增加 SA 迭代次数。

## 建议优先写进正文的结果

1. `single_p90_access_restore`：代表效率优先。
2. `weighted_gini_restore_l050`：代表公平性优先。
3. `single_triangle`：代表新增的 maximin equity 规则。
