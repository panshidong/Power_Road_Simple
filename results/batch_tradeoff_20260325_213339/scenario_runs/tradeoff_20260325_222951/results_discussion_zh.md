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
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 1925.452 | +0.0 | 0.231 | +0.0 | 0.384 | 2104.327 | 2104.327 | 0.284 | 否 |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 1855.934 | -3.6 | 0.211 | -8.6 | 0.300 | 1917.577 | 1917.577 | 0.287 | 否 |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 2470.813 | +28.3 | 0.174 | -24.8 | 0.274 | 2297.408 | 2297.408 | 0.235 | 否 |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 2278.121 | +18.3 | 0.242 | +4.5 | 0.437 | 2818.882 | 2818.882 | 0.419 | 否 |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 2044.511 | +6.2 | 0.159 | -31.3 | 0.222 | 1957.853 | 1957.853 | 0.207 | 是 |
| weighted_gini_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:gini_restore. | 1825.083 | -5.2 | 0.171 | -26.0 | 0.277 | 1762.055 | 1762.055 | 0.240 | 是 |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 2118.554 | +10.0 | 0.164 | -28.9 | 0.209 | 1923.105 | 1923.105 | 0.225 | 否 |
| weighted_maximin_time_avg_cri_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:maximin_time_avg_cri_loss. | 1925.452 | +0.0 | 0.231 | +0.0 | 0.384 | 2104.327 | 2104.327 | 0.284 | 否 |

## 表 2 代表性方案的敏感性范围

| 代表方案 | Gini restore 范围 | Min time-avg CRI 范围 | P90 restore 范围 | P90 access 范围 | Time-avg share access 范围 |
|---|---:|---:|---:|---:|---:|

## 结果讨论

- 效率最优方案是 `weighted_gini_restore_l050`，Triangle 为 `1825.083`，较基线改善 `-5.2%`。
- 公平性最优方案是 `single_p90_access_restore`，Gini restore 为 `0.159`，较基线改善 `-31.3%`。
- 新增的 maximin 视角下，`single_maximin_time_avg_cri` 取得最高的最小时间平均 CRI，`min_time_avg_cri=0.437`。
- `Triangle-Gini` 平面上的 Pareto 解为：`['single_p90_access_restore', 'weighted_gini_restore_l050']`。
- `share_access_final` 在当前假设下最终都接近完全恢复，更能区分方案的是 `p90_access_restore` 与 `time_avg_share_access`。基线的 `time_avg_share_access=0.284`。
- 如果后续需要进一步拉开 access 类目标，建议优先调整 critical locations 的布点或扩大受损资产集合；如果要加强搜索，则可以继续增加 SA 迭代次数。

## 建议优先写进正文的结果

1. `weighted_gini_restore_l050`：代表效率优先。
2. `single_p90_access_restore`：代表公平性优先。
3. `single_maximin_time_avg_cri`：代表新增的 maximin equity 规则。
