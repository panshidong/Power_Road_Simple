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
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 29373.055 | +0.0 | 0.202 | +0.0 | 0.476 | 33871.287 | 33685.718 | 0.360 | 否 |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 26852.222 | -8.6 | 0.168 | -17.2 | 0.491 | 30858.881 | 23028.731 | 0.454 | 否 |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 31948.962 | +8.8 | 0.172 | -15.0 | 0.360 | 33799.232 | 33619.240 | 0.264 | 否 |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 26880.386 | -8.5 | 0.182 | -10.1 | 0.492 | 30781.131 | 22950.981 | 0.454 | 否 |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 27956.272 | -4.8 | 0.173 | -14.3 | 0.484 | 33876.690 | 22702.233 | 0.453 | 否 |
| weighted_gini_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:gini_restore. | 27949.982 | -4.8 | 0.190 | -6.1 | 0.483 | 33876.690 | 22702.233 | 0.453 | 否 |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 26738.429 | -9.0 | 0.164 | -19.0 | 0.492 | 30587.050 | 22756.900 | 0.453 | 是 |
| weighted_maximin_time_avg_cri_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:maximin_time_avg_cri_loss. | 25327.035 | -13.8 | 0.207 | +2.4 | 0.490 | 30804.393 | 22974.243 | 0.469 | 是 |

## 表 2 代表性方案的敏感性范围

| 代表方案 | Gini restore 范围 | Min time-avg CRI 范围 | P90 restore 范围 | P90 access 范围 | Time-avg share access 范围 |
|---|---:|---:|---:|---:|---:|

## 结果讨论

- 效率最优方案是 `weighted_maximin_time_avg_cri_l050`，Triangle 为 `25327.035`，较基线改善 `-13.8%`。
- 公平性最优方案是 `guardrail_gini_restore`，Gini restore 为 `0.164`，较基线改善 `-19.0%`。
- 新增的 maximin 视角下，`guardrail_gini_restore` 取得最高的最小时间平均 CRI，`min_time_avg_cri=0.492`。
- `Triangle-Gini` 平面上的 Pareto 解为：`['guardrail_gini_restore', 'weighted_maximin_time_avg_cri_l050']`。
- `share_access_final` 在当前假设下最终都接近完全恢复，更能区分方案的是 `p90_access_restore` 与 `time_avg_share_access`。基线的 `time_avg_share_access=0.360`。
- 如果后续需要进一步拉开 access 类目标，建议优先调整 critical locations 的布点或扩大受损资产集合；如果要加强搜索，则可以继续增加 SA 迭代次数。

## 建议优先写进正文的结果

1. `weighted_maximin_time_avg_cri_l050`：代表效率优先。
2. `guardrail_gini_restore`：代表公平性优先。
3. `guardrail_gini_restore`：代表新增的 maximin equity 规则。
