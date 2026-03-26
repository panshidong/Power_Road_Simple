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
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 753.995 | +0.0 | 0.163 | +0.0 | 0.312 | 843.725 | 843.725 | 0.326 | 否 |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 608.931 | -19.2 | 0.218 | +33.9 | 0.446 | 823.724 | 823.724 | 0.341 | 否 |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 824.178 | +9.3 | 0.162 | -0.6 | 0.268 | 843.725 | 843.725 | 0.297 | 否 |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 795.564 | +5.5 | 0.243 | +49.7 | 0.396 | 976.531 | 976.531 | 0.320 | 否 |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 606.530 | -19.6 | 0.186 | +14.6 | 0.342 | 688.581 | 688.581 | 0.347 | 是 |
| weighted_gini_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:gini_restore. | 632.056 | -16.2 | 0.161 | -1.2 | 0.315 | 688.581 | 688.581 | 0.326 | 是 |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 755.388 | +0.2 | 0.156 | -4.3 | 0.293 | 823.724 | 823.724 | 0.316 | 是 |
| weighted_maximin_time_avg_cri_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:maximin_time_avg_cri_loss. | 596.828 | -20.8 | 0.235 | +44.5 | 0.467 | 831.381 | 831.381 | 0.346 | 是 |

## 表 2 代表性方案的敏感性范围

| 代表方案 | Gini restore 范围 | Min time-avg CRI 范围 | P90 restore 范围 | P90 access 范围 | Time-avg share access 范围 |
|---|---:|---:|---:|---:|---:|

## 结果讨论

- 效率最优方案是 `weighted_maximin_time_avg_cri_l050`，Triangle 为 `596.828`，较基线改善 `-20.8%`。
- 公平性最优方案是 `guardrail_gini_restore`，Gini restore 为 `0.156`，较基线改善 `-4.3%`。
- 新增的 maximin 视角下，`weighted_maximin_time_avg_cri_l050` 取得最高的最小时间平均 CRI，`min_time_avg_cri=0.467`。
- `Triangle-Gini` 平面上的 Pareto 解为：`['single_p90_access_restore', 'weighted_gini_restore_l050', 'guardrail_gini_restore', 'weighted_maximin_time_avg_cri_l050']`。
- `share_access_final` 在当前假设下最终都接近完全恢复，更能区分方案的是 `p90_access_restore` 与 `time_avg_share_access`。基线的 `time_avg_share_access=0.326`。
- 如果后续需要进一步拉开 access 类目标，建议优先调整 critical locations 的布点或扩大受损资产集合；如果要加强搜索，则可以继续增加 SA 迭代次数。

## 建议优先写进正文的结果

1. `weighted_maximin_time_avg_cri_l050`：代表效率优先。
2. `guardrail_gini_restore`：代表公平性优先。
3. `weighted_maximin_time_avg_cri_l050`：代表新增的 maximin equity 规则。
