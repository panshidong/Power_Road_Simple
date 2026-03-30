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
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 6531.351 | +0.0 | 0.143 | +0.0 | 0.105 | 5294.715 | 5507.901 | 0.218 | 否 |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 3106.464 | -52.4 | 0.090 | -36.8 | 0.232 | 3501.070 | 3465.421 | 0.141 | 是 |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 4450.718 | -31.9 | 0.073 | -49.0 | 0.335 | 4476.311 | 4476.311 | 0.208 | 否 |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 4370.719 | -33.1 | 0.239 | +67.4 | 0.416 | 5084.786 | 5084.786 | 0.330 | 否 |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 2474.309 | -62.1 | 0.244 | +70.6 | 0.164 | 2604.465 | 2604.465 | 0.229 | 否 |
| weighted_gini_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:gini_restore. | 4670.540 | -28.5 | 0.074 | -48.3 | 0.325 | 4602.282 | 4863.162 | 0.309 | 否 |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 2462.365 | -62.3 | 0.106 | -26.1 | 0.150 | 2755.036 | 2755.036 | 0.192 | 是 |
| weighted_maximin_time_avg_cri_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:maximin_time_avg_cri_loss. | 3882.080 | -40.6 | 0.057 | -59.9 | 0.411 | 4000.073 | 4000.073 | 0.277 | 是 |

## 表 2 代表性方案的敏感性范围

| 代表方案 | Gini restore 范围 | Min time-avg CRI 范围 | P90 restore 范围 | P90 access 范围 | Time-avg share access 范围 |
|---|---:|---:|---:|---:|---:|

## 结果讨论

- 效率最优方案是 `guardrail_gini_restore`，Triangle 为 `2462.365`，较基线改善 `-62.3%`。
- 公平性最优方案是 `weighted_maximin_time_avg_cri_l050`，Gini restore 为 `0.057`，较基线改善 `-59.9%`。
- 新增的 maximin 视角下，`single_maximin_time_avg_cri` 取得最高的最小时间平均 CRI，`min_time_avg_cri=0.416`。
- `Triangle-Gini` 平面上的 Pareto 解为：`['single_triangle', 'guardrail_gini_restore', 'weighted_maximin_time_avg_cri_l050']`。
- `share_access_final` 在当前假设下最终都接近完全恢复，更能区分方案的是 `p90_access_restore` 与 `time_avg_share_access`。基线的 `time_avg_share_access=0.218`。
- 如果后续需要进一步拉开 access 类目标，建议优先调整 critical locations 的布点或扩大受损资产集合；如果要加强搜索，则可以继续增加 SA 迭代次数。

## 建议优先写进正文的结果

1. `guardrail_gini_restore`：代表效率优先。
2. `weighted_maximin_time_avg_cri_l050`：代表公平性优先。
3. `single_maximin_time_avg_cri`：代表新增的 maximin equity 规则。
