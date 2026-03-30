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
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 2044.636 | +0.0 | 0.119 | +0.0 | 0.178 | 1906.293 | 1906.293 | 0.188 | 否 |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 819.288 | -59.9 | 0.128 | +8.2 | 0.512 | 784.871 | 784.871 | 0.491 | 是 |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 1813.834 | -11.3 | 0.074 | -37.2 | 0.374 | 1839.995 | 1839.995 | 0.295 | 否 |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 1415.025 | -30.8 | 0.221 | +86.3 | 0.486 | 1854.486 | 1854.486 | 0.289 | 否 |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 946.387 | -53.7 | 0.151 | +27.2 | 0.266 | 893.106 | 841.667 | 0.246 | 否 |
| weighted_gini_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:gini_restore. | 891.654 | -56.4 | 0.081 | -31.5 | 0.496 | 894.951 | 894.951 | 0.397 | 是 |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 856.286 | -58.1 | 0.091 | -23.4 | 0.543 | 842.132 | 842.132 | 0.492 | 是 |
| weighted_maximin_time_avg_cri_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:maximin_time_avg_cri_loss. | 1525.300 | -25.4 | 0.068 | -42.7 | 0.385 | 1536.764 | 1536.764 | 0.368 | 是 |

## 表 2 代表性方案的敏感性范围

| 代表方案 | Gini restore 范围 | Min time-avg CRI 范围 | P90 restore 范围 | P90 access 范围 | Time-avg share access 范围 |
|---|---:|---:|---:|---:|---:|

## 结果讨论

- 效率最优方案是 `single_triangle`，Triangle 为 `819.288`，较基线改善 `-59.9%`。
- 公平性最优方案是 `weighted_maximin_time_avg_cri_l050`，Gini restore 为 `0.068`，较基线改善 `-42.7%`。
- 新增的 maximin 视角下，`guardrail_gini_restore` 取得最高的最小时间平均 CRI，`min_time_avg_cri=0.543`。
- `Triangle-Gini` 平面上的 Pareto 解为：`['single_triangle', 'weighted_gini_restore_l050', 'guardrail_gini_restore', 'weighted_maximin_time_avg_cri_l050']`。
- `share_access_final` 在当前假设下最终都接近完全恢复，更能区分方案的是 `p90_access_restore` 与 `time_avg_share_access`。基线的 `time_avg_share_access=0.188`。
- 如果后续需要进一步拉开 access 类目标，建议优先调整 critical locations 的布点或扩大受损资产集合；如果要加强搜索，则可以继续增加 SA 迭代次数。

## 建议优先写进正文的结果

1. `single_triangle`：代表效率优先。
2. `weighted_maximin_time_avg_cri_l050`：代表公平性优先。
3. `guardrail_gini_restore`：代表新增的 maximin equity 规则。
