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
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 2421.806 | +0.0 | 0.204 | +0.0 | 0.195 | 2174.870 | 2174.870 | 0.296 | 否 |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 903.689 | -62.7 | 0.364 | +78.5 | 0.563 | 1179.265 | 1173.265 | 0.494 | 是 |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 1791.341 | -26.0 | 0.195 | -4.4 | 0.221 | 1606.530 | 1606.530 | 0.308 | 否 |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 1902.937 | -21.4 | 0.307 | +50.6 | 0.445 | 2451.315 | 2451.315 | 0.317 | 否 |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 1215.826 | -49.8 | 0.219 | +7.3 | 0.211 | 1091.770 | 1091.770 | 0.322 | 是 |
| weighted_gini_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:gini_restore. | 1232.990 | -49.1 | 0.256 | +25.7 | 0.375 | 1374.351 | 1374.351 | 0.380 | 否 |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 1241.211 | -48.7 | 0.191 | -6.3 | 0.223 | 1183.434 | 1183.434 | 0.292 | 是 |
| weighted_maximin_time_avg_cri_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:maximin_time_avg_cri_loss. | 1094.151 | -54.8 | 0.308 | +50.9 | 0.467 | 1576.810 | 1576.810 | 0.340 | 是 |

## 表 2 代表性方案的敏感性范围

| 代表方案 | Gini restore 范围 | Min time-avg CRI 范围 | P90 restore 范围 | P90 access 范围 | Time-avg share access 范围 |
|---|---:|---:|---:|---:|---:|

## 结果讨论

- 效率最优方案是 `single_triangle`，Triangle 为 `903.689`，较基线改善 `-62.7%`。
- 公平性最优方案是 `guardrail_gini_restore`，Gini restore 为 `0.191`，较基线改善 `-6.3%`。
- 新增的 maximin 视角下，`single_triangle` 取得最高的最小时间平均 CRI，`min_time_avg_cri=0.563`。
- `Triangle-Gini` 平面上的 Pareto 解为：`['single_triangle', 'single_p90_access_restore', 'guardrail_gini_restore', 'weighted_maximin_time_avg_cri_l050']`。
- `share_access_final` 在当前假设下最终都接近完全恢复，更能区分方案的是 `p90_access_restore` 与 `time_avg_share_access`。基线的 `time_avg_share_access=0.296`。
- 如果后续需要进一步拉开 access 类目标，建议优先调整 critical locations 的布点或扩大受损资产集合；如果要加强搜索，则可以继续增加 SA 迭代次数。

## 建议优先写进正文的结果

1. `single_triangle`：代表效率优先。
2. `guardrail_gini_restore`：代表公平性优先。
3. `single_triangle`：代表新增的 maximin equity 规则。
