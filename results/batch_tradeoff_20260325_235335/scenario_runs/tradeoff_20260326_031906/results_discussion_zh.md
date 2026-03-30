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
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 5888.215 | +0.0 | 0.246 | +0.0 | 0.176 | 6085.670 | 6085.670 | 0.299 | 否 |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 2331.465 | -60.4 | 0.241 | -2.2 | 0.227 | 2992.733 | 2992.733 | 0.310 | 否 |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 3781.606 | -35.8 | 0.194 | -21.4 | 0.157 | 3552.418 | 3552.418 | 0.305 | 否 |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 3731.526 | -36.6 | 0.246 | -0.2 | 0.468 | 4157.141 | 4437.721 | 0.486 | 否 |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 3132.962 | -46.8 | 0.189 | -23.4 | 0.474 | 2477.418 | 2477.418 | 0.557 | 否 |
| weighted_gini_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:gini_restore. | 2275.849 | -61.3 | 0.182 | -26.2 | 0.246 | 2085.570 | 2085.570 | 0.291 | 是 |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 2676.830 | -54.5 | 0.215 | -12.9 | 0.410 | 3154.723 | 3154.723 | 0.389 | 否 |
| weighted_maximin_time_avg_cri_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:maximin_time_avg_cri_loss. | 3507.159 | -40.4 | 0.268 | +8.6 | 0.388 | 5697.668 | 5728.182 | 0.347 | 否 |

## 表 2 代表性方案的敏感性范围

| 代表方案 | Gini restore 范围 | Min time-avg CRI 范围 | P90 restore 范围 | P90 access 范围 | Time-avg share access 范围 |
|---|---:|---:|---:|---:|---:|

## 结果讨论

- 效率最优方案是 `weighted_gini_restore_l050`，Triangle 为 `2275.849`，较基线改善 `-61.3%`。
- 公平性最优方案是 `weighted_gini_restore_l050`，Gini restore 为 `0.182`，较基线改善 `-26.2%`。
- 新增的 maximin 视角下，`single_p90_access_restore` 取得最高的最小时间平均 CRI，`min_time_avg_cri=0.474`。
- `Triangle-Gini` 平面上的 Pareto 解为：`['weighted_gini_restore_l050']`。
- `share_access_final` 在当前假设下最终都接近完全恢复，更能区分方案的是 `p90_access_restore` 与 `time_avg_share_access`。基线的 `time_avg_share_access=0.299`。
- 如果后续需要进一步拉开 access 类目标，建议优先调整 critical locations 的布点或扩大受损资产集合；如果要加强搜索，则可以继续增加 SA 迭代次数。

## 建议优先写进正文的结果

1. `weighted_gini_restore_l050`：代表效率优先。
2. `weighted_gini_restore_l050`：代表公平性优先。
3. `single_p90_access_restore`：代表新增的 maximin equity 规则。
