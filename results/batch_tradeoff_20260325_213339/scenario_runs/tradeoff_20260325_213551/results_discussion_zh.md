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
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 6856.396 | +0.0 | 0.173 | +0.0 | 0.303 | 6886.592 | 6886.592 | 0.295 | 否 |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 5645.090 | -17.7 | 0.151 | -12.9 | 0.375 | 7347.253 | 7347.253 | 0.288 | 是 |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 6438.715 | -6.1 | 0.135 | -22.2 | 0.334 | 7458.618 | 7458.618 | 0.252 | 是 |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 7255.360 | +5.8 | 0.244 | +41.0 | 0.462 | 7674.826 | 7674.826 | 0.450 | 否 |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 6856.396 | +0.0 | 0.173 | +0.0 | 0.303 | 6886.592 | 6886.592 | 0.295 | 否 |
| weighted_gini_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:gini_restore. | 6438.715 | -6.1 | 0.135 | -22.2 | 0.334 | 7458.618 | 7458.618 | 0.252 | 是 |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 5645.090 | -17.7 | 0.151 | -12.9 | 0.375 | 7347.253 | 7347.253 | 0.288 | 是 |
| weighted_maximin_time_avg_cri_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:maximin_time_avg_cri_loss. | 5466.537 | -20.3 | 0.219 | +26.7 | 0.340 | 6946.909 | 7296.723 | 0.360 | 是 |

## 表 2 代表性方案的敏感性范围

| 代表方案 | Gini restore 范围 | Min time-avg CRI 范围 | P90 restore 范围 | P90 access 范围 | Time-avg share access 范围 |
|---|---:|---:|---:|---:|---:|

## 结果讨论

- 效率最优方案是 `weighted_maximin_time_avg_cri_l050`，Triangle 为 `5466.537`，较基线改善 `-20.3%`。
- 公平性最优方案是 `single_gini_restore`，Gini restore 为 `0.135`，较基线改善 `-22.2%`。
- 新增的 maximin 视角下，`single_maximin_time_avg_cri` 取得最高的最小时间平均 CRI，`min_time_avg_cri=0.462`。
- `Triangle-Gini` 平面上的 Pareto 解为：`['single_triangle', 'single_gini_restore', 'weighted_gini_restore_l050', 'guardrail_gini_restore', 'weighted_maximin_time_avg_cri_l050']`。
- `share_access_final` 在当前假设下最终都接近完全恢复，更能区分方案的是 `p90_access_restore` 与 `time_avg_share_access`。基线的 `time_avg_share_access=0.295`。
- 如果后续需要进一步拉开 access 类目标，建议优先调整 critical locations 的布点或扩大受损资产集合；如果要加强搜索，则可以继续增加 SA 迭代次数。

## 建议优先写进正文的结果

1. `weighted_maximin_time_avg_cri_l050`：代表效率优先。
2. `single_gini_restore`：代表公平性优先。
3. `single_maximin_time_avg_cri`：代表新增的 maximin equity 规则。
