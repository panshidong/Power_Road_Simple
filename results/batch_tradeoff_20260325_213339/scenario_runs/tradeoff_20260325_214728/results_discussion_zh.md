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
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 36528.244 | +0.0 | 0.116 | +0.0 | 0.321 | 31728.644 | 31728.644 | 0.265 | 否 |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 34435.803 | -5.7 | 0.094 | -18.8 | 0.298 | 31641.289 | 31641.289 | 0.241 | 否 |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 36407.022 | -0.3 | 0.093 | -19.6 | 0.299 | 31673.869 | 31570.702 | 0.239 | 否 |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 36368.329 | -0.4 | 0.116 | +0.5 | 0.322 | 31697.860 | 31594.693 | 0.268 | 否 |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 36143.319 | -1.1 | 0.115 | -0.8 | 0.321 | 31324.822 | 31267.060 | 0.267 | 否 |
| weighted_gini_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:gini_restore. | 36343.456 | -0.5 | 0.116 | +0.4 | 0.322 | 31669.074 | 31565.907 | 0.267 | 否 |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 36308.096 | -0.6 | 0.095 | -18.3 | 0.298 | 31686.874 | 31686.874 | 0.240 | 否 |
| weighted_maximin_time_avg_cri_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:maximin_time_avg_cri_loss. | 34387.000 | -5.9 | 0.093 | -19.7 | 0.297 | 31350.665 | 31350.665 | 0.238 | 是 |

## 表 2 代表性方案的敏感性范围

| 代表方案 | Gini restore 范围 | Min time-avg CRI 范围 | P90 restore 范围 | P90 access 范围 | Time-avg share access 范围 |
|---|---:|---:|---:|---:|---:|

## 结果讨论

- 效率最优方案是 `weighted_maximin_time_avg_cri_l050`，Triangle 为 `34387.000`，较基线改善 `-5.9%`。
- 公平性最优方案是 `weighted_maximin_time_avg_cri_l050`，Gini restore 为 `0.093`，较基线改善 `-19.7%`。
- 新增的 maximin 视角下，`weighted_gini_restore_l050` 取得最高的最小时间平均 CRI，`min_time_avg_cri=0.322`。
- `Triangle-Gini` 平面上的 Pareto 解为：`['weighted_maximin_time_avg_cri_l050']`。
- `share_access_final` 在当前假设下最终都接近完全恢复，更能区分方案的是 `p90_access_restore` 与 `time_avg_share_access`。基线的 `time_avg_share_access=0.265`。
- 如果后续需要进一步拉开 access 类目标，建议优先调整 critical locations 的布点或扩大受损资产集合；如果要加强搜索，则可以继续增加 SA 迭代次数。

## 建议优先写进正文的结果

1. `weighted_maximin_time_avg_cri_l050`：代表效率优先。
2. `weighted_maximin_time_avg_cri_l050`：代表公平性优先。
3. `weighted_gini_restore_l050`：代表新增的 maximin equity 规则。
