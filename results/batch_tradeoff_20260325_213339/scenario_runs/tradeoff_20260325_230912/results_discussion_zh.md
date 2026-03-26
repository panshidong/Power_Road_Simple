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
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 8393.799 | +0.0 | 0.132 | +0.0 | 0.206 | 7460.177 | 7460.177 | 0.233 | 否 |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 6319.164 | -24.7 | 0.092 | -30.4 | 0.187 | 5381.587 | 5381.587 | 0.214 | 是 |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 8669.723 | +3.3 | 0.110 | -16.4 | 0.252 | 7066.897 | 7659.230 | 0.240 | 否 |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 10264.389 | +22.3 | 0.220 | +66.6 | 0.417 | 10196.018 | 10196.018 | 0.441 | 否 |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 7007.598 | -16.5 | 0.093 | -29.4 | 0.156 | 5672.909 | 5672.909 | 0.216 | 否 |
| weighted_gini_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:gini_restore. | 8393.799 | +0.0 | 0.132 | +0.0 | 0.206 | 7460.177 | 7460.177 | 0.233 | 否 |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 7393.009 | -11.9 | 0.118 | -10.6 | 0.142 | 6050.772 | 6050.772 | 0.220 | 否 |
| weighted_maximin_time_avg_cri_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:maximin_time_avg_cri_loss. | 8393.799 | +0.0 | 0.132 | +0.0 | 0.206 | 7460.177 | 7460.177 | 0.233 | 否 |

## 表 2 代表性方案的敏感性范围

| 代表方案 | Gini restore 范围 | Min time-avg CRI 范围 | P90 restore 范围 | P90 access 范围 | Time-avg share access 范围 |
|---|---:|---:|---:|---:|---:|

## 结果讨论

- 效率最优方案是 `single_triangle`，Triangle 为 `6319.164`，较基线改善 `-24.7%`。
- 公平性最优方案是 `single_triangle`，Gini restore 为 `0.092`，较基线改善 `-30.4%`。
- 新增的 maximin 视角下，`single_maximin_time_avg_cri` 取得最高的最小时间平均 CRI，`min_time_avg_cri=0.417`。
- `Triangle-Gini` 平面上的 Pareto 解为：`['single_triangle']`。
- `share_access_final` 在当前假设下最终都接近完全恢复，更能区分方案的是 `p90_access_restore` 与 `time_avg_share_access`。基线的 `time_avg_share_access=0.233`。
- 如果后续需要进一步拉开 access 类目标，建议优先调整 critical locations 的布点或扩大受损资产集合；如果要加强搜索，则可以继续增加 SA 迭代次数。

## 建议优先写进正文的结果

1. `single_triangle`：代表效率优先。
2. `single_triangle`：代表公平性优先。
3. `single_maximin_time_avg_cri`：代表新增的 maximin equity 规则。
