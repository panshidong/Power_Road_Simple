# Task B 正式结果表与讨论

## 数据来源
- 主结果文件：`tradeoff_summary.csv`
- 图：`tradeoff_scatter.png`
- 主设定：`w_e=0.5, w_a=0.5, cri_threshold=0.9, critical_access_threshold=0.9`
- 调度映射：`link_only`，使用 `new_bus_to_link.json`

## 表 1 主实验详细结果

| 实验 | 类型 | 描述 | Triangle | ΔTri % | Gini | ΔGini % | Min time-avg CRI | P90 restore | P90 access | Time-avg share access | Pareto |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| baseline_reference | baseline | Baseline repair sequence used as the reference point. | 310.830 | +0.0 | 0.351 | +0.0 | 0.525 | 521.950 | 629.845 | 0.481 | 否 |
| single_triangle | single | Pure efficiency objective that minimizes resilience-triangle loss. | 137.096 | -55.9 | 0.356 | +1.5 | 0.433 | 221.536 | 296.655 | 0.573 | 是 |
| single_var_restore | single | Equity objective that reduces dispersion in restoration times. | 137.876 | -55.6 | 0.356 | +1.5 | 0.430 | 221.536 | 296.655 | 0.573 | 否 |
| single_gini_restore | single | Equity objective that reduces restoration-time inequality using Gini. | 351.295 | +13.0 | 0.261 | -25.7 | 0.185 | 487.258 | 496.014 | 0.323 | 是 |
| single_p90_restore | single | Equity objective that shortens the time by which 90% of zones recover. | 137.876 | -55.6 | 0.356 | +1.5 | 0.430 | 221.536 | 296.655 | 0.573 | 否 |
| single_maximin_time_avg_cri | single | Equity via maximin: improve the worst-served zone by maximizing the minimum time-averaged CRI. | 309.409 | -0.5 | 0.284 | -19.2 | 0.575 | 394.491 | 479.939 | 0.472 | 否 |
| single_p90_access_restore | single | Critical-access objective that shortens the time by which 90% of zones regain access. | 162.114 | -47.8 | 0.339 | -3.4 | 0.282 | 270.110 | 256.110 | 0.391 | 否 |
| weighted_gini_restore_l025 | weighted_sum | Weighted trade-off objective. triangle + 0.25 * equity:gini_restore. | 187.331 | -39.7 | 0.310 | -11.7 | 0.485 | 301.449 | 379.325 | 0.488 | 否 |
| weighted_gini_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:gini_restore. | 201.552 | -35.2 | 0.282 | -19.5 | 0.430 | 265.541 | 303.187 | 0.443 | 是 |
| weighted_gini_restore_l100 | weighted_sum | Weighted trade-off objective. triangle + 1.0 * equity:gini_restore. | 201.552 | -35.2 | 0.282 | -19.5 | 0.430 | 265.541 | 303.187 | 0.443 | 是 |
| guardrail_gini_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:gini_restore. | 199.949 | -35.7 | 0.330 | -6.1 | 0.266 | 322.973 | 308.973 | 0.384 | 否 |
| weighted_p90_restore_l025 | weighted_sum | Weighted trade-off objective. triangle + 0.25 * equity:p90_restore. | 201.552 | -35.2 | 0.282 | -19.5 | 0.430 | 265.541 | 303.187 | 0.443 | 是 |
| weighted_p90_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:p90_restore. | 201.552 | -35.2 | 0.282 | -19.5 | 0.430 | 265.541 | 303.187 | 0.443 | 是 |
| weighted_p90_restore_l100 | weighted_sum | Weighted trade-off objective. triangle + 1.0 * equity:p90_restore. | 155.021 | -50.1 | 0.325 | -7.4 | 0.488 | 236.732 | 315.031 | 0.482 | 是 |
| guardrail_p90_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:p90_restore. | 160.214 | -48.5 | 0.376 | +7.2 | 0.318 | 314.723 | 300.723 | 0.413 | 否 |
| weighted_maximin_time_avg_cri_l025 | weighted_sum | Weighted trade-off objective. triangle + 0.25 * equity:maximin_time_avg_cri_loss. | 182.525 | -41.3 | 0.380 | +8.3 | 0.316 | 370.276 | 356.276 | 0.500 | 否 |
| weighted_maximin_time_avg_cri_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * equity:maximin_time_avg_cri_loss. | 185.462 | -40.3 | 0.295 | -15.9 | 0.460 | 301.017 | 379.315 | 0.468 | 是 |
| weighted_maximin_time_avg_cri_l100 | weighted_sum | Weighted trade-off objective. triangle + 1.0 * equity:maximin_time_avg_cri_loss. | 292.388 | -5.9 | 0.282 | -19.7 | 0.539 | 365.272 | 430.867 | 0.461 | 是 |
| guardrail_maximin_time_avg_cri | guardrail | Triangle minimization with a fairness guardrail. Guardrail on equity:maximin_time_avg_cri_loss. | 347.328 | +11.7 | 0.313 | -10.8 | 0.543 | 482.722 | 531.105 | 0.507 | 否 |
| weighted_p90_access_restore_l025 | weighted_sum | Weighted trade-off objective. triangle + 0.25 * critical_access:p90_access_restore. | 199.203 | -35.9 | 0.362 | +3.1 | 0.258 | 361.980 | 367.980 | 0.484 | 否 |
| weighted_p90_access_restore_l050 | weighted_sum | Weighted trade-off objective. triangle + 0.5 * critical_access:p90_access_restore. | 199.203 | -35.9 | 0.362 | +3.1 | 0.258 | 361.980 | 367.980 | 0.484 | 否 |
| weighted_p90_access_restore_l100 | weighted_sum | Weighted trade-off objective. triangle + 1.0 * critical_access:p90_access_restore. | 160.567 | -48.3 | 0.368 | +5.0 | 0.500 | 311.825 | 295.031 | 0.547 | 否 |
| guardrail_p90_access_restore | guardrail | Triangle minimization with a fairness guardrail. Guardrail on critical_access:p90_access_restore. | 137.096 | -55.9 | 0.356 | +1.5 | 0.433 | 221.536 | 296.655 | 0.573 | 是 |

## 表 2 代表性方案的敏感性范围

| 代表方案 | Gini restore 范围 | Min time-avg CRI 范围 | P90 restore 范围 | P90 access 范围 | Time-avg share access 范围 |
|---|---:|---:|---:|---:|---:|
| baseline_reference | 0.286 - 0.675 | 0.377 - 0.609 | 431.473 - 629.845 | 431.473 - 629.845 | 0.440 - 0.600 |
| single_triangle | 0.356 - 0.738 | 0.380 - 0.486 | 207.536 - 296.655 | 185.651 - 296.655 | 0.516 - 0.656 |
| single_gini_restore | 0.224 - 0.703 | 0.111 - 0.259 | 466.825 - 496.014 | 466.825 - 496.014 | 0.259 - 0.387 |
| single_maximin_time_avg_cri | 0.284 - 0.686 | 0.455 - 0.575 | 299.904 - 479.939 | 394.491 - 479.939 | 0.421 - 0.530 |
| guardrail_gini_restore | 0.259 - 0.745 | 0.184 - 0.348 | 308.973 - 328.973 | 264.013 - 328.973 | 0.327 - 0.490 |

## 结果讨论

- 效率最优方案是 `single_triangle`，Triangle 为 `137.096`，较基线改善 `-55.9%`。
- 公平性最优方案是 `single_gini_restore`，Gini restore 为 `0.261`，较基线改善 `-25.7%`。
- 新增的 maximin 视角下，`single_maximin_time_avg_cri` 取得最高的最小时间平均 CRI，`min_time_avg_cri=0.575`。
- `Triangle-Gini` 平面上的 Pareto 解为：`['single_triangle', 'single_gini_restore', 'weighted_gini_restore_l050', 'weighted_gini_restore_l100', 'weighted_p90_restore_l025', 'weighted_p90_restore_l050', 'weighted_p90_restore_l100', 'weighted_maximin_time_avg_cri_l050', 'weighted_maximin_time_avg_cri_l100', 'guardrail_p90_access_restore']`。
- `share_access_final` 在当前假设下最终都接近完全恢复，更能区分方案的是 `p90_access_restore` 与 `time_avg_share_access`。基线的 `time_avg_share_access=0.481`。
- 如果后续需要进一步拉开 access 类目标，建议优先调整 critical locations 的布点或扩大受损资产集合；如果要加强搜索，则可以继续增加 SA 迭代次数。

## 建议优先写进正文的结果

1. `single_triangle`：代表效率优先。
2. `single_gini_restore`：代表公平性优先。
3. `single_maximin_time_avg_cri`：代表新增的 maximin equity 规则。
