# Task B 正式结果表与讨论

## 数据来源
- 主结果文件：`tradeoff_summary.csv`
- 图：`tradeoff_scatter.png`
- 本文主表仅使用 `row_kind=optimized` 的 18 条主实验结果
- 主设定为 `w_e=0.5, w_a=0.5, cri_threshold=0.9, critical_access_threshold=0.9`

## 说明
- 表中的百分比变化均相对于 `baseline_reference`
- `Triangle` 越小越好
- `Gini restore` 越小越公平
- `P90 restore` 与 `P90 access restore` 越小越好
- `share_access_final` 在全部 18 条主实验中都等于 `0.9167`，区分度很弱，因此未放入主表

## 表 1 主实验详细结果

| 实验 | 类型 | Triangle | ΔTri % | Gini | ΔGini % | P90 restore | P90 access | Pareto |
|---|---|---:|---:|---:|---:|---:|---:|---|
| baseline_reference | baseline | 220.686 | 0.0 | 0.608 | 0.0 | 465.930 | 512.649 | 否 |
| single_triangle | single | 140.612 | -36.3 | 0.609 | +0.3 | 365.698 | 264.246 | 否 |
| single_var_restore | single | 228.381 | +3.5 | 0.486 | -20.0 | 351.839 | 426.037 | 否 |
| single_gini_restore | single | 217.043 | -1.7 | 0.432 | -29.0 | 394.281 | 320.448 | 否 |
| single_p90_restore | single | 147.169 | -33.3 | 0.572 | -5.9 | 364.669 | 287.864 | 否 |
| single_p90_access_restore | single | 140.612 | -36.3 | 0.609 | +0.3 | 365.698 | 264.246 | 否 |
| weighted_gini_restore_l025 | weighted_sum | 145.544 | -34.1 | 0.589 | -3.1 | 335.535 | 335.535 | 否 |
| weighted_gini_restore_l050 | weighted_sum | 192.816 | -12.6 | 0.446 | -26.6 | 267.010 | 305.933 | 否 |
| weighted_gini_restore_l100 | weighted_sum | 177.683 | -19.5 | 0.391 | -35.7 | 308.466 | 325.676 | 是 |
| guardrail_gini_restore | guardrail | 138.069 | -37.4 | 0.564 | -7.2 | 344.997 | 307.555 | 是 |
| weighted_p90_restore_l025 | weighted_sum | 149.468 | -32.3 | 0.586 | -3.6 | 294.748 | 333.671 | 否 |
| weighted_p90_restore_l050 | weighted_sum | 192.816 | -12.6 | 0.446 | -26.6 | 267.010 | 305.933 | 否 |
| weighted_p90_restore_l100 | weighted_sum | 142.011 | -35.7 | 0.613 | +0.9 | 344.438 | 319.790 | 否 |
| guardrail_p90_restore | guardrail | 140.612 | -36.3 | 0.609 | +0.3 | 365.698 | 264.246 | 否 |
| weighted_p90_access_restore_l025 | weighted_sum | 145.544 | -34.1 | 0.589 | -3.1 | 335.535 | 335.535 | 否 |
| weighted_p90_access_restore_l050 | weighted_sum | 173.165 | -21.5 | 0.599 | -1.5 | 435.625 | 327.490 | 否 |
| weighted_p90_access_restore_l100 | weighted_sum | 142.011 | -35.7 | 0.613 | +0.9 | 344.438 | 319.790 | 否 |
| guardrail_p90_access_restore | guardrail | 140.612 | -36.3 | 0.609 | +0.3 | 365.698 | 264.246 | 否 |

## 表 2 代表性方案的敏感性范围

下表使用同一条修复序列，在 27 组敏感性设定下重新评估得到的范围：

| 代表方案 | Gini restore 范围 | P90 restore 范围 | P90 access 范围 | 说明 |
|---|---:|---:|---:|---|
| baseline_reference | 0.531 - 0.656 | 410.385 - 512.649 | 386.381 - 512.649 | 基线本身对阈值设定较敏感 |
| single_triangle | 0.539 - 0.662 | 365.698 - 365.698 | 236.523 - 288.894 | 效率序列对 access 阈值较敏感，但对 P90 restore 稳定 |
| single_gini_restore | 0.364 - 0.653 | 394.281 - 394.281 | 320.448 - 320.448 | 该序列在主设定下公平性较好，但 Gini 对权重/阈值并不完全稳健 |
| weighted_gini_restore_l100 | 0.350 - 0.674 | 294.819 - 325.676 | 308.466 - 325.676 | Pareto 点之一，兼顾性最好，但公平性评价会随设定变化 |
| guardrail_gini_restore | 0.484 - 0.641 | 344.997 - 344.997 | 299.482 - 307.555 | Pareto 点之一，效率优势稳定，公平性改善中等 |

## 结果讨论

### 1. 纯效率目标确实能显著降低 Triangle，但不一定提升公平性

`single_triangle` 将 Triangle 从 `220.686` 降到 `140.612`，降幅达到 `36.3%`，说明如果只追求系统效率，当前模型可以很明显地压缩总体恢复损失。然而它的 `Gini restore` 从 `0.608` 微升到 `0.609`，几乎没有改善公平性，甚至略有恶化。这说明“恢复更快”并不自动等于“恢复更均衡”。

同样地，`single_p90_access_restore` 在当前 critical locations 设定下与 `single_triangle` 收敛到同一条序列，表明在目前的节点选择 `[1, 24]` 与阈值设定下，最早改善关键可达性时序的方案，恰好也与降低 Triangle 的方案高度一致。

### 2. 纯公平目标可以降低不平等，但代价不完全相同

`single_gini_restore` 将 `Gini restore` 从 `0.608` 降到 `0.432`，改善约 `29.0%`，而 `Triangle` 只下降 `1.7%`。这说明如果直接以 Gini 为目标，算法更倾向于重新分配恢复先后顺序，而不是最大化总体恢复速度。

`single_var_restore` 也带来明显的公平性改善，`Gini restore` 降低约 `20.0%`，但 `Triangle` 反而上升 `3.5%`。这说明不同公平性指标之间并不完全等价，且某些指标会把模型推向“更均匀但更慢”的恢复路径。

### 3. 加权和方法给出了最像“折中解”的方案

本次结果里，`weighted_gini_restore_l100` 是最值得关注的折中解之一。它把：
- `Triangle` 降低了 `19.5%`
- `Gini restore` 降低了 `35.7%`
- `P90 restore` 降低了 `33.8%`
- `P90 access restore` 降低了 `36.5%`

在 `Triangle` 与 `Gini restore` 的二维平面上，它是本次正式 run 的 Pareto 点之一。这条结果说明，当公平性项的权重足够高时，模型能找到一条既明显优于基线、又不会像纯公平目标那样严重牺牲效率的修复序列。

另一个值得注意的点是，`weighted_gini_restore_l050` 与 `weighted_p90_restore_l050` 收敛到了同一条序列，说明在中等权重下，Gini 与 P90 两类公平目标对调度的引导方向开始趋同。

### 4. Gini guardrail 是本次最“划算”的规则

`guardrail_gini_restore` 是另一个 Pareto 点，而且它在 18 条主实验中给出了**最小的 Triangle**：`138.069`，比基线降低 `37.4%`。更关键的是，它的 `Gini restore` 也从 `0.608` 降到 `0.564`，改善了 `7.2%`。

这意味着：与其完全把目标切换成“公平最优”，不如用一个温和的 Gini guardrail 去限制过度不均衡，反而可能同时拿到更好的效率结果。这在方法论上是很有价值的，因为它说明公平约束并不一定只是“加成本”；在某些情况下，它还能帮助搜索避免落到纯效率目标下的较差局部解。

### 5. Access 终值不太有用，Access 恢复时序更有用

所有 18 条主实验的 `share_access_final` 都是 `0.9167`，说明在本组设定下，最终是否达到 access 阈值几乎不能区分不同方案。真正有区分度的是**达到该阈值的时间**，也就是 `p90_access_restore`。

这也说明如果后续要写论文，“critical access status” 更适合用恢复时间统计来表达，而不适合只看最终覆盖率。

### 6. 当前结果已经出现了明显的目标对齐现象

18 条主实验只产生了 12 条不同序列，说明不少不同目标最终收敛到了相同的恢复顺序。最明显的例子是：

- `single_triangle`
- `single_p90_access_restore`
- `guardrail_p90_restore`
- `guardrail_p90_access_restore`

这四条实验都得到同一条序列：

`[(9, 10), 11, 15, (11, 14), 28, 32, 17]`

这说明在当前小规模案例、一个 crew、固定关键节点 `[1, 24]` 的设定下，系统效率和 access 时序之间存在明显一致性。因此，如果后续想要让 access 目标真正体现出独立性，可能需要：

- 更换 critical locations 的布点
- 增加 zone 划分差异
- 增加 crew 数量或扩大受损资产集合

### 7. 敏感性分析提示：结论方向稳定，但数值幅度会变化

从表 2 可以看到，同一条序列在不同 `w_e / w_a`、`CRI threshold`、`critical access threshold` 下，评价值会发生明显变化，尤其是 `Gini restore`。例如：

- `baseline_reference` 的 `Gini restore` 范围是 `0.531 - 0.656`
- `weighted_gini_restore_l100` 的 `Gini restore` 范围是 `0.350 - 0.674`

这表明“哪条序列更均衡”的方向性判断通常还在，但具体改善幅度会受参数设定影响。因此在论文里更稳妥的表述应是：

- 先报告主设定下的主结论
- 再说明 Pareto 点和代表性方案在敏感性设定下总体仍表现较强
- 避免过度强调某一个单独数字

## 建议写进正文的 3 个代表结果

如果需要在论文结果部分先挑少量方案重点展示，推荐这 3 条：

1. `single_triangle`
- 代表“纯效率优先”
- Triangle 改善最直观

2. `weighted_gini_restore_l100`
- 代表“效率-公平折中”
- 是本次最清晰的 Pareto 公平解

3. `guardrail_gini_restore`
- 代表“guardrail 规则”
- Triangle 最优，同时保持一定公平改善

这三条合在一起，已经足够支撑“效率、公平、折中”的完整叙述。
