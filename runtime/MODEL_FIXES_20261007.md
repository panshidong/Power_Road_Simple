# 2026-10-07 模型修正（v2 配置）

适用分支：`austin-model-fixes`。第一个提交 `Import the exact runtime code…` 原样导入了 2026-10-06 `core-unbounded-12` 结果包里的 `code_snapshot`，之后的提交都是相对这份代码的修改。

新模型只在 `configs/model-v2.toml` 中打开，通过 `extends = ["model-v2.toml"]` 被 `core-v2.toml`、`research-v2.toml` 引用。`research.toml` 和 `core-unbounded-12.toml` 保持 2026-10-06 的含义，旧配置和旧测试不受影响。

## 为什么要改：10-06 结果“反转”的原因

| 现象 | 原因 | 证据（10-06 结果） |
|---|---|---|
| JSH 均值比 CEN 差 51.5%，但 JSH 赢 17/24 场景 | 场景 15、6、16 三个场景决定了均值 | 去掉这三个场景，JSH 比 CEN 好 13.0% |
| 修好设备后供电反而下降 | 每个区域只有一个统一负荷倍率 | 场景 17：修好 p1uhs8 后 P1U 从 175.9 MW 降到 123.8 MW |
| 健康电网只供 26.3% 名义负荷 | TAMU 基准工况本身就有元件超过铭牌额定 | 健康满载时 P1U 有 20 台主变、618 条线路超铭牌，最高 4.26 倍 |
| 一个部分损坏的站让整区掉电、信号灯全灭 | 整区倍率降到 0.05，低于信号灯 0.10 的门槛 | 场景 15 的 p5uhs2：P5U 从 200 MW 降到 49 MW，320 个区可达性低于 0.9 |
| OD 策略比基础策略差 10–16% | OD 只有 3 对；分数为 0 的道路按名称字母排序 | 197 条受损道路中只有 3 条有 OD 分数；场景 15 的关键封路排在最后，约占 OD 劣势的 90% |
| 单次出行长达 831 分钟 | 静态高峰 UE 的 BPR 时间最高达自由流的 218 倍 | 自由流 25 分钟的行程在拥堵费用下要 822 分钟 |
| 表格噪声很大 | 每个站只在约 2 个构表场景中出现 | 场景间方差是 16 次排列抽样方差的 247 倍 |

结论：JSH−CEN 的损失差，与各区恢复到健康档位所需的时间差，相关系数为 0.94。10-06 的比较实际衡量的是：哪种排序碰巧先修好了卡住全区的那台变压器。

## 修改内容

### 1. 电力运行（`power.operator = "local"`）

旧规则 `uniform_grid`：整区所有负荷乘同一个倍率，按 1.0→0.01 逐档下降，直到所有元件都不超过铭牌额定。TAMU 健康满载工况本身已有元件超铭牌，所以健康电网被压到 20–40%。一台变压器受损或修复，会让整区跳一档；修复后供电下降也是这样产生的。

新规则：
- **正常额定**取 max(铭牌, 健康满载潮流)；健康电压低于 0.90 pu 的负荷节点，以健康电压作为下限。这些是合成数据本身的状况，不算灾害影响，所以健康电网供 100% 名义负荷。
- **灾后应急额定**：受损状态下，未受损元件按 `emergency_factor` × 正常额定检查（1.5，与 TAMU 全部 8,970 台变压器 EmergHKVA = 1.5 × kVA 的约定一致）。部分损坏的变压器只有 `factor` × 正常额定，不允许应急过载。
- **就地切负荷**：哪个元件越限，就只削减它下游（按求解出的有功方向）的负荷，削到额定的 98%，并重复求解。负荷节点低电压时，按 0.9 逐步削减所在馈线的负荷。
- **联络开关**：只能在不削减原有用户的前提下转供失电负荷。10-06 的 32,910 个区域状态里从未出现可用联络候选，这部分对结果没有实际影响。
- **兜底**：就地切负荷 10 轮仍不可行时，退回旧的整区倍率，并记录 `regional_fallback=true`。
- 结果新增字段：`shed_loads`、`shed_nominal_kw`、`shed_iterations`、`regional_fallback`、`preexisting_base_case`、`limit_reference`。`max_thermal_loading_ratio` 改为相对于实际执行的限值（应急额定或降额后额定），所以健康电网显示 0.667，不再是相对铭牌。
- 验收（`validation.py`）新增两项：local 模式下健康供电率必须 ≥ 0.999；全站故障和部分降额的供电都不能超过健康电网。

### 2. 修复时间（`repair_time_model = "severity"`）
- 变电站：完全失效 8 小时；部分损坏按损失容量在 2–6 小时之间插值（剩余 0.2–0.8 对应 2.8–5.2 小时）。
- 道路：封闭 4 小时；降容按损失容量在 1–3 小时之间插值。
- 这些是建模假设，不是校准过的恢复曲线，可在 `model-v2.toml` 中修改。

### 3. 抢修队数量（`crew_model = "per_damage"`）
- 每个工种每 3 个初始受损单元配 1 支队，至少 1 支、最多 6 支。主场景 8–15 个站对应 3–5 支电力队，3–13 个道路组对应 1–5 支道路队。
- 只取决于初始损坏，同一场景下各策略的队伍数相同。多支队同时从 depot 出发，按优先序依次领取任务。

### 4. 抢修队出行（`crew_travel_cost = "capped_congestion"`）
- 抢修队的路径费用取 min(UE 费用, 3 × 自由流时间)；封闭道路仍不可通行。
- IJSH 的“抢修可达性”项同样使用这个费用。
- 社区可达性（CRI）和 TSTT 仍用未封顶的 UE 费用。

### 5. Task C OD（`default_pairs = "zones_to_essential"`，`od_tie_break = "base_rule"`）
- OD 对：depot 到每个关键设施，加上每个有负荷的 TAZ 到其最近的关键设施（这正是 CRI 可达性项衡量的出行）。TAZ 只能经自己的连接弧进入路网，不会穿过其他 TAZ。
- 本机实测：597 对、1,791 条路径，51 秒，结果缓存，只算一次。有分数的道路组从 75 个增加到 2,586 个（31%）。24 个评估场景的 197 条受损道路中，有分数的从 3 条增加到 65 条。
- OD 分数相同（包括为 0）时，按对应基础策略的道路分数排序，不再按资产名称。

### 6. 构表（`core-v2.toml`）
- 24 个场景 × 16 次排列改为 72 个场景 × 4 次排列，并使用对偶排列：第 2k+1 次排列是第 2k 次的逆序。
- 原因：场景间方差远大于排列抽样方差，增加场景比增加排列更有效。
- 评估场景（24）、Task B（6 场景 × 12 次 SA）、核心阶段与 10-06 相同。

## 本机（源主机）已做的验证

本机没有运行 A/B/C 实验，只做了下面这些单区域或缩小规模的检查（同一 OpenDSSDirect 0.9.4）：

| 检查 | 结果 |
|---|---|
| P2U 健康 | 223.1/223.1 MW（100%），旧规则 40% |
| P5U 健康 | 998.4/998.4 MW（100%），旧规则 20%；1,405 条线路、29 台主变超铭牌，4,496 个负荷节点低于 0.9 pu 都作为基准状况接受 |
| 部分降额只影响本站 | P2U 6 个站、P5U 3 个站：其他站供电变化 ≤ 0.003 kW |
| 修复单调性 | P2U 12 次、P5U 8 次随机修复，供电从未下降 |
| 场景 17 的 P1U 复现 | 修好 p1uhs8 后从 584.9 MW 升到 617.9 MW（旧规则从 175.9 降到 123.8） |
| 每个受损状态的耗时（进程内，电路已载入） | P2U 2–3 s，P5U 8–15 s |
| 单元测试 | `tests/test_model_fixes.py` 11 项通过；仓库旧测试的失败与未修改快照相同（2 项，测试本身落后于快照代码） |
| 原生测试 `tests/test_power_local_native.py`（P2U） | 3 项通过，53 秒 |
| 缩小规模端到端（`Coupled`：全部 6 区 + TAP-B + 2 站 2 路修复完） | 健康 3,099.3/3,099.3 MW（100%），410 个信号灯全部有电，TAP-B gap 9.5e-5；修复时间 288/157/240/480 分钟与公式一致；出行 7–81 分钟；供电 0.9884→0.9936→1.0 单调；单进程 595 秒，峰值内存 2.9 GB |

## 仍需注意

- 修复时间、每队负责的单元数、应急倍率 1.5、拥堵上限 3 倍，都是可调的建模假设。灾害规模仍未校准（见 ASSUMPTIONS.md）。
- 部分损坏的站不能通过联络开关转供；数据里也没有可用的联络候选。
- 24 个评估场景仍然偏少：报告时请同时看中位数、胜负数和截尾均值。
- 道路封闭仍按原设计在 TAP-B 中设 9,999 分钟罚时。如果某条封闭道路是部分出行的必经之路，道路功能会大幅下降。缩小规模测试中，场景 15 的封闭路 6305-7252 让道路功能降到 0.54，贡献了 97% 的损失。电力损失改为按比例后，这一项在总损失中的比重会上升；建议先看 a-closure 的罚时敏感性，再决定是否改成“不可达出行另计”。
- 完整 core-v2 的运行时间没有实测。电力状态的求解次数应少于旧规则，但新增的 OD 表、更多构表场景和多支抢修队带来的事件数变化，都需要在目标机上观察。

## 在目标机上运行

目标机的 `runtime/` 有未提交的文件（运行代码、76 项测试、`prepared-unbounded-20261006/` 等）。先把代码和测试存进一个本地分支，再切换到本分支：

```bash
cd ~/Shidong/PowerRoad/Austin                       # 以实际路径为准
git status --short
git switch -c target-local-20261006
git add -u                                          # 只提交代码，不要加入 prepared*/、output/、cache/、build/
git add runtime/austin_runtime runtime/tests runtime/configs runtime/background_run.py
git commit -m "Target-machine runtime state of core-unbounded-12"
git fetch origin austin-model-fixes
# 下一行应无输出，说明本地代码与 10-06 快照一致；有输出就先人工核对差异
git diff 1098805 target-local-20261006 -- runtime/austin_runtime runtime/background_run.py
git switch -c austin-model-fixes --track origin/austin-model-fixes
git checkout target-local-20261006 -- runtime/tests  # 恢复目标机较新的测试；本分支新增的两个测试文件保留
cd runtime
../.venv-runtime/bin/python -m unittest discover -s tests -v
../.venv-runtime/bin/python -m austin_runtime.cli --config configs/core-v2.toml prepare
AUSTIN_NATIVE_TESTS=1 ../.venv-runtime/bin/python -m unittest tests.test_power_local_native -v   # 约几分钟
../.venv-runtime/bin/python background_run.py prepare --config configs/core-v2.toml --profile core --workers 12 --output output/core-v2-20261008
../.venv-runtime/bin/python background_run.py start output/core-v2-20261008
../.venv-runtime/bin/python background_run.py status output/core-v2-20261008
```

`1098805` 是本分支导入快照的第一个提交。

运行开始时会自动做物理验收（validation.json），新增的两项检查不通过就会停止。不要把 core-v2 的结果和 10-06 的结果放进同一输出目录；指纹不同，程序也会拒绝。
