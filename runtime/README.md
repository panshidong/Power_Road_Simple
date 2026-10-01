# Austin A/B/C 多核恢复实验运行时

本目录把三个任务迁移到 Austin_sdb + TAMU Austin 六个完整区域 OpenDSS 电网。它独立于旧文章项目，不导入旧实验结果，也不需要 `TaskA_v4/`、`TaskB_powerfix/`、`OD/` 或 `power_sys_test/` 才能运行。

**交付状态：源码、配置、安装/运行入口、测试与分析程序已写入；按用户要求，本轮没有安装依赖、编译、运行测试或启动任何求解器。可运行性、数值可行性和性能必须在目标机器验收。先前 Austin 数据检查不等于这个新运行时已通过验证。** 机器可读状态见 [IMPLEMENTATION_STATUS.json](IMPLEMENTATION_STATUS.json)。

## 在新机器开始

目标环境为 Linux，Python >= 3.11，GCC、make、curl、Python venv、系统 C/C++ 运行库。可在 Debian/Ubuntu 上自行安装 `python3-venv python3-dev build-essential curl libgomp1`。本目录不提供 Windows 原生进程锁/取消实现；WSL2 Linux 可以作为目标环境，但需检查可用内存。

在包含完整 `Austin/` 的仓库 clone 中执行：

```bash
cd Austin
bash runtime/bootstrap.sh
bash runtime/run_all.sh smoke 2
bash runtime/run_all.sh research
```

如果仓库本身就是 Austin 独立仓库，进入仓库根目录即可，不再重复 `cd Austin`。以上命令**仅供用户在新机器执行**，没有在源主机执行过。

`bootstrap.sh` 建立 `.venv-runtime/`，安装本项目，在 `runtime/build/tap-b/` 用 **make parallel** 编译含 `PARALLELISM=1` 的 TAP-B，下载并校验缺失原始文件，运行单元测试，解包完整区域模型、建立资产/TAZ/中心性表，最后打印环境和任务计划。它不会开始 A/B/C 批量实验。默认编译并发 4，可用 `AUSTIN_BUILD_JOBS` 调整；可用 `AUSTIN_PYTHON` 指定解释器。

示例给 smoke 显式使用 2 个 worker，以覆盖并行路径；资源预算不足时先调低线程/内存预算或换更充足的机器。`smoke` 保留全部六个区域和完整路网，只缩减场景数、Shapley 排列数和 SA 迭代数。它仍可能较慢，不能把它当作小电网。`research` 使用完整实验矩阵。两者输出目录和数值指纹不同，不混用训练表。

完整任务入口会先执行验收：六个区域健康工况、一个关联信号灯的变电站完整故障、部分降额、重新编译后的修复复原、TAP-B 逐弧身份/收敛/流量守恒检查。验收不通过就停止，保留 `validation.json` 和异常；不会用简化供电计数或放松 TAP-B gap 补成成功。

## 克隆需要包含什么

`Austin/.gitignore` 已调整为保留冻结的 `data/processed/`、道路/坐标/信号原始快照、来源锁及 TAP-B 源代码。单个 320 MB 官方电网 ZIP 被排除，`scripts/fetch_sources.py --download-missing` 按 `sources.lock.json` 的公开来源及 SHA256 取回。信号灯必须保留原始 JSON 快照；若只下载后来的实时清单，哈希不同会停止，不能称为同一数据版本。

运行时需要以下项目：

- `runtime/` 全部源码、配置、测试与 shell 入口。
- `data/processed/` 全部冻结表、`data/raw/` 中的冻结小文件、`sources.lock.json`。
- `vendor/tap-b/{src,include,Makefile,LICENSE.md,SOURCE.json}` 及已有依赖文件。
- `scripts/fetch_sources.py`、数据说明和来源引用。

不要从旧的 `Austin_v0_prepared.zip` 或 `Austin_v1_georeferenced.zip` 推断本运行时已经包含在里面；这些是较早的数据包，本轮没有重新打包它们。不要把 `.venv*`、求解结果或旧主机二进制当作安装依赖复制。

本项目发布到 `panshidong/Power_Road_Simple` 的独立分支 `austin-runtime`，该分支根目录直接是当前 Austin 项目内容，采用独立目录布局，旧文章分支不变。冻结数据和运行时一起纳入版本控制，第三方许可证与来源保留。请按根目录 [RUN_ON_NEW_MACHINE.md](../RUN_ON_NEW_MACHINE.md) 克隆和启动；不要在 clone 后再进入一层不存在的 `Austin/` 子目录。

## 并行与资源

实验以“一个场景 + 一个策略”为进程任务；构表以一个场景为任务。每个进程有独立 OpenDSS context，每个原生 TAP-B 调用有独立工作目录，避免 `s.txt` 等固定文件覆盖。进程采用 spawn；默认每批最多安排 workers×4 个任务，批内动态派发，批末整体更换进程池以释放原生内存。Python 较早版本的单 worker 自动回收存在挂起问题，因此不使用 `max_tasks_per_child`；参见 [Python 官方说明](https://docs.python.org/3/library/concurrent.futures.html#concurrent.futures.ProcessPoolExecutor)。批边界需要等待该批最慢任务。

默认每个 TAP-B 使用 4 线程，BLAS/OMP 限制为 1 线程。自动 worker 数取 CPU 预算和可用内存预算的较小值：CPU affinity/cgroup v2 quota 减保留的 2 核，再除以 TAP-B 线程数；内存按每 worker 16 GiB 估计，并考虑可读取的 cgroup v2 限制。**16 GiB 尚未实测，不是内存上界。** 先查看 smoke 结果的 `resources`，再调整，低于预算时会拒绝启动。

```bash
# 已准备数据后，在新机器显式使用 16 个实验进程。
bash runtime/run_all.sh research 16

# 仅列出实验，不导入求解器。
.venv-runtime/bin/python runtime/run.py plan

# 重新汇总已完成结果。
.venv-runtime/bin/python runtime/run.py analyze
```

`resources.worker_lifetime_peak_rss_mib` 是当前 worker 生命周期内最大 RSS；`native_children_peak_rss_mib` 是其原生子进程统计。不是每个任务的独立峰值，也不包含所有并发进程总和。默认预算是准入估计，没有运行中强制内存限额。

全量构表本身有 400 × 120 = 48,000 条排列路径，每条需要多个电网/交通状态。缓存可减少重复计算，但本轮没有测得耗时或空间上界。将 `cache` 和 `output` 配到足够大的本地 SSD；完整研究可能产生大量文件，不承诺小时级完成。OpenDSS 与 TAP-B 的结果缓存均使用进程锁和原子写入。缓存不会自动删除，避免后台清理破坏复现。

## 阶段、配置与续跑

配置以 [configs/research.toml](configs/research.toml) 为基础，`--config` 指定的 TOML 递归覆盖。路径 `output/cache/prepared` 相对于 `runtime/`，也可以用绝对路径。全局选项放在命令前面：

```bash
.venv-runtime/bin/python runtime/run.py --config runtime/configs/research.toml prepare
.venv-runtime/bin/python runtime/run.py --config runtime/configs/research.toml validate
.venv-runtime/bin/python runtime/run.py run --stage construct --workers 16
.venv-runtime/bin/python runtime/run.py run --stage tables
.venv-runtime/bin/python runtime/run.py run --stage a-main
.venv-runtime/bin/python runtime/run.py run --stage b-main
.venv-runtime/bin/python runtime/run.py run --stage c-main
.venv-runtime/bin/python runtime/run.py analyze
```

完整顺序和内容见 [TASK_MATRIX.md](TASK_MATRIX.md)。`a-cases` 依赖完整 `a-main`；`b-sensitivity` 的固定序列复评依赖 `b-main`；A/C 排名评估依赖 `construct` 和 `tables`。`run_all.sh` 自动按依赖顺序执行全部阶段。

同一命令再次执行会跳过签名相同的成功任务，从构表排列或 SA 迭代检查点继续。普通任务只从任务边界续跑，但状态缓存可以复用。相同输出目录有运行锁，防止两套调度器同时写。Ctrl+C 会取消本次调度器拥有的 worker 进程组及其 TAP-B 子进程，已写入的结果保留。通过 SIGKILL 强行杀父进程不保证清理子进程；优先 Ctrl+C。

输入、代码、Python/依赖版本、TAP-B 二进制、TAP-B 线程数或数值配置改变会导致指纹变化。不能覆盖旧输出继续混跑，需设置新 `output`。只改变 worker 数/资源预算或路径不改变数值指纹。不要在运行期间执行安装、编译、prepare 或编辑同一套输入。

道路参数、站点、抢修时间等影响准备目录配置指纹；改动它们后需要重新 `prepare`，并使用新实验输出。`prepare` 会写派生副本，不改原始 ZIP 或冻结处理表。正常续跑不需要再执行 `bootstrap.sh`，避免无意升级依赖或重新编译。

## 输出与判读

`runtime/output/research/` 默认包含：

| 文件/目录 | 用途 |
|---|---|
| `run_manifest.json`、`validation.json` | 环境、配置指纹、实际验收证据 |
| `plans/` | 各阶段确切任务及失败清单 |
| `results/<stage>/*.json` | 场景、完整恢复事件、派工、每区域物理运行诊断、求解质量、指标、资源占用 |
| `results/<stage>/*.error.json` | 失败原因和堆栈，不作为零损失结果 |
| `checkpoints/` | Shapley 与 SA 的可续跑状态 |
| `tables.json`、`case_selection.json` | 多个构表检查点及覆盖率、事先冻结的代表案例 |
| `analysis/scenario_rows.csv`、`aggregate.csv` | 单场景与汇总指标、截尾与补位计数 |
| `analysis/paired_statistics.json` | 配对均值/中位数/截尾均值、bootstrap 区间、胜负与精确符号检验 |
| `analysis/table_coverage.json`、`table_stability.csv` | 采样覆盖、Spearman/Kendall 稳定性 |
| `analysis/sequence_agreement.csv`、`od_coverage.csv`、`od_asset_scores.csv` | 专业内次序一致率、静态 OD 路径覆盖与分数；OD_X 与 X 的电力顺序一致性检查 |
| `analysis/task_b_pareto.csv` | 同一场景内效率—公平 Pareto 标签 |
| `analysis/figures/` | 主结果、变体、敏感性、代表恢复曲线和 CRI 热图，PNG/SVG |
| `analysis/status.json`、`SUMMARY.md` | 已计划阶段的完成/缺失/失败情况；未计划阶段不会被假称完成 |

原始事件保留每个有效 TAZ 的供电/可达性数组，足以重新计算 CRI、恢复时间与替代权重。分析时不把全研究事件数组同时载入内存。配对统计只使用相同 scenario seed 的交集，并报告未配对数。数值字段 `mean_reduction` 等表示正向改善；对 `min_time_avg_cri` 使用候选减基线，并记录 `higher_is_better=true`。

恢复未完成（`complete=false`、到达 1440 分钟或无可达工作）是显式截尾轨迹，与求解失败区分。统计包含这些轨迹，需一起报告截尾数量。只有数值成功的结果才进入分析；失败不被填零或偷偷剔除后宣称全样本完成。

## 电网范围和关键假设

使用提供方的原始六个完整区域电网，包含区域内上级网络、448 馈线、跨馈线控制、变压器、调压器、电容器和负荷。AC 求解由 [OpenDSSDirect.py / DSS-Extensions](https://dss-extensions.org/OpenDSSDirect.py/) 完成，数据来自 [TAMU combined T&D synthetic dataset](https://electricgrids.engr.tamu.edu/combined-td-synthetic-dataset/)。

运行方式是离散统一减载搜索加有限的原有联络开关闭合搜索；每个候选都要通过 AC 收敛、电压、线路电流、变压器 kVA 和有功平衡检查。它是物理约束下的可行启发式，没有复刻 IEEE33 的 LinDistFlow MILP，也不声称全局最优调度。

**六区域仍各自采用原始 230 kV 理想边界源，没有把另一个版本的 Travis150 AUX 拼成全 Austin 输配联合 AC 求解。** 相邻区域的输电拥塞、发电机出力和区域间联络约束没有建模。详细迁移差异、信号候选关系、TAZ 分配和故障定义均见 [ASSUMPTIONS.md](ASSUMPTIONS.md)。

现阶段尤其要关注表覆盖：400 个构表场景对 Austin 众多道路并不充足。默认对缺分资产显式用 CEN 补位，并在每个结果列出资产；这类结果应标为含补位的迁移试验。使用 [configs/strict_scores.toml](configs/strict_scores.toml) 可拒绝任何缺分；如需正式研究，可增加独立构表样本后检查覆盖与稳定性，不能把未采样道路的 CEN 值称作已估计的 Shapley 值。
