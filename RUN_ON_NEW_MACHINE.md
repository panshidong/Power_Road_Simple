# 在新机器运行 Austin A/B/C

分支：[`austin-runtime`](https://github.com/panshidong/Power_Road_Simple/tree/austin-runtime)。这个分支的仓库根目录直接是 Austin 项目，采用独立目录布局，不依赖旧文章目录。运行代码尚未在源机器构建、测试或执行；以下验收步骤必须在目标机器完成。

已有 clone 并准备在本机开启 Codex 会话时，先读 [AGENTS.md](AGENTS.md) 和 [CODEX_HANDOFF.txt](CODEX_HANDOFF.txt)，其中记录了当前进度及尚未解决的问题。下文安装步骤面向全新环境；已有安装先检查日志和运行状态，不必重装。

## 1. 克隆与系统依赖

以下适用于 Ubuntu/Debian Linux。其他 Linux 发行版安装同等依赖即可。

```bash
git clone --single-branch --branch austin-runtime --depth 1 \
  https://github.com/panshidong/Power_Road_Simple.git Austin
cd Austin
sudo apt-get update
sudo apt-get install -y python3 python3-venv python3-dev build-essential curl libgomp1
```

要求 Python >= 3.11。若系统默认版本较旧，请安装新版，再用例如 `AUSTIN_PYTHON=python3.12 bash runtime/bootstrap.sh`。不需要原机器上的 `.venv`、绝对路径、已编译二进制或旧论文文件。

## 2. 准备环境和数据

```bash
bash runtime/bootstrap.sh
```

该命令安装依赖，单独编译并行 TAP-B，下载缺失的 320 MB TAMU 原始电网 ZIP 并校验 SHA256，运行单元测试，准备原始区域电网和资产/空间/中心性表，打印任务计划。准备中心性表会使用计算资源，但不会开始 A/B/C 正式实验。若失败，先解决终端报错，不跳过验收。

冻结路网、节点坐标、处理后的电网表及信号灯快照已纳入分支。若官方 ZIP 下载失败，可从原主机复制 `Austin/data/raw/power/syn-Austin-TDgrid-v03.zip` 到新机器同一相对路径，再重试；不得拿另一个发布版本替换同名文件。准备脚本会检查来源锁中的哈希。

## 3. 完成物理 smoke 验收

```bash
bash runtime/run_all.sh smoke 2
```

只验证六区域 AC 健康/故障/修复状态及 TAP-B 质量，通过后结束，不执行 A/B/C 实验矩阵。实际 AC 并发由 CPU/内存预算确定，最多六个；这里的 `2` 保留原有预算参数兼容性，不代表只运行两个 AC 区域进程。完整电网的首次物理验收仍可能较慢，目前没有实测耗时保证。

检查 `runtime/output/smoke/validation.json` 与 `smoke_summary.json` 的 `passed=true`。smoke 不生成分析结果，也不证明完整研究流程已通过。旧版本 smoke 包含 225 个实验任务，现已改名为 `mini-research`；需要检查小样本 A/B/C 全流程时才显式执行 `bash runtime/run_all.sh mini-research 2`。先根据日志评估物理求解耗时，再决定是否启动完整矩阵。

如果正运行旧 smoke，先在其运行终端用 Ctrl+C 停止并等待退出，再 `git pull --ff-only origin austin-runtime`。已经使用 `flows.txt` 修复版的用户可沿用原输出目录，例如：

```bash
bash runtime/run_all.sh smoke 2 output/smoke-v2
```

匹配且已通过的 `validation.json` 会复用并直接结束。入口拆分不改变物理代码、数值配置和现有缓存指纹；无需重新安装、编译或 prepare。已有实验结果保留。

## 4. 后台运行完整实验

```bash
mkdir -p runtime/output
nohup bash runtime/run_all.sh research >> runtime/output/research.log 2>&1 &
tail -f runtime/output/research.log
```

默认 worker 数根据 CPU 和可用内存自动决定，每个 TAP-B 使用 4 线程；不需要填满所有逻辑核。若已根据验收结果确定预算，可显式指定进程数，例如：

```bash
nohup bash runtime/run_all.sh research 16 >> runtime/output/research.log 2>&1 &
```

这两个启动命令任选其一；不要同时启动。Ctrl+C 退出 `tail -f` 只关闭日志查看，不停止后台实验。若需便于交互中断，也可在 tmux 中直接运行 `bash runtime/run_all.sh research`。

主要结果在 `runtime/output/research/analysis/`，原始恢复事件及物理诊断在 `runtime/output/research/results/`。正式任务有大量潮流和交通分配调用，运行时间和缓存规模尚无实测上界，请把缓存/结果放在空间充足的 SSD。

## 5. 续跑与单独汇总

程序退出或机器重启后，再执行同一 research 命令即可跳过成功任务并恢复检查点。不要重新 bootstrap、改依赖或重新编译后混用已有结果。相同输出目录有运行锁；仍有进程时第二个调度器会拒绝启动。

```bash
# 同配置续跑
nohup bash runtime/run_all.sh research >> runtime/output/research.log 2>&1 &

# 仅对已完成结果重新汇总，不启动求解器
.venv-runtime/bin/python runtime/run.py analyze
```

更改数值配置、代码、依赖或 TAP-B 二进制需要新输出目录；只改变 worker 数可以继续原实验。高级配置、分阶段运行、失败诊断见 [runtime/README.md](runtime/README.md)，完整实验清单见 [runtime/TASK_MATRIX.md](runtime/TASK_MATRIX.md)。

## 解释结果前

这是 TAMU 六个区域原始模型的 AC 运行与道路耦合，保留各区域理想上级电源，没有联合求解完整上级输电网。空间供电关系、关键设施和抢修时间等假设见 [runtime/ASSUMPTIONS.md](runtime/ASSUMPTIONS.md)。

默认 400 个构表场景可能覆盖不了 Austin 的全部道路。结果明确记录未采样资产的中心性补位；正式解释 JSH/IJSH 前应检查 `table_coverage.json`、补位数及截尾恢复数。严格拒绝补位的配置见 `runtime/configs/strict_scores.toml`。

## 已安装旧版且 smoke 报 s.txt 不存在时

该错误由原生 TAP-B 输出 `flows.txt`、旧接口却读取 `s.txt` 引起。本次同时加入启动电网验收的区域进程池、每 30 秒进度信息和失败现场保留。停止旧进程后，在 Austin 根目录依次运行：

```bash
git pull --ff-only origin austin-runtime
```

```bash
.venv-runtime/bin/python -m unittest discover -s runtime/tests -v
```

测试通过后运行：

```bash
bash runtime/run_all.sh smoke 2 output/smoke-v2
```

不需要重新 bootstrap 或编译。新输出位于 `runtime/output/smoke-v2/`，避免与旧代码指纹混用；旧输出和缓存保留。观察 `AC regional pool: N processes` 和各区域进度；验收区域进程数会根据实际可用 CPU/内存自动选择，最多六个，与实验的两个 worker 分别计数。源机器没有执行新增测试或求解器，因此性能和完整验收仍待目标机器确认。
