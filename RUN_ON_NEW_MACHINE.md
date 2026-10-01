# 在新机器运行 Austin A/B/C

分支：[`austin-runtime`](https://github.com/panshidong/Power_Road_Simple/tree/austin-runtime)。这个分支的仓库根目录直接是 Austin 项目，采用独立目录布局，不依赖旧文章目录。运行代码尚未在源机器构建、测试或执行；以下验收步骤必须在目标机器完成。

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

## 3. 用两进程完成小批量验收

```bash
bash runtime/run_all.sh smoke 2
```

先验证六区域 AC 健康/故障/修复状态及 TAP-B 质量，再执行缩小样本的 A/B/C 全流程。`smoke` 仍包含整个 Austin 路网和六个区域电网，因此不是瞬间结束的小算例。默认两进程至少需要满足约 32 GiB 可用内存预算；该预算是尚未实测的估计，并非实际峰值保证。

检查 `runtime/output/smoke/validation.json` 的 `passed`，以及 `runtime/output/smoke/analysis/status.json` 的 `complete`。不要把求解失败或未完成分析当作验收通过。每条结果中的 `resources` 可用于估计目标机器实际内存需求。

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
