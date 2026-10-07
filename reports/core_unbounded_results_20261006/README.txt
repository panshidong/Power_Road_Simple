Austin 无观察窗核心版完整结果报告
生成：2026-10-06 21:18 CDT
运行：/home/acc/Shidong/PowerRoad/Austin/runtime/output/core-unbounded-12-20261006
指纹：2928ebcf3753c3f66e5be082ef1c8eef411c0ca2117b72f342f49c6b80a07c63

本版：12 worker / TAP-B 2线程 / 无观察窗 / 288任务已完成 / 264条评估全部修复。
与旧版限窗报告使用独立文件夹，不要混淆。两次P3U原生故障经恢复和独立确认通过。

1. Austin_unbounded_report.pdf：14页中文报告，包含全部17组结果及18条超过24h的完成时间。
2. Austin_unbounded_report.html：完整中文报告、264结果筛选、全部90配对统计及下载链接。
3. Austin_unbounded_all_results.xlsx：30张表，完整22指标及逐场景、派工、事件、区域等数据。
4. data/：UTF-8 BOM CSV。嵌套对象以JSON保存。Excel超过32767字符的单元格提示参看CSV/原JSON。
5. figures/：四张PNG及矢量PDF图。qa/：PDF每页预览。
6. Austin_unbounded_complete_results_20261006.zip：以上文件及288份原始任务JSON、原分析、计划、
   检查点、代码/配置、终审/物理验收、两次原生恢复的失败和独立确认记录。

ZIP original_run/results/保留完整逐事件区域电力/可达性向量与AC物理字段。
6754个新版缓存、旧缓存/结果和大体量原模型输入留在原项目，未删除或重复归档。
这是一份全部任务结果包，不是完整运行环境镜像，也不是完整研究矩阵结果。

old_vs_new_*为新旧版本比较；new_trajectory_24h_decomposition是本版轨迹事后分段积分。
旧/新TAP-B线程数不同；新旧差异不能全部归于去掉观察窗。
SHA256SUMS.txt位于ZIP内，校验每个成员；package_verification.json校验ZIP本身。
verification.json记录原始结果与前版结果SHA256、22指标重算及原物理门槛检查。

build_report.py可在原项目中用.venv-runtime/bin/python重新离线导出，不启动求解器。
程序验收通过不代表真实灾害、健康降载、交通基准和修复时间假设已获现实校准。
