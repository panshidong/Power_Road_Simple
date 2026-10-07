Austin 无观察窗核心版：GitHub 结果交付

本目录对应2026-10-06完成的 core-unbounded-12-20261006 运行。
288/288任务成功；264/264评估轨迹全部完成修复，unfinished repairs=0。
12 worker，每TAP-B调用2线程；机器计算约11小时41分。
18条模型恢复轨迹超过原24h窗口，最长2357.196分钟（约39小时17分）。

文件入口

完整结果包（约30.5 MB）：
https://github.com/panshidong/Power_Road_Simple/raw/refs/heads/austin-runtime/reports/core_unbounded_github_20261006/Austin_unbounded_complete_results_20261006.tar.xz

中文PDF报告（14页）：
https://github.com/panshidong/Power_Road_Simple/blob/austin-runtime/reports/core_unbounded_results_20261006/Austin_unbounded_report.pdf

全部结果Excel（30张表）：
https://github.com/panshidong/Power_Road_Simple/raw/refs/heads/austin-runtime/reports/core_unbounded_results_20261006/Austin_unbounded_all_results.xlsx

逐项CSV、HTML、图表、生成脚本及核查记录：
https://github.com/panshidong/Power_Road_Simple/tree/austin-runtime/reports/core_unbounded_results_20261006

压缩格式与完整性

原桌面ZIP为79,296,236字节，超过GitHub连接接口的单次传输上限。
本次将同样556个包内文件无损重打包为.tar.xz，逐文件SHA256核对完全一致；
没有缩减结果、重算指标或改写原始JSON。原ZIP仍保留在本机及Windows桌面。
原报告中的ZIP下载链接对应本地桌面版；从GitHub获取完整结果请使用上方.tar.xz链接。

.tar.xz SHA256：
801e2841af22628aa88f91a425f43108098edde525533f3d271f6d1b4e8b08a1
原ZIP SHA256：
610c354769782b4a767cd1bbbbf4f523b89ccd513971f937e8d9086126805dbe

在WSL/Linux中解压：
mkdir Austin_results
tar -xJf Austin_unbounded_complete_results_20261006.tar.xz -C Austin_results

解压后：original_run/results/包含全部288份任务JSON；data/提供30份CSV；
code_snapshot/保留与结果对应的冻结运行源码和配置；original_run/control/包含验收/终审；
original_run/scratch/保留两次P3U原生故障、受限重算与独立确认的输入输出。
SHA256SUMS.txt可核查全部原包成员；archive_verification.json记录无损重打包验证。
publish_manifest.json记录此次GitHub提交文件的SHA256和Git blob SHA。

结果解释边界

这是当前核心288任务，不是完整研究矩阵。两次P3U原生SIGSEGV已恢复并独立确认，
但不证明原生根因已解决。B公平性guardrail仍5/6未达标；健康基准实际可行供电为
名义需求的26.26915%；道路贡献分数稀疏及灾害/交通/设施/维修时长校准问题仍保留。
新旧版TAP-B线程数也发生变化，因此不能把所有版本差异都归为去掉观察窗。

本次用户明确授权将运行结果推送GitHub。提交仅涉及结果交付目录；没有纳入工作区
尚未提交的运行代码改动、旧论文项目、虚拟环境或全部物理缓存。
需要核对生成本次结果的代码时，应使用包内code_snapshot和原运行配置，
不能把仓库其他尚未同步的运行时文件当作本次结果的精确生成版本。
