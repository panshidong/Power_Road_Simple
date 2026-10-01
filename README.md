# Austin 大算例数据 v1：带来源坐标的 Austin_sdb

独立目录 `/home/workenv/Austin/`，2026-09-30（America/Chicago）准备。
现有 Sioux Falls + IEEE33、`power_sys_test`、A/B/C 实验和历史结果没有改动。

**结论：已找到与现有 Austin_sdb 完全对应的 7,466 个节点坐标，无需替换路网或重造 OD。完整 TAP-B 基准已运行通过，并增加信号失电容量场景。路网、电网、信号灯已可做空间候选关联；人工关联仍待审核；新增 A/B/C 运行时见下节，按本轮要求尚未执行。**

新机器从 [RUN_ON_NEW_MACHINE.md](RUN_ON_NEW_MACHINE.md) 开始：克隆 `austin-runtime` 分支后，仓库根目录就是 Austin 项目。

## 新增 A/B/C 独立运行时（尚未执行）

`runtime/` 提供多进程运行入口、TAMU 区域 AC 潮流、TAP-B 并行构建、任务 A/B/C 实验矩阵、断点续跑和分析输出。请从 [runtime/README.md](runtime/README.md) 开始；假设与迁移差异在 [runtime/ASSUMPTIONS.md](runtime/ASSUMPTIONS.md)。

按本轮要求，此运行时没有安装、构建、测试或启动；下文的检查结果仅针对先前数据版本，不能作为新运行时已验收的证据。现有 v0/v1 ZIP 也不包含这个后来新增的运行时。

## 先看这些文件

- `reports/geographic_preview.html`：离线交互地图，显示道路、节点、信号灯、变电站；可点选信号灯查看道路和供电候选，输入 `N:4089` 可定位道路节点。路段用来源端点直线示意。
- `reports/dataset_validation.json`：实际数据规模、拓扑检查、匹配统计和缺口。
- `reports/tapb_validation.json`：完整路网的 TAP-B 收敛、逐节点流量守恒与平行边检查。
- `reports/road_coordinate_audit.json`：两个公开来源的全表坐标、路段与 OD 一致性证据。
- `reports/road_geography_validation.json`：空间候选规则、覆盖率、歧义与距离统计。
- `reports/tapb_signal_outage_validation.json`：信号容量下降场景的 TAP-B 运行证据。
- `reports/opendss_validation.json`、`reports/power_exception_audit.json`：全量馈线检查及原始模型例外。
- `data/processed/power/distribution_v03/solver_quality.csv`：每条馈线应使用的模型入口。
- `sources.lock.json`：原始来源、版本、字节数与 SHA256。当前信号清单会更新，务必保留原始快照。
- `DATA_SCHEMA.md`：字段、单位、关系与使用边界。

已提供 `Austin_v1_georeferenced.zip`（上一版 `Austin_v0_prepared.zip` 保留作历史快照），包含标准表、说明、验证报告、脚本及小体积原始快照；320 MB 的官方 T+D ZIP 单独保存在 `data/raw/power/`，不重复装入数据包。该官方 ZIP 可凭固定来源链接和 SHA256 重取。

## 已生成的数据

| 内容 | 实际规模 | 位置 |
|---|---:|---|
| Austin_sdb 路网 | 7,466 个有坐标节点；1,117 分区；18,710 有向链接 | `data/processed/road/` |
| OD | 234,729 个正需求条目；总需求 695,013 | `road/od.csv.gz` |
| TAMU v03 配电馈线 | 448 馈线；128 配电变电站；307,236 负荷 | `power/distribution_v03/` |
| 配电详细表 | 853,449 个馈线内母线坐标记录；721,341 线路；132,427 变压器记录 | 同上，压缩 CSV |
| Travis150 电力—燃气版本的电力 AUX | 173 母线；140 变电站；308 二端支路；5 三绕组变压器；39 发电机 | `power/transmission_travis150_aux/` |
| 市政府信号清单 | 原始 1,345 条；筛出 965 个已启用主交通信号灯 | `signals/inventory.csv` |
| 信号灯→合成负荷候选 | 889 个最近点在 200 m 内；76 个超距；每灯保留三个候选 | `coupling/signal_power_candidates.csv` |
| 信号灯→道路路口候选 | 428 个暂选路口；其余 537 个保留未选候选 | `coupling/signal_road_candidates.csv` |
| 道路进口与联合候选 | 1,427 个有向进口；410 个灯具同时有道路和供电候选 | `coupling/signal_approach_links.csv`、`signal_road_power_provisional.csv` |
| 配电站→道路出入候选 | 128 个站各取 3 点；82 个最近点在 500 m 内 | `coupling/substation_road_access_candidates.csv` |
| 配电站→输电母线候选 | 128 个站各有一个同位置、同电压候选 | `coupling/distribution_transmission_candidates.csv` |

这里的道路版本是项目已有的 `spartalab/tap-b/net/Austin_sdb_*`，固定在子模块提交 `5cde0110979857dbe3d3df0be45c811a8352bf78`。本地原始文件与该提交公开下载的 SHA256 一致。它与 TransportationNetworks 的 **7,388 节点 Austin** 不是同一个版本，不能混用编号或 OD。

## 坐标来源与版本核对

坐标主来源为 [SPARTA TAP_Demand 的 Austin_sdb_node.txt](https://github.com/spartalab/TAP_Demand/blob/9e7463fbe2911924767b08a7e7429178efdb1e71/Austin_sdb/Austin_sdb_node.txt)，以 [DSTAP 的 Austin_sdb_node.txt](https://github.com/venktesh22/DSTAP/blob/2e2d95234f3c580e26683231832866678f31385b/Networks/Austin_sdb/Austin_sdb_node.txt) 交叉核对。实际检查：

- 两表包含完全相同的 7,466 个节点编号，经纬度逐点相同，没有空值或零值占位。
- TAP_Demand 的 18,710 条链接与本地 TAP-B 的全部十个数值字段、行序完全一致；234,729 个 OD 条目及需求也完全一致。没有按“名字相同”直接拼表。
- DSTAP 拓扑编号、链接行序和 OD 一致，但 286 条链接 FFT 被改过。因此只用它核对坐标，没有采用它的道路时间参数。差异清单为 `reports/coordinate_source_road_differences.csv`。
- 经纬度范围约为 longitude [-98.270406, -97.162178]、latitude [29.773576, 30.905132]。X/Y 是经纬度十进制度，但来源未注明测地基准；按兼容 WGS84 绘图和球面距离处理是显式假设，不能宣称测绘精度已验证。
- `nodes.geojson` 为来源点位；`links.geojson` 是链接端点之间的直线，不是实测中心线。推测原长度为英尺时，42 条物理链接的端点距离超过该长度的 1.2 倍，详见 `geometry_audit.csv`；没有用直线距离覆盖 TAP-B 原长度/时间。

Austin 确实是公开交通分配基准案例，但要区分 TransportationNetworks 的 7,388 节点版和这里的 7,466 节点 Austin_sdb。当前数据沿用后者，保留全部原始容量、时间、OD 和质心规则。

TAMU 官网的“150-bus”标签、论文 v03 的“160 transmission buses”和 AUX 中的 173 条母线记录来自不同描述或发布文件。此处以实际原始文件为准，保留辅助母线、支路 circuit ID 和三绕组对象，不裁成 150，也不把另一发布版的 AUX 宣称为 v03 PWB 的无损导出。

## 已完成的运行检查

Austin_sdb 完整路网和归一化后的输入均运行成功。转换后 TAP-B 用 13 次迭代达到 relative gap **9.4263439647e-5**，TSTT 为 **20,579,939.026198（原始数据单位）**。全部 18,710 行流量都有输出。所有正的区际 OD 可达，检查时没有允许车辆穿过其他分区质心；区内需求 19,873 原样保留，不承担路网路径流量。

信号失电试验也已实际跑通：暂选 `p5uhs3_1247` 关联的 45 个灯，将 147 条有向进口容量减半；13 次迭代达到 gap **9.3803883529e-5**，TSTT **20,692,218.992426**，相对基准增加约 **0.546%**。两次运行的最大节点流量不平衡均约 **2e-6**（输出六位小数的舍入量级）。变化幅度仅是该假设场景的数值结果，不是实际停电影响预测。

全部 448 条馈线都做了独立 OpenDSS 检查：**446 条原样独立收敛且全部负荷端带电**，读入负荷数及标称 P/Q 与转换表一致。两条例外不是丢失数据：

- `p3uhs2_1247--p3udt31865` 的电容控制引用同一站另一馈线的线路。
- `p5uhs1_1247--p5udt4629` 的电容控制引用 P5U 另一站馈线的线路。

P3U 与 P5U 两个原始区域模型均已验证收敛，全部负荷端带电，见 `power_exception_audit.json`；不能删除控制器来假装原始独立馈线已通过。原始独立检查的 `passed=false` 保留。

11 条独立馈线的 `Set Voltagebases` 漏掉源侧 12.47 kV，会把部分中压节点显示为约 1.78 pu。已单独生成 `Master_voltagebases_fixed.dss`，只补电压基准；逐条对比物理电压，最大变化为 **0 V**。ZIP 原始数据未修改。无负荷断开支路端点的零电压另外记录，不等同于负荷失电。

P3U/P5U 区域模型的最低带电节点电压分别约为 0.902/0.858 pu，后者提示需要后续电压限值分析。区域模型使用自身理想边界电源，不能等同于全网联合潮流。

这些检查证明数据能导入且基态能求解；**不等于完成热限值、保护配合、故障恢复、全网 T+D 联合仿真或 A/B/C 策略实验**。输电 AUX 目前仅转换为完整对象表，没有冒称通过输电潮流。

## 如何运行

依赖已装在本目录 `.venv/`，本地 TAP-B 的源代码、MIT 许可证及同一可执行文件已复制到 `vendor/tap-b/`。无需依赖旧项目的运行目录。

```bash
cd /home/workenv/Austin
make sources        # 验证冻结原始数据的 SHA256
make build          # 从原始包重建标准表，含拓扑、需求与引用检查
make preview        # 重建离线地理检查页
make validate-road  # 独立目录运行 TAP-B，目标 gap 1e-4
make validate-signal-scenario # 暂定供电关联下的灯控进口容量变化试验
make validate-power # 448 馈线检查 + 已知例外检查
make verify         # 核查交付文件、统计与输入快照
```

新环境先执行 `make setup`。缺失原始文件可执行 `.venv/bin/python scripts/fetch_sources.py --download-missing`；脚本遇到内容变化会停止，不会覆盖冻结版本。特别是市政府实时清单，日后重新下载不保证得到相同快照。

`runs/` 包含验证时从 ZIP 解出的 OpenDSS 模型及 TAP-B 输出。每种运行使用独立目录。11 条修正后的入口也在其中，可由 `make validate-power` 重建。运行原始模型前检查 `solver_quality.csv`；区域模型与独立馈线的边界电源不同，结果不能不加区分地合并成同一种基态。

## 信号灯怎么补，以及当前缺什么

第一版使用 Austin 市政府的实际信号位置，筛选 `TRAFFIC + TURNED_ON + PRIMARY`；保留其他记录供追溯。信号清单**没有完整周期、绿信比或相位表**，这些字段留空。现有小算例用停电导致容量下降来表达信号影响，这一聚合模型本身不要求完整配时；若后续要研究信号控制，仍需另外取得配时数据。

供电候选依据官方电网包中合成用户负荷的经纬度，用球面距离取最近三个，暂用 200 m 阈值。它是**建模假设，不是真实 Austin Energy 接线记录**；TAMU 电网本身也是合成电网。200 m 阈值和停电容量系数 0.5 都需要敏感性分析。超距记录不强行分配，当前信号设施与历史交通需求也不是同年份观测。

道路坐标已补齐。路口候选仅考虑至少有三个不同物理邻居、同时有进入和离开物理边的节点；排除分区质心、死路与两邻居几何点。每个灯保留最近三个路口。默认暂选要求最近点 ≤50 m、第二近点至少远 20 m，并且不能有两个灯同时选中同一节点。428 个灯满足要求；518 个超距、17 个近邻歧义、2 个同节点冲突均留作待核对。当前设施与历史抽象路网的差异会影响覆盖率，不应通过无上限最近邻强行填满。

`signal_approach_links.csv` 对暂选节点枚举所有进入它的物理链接，以 `link_id` 保留平行边，排除质心连接边。灯控进口方向、转向和匝道仍需人工核对；“全部进入边受影响”仅作为试验假设。410 个灯同时满足道路候选和合成负荷供电候选，见 `signal_road_power_provisional.csv`。

`substation_road_access_candidates.csv` 提供每个配电站最近的三个有进入/离开边的物理道路点，82 个站最近点在 500 m 内。距离不是实际车行入口，尚不代表可达维修车道或已完成调度映射。

`data/scenarios/signal_outage_demo/` 选择联合候选灯数量最多的合成配电站，假定其候选灯失电，将相应物理进口容量乘 0.5，然后重新求解同一个 OD。输入文件、灯表、受影响 link_id、求解日志与逐条流量均保存。这只验证数据关联和 TAP-B 容量修改可以运行，不等于已做电网故障潮流、停电传播或恢复优化。

已有代码还需后续适配：

1. `bus_to_link.json` 当前是单母线→单道路对；Austin 需要一对多供电与信号进口关联。
2. 当前道路容量调整会同时影响正反向道路；实际灯控应核对进口方向、匝道和分隔道路。
3. Austin_sdb 有 **7 组重复有向端点对**，数据已保留。旧代码用 `(u,v)` 作唯一键会丢失平行边，需改用 `link_id`。
4. `FIRST THRU NODE=1118`，质心连接边应与物理道路区分；旧的简易最短路实现不自动遵守这条限制。
5. `power_system.py` 是 IEEE33 径向潮流模型；TAMU 网状输电及不平衡配电不能直接换数组套入。
6. v03 配电标称负荷合计约 3,099.278 MW，Travis150 AUX 有功负荷约 3,254.220 MW，相差约 5%。地理候选对齐不代表工况一致，联合仿真前须协调负荷缩放、损耗与控制设置。
7. 小算例的固定 33 母线计数、Sioux Falls 基准 TSTT、损坏集合、分区和关键目的地都需重新定义。

## 原始来源与引用

- [SPARTA TAP-B](https://github.com/spartalab/tap-b)：算法及 Austin_sdb 文件；固定版本见 `sources.lock.json`。
- [SPARTA TAP_Demand](https://github.com/spartalab/TAP_Demand)：节点坐标主来源及完全一致的道路/OD 佐证。使用其代码或数据时，发布方要求引用 Jake Robbennolt, Dale Robbennolt, Stephen D. Boyles, *A Relaxed Singly Constrained Static Traffic Assignment Model with Elastic Demand: Application to Telework and Urban Development Scenarios in Austin, Texas*；原 README 已保留。
- [Venktesh Pandey 的 DSTAP](https://github.com/venktesh22/DSTAP)：节点坐标交叉验证；相关论文 Jafari, Pandey, Boyles (2017), *A decomposition approach to the static traffic assignment problem*, DOI [10.1016/j.trb.2017.09.011](https://doi.org/10.1016/j.trb.2017.09.011)。
- [TNTP Austin 对照说明](https://github.com/bstabler/TransportationNetworks/tree/master/Austin)：仅用于辨别版本，没有混入本数据。
- [TAMU Combined T+D Synthetic Dataset](https://electricgrids.engr.tamu.edu/combined-td-synthetic-dataset/)：`syn-Austin-TDgrid-v03.zip`，320,294,926 字节，含 OpenDSS、坐标、GIS 和 PWB。公开下载链接从官方公开页面元数据取得，未提交个人信息或下载登记。
- [TAMU Travis150 Electric-Gas](https://electricgrids.engr.tamu.edu/synthetic-gas-electric-test-case-for-the-travis-150-system/)：独立的 `.aux/.pwb/.pwd/.xlsx` 包，544,266 字节。AUX 是可解析文本。
- [Austin Traffic Signals and Pedestrian Signals](https://data.austintexas.gov/d/p53x-x73x)：市政府清单；本地保留下载快照。
- Li et al. (2020), *Building Highly Detailed Synthetic Electric Grid Data Sets for Combined Transmission and Distribution Systems*, DOI [10.1109/OAJPE.2020.3029278](https://doi.org/10.1109/OAJPE.2020.3029278)。
- [OpenDSSDirect.py 官方文档](https://dss-extensions.org/OpenDSSDirect.py/notebooks/GettingStarted.html)：Linux 可用的 DSS 引擎接口。

TAMU 数据是供研究使用的合成系统，不代表当地真实电网。源数据和软件的使用条件按各发布方说明保留；此目录不额外宣称统一的开放许可证。
