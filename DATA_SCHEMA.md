# 数据表与使用约定

CSV 为 UTF-8、首行列名；大表用 gzip 压缩。原始 ZIP 和来源快照在 `data/raw/`，标准表在 `data/processed/`。缺失坐标/编号/配时使用空值，不用 0 代替。原始标识符原样保留，DSS 标识符按其不区分大小写的规则转小写。

## road

- `nodes.csv`：`node_id` 为 1–7466。1–1117 是 centroid；1118–7466 是 physical。`longitude/latitude` 取自 TAP_Demand，并与 DSTAP 全表逐点核对；`coordinate_status=source_crosschecked_datum_undeclared`。来源给出经纬度，未声明测地基准。
- `links.csv`：`link_id` 为原始文件顺序的稳定、从 1 开始的唯一标识。`from_node/to_node` 是有向端点，不能替代 link_id。`centroid_connector=1` 表示至少一个端点小于 FIRST THRU NODE。`capacity_source/length_source/free_flow_time_source/speed_source/toll_source` 保留源数值；`bpr_alpha/bpr_beta/link_type` 保留源参数。
- `od.csv.gz`：`origin/destination/demand_source/intrazonal`。包括区内需求，未缩放到每小时或车辆数。总量与源头 `TOTAL OD FLOW` 核对。
- `Austin_net.tntp`、`Austin_trips.tntp`：直接供本项目 vendor TAP-B 使用。路网头归一化成八行，所有数值、顺序、平行边保留。

`Austin_node.tntp` 为 `Node X Y` 三列坐标；X 经度，Y 纬度。`nodes.geojson` 为点，`links.geojson` 为端点直线；保留独立 link_id，平行边可能重叠显示。GeoJSON 按 WGS84 兼容坐标解释是绘图假设。`geometry_audit.csv` 记录球面端点距离，以及原长度按英尺解释的比值；没有用于改写模型长度。

来源没有足够明确的时段和单位说明，本版不为字段强加已验证的小时/分钟标签。长度、速度和时间之间呈现的换算关系可作线索，不能代替版本文档。

平行边为 `(6018,6016)`、`(6703,4574)`、`(4574,6703)`、`(6010,6009)`、`(2531,2541)`、`(4335,4334)`、`(4334,4335)`，每对两条；没有合并不同容量或阻抗的边。

## power/distribution_v03

主键应保留 feeder 范围。馈线名称形式为 `p1rhs0_1247--p1rdt5663`，`--` 前的部分为配电站 ID。馈线不是输电母线。

- `feeders.csv`：每行一个馈线，包含 substation、region、source_bus、source_kv、负荷数、线路数、变压器数和总 kW/kvar；`master_in_zip` 指向原始包的入口。`all_loads_topologically_connected` 是带电拓扑的粗检查，不替代多相潮流。
- `substations.csv`：原始 `all_substations.csv`。`Total Real Load` 单位 W；`Total Reactive Load` 为 var；`kV` 为上游接入额定电压；`Latitude/Longitude` 来自原始包。
- `buses.csv.gz`：`feeder_id/bus_id/longitude/latitude`。853,449 是馈线内坐标记录数，不是独立相节点数，也未宣称跨馈线母线完全不重复。
- `loads.csv.gz`：`load_id`、所属馈线/站、母线、完整端子 `bus_terminal`、相数、连接方式、额定 kV、标称 kW/kvar、经纬度。kV 保持 OpenDSS 该负荷连接方式下的语义，不能不区分相间/相地直接换算。
- `lines.csv.gz`：线路的母线端点和完整相位端子、长度与原始单位、linecode、switch、enabled。保留禁用线段。不能仅凭 `switch=y` 判断闭合；当前原始 `enabled` 状态需结合模型读取。
- `transformers.csv.gz`：绕组母线数组及 `raw_dss_command`。保留多绕组、各绕组 kV/kVA/%r、XHL/XHT/XLT 等，不将三绕组简化为普通二端边。
- `linecodes.csv.gz`：完整 DSS 阻抗/电纳矩阵及 normamps 命令；按 feeder+linecode 查找，矩阵单位以其 `units` 为准。
- `controls.csv.gz`：保险丝、电容及电容控制、调压控制原始命令和所在文件。它们可能有跨馈线引用。
- `solver_quality.csv`：由求解校验生成。`eligibility=import_solve_ready` 仅指已有可导入/求解的推荐入口，不代表电压或电流限制合格。`requires_regional_model=1` 的记录应使用 regional 入口，不能把区域结果当成独立馈线结果。

上述详细表针对叶馈线。更上游的变电站变压器/调压器及区域 69/230 kV 元件保存在原始 ZIP，区域 OpenDSS 校验读取完整目录。没有用表格静默代替完整电气模型。

## power/transmission_travis150_aux

`.aux` 文件原样文本保留；每种 `DATA(object, [fields])` 对象导出同名 CSV。字段采用 PowerWorld 原名，单元格保留源字符串，包括空值、开关状态和 circuit ID。

`Branch.csv` 合并二端线路及二绕组变压器的字段并集；`3WXFormer.csv` 独立保留三绕组模型。没有把高压线路强行树化，没有丢掉辅助母线。`Bus.csv` 的 173 行与标签“150”不同是原始数据事实。

`bus_geography.csv` 在 Bus 缺经纬度时按 `SubNum` 引用原始 Substation 坐标，并标明 `coordinate_source`；无坐标仍留空。它不是道路维修点映射。

本版不提供未经检验的 MATPOWER case，尤其不能忽略 AUX 的三绕组设备、变压器参数基准、并联支路和控制设置。

## signals 与 coupling

`signals/inventory.csv` 保留全量 1,345 条设施，含类别、状态、主/次控制、启用日期、修改日期和坐标。`use_in_v0=1` 仅限 `TRAFFIC + TURNED_ON + PRIMARY + 有坐标`，共 965 条。PHB、次控制器和非在役记录仍保留但不进入默认候选统计。

`cycle_seconds/green_ratio` 空，`timing_status=not_in_inventory`。`road_node_id` 仅对 428 个暂选道路候选填入；`road_match_status` 区分 `provisional_spatial_candidate`、`over_50m`、`ambiguous_nearby_nodes`、`multiple_signals_same_node` 和未进入主信号筛选的记录。字段 `use_in_v0` 为兼容保留，v1 仍使用同一筛选条件。

`signal_power_candidates.csv` 为每个入选信号灯保留最近三个合成用户负荷：

- `candidate_rank`：按球面距离排序，距离单位 m。
- `candidate_within_200m`：是否小于等于 200 m。
- `selected_provisional`：仅最近一个且在 200 m 内为 1；仍然是待审核假设。
- `mapping_kind=synthetic_nearest_load_assumption`；`real_utility_connection_verified=0`。
- `distribution_bus_id/feeder_id/substation_id/load_id`：供追溯候选电气位置，不代表新增电力负荷。信号灯的实际瓦数没有编造。
- 未匹配的 76 个信号仍保留最近候选作为检查信息，不被强制连上电网。

距离用原始经纬度计算的球面距离（地球平均半径 6,371,008.8 m）。这些距离是点间球面距离，不是沿街道/导线的距离，也不表示实际接线。

`distribution_transmission_candidates.csv` 使用 10 m 内同额定 kV 候选。128 个站都只有一个候选，但 `cross_release_verified=0`，因为地理位置一致不能证明两发布版电气工况相同。

`signal_road_candidates.csv`：965 个灯各保留三个合格物理路口；包括 road_node_id、路口经纬度、距离、次近点与最近点距离差、物理邻居数、进入边数。`selected_provisional=1` 需满足最近点 ≤50 m、间隔 ≥20 m、且无多个灯选中同一节点；共 428 个。`manual_verified=0` 表明都未人工确认。

`signal_approach_links.csv`：428 个灯的 1,427 个有向物理进口，主键为 signal_id + link_id；用进入该节点的全部物理链接作情景假设，剔除质心连接边。`outage_capacity_factor=0.5` 是未校准假设。

`signal_road_power_provisional.csv`：410 个道路与供电候选的交集。保留 signal_id、road_node_id、feeder_id、substation_id、load_id、distribution_bus_id 和两类距离。不是电网真实接线；也未增加未核实的信号灯电力负荷。

`substation_road_access_candidates.csv`：128 个配电站各三个候选道路节点，必须有物理进入/离开边；candidate_within_500m 仅用于筛选。actual_access_verified 始终为 0。

`data/scenarios/signal_outage_demo/` 保存可直接给 TAP-B 的容量变化网络、信号/路段清单及规则。`runs/tapb_baseline/link_flows.csv.gz` 和 `runs/tapb_signal_outage/link_flows.csv.gz` 逐条保留稳定 link_id；不能以端点对合并平行边。日志/报告记录原始单位 TSTT 和收敛 gap。

## 缺口的机器可读表示

`reports/dataset_validation.json` 中 `full_geographic_road_power_coupling=false` 和 `full_abc_restoration_experiments=false`。没有生成假 `bus_location.json` 或 `bus_to_link.json` 去绕过缺失输入。

节点编号、跨仓库坐标、配套路段和 OD 一致性已验证；所有节点有来源坐标。仍缺测地基准声明、道路实际中心线、人工确认的灯控进口/电力接线、维修实际出入口、场景与恢复模型。`full_geographic_road_power_coupling=false` 指实际关联尚未验证，不再表示缺节点坐标。
