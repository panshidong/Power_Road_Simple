# 三任务覆盖表

本表对应新运行时的最终实验矩阵，迁移研究问题和分析内容，不承诺重现旧网络数值。原代码仅用于只读参考，运行时不依赖那些目录。

参考实现：Task A `TaskA_v4`、封闭罚时/案例 `TaskA_v5`、在线贪心 `TaskA_v6`；Task B `TaskB_powerfix`；Task C `OD/taskc_offline`；物理供电思想参考 `power_sys_test/power_sys` 并采用 TAMU 原始区域模型替代 IEEE33 参数。

| 阶段 | research 配置内容 | 主要产物 |
|---|---|---|
| `construct` | 400 场景 × 120 排列；同时构建 JSH 与 alpha=0,.5,1,2 的 IJSH | 场景边际均值、Monte Carlo 标准误、排列检查点 |
| `tables` | 50/100/200/400 场景累计表；CEN 全表 | 分资产采样计数、覆盖表 |
| `a-main` | 300 独立场景 × CEN/JSH/IJSH | H1/H2/H3 配对效率差异 |
| `a-shift` | small/large/light/clustered，各 100 × 三策略 | 分布偏移稳健性 |
| `a-alpha` | 300 × alpha=0,.5,1,2，共享联盟构表 | 可达性权重敏感性 |
| `a-closure` | 300 × 三策略，封闭罚时 1000；主实验为 9999 | 全样本罚时复评 |
| `a-greedy` | 300 × CEN/JSH/IJSH/GREEDY | 离线表对在线立即增益规则；GREEDY 决策增益日志 |
| `a-cases` | 主结果选中位/最大/最小改善案例 × 三策略 × 1000/9999/100000 | 配对代表轨迹与派工事件 |
| `b-main` | 100 × 下列 8 组；每组 SA 80 次（基线不优化） | 效率、公平性、Pareto、优化轨迹 |
| `b-sensitivity` | 第 65 场景：8 组固定最优序列换 CRI 权重两组、阈值两组；另 4 指标 × lambda=.25/.75 重新优化 | 度量敏感性与权衡系数敏感性分开记录 |
| `c-main` | 300 × CEN/JSH/IJSH/OD_CEN/OD_JSH/OD_IJSH；k=3 | 三个固定电力次序的道路排序对照 |
| `c-k` | 前 100 × k=1/2/5 × CEN/IJSH/OD_CEN/OD_IJSH | 路径数敏感性，明确保留共同基线 |
| `c-decomp` | 300 × 三基础策略及 IJSH路+JSH电、JSH路+IJSH电 | 排序来源分解 |
| `c-shift` | 四种偏移各 100 × 六策略 | OD 策略分布外检验 |
| `c-alternatives` | 下列十种关键点集合；含 allpairs、central、roadcentral 偏移，以及 allpairs 的 k 敏感性 | 关键点选择与 OD 配对模式敏感性 |
| `analyze` | 所有已完成、有同一指纹的结果 | CSV/JSON、配对 bootstrap/sign test、表稳定性、PNG/SVG |

Task B 八组：`baseline_reference`、`single_triangle`、`single_gini_restore`、`single_maximin_time_avg_cri`、`single_p90_access_restore`、`weighted_gini_restore_l050`、`weighted_maximin_time_avg_cri_l050`、`guardrail_gini_restore`。

Task B 补充重新优化的四个指标：Gini 恢复时间、P90 CRI 恢复时间、1−最小时间平均 CRI、P90 可达性恢复时间。固定序列复评用电/可达性权重 (.3,.7)、(.7,.3)，CRI 和 access 阈值 .85/.95。优化中保留每次交换、接受与温度的检查点，不把复评误称为再次优化。

Task C 十种替代集合：`allpairs`（四个默认关键点间所有有向配对）、`central`（两个高负荷站）、`roadcentral`（三个高道路中心性节点）、`randA/randB/randC`（两站，固定不同种子）、`size4/size7/size11`（相应站数）、`allnodes`（所有候选物理道路端点）。每种均做 300 × 六策略；allpairs、central、roadcentral 另做四种偏移 × 100 × 六策略；allpairs 另做 k=1/2/5 × 前100 × 四策略。无法到达的 OD 对单列，不能静默当作已覆盖。

共同指标包括标准/关键加权损失、供电与道路分项损失、Gini/方差/P90/P95 恢复时间、最小/平均时间平均 CRI、P90/P95 可达性恢复、最终 CRI 分布、shelter 初末可达性、截尾数、未采样补位数和求解诊断。图表包含主比较、偏移/k/alpha/点集、表稳定性、代表案例功能曲线、Task B 代表场景 CRI 热图及效率—公平图。所有原始指标和事件可再作图，不复用旧文章配图。

相比旧任务，Austin 必需且已记录的变化包括：负荷 kW 供电指标、物理变压器故障、地理簇、最近 TAZ、假设设施、采样介数、训练覆盖不足时显式 CEN 补位、严格 UE 收敛、封闭路段不可抢修穿行和 horizon 截尾。完整定义见 [ASSUMPTIONS.md](ASSUMPTIONS.md)。
