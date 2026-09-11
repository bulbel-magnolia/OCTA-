你现在负责执行 OCTA 项目下一轮“SV + 师兄拖尾指标定义迁移分析”。

本任务是一个方法学交叉分析。目标是在保持正式 SV 算法和当前数据完全不变的条件下，把师兄 OMAG 拖尾分析中的 ROI、背景、tail AUC、denominator、normalized endpoint，以及在第二阶段可明确复现的定位/assessability 规则迁移到 SV 数据，观察其直径趋势是否与师兄结果呈现相似形态。

本任务不是为了调参得到和师兄一样的结果。任何趋势方向都必须原样报告。

==================================================
一、工作区、冻结基线与分支
==================================================

仓库：
bulbel-magnolia/OCTA-

本地工作树：
C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2

当前已冻结成果分支：
analysis/sv-diameter-stage2

当前结果解释及比较定位：
analysis/sv_fixed_diameter_selective_v2/RESULT_INTERPRETATION_20260911.md

当前固定直径正式结果：
analysis/sv_fixed_diameter_selective_v2/

正式固定直径输入：
analysis/sv_physical_diameter_geometry_v2/

有效政策：
analysis/ACTIVE_ANALYSIS_POLICY.json

当前远端基线至少应包含 commit：
2fd6cb6e8e0cc23b293e77f312983c058df0564e

先核对当前 HEAD 和远端分支实际状态。
不要假定本地一定已经同步到该 commit。

为本任务新建独立分支：
analysis/sv-senior-metric-transfer-v1

从当前 analysis/sv-diameter-stage2 最新有效 commit 建立。

新结果目录固定为：
analysis/sv_senior_metric_transfer_v1/

禁止修改或覆盖：
analysis/sv_fixed_diameter_selective_v2/
analysis/sv_physical_diameter_geometry_v2/
任何历史结果目录。

==================================================
二、科学问题
==================================================

本轮只回答四个问题：

1. 对相同 SV_raw 数据采用师兄的 absolute-tail 定义后，
   是否出现类似师兄的“中间直径 absolute AUC 较高”的形态？

   重点检查：
   D128 < D235 > D285
   是否在匹配流速中出现。

   注意：
   师兄原数据的小直径是 D185，
   我们的数据是 D128。
   所以这里只比较“形态是否类似”，
   不能称为对师兄结果的直接复现。

2. 采用师兄的 vessel-core denominator 定义后，
   denominator 是否随直径增加，例如出现：

   D128 < D235 < D285

3. 使用师兄的 normalized endpoint：

   normalized_AUC = raw_AUC / vessel_core_P95

   后，是否出现：

   D128 > D235 > D285

4. 如果结果与当前正式 SV v2 的
   tail mean / RI 结果不同，
   差异主要在“指标定义迁移”阶段已经出现，
   还是只有进一步迁移师兄的定位/assessability规则后才出现？

==================================================
三、明确保持不变的内容
==================================================

以下内容绝对不变：

1. 正式信号仍然是：

   SV_raw(x,z) = Var_t(|E(x,z,t)|)

   MATLAB 对应：
   var(abs(IMG),1,3)

   方差分母 N。

2. 不换成 OMAG。

3. 不换数据。

4. 主队列仍为：

   D128 × 1/3/5/7/10 mm/s
   D235 × 1/3/5/7/10 mm/s
   D285 × 1/3/5/7/10 mm/s

   共 15 个 scan volumes。

5. D500 不进入本轮任何主分析、敏感性分析或趋势比较。

6. 不重新采集数据。

7. 不默认重新从 .oct 重建 SV。

   优先复用已经过身份核验的 retained SV_raw arrays
   及现有 manifests。

8. dx = 12.7 μm/pixel
   dz = 6.7 μm/pixel

9. B-scan/frame 是 slow-axis 空间位置，
   不是时间重复或独立实验重复。

10. 每个 diameter×flow 只有一个 scan volume。

不做 frame-level p-value，
不把帧数作为独立实验 n，
不进行复杂推断模型。

==================================================
四、第一步必须先做：师兄代码定义审计
==================================================

在写新分析代码之前，先定位并读取师兄原始交付代码。

优先搜索以下文件名或对应内容：

export_p_bld_ed_volume.m
metrics.py
formal_9scan_freeze.py
METHOD_PARAMETERS.md

以及：
Bscan_tail_auto_quantification_for_junior_20260903
相关目录、ZIP解压目录或历史副本。

重点提取并写入：

analysis/sv_senior_metric_transfer_v1/senior_method_contract.json

必须明确记录：

- central ROI 的精确像素构造方式；
- 0.4D lateral width 的 rounding/indexing 规则；
- 左右 background ROI 的宽度；
- vessel edge 到 background ROI 的 gap；
- 左右 background 如何合并；
- depth-wise median 的实现；
- excess 的具体公式；
- positive truncation 的实现；
- operational lower boundary；
- guard；
- 200 μm AUC 的端点及像素数；
- AUC 是 sum×dz、trapz 还是其他实现；
- vessel-core 0.10D–0.45D 的具体端点规则；
- P95 的 percentile 实现；
- denominator 无效规则；
- frame assessability 规则；
- scan-level pooling/median 规则；
- 师兄定位/tracking 规则；
- 所有依赖 OMAG 信号尺度的阈值或参数。

必须使用师兄代码的真实实现。
README/文档与代码不一致时，
优先记录代码实际执行定义，同时在 audit 中写明差异。

严禁根据我们此前的口头总结重新猜实现细节。

如果师兄原始代码在本机找不到，
第一阶段仍可按下文已经冻结的核心定义执行，
第二阶段涉及 tracking/assessability 的部分标记 BLOCKED，
不得自行发明规则。

==================================================
五、第一阶段：只迁移师兄“指标定义”
==================================================

阶段名称：
STAGE_A_METRIC_TRANSFER

这是本轮最重要的主比较。

保持当前已经冻结的：

- X4
- z_top
- 原 geometry valid 身份
- 原始 frame_index 0–499 坐标

不重新定位血管。
不重新定义原 geometry validity。

新增的师兄式 metric assessability
单独保存为新字段，
不得覆盖原 valid flag。

------------------------------------------
A1. 师兄式 central vessel profile
------------------------------------------

对于每个原 geometry-valid frame：

使用真实直径 D。

横向 central ROI：

宽度 = 0.40 × D

中心仍使用冻结的 X4。

必须按照师兄代码的实际 pixel rounding/indexing
精确复现。

不要使用我们当前 full-D tail width 作为 central profile。

在每一个轴向深度 z：

P(z) =
central ROI 内 SV_raw 的横向 median

注意：

这里不是 mean，
不是 fractional weighted source mean，
而是师兄式横向 median profile。

------------------------------------------
A2. 师兄式左右 background
------------------------------------------

按照师兄代码精确建立左、右背景 ROI。

根据当前已知冻结定义，预期约为：

background lateral width ≈ D/3
vessel edge 到 background 的 gap ≈ D/4

但实际实现必须以师兄代码为准。

左右背景在同一 depth z 取值，
然后按照师兄代码的方式合并。

如果原方法是：

左右背景像素合并后取 median

则必须这样实现。

不得：
- 挑较暗的一侧；
- 挑较稳定的一侧；
- 根据结果换背景；
- 使用当前 v2 的 ±1.5D pseudo ROI 代替师兄背景。

得到：

B(z)

------------------------------------------
A3. 师兄式 excess
------------------------------------------

精确使用：

E(z) = max(P(z) - B(z), 0)

如果师兄代码存在更细节的有限值处理，
按代码执行。

这里允许 background subtraction
和 positive truncation，
因为本轮就是迁移师兄的指标定义。

但必须明确标记：

这是 secondary method-transfer analysis，
不替代当前正式 raw-SV 主分析。

------------------------------------------
A4. 操作性血管下边界
------------------------------------------

按照师兄定义：

z_lower =
z_top + round(D / dz)

其中：

dz = 6.7 μm/pixel

若师兄代码存在 MATLAB/Python 索引偏移或
inclusive/exclusive 差异，
严格复现代码。

该下边界是几何先验，
不是检测到的真实下壁。

------------------------------------------
A5. Tail start
------------------------------------------

师兄主 guard：

2 pixels

因此：

tail_start =
z_lower + 2 pixels

也就是约：

13.4 μm

不得改成当前正式分析的 guard=0。

------------------------------------------
A6. 主 absolute tail endpoint
------------------------------------------

从 tail_start 开始，
向下固定 200 μm。

严格使用师兄代码的 endpoint/indexing 方式。

主指标：

raw_AUC_200

根据当前已知定义应为类似：

raw_AUC_200 =
dz × Σ E(z)

最终实现以师兄代码为准。

明确保存单位：

SV instrument units × μm

不得把它称为当前正式的 tail mean。

同时保存：

- 有效深度 pixel 数；
- 实际积分深度；
- FOV 是否完整；
- AUC 是否可定义；
- 无效原因。

------------------------------------------
A7. 师兄式 vessel-core denominator
------------------------------------------

在血管内部：

从 z_top 向下：

0.10D 到 0.45D

的核心范围内，
使用同一个 excess profile E(z)。

主 denominator：

vessel_core_P95

即该范围内 E(z) 的 P95。

percentile 的具体算法及 endpoint
必须复现师兄代码。

不得替换成：
- whole-source mean；
- median；
- max；
- 当前 v2 source mean。

如果 denominator <= 0、
无有限样本、
范围不完整或原师兄规则判定无效，
输出 NaN 和明确 reason。

------------------------------------------
A8. 师兄式 normalized endpoint
------------------------------------------

定义：

normalized_AUC =
raw_AUC_200 / vessel_core_P95

只在 numerator 和 denominator 都有效，
且 denominator 满足师兄合法性规则时计算。

不要把这个量标为“无量纲”，
除非师兄代码/单位定义能够证明如此。

由于 numerator 是深度 AUC、
denominator 是信号幅度，
按当前数学定义其量纲应保留长度成分。

输出中建议名称明确使用：

normalized_AUC_um

或：
rAUC_senior_definition

不得与当前正式 RI 混名。

------------------------------------------
A9. Stage A 的 frame-level 输出
------------------------------------------

保留全部 7500 个 nominal positions。

至少输出：

frame_index
scan_id
diameter_um
flow_mm_s
original_geometry_valid
x4_frozen
z_top_frozen
senior_metric_assessable
assessability_reason
central_x_start/end
background_L_start/end
background_R_start/end
z_lower
tail_start
tail_end
raw_AUC_200
vessel_core_P95
normalized_AUC
以及必要 QC 字段。

原 482 个 geometry invalid 位置仍为缺失，
不得压缩掉。

新增 metric invalid 只能作为 secondary mask，
不能更改原 geometry validity。

------------------------------------------
A10. Stage A 的 volume-level 主值
------------------------------------------

按照师兄 scan-level aggregation 定义，
对每个 scan volume 汇总。

如果师兄正式规则是：
在符合 assessability 且 metric 有效的 frames 中
取 pooled median，

则严格如此。

每卷至少输出：

n_nominal
n_geometry_valid
n_metric_assessable
n_auc_defined
n_denominator_defined
n_normalized_defined

median_raw_AUC_200
median_vessel_core_P95
median_normalized_AUC

以及 min/max 或必要分布摘要。

15 卷逐卷结果全部保留。

==================================================
六、Stage A 的主要趋势判断：执行前冻结
==================================================

不要运行后再定义“像不像师兄”。

运行之前固定以下三个形态诊断：

A. absolute AUC 师兄式形态：

D128 < D235 > D285

对每个匹配 flow 单独判断。

输出：
5 个 flow 中有多少个满足。

B. denominator 师兄式形态：

D128 < D235 < D285

对每个匹配 flow 单独判断。

输出：
5 个 flow 中有多少个满足。

C. normalized endpoint 师兄式形态：

D128 > D235 > D285

对每个匹配 flow 单独判断。

输出：
5 个 flow 中有多少个满足。

同时输出每个 flow：

D128→D235
D235→D285
D128→D285

的：
- absolute difference
- percent difference

不做 p-value。

不能因为 4/5 比 5/5 不漂亮就改定义。

==================================================
七、第二阶段：迁移师兄 tracking / assessability
==================================================

阶段名称：
STAGE_B_FULL_SENIOR_PROCESS_TRANSFER

只有在师兄原始 tracking 和 assessability
实现能够明确恢复时执行。

Stage B 不替代 Stage A，
它是额外敏感性分析。

目的：

判断如果连定位、tracking、
frame assessability 也采用师兄规则，
结果是否进一步向师兄图形变化。

------------------------------------------
B1. 不覆盖冻结坐标
------------------------------------------

保留：

x4_frozen
z_top_frozen
original_geometry_valid

新增：

x_senior
z_top_senior
senior_tracking_status
senior_full_assessable

不能覆盖原字段。

------------------------------------------
B2. 参数迁移原则
------------------------------------------

师兄 tracking / assessability
中的参数全部按原代码冻结。

不得为了适应 SV 而人工调参。

但必须把参数分成两类：

1. 几何、比例、相对统计或无量纲参数
   → 可直接迁移。

2. 明确依赖 OMAG 绝对信号单位/绝对幅度尺度的阈值
   → 不得悄悄套到 SV 后再调。

如果原方法中存在这类 OMAG-specific absolute threshold：

- 在 senior_method_contract.json 中标出；
- 判断师兄代码是否本身有尺度标准化；
- 如果没有客观映射规则，
  对该 tracking/assessability 子步骤标记 BLOCKED；
- 不自行寻找让 SV 工作“最好”的阈值。

如果 tracking 方法完全可以按原代码直接在 SV 上运行，
则原参数不动执行。

------------------------------------------
B3. Stage B 指标
------------------------------------------

对 Stage B 通过师兄完整流程的 frames，
重新计算完全相同的：

raw_AUC_200
vessel_core_P95
normalized_AUC

输出同样的 volume median 和趋势判断。

必须同时做两种比较：

1. Stage A vs Stage B 各自原 support；
2. Stage A vs Stage B common-support 比较。

这样才能区分：
- 数值变化来自定位；
- 还是来自 frame selection/coverage 改变。

==================================================
八、和当前正式 SV v2 的比较
==================================================

当前正式 v2 结果保持只读。

读取：

analysis/sv_fixed_diameter_selective_v2/

至少比较：

正式：
tail100
tail500
RI100
RI500
whole-source mean

与新方法：
raw_AUC_200
vessel_core_P95
normalized_AUC

强调：
这些量定义不同、单位不同，
不得直接比较绝对数值大小。

只比较：

- 直径方向；
- 匹配 flow 下的排序；
- 定义改变后趋势是否改变。

建立：

metric_definition_direction_ledger.csv

每个 flow 一行或长表，
清楚记录：

current_tail100_direction
current_tail500_direction
current_RI100_direction
current_RI500_direction
stageA_AUC_direction
stageA_denominator_direction
stageA_normalized_direction
stageB_AUC_direction（若可用）
stageB_denominator_direction
stageB_normalized_direction

==================================================
九、不要在本轮重做的内容
==================================================

本轮不要重新做：

- 6×20 coupling 主分析；
- C16；
- whole-source coupling；
- lag；
- adjacent difference；
- current ±1 pixel sensitivity；
- current pseudo-control coupling；
- current structural coupling；
- D500；
- normalized-tail 0–1D；
- 4/8 source bins；
- frame-level hypothesis testing；
- pooled-frame p-value；
- 复杂回归/混合模型。

本轮核心就是：
把师兄的“拖尾量化定义”迁移到 SV，
并判断曲线方向。

==================================================
十、最少必要 QC
==================================================

必须完成以下 QC：

1. 输入身份
   - SV_raw 数组 SHA/manifest 与固定 v2 输入一致；
   - 不重新生成不同版本 SV。

2. 几何核查
   - D128/D235/D285 的 central ROI 宽度正确；
   - background width/gap 按师兄代码；
   - tail_start 和 200 μm window 正确；
   - core 0.10D–0.45D 正确。

3. 公式核查
   对固定测试帧独立验证：

   E(z) = max(P-B,0)

   raw_AUC
   P95
   normalized_AUC = AUC/P95

4. scan-level median
   从 framewise 表独立重算每卷 median。

5. 端点/rounding
   专门测试：
   - D/dx 非整数；
   - D/dz 非整数；
   - 0.4D；
   - D/3；
   - D/4；
   - 0.10D；
   - 0.45D；
   - 200 μm；
   的边界。

6. NaN/invalid
   不填零。
   每一种 NaN 都有 reason code。

7. FOV
   background 或 tail window 越界时，
   不平移 ROI 去“救”数据。
   按师兄规则判定 assessability。

8. 结果完整性
   不因低 AUC、
   denominator 小、
   normalized 值异常、
   趋势不符合预期而筛帧。

==================================================
十一、结果图
==================================================

主图必须仿照师兄三联图的逻辑，
但明确标注为 SV：

figure:
senior_style_sv_metrics.png
以及 PDF。

横轴：

128, 235, 285 μm

五条线：

1, 3, 5, 7, 10 mm/s

Panel A：
SV with senior definition:
raw tail AUC 0–200 μm

Panel B：
senior-definition vessel-core P95 denominator

Panel C：
senior-definition normalized AUC

图标题必须说明：

SV signal + senior metric definitions

不要写成 OMAG。

如果 Stage B 成功，
另外生成：

senior_style_sv_metrics_stageB.png

不要把 Stage A 和 Stage B 混在同一曲线里。

再生成一张：

definition_direction_comparison.png

只展示趋势方向，
比较：

current SV v2
vs
Stage A
vs
Stage B（若可用）

不比较不同单位的 y 值。

==================================================
十二、结果解释模板
==================================================

README 必须直接回答：

1. Stage A 的 raw AUC 是否呈
   D128 < D235 > D285？

   满足几种 flow？
   不满足的是哪些 flow？

2. Stage A denominator 是否呈
   D128 < D235 < D285？

3. Stage A normalized AUC 是否呈
   D128 > D235 > D285？

4. 和当前正式 tail/RI 相比，
   哪些直径排序保持，
   哪些发生改变？

5. Stage B 相比 Stage A 是否改变趋势？

6. 如果 Stage A 已经接近师兄趋势：
   结论写成：
   “趋势差异对 metric definition 敏感。”

   不写成：
   “已经证明指标定义就是原因。”

7. 如果只有 Stage B 接近师兄：
   说明 localization/assessability
   对趋势有明显影响。

8. 如果 Stage A/B 都不接近师兄：
   说明在相同数据和 SV signal 下，
   迁移师兄定义仍不能复现师兄曲线形态；
   这会增加 signal algorithm / dataset
   对差异贡献的可能性，
   但不能单独区分 OMAG algorithm 与两套原始数据。

9. normalized endpoint 即使和师兄同方向，
   也要同时报告 raw AUC 和 denominator，
   不能只展示归一化后的漂亮结果。

==================================================
十三、输出目录结构
==================================================

analysis/sv_senior_metric_transfer_v1/

至少包含：

README.md
analysis_plan.json
senior_method_contract.json
provenance.json
run_status.json
input_manifest.csv
validation.json

stageA/
  framewise_metrics_*.csv.gz
  volume_metrics.csv
  diameter_contrasts_by_flow.csv
  trend_ledger.csv
  coverage.csv

stageB/
  framewise_metrics_*.csv.gz
  volume_metrics.csv
  trend_ledger.csv
  coverage.csv
  或 BLOCKED.md

comparison/
  metric_definition_direction_ledger.csv
  stageA_stageB_common_support.csv
  current_vs_senior_definition.csv

figures/
  senior_style_sv_metrics.png
  senior_style_sv_metrics.pdf
  definition_direction_comparison.png
  definition_direction_comparison.pdf
  stageB 图（若执行成功）

validation/
  fixed_frame_checks.csv
  rounding_endpoint_checks.csv
  scan_median_reconstruction.csv

以及：
output_sha256.csv

==================================================
十四、分析身份和解释边界
==================================================

README 必须明确：

这是：
method-definition transfer / sensitivity analysis

不是新的独立实验验证。

使用的是与当前 SV v2 相同数据。

历史师兄结果在分析前已知，
因此不能称为盲验证或预注册。

D128 与师兄 D185 不同，
所以“与师兄相同趋势”
只能表示 qualitative shape similarity，
不能称 exact replication。

不同 flow 不是同条件重复。

15 volumes 是 15 个 diameter×flow 条件扫描，
不是每个直径有 5 次独立重复。

不进行机制归因。

==================================================
十五、Git 与交付
==================================================

完成全部可执行阶段并通过验证后：

1. 保持现有冻结结果不变。
2. 只提交：
   analysis/sv_senior_metric_transfer_v1/
   以及本轮直接需要的新脚本。
3. 不提交原始 .oct、大型重复 raw arrays 或环境目录。
4. commit message：

   Add SV senior-metric transfer analysis

5. push 到：

   origin/analysis/sv-senior-metric-transfer-v1

6. 不 merge 到 main/master。
7. 不 force push。
8. 不改写历史。

最后向用户报告：

- Stage A 是否成功；
- Stage B 是否成功或为何 BLOCKED；
- 三个主趋势各满足多少个 flow；
- 与当前正式 v2 相比最重要的方向变化；
- 输出目录；
- validation 状态；
- commit SHA；
- GitHub commit URL；
- 是否存在未提交文件；
- 是否有任何 BLOCKED / MISSING 项。

执行完成后停止。
不要继续自行扩展分析。