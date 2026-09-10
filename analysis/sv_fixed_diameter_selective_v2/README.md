# 固定真实直径 ROI v2：三个核心问题的新答案

本轮已完成三直径、五个匹配流速共 15 卷的新版描述、空间关联、匹配背景、负对照与限定敏感性分析。下文的结论均对应新版固定 D 几何；保持、减弱、反向和不可定义的结果分别列出。完整逐帧、逐 cell、逐 lag 结果随本目录交付。

验收：[validation.json](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/validation.json>) 中 102 项最终数值、覆盖和身份检查通过。独立实现核查是数值验证；本批数据没有独立实验验证。

## 阅读约定与分析身份

正式信号为 SV_raw=var(abs(IMG),1,3)，方差分母 N，保持 linear raw SV。Q=sum(SV_raw×w)×12.7×6.7，mean=Q/area。Source 是真实 D×D 圆形几何，保留冻结的 X4、z_top；source 的 6 个归一化层与 tail 的 20 个 25 μm bands 均固定。下界 z_top+D/6.7 来自几何先验。

每卷保留原 frame_index=0…499，总计 7500 个名义位置；D128/D235/D285 分别有 2422/2225/2371 个原有效帧，482 个原无效位置保持缺失。背景与四方向扰动没有新增原 geometry 无效帧。所有相关在单卷原坐标上计算，不合并不同 scan 的帧。

下文 median [min,max] 是同一直径五个不同流速卷的描述，不是同条件五次重复，也不是置信区间。五个原值保存在各 across_flow 表。跨 scan 采集设置一致性未核实，绝对 SV 采用记录仪器单位；RI 无量纲。空间 lag 采用帧单位，不换算时间或物理距离。

合同在科学统计前冻结：[analysis_plan.json](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/analysis_plan.json>)、[selection_log.md](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/selection_log.md>)、[column_mapping.json](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/column_mapping.json>)。历史同数据结果参与了 C16 和诊断范围选择，informed_by_historical_results=true、independent_validation=false。D500 的 8 卷按既有事后排除决定仅保留 [d500_exclusion_record.csv](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/d500_exclusion_record.csv>)，本轮读取 D500 数组数为 0。

## 问题一：匹配流速下的 absolute source/tail 与相对 tail/source

五个匹配流速下，tail100、tail500、RI100 和 RI500 均满足 D128 > D235 > D285。Source mean 的关系不同：D128 低于另外两直径；D285 相对 D235 在 1/3/5/7 mm/s 较高，在 10 mm/s 较低。不能用一条“直径越大、所有信号都越低”的结论概括。证据：[Q1-01](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/volume_metrics.csv>)、[Q1-02](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/diameter_contrasts_by_flow.csv>)。

下表列出全部 15 卷。Source 按 10⁸、tail 按 10⁷ 仪器单位展示；RI 是先逐帧取比值再取卷内中位数。

| D (μm) | 流速 | source ×10⁸ | tail100 ×10⁷ | tail500 ×10⁷ | RI100 | RI500 |
| --- | --- | --- | --- | --- | --- | --- |
| 128 | 1 | 1.849 | 8.209 | 4.331 | 0.450 | 0.235 |
| 128 | 3 | 1.887 | 8.098 | 4.234 | 0.434 | 0.224 |
| 128 | 5 | 1.921 | 7.941 | 4.209 | 0.421 | 0.223 |
| 128 | 7 | 1.453 | 8.010 | 4.317 | 0.553 | 0.297 |
| 128 | 10 | 1.938 | 8.178 | 4.308 | 0.430 | 0.223 |
| 235 | 1 | 2.067 | 4.902 | 3.037 | 0.238 | 0.150 |
| 235 | 3 | 2.153 | 4.996 | 3.111 | 0.237 | 0.147 |
| 235 | 5 | 2.148 | 4.863 | 3.059 | 0.228 | 0.142 |
| 235 | 7 | 2.149 | 4.957 | 3.043 | 0.234 | 0.143 |
| 235 | 10 | 2.200 | 4.998 | 3.071 | 0.229 | 0.140 |
| 285 | 1 | 2.176 | 4.522 | 2.886 | 0.205 | 0.133 |
| 285 | 3 | 2.238 | 4.534 | 2.927 | 0.202 | 0.131 |
| 285 | 5 | 2.205 | 4.555 | 2.909 | 0.204 | 0.131 |
| 285 | 7 | 2.164 | 4.549 | 2.952 | 0.208 | 0.137 |
| 285 | 10 | 2.114 | 4.347 | 2.809 | 0.202 | 0.131 |

RI500 的逐流速百分比变化如下；差值和其余指标的同类比较均保存在对比长表。

| 参考 D | 比较 D | 1 mm/s (%) | 3 mm/s (%) | 5 mm/s (%) | 7 mm/s (%) | 10 mm/s (%) |
| --- | --- | --- | --- | --- | --- | --- |
| 128 | 235 | -36.469 | -34.283 | -36.222 | -51.914 | -37.237 |
| 128 | 285 | -43.601 | -41.590 | -41.360 | -53.928 | -41.110 |
| 235 | 285 | -11.226 | -11.119 | -8.056 | -4.189 | -6.172 |

全部 20-band mean 与逐帧 RI 的卷内摘要见 [Q1-03](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/raw_depth_profiles.csv>)。深度坐标使用 band 中点 12.5、37.5、…、487.5 μm；末 band 是 475–500 μm 的平均，不能解释为在 500 μm 的点测量。固定空间块 0–99、100–199、…、400–499 的五项指标均保留，作为覆盖与异质性描述，不进行跨卷解剖同位配对：[Q1-04](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/spatial_block_summary.csv>)。

Q 和实际面积的描述见 [absolute_integral_metrics.csv](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/absolute_integral_metrics.csv>)，解析圆面积、离散面积偏差见 [geometry_area.csv.gz](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/geometry_area.csv.gz>)。未将栅格圆面积近似要求为解析面积的机器精度恒等。

| 状态 | 新版判断 |
| --- | --- |
| 保持 | 三直径 RI 的匹配流速排序保持；绝对 tail 的排序也逐流速成立。 |
| 减弱／改变 | 新旧几何下 RI 的具体数值改变，必须使用新版数值。不能沿用旧 RI 的幅度；同身份逐帧对照已单独保存。 |
| 反向／例外 | 10 mm/s 的 D285 source mean 低于 D235，与其余四个匹配流速的方向相反。 |
| 不可定义／未核实 | 482 个原 geometry 无效位置的逐帧 RI 为缺失；本批有效帧未出现额外 source 零分母。跨 scan 设置一致性及经校准的绝对强度可比性未得到独立核实。 |

新旧 RI 对照证据：[Q1-05](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/geometry_version_direction_ledger.csv>)。逐帧 RI_new/RI_old 恒等式已核对；未把该恒等式套用到各自的卷中位数，也未计算“分母解释百分比”。

## 问题二：同卷关联的深度分布与零位对应

C16 固定为 S3–S6 × T1–T4。其定义是全部 16 个有符号 lag0 Spearman rho 的中位数；任一 cell 不可定义时主值为 NaN。C100 则是按 Q/area 合并 S3–S6 后与 pooled tail100 的标量相关，两者不是同一指标。

本轮 raw 和 51-frame 地图各有完整 1800 个 cell。51-frame 下 C16 的全部 240 个 cell 都为正，但完整地图保留了 320 个负值；raw 地图保留 166 个负值。两种主地图均无不可定义 cell，也无恰好为零的 cell。方向一致的候选区域并不代表全深度都正相关。证据：[Q2-01](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/fixed_region_summary.csv>)、[Q2-02](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/depth_maps_raw51.csv>)。

| D | 模式 | 流速 1 | 3 | 5 | 7 | 10 | median [min,max] | 正/负/零/可定义卷 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 128 | raw | 0.190 | 0.177 | 0.186 | 0.165 | 0.154 | 0.177 [0.154, 0.190] | 5/0/0/5 |
| 128 | detrended_51 | 0.164 | 0.151 | 0.152 | 0.142 | 0.153 | 0.152 [0.142, 0.164] | 5/0/0/5 |
| 235 | raw | 0.268 | 0.171 | 0.203 | 0.196 | 0.279 | 0.203 [0.171, 0.279] | 5/0/0/5 |
| 235 | detrended_51 | 0.182 | 0.122 | 0.168 | 0.134 | 0.221 | 0.168 [0.122, 0.221] | 5/0/0/5 |
| 285 | raw | 0.184 | 0.225 | 0.245 | 0.401 | 0.267 | 0.245 [0.184, 0.401] | 5/0/0/5 |
| 285 | detrended_51 | 0.137 | 0.084 | 0.134 | 0.110 | 0.172 | 0.134 [0.084, 0.172] | 5/0/0/5 |

51-frame 的 C16 范围为 0.084–0.221，说明该固定区域的关联幅度有限。本轮不设强／弱阈值。六个 source 层各自对 T1–T4 的四项中位数均完整保留；没有重选热点或用最大 cell 替代固定区域。

三个代表量的 51-frame 结果如下。Specificity 采用固定 42 个远端 lag；matched specificity 用每个远端 lag 的共同配对坐标核查组成差异。正 specificity 不等同于峰必在零位。证据：[Q2-03](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/whole_coupling.csv>)。

| D | 配对 | rho0 median [min,max] | specificity | matched specificity | 峰在 0 的卷数 | 峰 lag 范围 |
| --- | --- | --- | --- | --- | --- | --- |
| 128 | C100 | 0.341 [0.326, 0.363] | 0.345 [0.294, 0.365] | 0.331 [0.282, 0.362] | 5 | 0…0 |
| 128 | W100 | 0.370 [0.321, 0.400] | 0.366 [0.316, 0.404] | 0.368 [0.314, 0.386] | 5 | 0…0 |
| 128 | W500 | 0.434 [0.349, 0.441] | 0.439 [0.340, 0.446] | 0.426 [0.333, 0.454] | 5 | 0…0 |
| 235 | C100 | 0.324 [0.304, 0.486] | 0.328 [0.277, 0.392] | 0.317 [0.272, 0.380] | 4 | -2…0 |
| 235 | W100 | 0.269 [0.246, 0.374] | 0.249 [0.239, 0.284] | 0.249 [0.234, 0.288] | 3 | -2…0 |
| 235 | W500 | 0.340 [0.321, 0.407] | 0.314 [0.294, 0.356] | 0.311 [0.304, 0.369] | 3 | -1…0 |
| 285 | C100 | 0.336 [0.222, 0.421] | 0.296 [0.257, 0.418] | 0.302 [0.265, 0.401] | 4 | 0…1 |
| 285 | W100 | 0.229 [0.200, 0.317] | 0.217 [0.212, 0.244] | 0.222 [0.215, 0.252] | 3 | 0…1 |
| 285 | W500 | 0.287 [0.240, 0.359] | 0.280 [0.217, 0.307] | 0.296 [0.222, 0.320] | 5 | 0…0 |

原坐标相邻差分的关系总体更小，D285 的 W100 在 3 mm/s 为 -0.022256，为反向结果；D235 的 W500 在 10 mm/s 为 0.003177，接近零。其余全部卷和配对没有省略：[Q2-04](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/adjacent_difference.csv>)。

31/101-frame 敏感性仅覆盖固定 16 cells 与三个代表量。相对 51-frame，代表 rho0 的变化为 -0.161 至 0.151，specificity 变化为 -0.154 至 0.336；需按卷与窗口读取，不能将主窗口的幅度推广为尺度不变。[Q2-05](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/detrend_sensitivity_changes.csv>)。

同身份、同一固定区域定义的新旧 C16 对比如下。这里先在各卷取 C16，再描述五流速；没有与“先对每个 cell 做跨流速共识再取区域中位数”混用。

| D | 旧版同定义 C16 | 新版同定义 C16 | 减弱卷数 | 增加卷数 |
| --- | --- | --- | --- | --- |
| 128 | 0.182 [0.171, 0.221] | 0.152 [0.142, 0.164] | 5 | 0 |
| 235 | 0.153 [0.126, 0.202] | 0.168 [0.122, 0.221] | 1 | 4 |
| 285 | 0.131 [0.092, 0.185] | 0.134 [0.084, 0.172] | 2 | 3 |

这个对照显示 D128 的 C16 在五卷均减弱，D235 为四卷增加、一卷减弱，D285 方向混合。按逐卷 C16 的五流速中位数，旧版 D128 高于 D235，新版变为 D235 高于 D128；这不是 RI 排序的变化。[Q2-06](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/geometry_version_direction_ledger.csv>)。

| 状态 | 新版判断 |
| --- | --- |
| 保持 | 固定 C16 的正方向在三个直径、全部流速保留；三个代表量 51-frame 的 specificity 均为正。 |
| 减弱 | D128 同定义 C16 相对旧几何全五卷降低；D285 raw 的广泛正关联在去趋势后明显收缩；原始相邻差分中存在接近零的结果。 |
| 反向／例外 | 全深度地图有负相关；D285、3 mm/s、W100 的相邻差分为负；部分代表量的峰在 −2、−1 或 +1 帧；旧 C16 的 D128/D235 汇总排序反转。 |
| 不可定义 | 主 6×20 地图、C16 和三个代表量的主摘要本批均可定义，NaN 数为 0；不足 3 对、常量或任一必需 lag 缺失的规则已经测试，不以 0 填补。lag 的 μm 或时间含义不可定义。 |

## 问题三：背景、负对照和定位敏感性如何约束解释

左右背景中心均按 X4±1.5D/dx 平移，使用相同形状和深度；7018 个原有效位置均满足双侧完整 FOV。real/L/R 先取共同掩码，再分别去趋势，因此本批 real-matched 与 real-all 的支持恰好相同。Combined BG 由左右像素及 fractional weights 拼接，局部 Z 使用加权中位数和 MAD，不是两侧 Z 的平均。覆盖见 [background_coverage.csv](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/background_coverage.csv>)。

15 卷的 C16_excess 均为正，范围 0.038–0.158；固定区域的 real-minus-L、real-minus-R 也均为正。但全图 1800 个 cell 有 540 个 excess≤0，其中固定 C16 内有 12/240 个。不能将区域中位数为正改写成“所有深度都超过背景”。[Q3-01](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/fixed_region_control_excess.csv>)、[Q3-02](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/fixed_region_control_cells.csv>)。

| D | 流速 | real C16 | L C16 | R C16 | C16 excess | median(real−L) | median(real−R) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 128 | 1 | 0.164 | 0.037 | 0.031 | 0.132 | 0.126 | 0.127 |
| 128 | 3 | 0.151 | 0.081 | 0.101 | 0.066 | 0.074 | 0.083 |
| 128 | 5 | 0.152 | -0.002 | 0.046 | 0.125 | 0.169 | 0.108 |
| 128 | 7 | 0.142 | 0.045 | 0.029 | 0.141 | 0.103 | 0.111 |
| 128 | 10 | 0.153 | -0.005 | 0.022 | 0.120 | 0.149 | 0.106 |
| 235 | 1 | 0.182 | 0.024 | 0.001 | 0.149 | 0.154 | 0.162 |
| 235 | 3 | 0.122 | 0.025 | 0.029 | 0.097 | 0.106 | 0.101 |
| 235 | 5 | 0.168 | 0.002 | 0.015 | 0.158 | 0.166 | 0.131 |
| 235 | 7 | 0.134 | 0.003 | 0.014 | 0.117 | 0.131 | 0.125 |
| 235 | 10 | 0.221 | 0.072 | 0.071 | 0.132 | 0.135 | 0.128 |
| 285 | 1 | 0.137 | 0.066 | 0.082 | 0.079 | 0.078 | 0.087 |
| 285 | 3 | 0.084 | 0.050 | 0.048 | 0.038 | 0.035 | 0.033 |
| 285 | 5 | 0.134 | 0.049 | 0.018 | 0.081 | 0.083 | 0.104 |
| 285 | 7 | 0.110 | 0.064 | 0.034 | 0.070 | 0.078 | 0.059 |
| 285 | 10 | 0.172 | 0.043 | 0.072 | 0.126 | 0.132 | 0.108 |

C16 excess 是逐 cell 相减后取 16 项中位数；表中左右各自的区域中位数不能相减来重构该列。系数差只是一项描述性诊断，不是背景校正后的正式 SV/coupling，也不是无偏因果效应。

左右 pseudo 的 90 个代表量相关和 specificity 本批也都为正，rho0 为 0.012–0.270，说明正的局部关联或正 specificity 本身并非 real ROI 独有。D128、1 mm/s、左侧 W500 的峰位于 +50 帧边界，rho0=0.097、specificity=0.120；原样保留边界峰，未扩大 lag 范围或推断范围外峰位。完整 L/R lag 与 matched specificity：[Q3-13](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/pseudo_whole_coupling.csv>)。

背景对比随深度减小，并在部分深 band 反向。以下为五流速卷内 Z 中位数的 median [min,max]；没有设置 Z>1/2/3 阈值。负 Z 表示该 ROI 像素中位数低于局部背景中位数，不能解释为负的正式 SV。[Q3-03](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/local_contrast_profiles.csv>)。

| D | source | tail100 | tail500 | t10 | t20 |
| --- | --- | --- | --- | --- | --- |
| 128 | 5.433 [3.941, 5.539] | 1.960 [1.900, 2.130] | 0.512 [0.501, 0.542] | 0.414 [0.398, 0.440] | 0.059 [0.039, 0.067] |
| 235 | 5.810 [5.609, 5.972] | 0.928 [0.897, 0.951] | 0.156 [0.152, 0.168] | 0.034 [0.028, 0.054] | -0.069 [-0.078, -0.059] |
| 285 | 5.777 [4.374, 5.781] | 0.853 [0.667, 0.889] | 0.175 [0.122, 0.188] | 0.060 [0.040, 0.077] | -0.042 [-0.046, -0.020] |

D235 和 D285 的 T20（475–500 μm）局部 Z 在五个流速下均为负；D128 的末 band 局部对比仍为正但接近零。本轮没有据此生成仪器噪声底或“检测深度”。同深度 real-tail 与 combined BG-tail 的 51-frame rho 在全 20-band、全流速中的范围为 D128 −0.112 至 0.118、D235 −0.152 至 0.117、D285 −0.141 至 0.099，幅度接近零且包含反向。Crossed controls 同样存在正负方向：例如 D235 的 W500，左侧五卷均为负（−0.087 至 −0.029），右侧五卷均为正（0.014 至 0.110）。全部 raw/51 和 L/R 值见 [Q3-04](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/background_real_tail_covariation.csv>)、[Q3-05](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/crossed_controls.csv>)。

背景本身沿空间坐标变化。Combined 背景的 source、tail100、tail500 固定末端块中位数相对起始块均下降。三者相对差范围分别为 -15.5% 至 -2.9%、-22.3% 至 -1.0%、-19.9% 至 -2.8%。分母为各自整条背景中位数曲线的中位数；原端点范围为 0–50 与 449–499（均含端点），实际支持数随表提供。[Q3-06](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/background_drift.csv>)。

全部四个单轴平移都已计算：横向 ±1 像素（±12.7 μm），轴向 ±1 像素（±6.7 μm）。使用 baseline 与四偏移共同支持后，RI500 的卷内中位数变化为 -5.47% 至 6.93%；RI100 为 -7.90% 至 9.90%。C16 的变化范围为 -0.032 至 0.022，四方向各卷 C16 仍为正。[Q3-07](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/position_sensitivity_metric_changes.csv>)、[Q3-08](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/position_sensitivity_pairs.csv>)。这是小幅系统位置偏移的敏感性，不是定位误差真值，也未验证逐帧随机抖动。

在同一固定平移方向下继续匹配流速，四个方向的全部五个流速均保留 RI500 的 D128 > D235 > D285 排序；幅度变化与排序保持分别记录。[Q3-14](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/position_sensitivity_metrics.csv>)。

几何共享像素检查发现每卷 source 的最下层与 T1 存在共享像素，故按预定触发 A 补充 S3–S6×T2–T4 的 C12，以及 pooled S3–S6 对 tail25–100 的 51-frame 相关。C12 为 0.089–0.208，pooled 相关为 0.205–0.425，全部为正。该结果不等于排除了 PSF 串扰。[Q3-09](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/conditional_A_summary.csv>)。Neff、非零像素数、sum(w²)、层轴向像素跨度和解析面积误差均已保存；Neff 不是独立散斑数。

45 帧结构/SV 并列图均已审核。结构图存在轴向亮度梯度，完整 source 轮廓和独立下壁不足以判定，因此按合同保留“图像不足以判断”。D285 的 1、3 mm/s 首帧中 source 下缘附近的结构亮区触发了 B，补看最大相邻 ΔX4 与 Δz_top 的四组帧对；这些图仍未提供下壁真值，未修改定位或筛帧。[Q3-11](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/localization_review.csv>)、[Q3-12](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/conditional_B_review.csv>)。

触发 D 后，仅补算 real/L/R 的 source、tail100、tail500 structural mean。Real ROI 的同区域 structural–SV 51-frame rho 范围为 0.543–0.803，表明所测 SV 与结构变化共同出现；不能将这种共变量化为机制贡献比例。左右 pseudo 的对应值和结构 source–tail 相关也全部交付。[Q3-10](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/conditional_D_structural_correlations.csv>)。

C 的两侧反向判断条件未触发，图像亦未发现足以触发的单侧背景结构失配，因此未增加 ±2.5D 对照。报告不依赖精确单一 source bin 的位置，E 未触发。各项 trigger、帧 ID、范围与完成状态见 [conditional_extension_log.csv](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/conditional_extension_log.csv>)。

| 状态 | 新版判断 |
| --- | --- |
| 保持 | 主双侧匹配背景完整；C16 的区域中位超额在全部 15 卷为正，四个固定方向的位置扰动下 C16 也为正。 |
| 减弱 | 单个 cell 并非普遍超出对照；深部局部对比趋近零；窗口和位置会改变幅度，背景漂移与结构共变限制了仅凭正相关作出的解释。 |
| 反向 | 540 个全图 excess≤0，其中 12 个位于 C16；D235/D285 最后 band 的局部 Z 在五流速均为负。 |
| 不可定义 | 原 geometry 无效位置保留缺失；本批有效像素区域未出现额外背景 MAD=0。独立下壁、PSF、噪声底、独立散斑数、校准检测深度及因果效应均不能由当前数据定义。 |

## 追溯、验证与完整交付

每个正文证据编号对应 [answer_evidence.csv](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/answer_evidence.csv>) 的表名、SHA256、筛选条件及 CSV 行号；对应原行的完整值另存 [answer_evidence_rows.csv.gz](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/answer_evidence_rows.csv.gz>)。全覆盖矩阵见 [coverage_ledger.csv](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/coverage_ledger.csv>)，三个问题的完成对照见 [core_question_coverage.csv](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/core_question_coverage.csv>)。

原输入和历史目录的前后 SHA 均一致：[input_hashes_after.csv](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/input_hashes_after.csv>)。正式 v2 输入身份、运行环境和起始 HEAD 见 [provenance.json](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/provenance.json>)；本轮没有 checkout、pull、reset、commit、push 或远端发布。

数值核查包括同 16×16 定义的 45 帧 real/pseudo/四偏移积分、分层重构、weighted median/MAD、原坐标 lag 与有效配对、常量和 ties、已知 shift、正比例缩放、共同掩码、相邻差分和 matched specificity。非零积分相对误差阈值 1e-10，rho 差阈值 1e-12，零参考单列。[independent_integration_checks.csv.gz](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/independent_integration_checks.csv.gz>)、[statistical_independent_checks.csv.gz](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/statistical_independent_checks.csv.gz>)、[control_position_independent_checks.csv.gz](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/control_position_independent_checks.csv.gz>)。

过程中发现并修复了加权半质量边界的浮点累加差异：边界情况用精确二进制权重算术判断，独立 Fraction 实现复核。首次失败的本轮检查点保留在 failed_weighted_roundoff_checkpoints/，不属于有效科学交付；修正后全量像素统计和 45 帧检查通过。没有改变信号、阈值或分析合同。记录见 [weighted_numerical_fix.json](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/weighted_numerical_fix.json>)。

## 图表索引

- D128：[全部五卷 raw/51 地图](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/figures/depth_atlas_D128.png>)；[real-matched／L／R／excess 全图](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/figures/control_atlas_D128.png>)；[全部五卷三代表 lag](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/figures/lag_atlas_D128.png>)。
- D235：[全部五卷 raw/51 地图](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/figures/depth_atlas_D235.png>)；[real-matched／L／R／excess 全图](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/figures/control_atlas_D235.png>)；[全部五卷三代表 lag](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/figures/lag_atlas_D235.png>)。
- D285：[全部五卷 raw/51 地图](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/figures/depth_atlas_D285.png>)；[real-matched／L／R／excess 全图](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/figures/control_atlas_D285.png>)；[全部五卷三代表 lag](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/figures/lag_atlas_D285.png>)。

[匹配流速强度和 RI](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/figures/matched_flow_metrics.png>)；[全 20-band 强度／RI／局部 Z](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/figures/depth_profiles.png>)；[五流速系数中位共识图](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/figures/consensus.png>)；[固定区域对照](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/figures/fixed_region_controls.png>)；[四方向敏感性](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/figures/position_sensitivity.png>)。所有统计图另有同名 PDF；地图相关色轴固定 [−1,1]，excess 为 [−2,2]。

原始曲线、各模式实际窗口支持数、完整逐 cell 值、可定义原因、负值表及所有单侧背景均在对应 CSV 中。[figure_manifest.csv](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/figure_manifest.csv>) 提供图与输入表映射。完整输出文件身份见 [output_sha256.csv](<C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2/analysis/sv_fixed_diameter_selective_v2/output_sha256.csv>)。

复现顺序：initialize.py（仅空的新输出目录）→ verify_statistics.py → analyze_tables.py → extract_pixels.py → render_localization.py → analyze_controls.py → 对 45 帧完成并记录图像审核 → conditional_structural.py → 完成触发 B 的图像审核 → summarize_results.py → make_figures.py → validate_delivery.py → build_report.py → finalize_delivery.py。record_visual_review.py 和 finish_visual_review.py 保存的是本次实际图像审核记录；更换输入后必须重新看图，不能把这些标签当成自动定位判定。现有成功输出不应自动覆盖；需要复跑时使用另一个明确命名的输出副本。
