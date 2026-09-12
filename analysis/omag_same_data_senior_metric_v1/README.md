# 同数据 OMAG × 师兄指标

共同 preprocessing identity gate：**PASS**。15 卷 raw OCT 的前后 SHA 与正式来源一致；45 帧重建 SV 与 retained SV 逐像素完全相同（最大绝对误差 0）。同一 complex IMG 的独立读取检查及 raw OMAG 固定帧检查通过。

O-A absolute AUC 与 SV-A 的相邻直径方向在 **0/5** 个 flow 不同，以下为完全相同 geometry、endpoint-specific common support 的主比较。

| Flow (mm/s) | SV-A raw AUC | OMAG O-A raw AUC | SV-A P95 | O-A P95 | SV-A normalized | O-A normalized |
|---|---|---|---|---|---|---|
| 1 | D128>D235>D285 | D128>D235>D285 | D128<D235<D285 | D128<D235<D285 | D128>D235>D285 | D128>D235>D285 |
| 3 | D128>D235>D285 | D128>D235>D285 | D128<D235<D285 | D128<D235<D285 | D128>D235>D285 | D128>D235>D285 |
| 5 | D128>D235>D285 | D128>D235>D285 | D128<D235<D285 | D128<D235<D285 | D128>D235>D285 | D128>D235>D285 |
| 7 | D128>D235>D285 | D128>D235>D285 | D128<D235>D285 | D128<D235<D285 | D128>D235>D285 | D128>D235>D285 |
| 10 | D128>D235>D285 | D128>D235>D285 | D128<D235>D285 | D128<D235<D285 | D128>D235>D285 | D128>D235>D285 |

预先冻结趋势的满足流速数：

| 预定义形态 | O-A | O-B |
|---|---|---|
| raw AUC: D128<D235>D285 | 0/5 | 0/5 |
| P95: D128<D235<D285 | 5/5 | 5/5 |
| normalized: D128>D235>D285 | 5/5 | 5/5 |

以上统计采用各阶段 native support；完整排序（包括非相邻直径的次序）、native/common support 和逐算法直径百分比变化保存在 comparison/。

## 判读

OMAG O-A 的 absolute AUC 也在全部五个流速呈 D128>D235>D285。师兄历史中间峰不能由 OMAG algorithm alone 在当前数据中重现。D128 与历史 D185 的直径范围、原始采集及其他 historical dataset differences 仍需区分；本轮未检验这些因素。
O-A → O-B 在 common support 下改变 0/5 个 flow 的 absolute AUC 相邻直径方向。
O-A → O-B 在 common support 下改变 0/5 个 flow 的 P95 denominator 相邻直径方向。
O-A → O-B 在 common support 下改变 0/5 个 flow 的 normalized AUC 相邻直径方向。

O-B 与 SV-B 的比较包含信号算法、算法依赖的定位和纳入变化，属于 full-pipeline comparison。O-A/O-B 的 common-support 结果用于检查定位对 denominator 排序的影响，不把覆盖变化和算法差异混成纯 signal algorithm effect。

Normalized 方向相同不能证明算法等效；应同时读取 numerator 与 P95 的方向及各自的直径百分比变化。本轮不计算 OMAG AUC/SV AUC，也不以跨算法绝对单位的比值描述“高多少倍”。

## 四条 P95 轨迹与定位

| Flow | SV-A P95 | O-A P95 | SV-B P95 | O-B P95 | O-A common P95 | O-B common P95 |
|---|---|---|---|---|---|---|
| 1 | D128<D235<D285 | D128<D235<D285 | D128<D235<D285 | D128<D235<D285 | D128<D235<D285 | D128<D235<D285 |
| 3 | D128<D235<D285 | D128<D235<D285 | D128<D235<D285 | D128<D235<D285 | D128<D235<D285 | D128<D235<D285 |
| 5 | D128<D235<D285 | D128<D235<D285 | D128<D235<D285 | D128<D235<D285 | D128<D235<D285 | D128<D235<D285 |
| 7 | D128<D235>D285 | D128<D235<D285 | D128<D235<D285 | D128<D235<D285 | D128<D235<D285 | D128<D235<D285 |
| 10 | D128<D235>D285 | D128<D235<D285 | D128<D235<D285 | D128<D235<D285 | D128<D235<D285 | D128<D235<D285 |

## 数据、计算与核验

O-A 保留 7500 nominal positions、原 geometry-valid 7018、原 invalid 482。有效 AUC/P95/normalized pooled 帧数为 7018/7018/7018。O-B 对应为 7500/7500/7500。0 AUC 保留，P95 必须 finite 且 >1e-12；不因低信号或趋势新增排除。

15 卷 OMAG 采用已验证 retained export 的 raw 第二输出无损复用：每个来源文件 SHA、frame index、OMAG NumEV、共同预处理元数据均核验，并由每卷三帧原始 OCT 独立重建校验。大型中间数组为 float64 NPY，逻辑顺序 z,x,frame，frame=0…499；SHA 和路径见 omag_reconstruction/。它们未进入 Git。没有声称本轮从 OCT 重新导出了全部 7500 帧。

O-A 直接调用上一轮冻结 metric；O-B 使用原 senior tracking/X4/body assessability（alpha=.15、阈值≥.60、uncertain exclusion、原 numerical floors），定位工作数组按原 loader 转 float32，无 SV-specific scale mapping。指标保留原 raw float64，与已有 SV-A/SV-B adapter 的精度一致。原始 senior exporter 的 single 存储转换没有施加到无损 retained raw 上；这一点在合同中明确，避免额外量化成为主比较的变化项。

中央/单侧背景/gap 宽度分别为 4/8/9、3/6/7、2/5/6 px。先 round(D/dx)，z_top 按 Python ties-to-even；D/dz 为19/35/43，core offsets [2,9)/[4,16)/[5,20)，P95 为 linear quantile。P、合并双侧 B 按深度取 median，E=max(P−B,0)；AUC=[tail_start+1,tail_start+30)，29×6.7=194.3 μm。Normalized=AUC/P95，单位 μm。

核验通过：45 帧共同预处理、45 帧 OMAG 独立重建、O-A/O-B 指标独立复算 45/45 帧、90 个 volume endpoint medians、135 个 endpoint-specific common-support 组合；15 卷原始 OCT 前后 SHA 不变，495 个冻结文件前后 SHA 不变。所有 7500 个 OMAG frame finite。

这是同数据、同共同 preprocessing、同 senior metric 的 algorithm comparison within this dataset。每个 diameter×flow 只有一卷，空间帧不作为独立重复；不做帧级显著性检验或因果机制归因。历史最小直径为 D185，本轮为 D128，仅可比较 qualitative middle-diameter-peak shape，不能称 exact replication。

## 复现

在本目录依次运行 `python -B prepare.py`；MATLAB R2023a `rebuild_fixed`；`python -B validate_gate.py`；`python -B run_omag.py`；`python -B compare_report.py`；`python -B finalize.py`。原始路径和 retained 路径记录于 manifests；需在新的空输出副本复现，脚本拒绝覆盖完成的中间数组/汇总。依赖为 Python/NumPy/pandas/SciPy/Matplotlib 与 MATLAB 原函数；版本及源文件 SHA 见 provenance/contract。

BLOCKED/MISSING：无。
