# SV-OCTA 数据质量 / 背景可检测性 / 负对照 QC 审计 v1

本目录为独立诊断层。正式 raw SV、continuity-first v2.1 geometry 及两个冻结 coupling 目录保持原样；没有背景相减、重算 formal coupling、自动排除或物理机制解释。独立实验单位为 diameter × flow × scan volume，B-scans 是卷内空间位置。

Starting HEAD：`cfac5eb929208195c3a6a65ce79c251ce1e59641`；branch：`analysis/sv-diameter-stage2`。已执行 fetch、checkout、pull --ff-only；未 merge main。最终提交 SHA 由 Git commit 与交付回复标识，避免自引用。

## 固定方法与统计约定

SV = var(abs(IMG),1,3)，分母 N，全部正式 mean/Q 使用 linear raw SV。Source 是冻结 X4 中心的 X1 × physical-D ellipse，dx=12.7 μm、dz=6.7 μm、16×16 supersample；tail 从 true physical bottom 起，guard=0，X1 宽，无 cone/spreading。左右中心只移动 ±1.5 X1_px，边到边间距 0.5 X1。连续几何严格平移，仍在固定图像子像素格点积分；分数像素平移会使 rasterized 有效面积略有变化，不能要求左右离散权重逐元素相同。每侧需完整 source 与 0–500 μm tail 在 FOV 内才有效，不挪回、不缩小、不因高 SV 删除。

ROI mean、median、MAD 与 population SD 均采用 fractional weights；weighted median 是累计质量过半处，恰好一半时取相邻值均值。Combined BG 拼接有效侧像素和权重，按总像素权重计算，非左右 Z 的平均。Z=(weighted real median−weighted BG median)/(1.4826 weighted BG MAD)；CNR-like=(mean difference)/weighted BG SD。SBR=real mean/BG mean，SV SBR dB=10log10；contrast=(real−BG)/(real+BG)。Z/CNR 零分母记 NaN，有计数；不加 epsilon。Asymmetry 保留 signed log ratio、absolute log ratio、absolute difference 及 2|L−R|/(L+R)，log epsilon 固定 float64 tiny。Z>0/1/2/3 的 fraction 同时给出可评价帧分母和全部有效 geometry 帧分母。Kish Neff=(sum w)^2/sum(w²)，所有 region/side 都保留。

6 source bins 严格沿用 frozen normalized bin boundaries；20 tail bands 为 25 μm 等宽。Pooled 0–100 μm 直接整块 ROI 重算 median/MAD/Z，不平均四个 Z。背景稳定性分别保存每个 ROI 的 left/right/combined、pixel mean/median 的 slow-axis curve。CV 用 population SD/mean；RCV=1.4826 MAD/median；dynamic range=(p95−p5)/median。Drift 为固定原始 0–50 与 449–499 两块的 median 差，另有相对值、每原始 frame slope、51-frame trend max−min 及 residual MAD，无 p-value。

所有 detrend 在原始 0…499 网格上做 centered 51-frame moving median，端点截短，min_periods=1，忽略邻域 NaN，缺失中心保持缺失。Pseudo mean–mean 6×20 的 Spearman、lag ±50、signed peak/tie、specificity 直接调用冻结数值函数，不运行其正式分析入口。Lag=Δ 为 Source(i) vs Tail(i+Δ)，不循环。Specificity=rho(0)−median(rho at 30≤|lag|≤50)。High region 与上一轮相同：正 rho cells 中最高 ceil(10% positive cells)，截止 ties 全纳入。Pseudo 自身 high region 与在 frozen real high region 上的值分别保存；real/control excess 比较使用同一 real region，避免不同 cells 的混淆。左右 rho 的 cellwise median 形成 control rho，Δrho 仅 diagnostic。

Geometry jitter 只在原始 indices 相差 1 时计算；|Δg| > median+5×1.4826 MAD 仅 flag。每直径五个 common flows 1/3/5/7/10 独立卷；D500 2/9/12 单独 extension，不进入 common 汇总。无 frame-level tests、p-values、mixed models 或复杂预测模型。

## 输入与验证

| metric | max_absolute_error | max_relative_error |
| --- | --- | --- |
| source_mean | 5.960464477539063e-08 | 3.031258819239973e-16 |
| source_q | 0.0078125 | 2.644866984456457e-16 |
| tail100_mean | 5.960464477539063e-08 | 3.8270081669490733e-16 |
| tail100_q | 0.001953125 | 2.648443246066629e-16 |
| tail500_mean | 2.9802322387695312e-08 | 4.958140143760652e-16 |
| tail500_q | 0.00390625 | 2.6413540709509414e-16 |

10,931 valid arrays（D128 NPZ=2,422，MAT=8,509）与上一轮 identity 完全一致；25 release ZIP SHA 全部核对。另审计 569 个 invalid-geometry raw arrays：D128 使用冻结 release NPZ SHA，其他 MAT 保存新 SHA receipt（这些无效帧不在既有 coupling identity 表中，未伪称其 SHA 为旧表验证）。全部 11,500 arrays 可用，共 2,018,250,000 pixels，finite=100%，NaN/Inf/negative=0。逐数组 shape、zero fraction、min、median、p95/p99/p99.9/max 保存在 raw audit CSV。

背景左右 coverage 均 10,931/10,931；最大 overlap weight=0.0；连续 translation 参数最大误差=4.55e-13；full-image ellipse spot-check 权重误差=0。独立 scalar lag spot-check=9292 条，最大 rho 差=4.44e-16，原始 index 配对计数通过。完整冻结文件 SHA 检查 135 项通过。

## Common-grid 数字汇总

所有单元格均为五卷值的 median [min,max]；五个原值另保存在 CSV。

| D | source Z | source SBR | BG RCV | BG relative drift | valid fraction |
| --- | --- | --- | --- | --- | --- |
| 128 | 3.716 [2.906, 3.741] | 7.039 [5.597, 7.418] | 0.084 [0.083, 0.093] | -0.125 [-0.146, -0.072] | 0.972 [0.936, 0.986] |
| 235 | 4.841 [4.656, 4.887] | 8.297 [8.068, 8.504] | 0.050 [0.044, 0.060] | -0.064 [-0.084, -0.033] | 0.896 [0.848, 0.922] |
| 285 | 4.504 [3.730, 4.554] | 8.254 [6.782, 8.343] | 0.059 [0.058, 0.063] | -0.096 [-0.121, -0.064] | 0.948 [0.906, 0.988] |
| 500 | 1.516 [1.330, 1.591] | 3.921 [3.635, 4.225] | 0.060 [0.057, 0.073] | -0.151 [-0.203, -0.145] | 0.980 [0.942, 1.000] |

Proximal-tail robust Z：

| D | 0–25 | 25–50 | 50–75 | 75–100 | pooled 0–100 |
| --- | --- | --- | --- | --- | --- |
| 128 | 2.333 [2.206, 2.537] | 1.833 [1.799, 2.018] | 1.554 [1.473, 1.689] | 1.301 [1.230, 1.365] | 1.699 [1.638, 1.852] |
| 235 | 1.596 [1.509, 1.611] | 1.268 [1.185, 1.303] | 0.992 [0.924, 0.999] | 0.766 [0.702, 0.776] | 1.119 [1.042, 1.124] |
| 285 | 1.436 [1.139, 1.500] | 1.126 [0.885, 1.140] | 0.872 [0.674, 0.890] | 0.683 [0.544, 0.707] | 0.997 [0.797, 1.030] |
| 500 | 0.397 [0.360, 0.422] | 0.293 [0.269, 0.306] | 0.206 [0.194, 0.217] | 0.152 [0.129, 0.180] | 0.258 [0.230, 0.272] |

Source-bin robust Z：

| D | S1 | S2 | S3 | S4 | S5 | S6 |
| --- | --- | --- | --- | --- | --- | --- |
| 128 | 1.787 [1.726, 1.935] | 4.023 [2.683, 4.226] | 3.828 [2.434, 3.915] | 3.835 [3.198, 4.052] | 4.074 [3.698, 4.237] | 4.566 [4.393, 4.779] |
| 235 | 6.347 [6.210, 6.904] | 10.222 [9.524, 10.564] | 6.531 [6.402, 6.619] | 4.221 [4.050, 4.358] | 3.129 [3.047, 3.211] | 2.178 [2.061, 2.252] |
| 285 | 6.532 [6.044, 8.124] | 10.401 [7.982, 10.480] | 6.561 [5.200, 6.578] | 3.999 [3.235, 4.082] | 2.799 [2.177, 2.964] | 1.935 [1.483, 2.004] |
| 500 | 6.428 [6.301, 8.819] | 3.712 [3.370, 3.916] | 1.418 [1.234, 1.484] | 0.945 [0.845, 0.989] | 0.763 [0.679, 0.804] | 0.612 [0.526, 0.655] |

Geometry：

| D | geometry | median |delta| pixel | p95 |delta| pixel | median |delta| um | flag fraction |
| --- | --- | --- | --- | --- | --- |
| 128 | X4 | 0.537 [0.512, 0.558] | 1.632 [1.514, 1.775] | 6.820 [6.507, 7.086] | 0.000 [0.000, 0.002] |
| 128 | X1 | 2.000 [2.000, 2.000] | 6.000 [5.000, 6.000] | 25.400 [25.400, 25.400] | 0.000 [0.000, 0.002] |
| 128 | z_top | 0.000 [0.000, 0.000] | 1.000 [1.000, 1.000] | 0.000 [0.000, 0.000] | 0.072 [0.071, 0.082] |
| 235 | X4 | 0.360 [0.331, 0.427] | 1.058 [0.930, 1.100] | 4.578 [4.205, 5.428] | 0.000 [0.000, 0.003] |
| 235 | X1 | 2.000 [2.000, 2.000] | 6.000 [6.000, 6.000] | 25.400 [25.400, 25.400] | 0.002 [0.000, 0.005] |
| 235 | z_top | 0.000 [0.000, 0.000] | 1.000 [1.000, 1.000] | 0.000 [0.000, 0.000] | 0.075 [0.056, 0.089] |
| 285 | X4 | 0.412 [0.391, 0.428] | 1.214 [1.170, 1.281] | 5.230 [4.961, 5.439] | 0.002 [0.000, 0.007] |
| 285 | X1 | 2.000 [2.000, 2.000] | 6.000 [5.000, 7.000] | 25.400 [25.400, 25.400] | 0.004 [0.002, 0.007] |
| 285 | z_top | 0.000 [0.000, 0.000] | 1.000 [1.000, 1.000] | 0.000 [0.000, 0.000] | 0.071 [0.053, 0.089] |
| 500 | X4 | 0.770 [0.744, 0.863] | 2.235 [2.031, 2.407] | 9.782 [9.446, 10.963] | 0.000 [0.000, 0.002] |
| 500 | X1 | 2.000 [2.000, 2.000] | 6.000 [6.000, 7.000] | 25.400 [25.400, 25.400] | 0.006 [0.000, 0.012] |
| 500 | z_top | 0.000 [0.000, 0.000] | 1.000 [1.000, 1.000] | 0.000 [0.000, 0.000] | 0.066 [0.062, 0.076] |

Negative controls（51-frame detrended）：

| D | real high | L at real high | R at real high | cellwise excess high | real > both fraction | BG–tail rho |
| --- | --- | --- | --- | --- | --- | --- |
| 128 | 0.237 [0.236, 0.286] | 0.063 [0.011, 0.152] | 0.050 [0.003, 0.117] | 0.204 [0.096, 0.272] | 1.000 [0.909, 1.000] | 0.027 [0.012, 0.031] |
| 235 | 0.210 [0.169, 0.245] | 0.027 [0.020, 0.151] | 0.045 [-0.010, 0.109] | 0.156 [0.121, 0.187] | 1.000 [0.875, 1.000] | -0.014 [-0.017, 0.021] |
| 285 | 0.167 [0.156, 0.224] | 0.040 [-0.005, 0.060] | 0.055 [0.049, 0.083] | 0.148 [0.124, 0.165] | 1.000 [1.000, 1.000] | -0.006 [-0.019, -0.001] |
| 500 | 0.117 [0.116, 0.148] | 0.034 [0.017, 0.071] | 0.044 [0.025, 0.065] | 0.098 [0.073, 0.115] | 1.000 [0.857, 1.000] | 0.007 [-0.019, 0.015] |

QC → real high-region rho 的 descriptive Spearman：

| predictor | D128 | D235 | D285 | D500 |
| --- | --- | --- | --- | --- |
| source_z | -0.200 | 0.700 | 0.500 | 0.400 |
| proximal_z | -0.100 | 0.600 | 0.600 | 0.000 |
| background_rcv | -0.900 | 0.200 | 0.000 | 0.500 |
| z_top_jitter | NaN | NaN | NaN | NaN |
| X1_jitter | NaN | NaN | NaN | NaN |
| valid_fraction | 0.500 | 0.975 | 0.100 | 0.462 |
| negative_control_rho | 0.000 | 0.800 | 0.500 | 0.500 |

## Q1–Q15：纯数据答案

**Q1**　Whole-source robust Z 的五卷 median [min,max]：D128 3.716 [2.906, 3.741]；D235 4.841 [4.656, 4.887]；D285 4.504 [3.730, 4.554]；D500 1.516 [1.330, 1.591]。

**Q2**　D500 与 D128/D235/D285 的五个匹配 flow 比较，source Z 均更低（各 5/5）；具体值见 Q1。

**Q3**　直接从完整 0–100 μm ROI 计算的 pooled tail robust Z：D128 1.699 [1.638, 1.852]；D235 1.119 [1.042, 1.124]；D285 0.997 [0.797, 1.030]；D500 0.258 [0.230, 0.272]。四个 25 μm band 见下表。

**Q4**　将既有每卷 real high-region cells 对应的 tail-band volume Z 取中位数（重复 tail bin 按入选 cell 保留），结果为 D128 1.833 [1.636, 2.018]；D235 1.239 [0.998, 1.457]；D285 1.126 [0.885, 1.140]；D500 0.165 [0.068, 0.195]。这是相对局部背景的描述性分离度，Z>3 未经校准，不能据此称为“可靠检测”或统计显著。

**Q5**　T1（0–25 μm）为 D128 2.333 [2.206, 2.537]；D235 1.596 [1.509, 1.611]；D285 1.436 [1.139, 1.500]；D500 0.397 [0.360, 0.422]；T20（475–500 μm）为 D128 0.064 [0.048, 0.085]；D235 -0.026 [-0.032, -0.014]；D285 -0.011 [-0.032, 0.002]；D500 -0.064 [-0.069, -0.060]。完整 20 层见 Figure 2；不假定逐层严格单调。

**Q6**　S1–S6 的五卷 Z 见下表；各直径最高 median source-bin 为 D128 S6；D235 S2；D285 S2；D500 S1。

**Q7**　每卷最大 coupling cell 所在 source bin 与最大 source detectability bin 完全一致的卷数：D128 0/5；D235 0/5；D285 0/5；D500 0/5。两种位置不作等同或因果解释。

**Q8**　Source matched BG median curve 的 RCV：D128 0.084 [0.083, 0.093]；D235 0.050 [0.044, 0.060]；D285 0.059 [0.058, 0.063]；D500 0.060 [0.057, 0.073]；首尾 51 原始帧块 median 差 / 全卷 median：D128 -0.125 [-0.146, -0.072]；D235 -0.064 [-0.084, -0.033]；D285 -0.096 [-0.121, -0.064]；D500 -0.151 [-0.203, -0.145]。

**Q9**　Valid geometry fraction：D128 0.972 [0.936, 0.986]；D235 0.896 [0.848, 0.922]；D285 0.948 [0.906, 0.988]；D500 0.980 [0.942, 1.000]。三项 jitter 的 median 与 p95 见表，仅使用真正相邻原始 indices；没有以 QC 删除帧或卷。

**Q10**　每卷 20 个 real-tail / combined-BG-tail detrended rho 的中位数，再汇总五卷：D128 0.027 [0.012, 0.031]；D235 -0.014 [-0.017, 0.021]；D285 -0.006 [-0.019, -0.001]；D500 0.007 [-0.019, 0.015]。原始和 detrended 各 band/侧结果均保留。

**Q11**　Pseudo whole-map median rho：L D128 0.033 [0.028, 0.065]；D235 0.024 [0.013, 0.039]；D285 0.030 [0.012, 0.038]；D500 0.022 [0.017, 0.051]；R D128 0.024 [0.015, 0.067]；D235 0.025 [0.012, 0.051]；D285 0.028 [0.019, 0.036]；D500 0.022 [-0.003, 0.044]。Real-vs-L map Spearman similarity：D128 0.126 [-0.054, 0.492]；D235 0.126 [0.053, 0.397]；D285 0.037 [-0.071, 0.221]；D500 -0.066 [-0.100, 0.023]；Real-vs-R：D128 0.314 [0.186, 0.459]；D235 0.021 [-0.143, 0.341]；D285 0.157 [0.054, 0.215]；D500 0.050 [-0.008, 0.152]。

**Q12**　Real high region 的 median rho：D128 0.237 [0.236, 0.286]；D235 0.210 [0.169, 0.245]；D285 0.167 [0.156, 0.224]；D500 0.117 [0.116, 0.148]；同一组 cells 的逐 cell real-minus-control，再取 median：D128 0.204 [0.096, 0.272]；D235 0.156 [0.121, 0.187]；D285 0.148 [0.124, 0.165]；D500 0.098 [0.073, 0.115]。Real 大于双侧 controls 的 cell fraction：D128 1.000 [0.909, 1.000]；D235 1.000 [0.875, 1.000]；D285 1.000 [1.000, 1.000]；D500 1.000 [0.857, 1.000]。这些不是 background-corrected formal coupling。

**Q13**　七项 QC 对 real high-region rho 的 within-diameter descriptive Spearman 见表，每项 n=5；常量指标的 rho 为 NaN。其余三个 coupling outputs 和 pooled-20 exploratory coefficients 保存在 qc_vs_coupling_summary.csv；pooled 结果可能受 diameter 混杂。

**Q14**　Mixed evidence。D500 source 与 proximal-tail detectability 在五个 matched common flows 中均低于 D235/D285；geometry valid fraction、jitter、background stability 和 negative controls 需分别按表读取，不能合并成一个 SNR。较低局部可检测性是较弱观测 coupling 的可能贡献因素；本轮描述性数据不能量化其解释比例，也不能断言完全解释或无关。

**Q15**　本轮不建立 post-hoc exclusion threshold。Z>0/1/2/3 与 jitter robust flag 仅为预先固定的描述性标记，全部 23 卷保留。

## Structural OCT 与噪声参考

23 卷 retained arrays 均含共配准、同形状、全有限 stru_amp。该量按代码为 repeat amplitude mean，只报告 structural source-to-local-background contrast、local SBR 和 BG variability；structural ratio 未转换为 dB。没有在 retained cohort 的 manifest、README、signal implementation 或 metadata 中识别到有明确 provenance 的独立 dark/air/noise-floor reference。No calibrated instrument-noise-referenced OCT SNR could be computed from the available retained data。结构 QC 与 SV/coupling 关系仅为 same-cohort volume-level descriptive Spearman。

## 异常、缺失与范围边界

569 个 geometry 缺失在名义 11,500 行中保留；没有 QC 新增 exclusion。无 FOV coverage 缺失，无 raw nonfinite/negative。左右不对称按事实保留，没有按值删除“疑似其他血管”；asymmetry 本身不标定污染来源。零分母计数、constant-input NaN 和 Neff 全保留 CSV。仪器 noise reference 不可用，未计算 true instrument SNR；本轮不能从描述性 n=5 关系确定 measurement detectability 对 coupling 降低的因果解释比例。全部任务验证以 validation.json 为准。

## 文件与复现

逐帧文件按 D128/D235/D285/D500 与 D500_extension 分组 gzip；含名义帧、valid flag、全 source/bin/tail 的 real/L/R/BG mean/median/MAD/SD/Q/area/Neff、contrast、阈值、structural 指标。volume_region_detectability_summary 为所有 ROI 的汇总；用户指定的 source/tail/proximal/geometry/background/negative-control/relationship CSV 均单独输出。完整 lag 曲线按 scan×side×mode gzip；slow-axis raw/trend/residual 另有 gzip。Figure 1–10 同时提供 PNG/PDF，mean/median 统计不将 B-scans 当误差重复；所有 flow 全展示，coupling 色轴 [-1,1]，excess [-2,2]，同 metric 共享含零尺度。

在仓库根目录执行（需 input_manifest 内相同 SHA retained arrays；文件路径可按已有 identities 解析）：

```powershell
python analysis/sv_data_quality_qc_v1/derive_qc.py
python analysis/sv_data_quality_qc_v1/audit_additional_arrays.py
python analysis/sv_data_quality_qc_v1/analyze_qc.py
python analysis/sv_data_quality_qc_v1/make_figures.py
python analysis/sv_data_quality_qc_v1/verify_qc.py
python analysis/sv_data_quality_qc_v1/build_report.py
```

## 最重要的 12 条数据事实

1. Whole-source robust Z 的五卷 median [min,max]：D128 3.716 [2.906, 3.741]；D235 4.841 [4.656, 4.887]；D285 4.504 [3.730, 4.554]；D500 1.516 [1.330, 1.591]。

2. 直接从完整 0–100 μm ROI 计算的 pooled tail robust Z：D128 1.699 [1.638, 1.852]；D235 1.119 [1.042, 1.124]；D285 0.997 [0.797, 1.030]；D500 0.258 [0.230, 0.272]。四个 25 μm band 见下表。

3. Source matched BG median curve 的 RCV：D128 0.084 [0.083, 0.093]；D235 0.050 [0.044, 0.060]；D285 0.059 [0.058, 0.063]；D500 0.060 [0.057, 0.073]；首尾 51 原始帧块 median 差 / 全卷 median：D128 -0.125 [-0.146, -0.072]；D235 -0.064 [-0.084, -0.033]；D285 -0.096 [-0.121, -0.064]；D500 -0.151 [-0.203, -0.145]。

4. Valid geometry fraction：D128 0.972 [0.936, 0.986]；D235 0.896 [0.848, 0.922]；D285 0.948 [0.906, 0.988]；D500 0.980 [0.942, 1.000]。三项 jitter 的 median 与 p95 见表，仅使用真正相邻原始 indices；没有以 QC 删除帧或卷。

5. 每卷 20 个 real-tail / combined-BG-tail detrended rho 的中位数，再汇总五卷：D128 0.027 [0.012, 0.031]；D235 -0.014 [-0.017, 0.021]；D285 -0.006 [-0.019, -0.001]；D500 0.007 [-0.019, 0.015]。原始和 detrended 各 band/侧结果均保留。

6. Pseudo whole-map median rho：L D128 0.033 [0.028, 0.065]；D235 0.024 [0.013, 0.039]；D285 0.030 [0.012, 0.038]；D500 0.022 [0.017, 0.051]；R D128 0.024 [0.015, 0.067]；D235 0.025 [0.012, 0.051]；D285 0.028 [0.019, 0.036]；D500 0.022 [-0.003, 0.044]。Real-vs-L map Spearman similarity：D128 0.126 [-0.054, 0.492]；D235 0.126 [0.053, 0.397]；D285 0.037 [-0.071, 0.221]；D500 -0.066 [-0.100, 0.023]；Real-vs-R：D128 0.314 [0.186, 0.459]；D235 0.021 [-0.143, 0.341]；D285 0.157 [0.054, 0.215]；D500 0.050 [-0.008, 0.152]。

7. Real high region 的 median rho：D128 0.237 [0.236, 0.286]；D235 0.210 [0.169, 0.245]；D285 0.167 [0.156, 0.224]；D500 0.117 [0.116, 0.148]；同一组 cells 的逐 cell real-minus-control，再取 median：D128 0.204 [0.096, 0.272]；D235 0.156 [0.121, 0.187]；D285 0.148 [0.124, 0.165]；D500 0.098 [0.073, 0.115]。Real 大于双侧 controls 的 cell fraction：D128 1.000 [0.909, 1.000]；D235 1.000 [0.875, 1.000]；D285 1.000 [1.000, 1.000]；D500 1.000 [0.857, 1.000]。这些不是 background-corrected formal coupling。

8. 每卷最大 coupling cell 所在 source bin 与最大 source detectability bin 完全一致的卷数：D128 0/5；D235 0/5；D285 0/5；D500 0/5。两种位置不作等同或因果解释。

9. 20 common-grid volumes / 9,456 valid frames；3 extension volumes / 1,475 valid frames；共 11,500 原始位置中 569 个 geometry 缺失。

10. 全部 11,500 retained raw SV arrays 为 351×500，2,018,250,000 pixels；NaN=0，Inf=0，负值=0。

11. 10,931 个正式有效帧逐一 SHA256 通过；左右 matched ROI 在全部有效帧均完整覆盖，无 source/tail 重叠。

12. 23 卷 10,931 个有效帧均有同形状、有限的 retained stru_amp；只计算 structural local contrast/SBR，没有 calibrated noise-floor reference。

## 背景不对称与结构对照补充

| D | source BG pixel median | source BG pixel MAD | BG absolute log asymmetry | structural source local SBR |
| --- | --- | --- | --- | --- |
| 128 | 14695428.859 [14490209.730, 15311854.638] | 10270452.772 [10025416.354, 10548551.043] | 0.083 [0.076, 0.090] | 0.968 [0.917, 0.972] |
| 235 | 15129131.276 [14730006.674, 15557117.168] | 10566763.486 [10267577.844, 10852545.210] | 0.046 [0.044, 0.057] | 0.869 [0.869, 0.889] |
| 285 | 15265245.791 [14657857.786, 18566692.801] | 10693966.493 [10255687.291, 13084564.765] | 0.049 [0.046, 0.055] | 0.828 [0.819, 0.835] |
| 500 | 16911997.495 [15761031.191, 18465438.182] | 11880669.265 [11101613.998, 13008776.506] | 0.079 [0.064, 0.083] | 0.445 [0.435, 0.447] |

左右背景不对称完整保留，包括每卷 median/p90/p95；上述数值没有定位其来源，也没有据此认定或移除其他血管。D500 structural source/local-background SBR 低于其他三直径，且与更低 SV detectability 同时出现；这些是 depth-matched local contrast，不是噪声参考的仪器 SNR。noise-reference 的 metadata 检索范围是预先固定每卷最接近 index 249 的有效帧，共 23 个记录，加上 cohort manifests、README 和 source implementation；完整检索记录另存 CSV，不声称穷尽每个原始 acquisition metadata。

Common-grid 的 X1 median absolute jitter 全为 2 pixels，z_top median absolute jitter 全为 0；两项作为候选 predictor 在每直径五卷内恒定，相关系数不可定义（NaN），不是漏算。X1/z_top 的 p90/p95/p99/max/MAD 和 flag fraction 仍完整保留；z_top p95 为 1 pixel。零 MAD 时的 jitter flag 按原先固定公式执行，没有额外放宽。

最后执行 `python analysis/sv_data_quality_qc_v1/finalize_delivery.py` 核对交付范围及输出身份，生成 output_sha256.csv（排除其自身）；intermediate checkpoint 和缓存不提交。所有文本固定 UTF-8/LF，确保 Git checkout 后 SHA 一致。
