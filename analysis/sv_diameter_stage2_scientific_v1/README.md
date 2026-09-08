# SV 人工血管拖尾直径分析 — Stage 2 scientific analysis v1

本目录只分析已经冻结的 no-background raw-SV 结果，不修改 Stage 1 定量定义、连续性优先 v2.1 几何或 D128 冻结文件。独立实验单位仍然是 **scan volume**；B-scan 只表示同一 volume 内的 slow-axis 空间位置。

## 1. 最主要的现象

在 1/3/5/7/10 mm/s 五个共同流速下，volume-median RI_tail 再取中位数依次为：D128 **0.259467**、D235 **0.175291**、D285 **0.160953**、D500 **0.257642**。

也就是说，128→235 um 时 RI_tail 下降 **32.44%**，235→285 um 再下降 **8.18%**。第一段下降明显更大。D500 相对 D285 又回升 **60.07%**。这个结果仍然是描述性扫描卷比较，不把它写成直径的因果效应。

## 2. D500 为什么 RI 又高起来

D500 的 vessel raw SV 相比 D285 低 **45.60%**，500 um absolute tail raw SV 也低 **11.92%**。因此 D500 并不是“绝对拖尾变强”。它的分母下降得更明显，使 tail/source 这个比值重新升高。

逐个共同流速直接比较时，D500 相对 D285 的 source 变化中位数为 **-45.29%**，tail 变化中位数为 **-11.22%**，RI 变化中位数为 **+60.89%**。这说明 D500 的 RI 回升主要要从低 vessel denominator 来理解，同时 numerator 本身没有升高。

## 3. 这种直径排序是不是只出现在 volume 的某一小段

把每个 volume 的 500 个 slow-axis 位置固定分成五段后，共得到 **25** 个“共同流速 × 空间段”比较。D128 > D235 > D285 的 RI 排序出现在 **25/25** 段；其中 128→235 的下降幅度大于 235→285 的段数为 **25/25**。

D500 RI 高于 D285 出现在 **25/25** 段；D500 vessel source 低于 D285 出现在 **25/25** 段。这里的‘段’仍然只是空间稳健性检查，不是新的独立重复。

## 4. 限制几何和定位来源以后，主现象还在不在

预先定义的 6 个子集（全部有效帧、direct-only、X1 中央 80%、z_top 中央 80%、X1+z_top 中央 80%、direct+central geometry）中，128→235→285 的 RI 递减保留在 **6/6** 个子集；D500 相对 D285 的 RI 回升保留在 **6/6** 个子集。

这一步的含义很直接：如果删掉几何极端位置、只保留直接候选支持帧后排序仍然存在，那么主要现象就不容易用‘某个直径组恰好收进了更多异常定位帧’来解释。具体每个子集的数值见 `restricted_subset_summary.csv` 和 `restricted_subset_contrasts.csv`。

## 5. 直径差异随拖尾深度怎么变化

对每个固定流速先比较 RI(r)，再对五个共同流速的百分比差异取中位数：在 vessel bottom 附近（0 um），128→235 为 **-40.27%**、235→285 为 **-10.55%**；到 500 um 深度分别为 **-17.93%** 和 **-5.53%**。

D500 相对 D285 的 RI(r) 差异在 0 um 为 **+17.20%**，在 500 um 为 **+93.60%**。深部仍然存在相对差异时，应理解为 raw-SV floor 与 vessel denominator 共同决定的 normalized floor；这里没有做背景扣除，所以不能把深部非零 RI(r) 当成拖尾终点或 detection depth。

## 6. Flow 的描述性结构

- D128: 五个共同流速的 RI_tail 范围为 **0.256329–0.323180**，相对于该直径跨流速中位数的范围约 **25.8%**；最低点在 3 mm/s，最高点在 7 mm/s。
- D235: 五个共同流速的 RI_tail 范围为 **0.169371–0.177773**，相对于该直径跨流速中位数的范围约 **4.8%**；最低点在 10 mm/s，最高点在 1 mm/s。
- D285: 五个共同流速的 RI_tail 范围为 **0.158516–0.162519**，相对于该直径跨流速中位数的范围约 **2.5%**；最低点在 5 mm/s，最高点在 1 mm/s。
- D500: 五个共同流速的 RI_tail 范围为 **0.234979–0.265141**，相对于该直径跨流速中位数的范围约 **11.7%**；最低点在 1 mm/s，最高点在 5 mm/s。

不同直径的最高/最低流速位置并不统一，因此当前数据没有显示一个可以跨直径直接概括成‘随 flow 单调增加/降低’的共同结构。每个 diameter×flow 只有一个 independent volume，这里不进行 flow 的显著性检验或因果解释。

## 7. D128 100/200/300 um 窗口如何补齐

D128 历史 no-background 主表只保存了 500 um scalar，但冻结 release 保存了逐帧二维 `sv_raw`。本 Stage 2 仅在新目录中重新读取这些数组，用同一 X1/X4/z_top、真实 128 um 直径和同一 `quantify_raw_sv` 计算 100/200/300/500 um；不覆盖任何 D128 文件。

2422/2422 帧完成只读 replay。与冻结 D128 500 um 结果相比，source/tail/RI 最大相对误差分别为 **4.859e-16 / 4.119e-16 / 1.304e-15**，均低于 1e-10。

## 8. 当前可以支持的 Stage 2 判断

1. **128→235→285 的相对拖尾下降是一个跨共同流速、跨 slow-axis 空间段并对几何/QC 限制具有稳健性的描述性模式。**
2. **主要下降集中在 128→235；235→285 是较小的进一步下降。**
3. **D500 的 RI 回升不代表 absolute tail raw SV 更高；低 vessel raw SV denominator 是核心组成因素。**
4. **直径差异在靠近 vessel 的区域更容易拉开，深部逐渐受到 raw-SV floor / vessel denominator 的影响。**
5. **不同直径下的 flow 变化形态不完全一致，当前不建立统一 Flow×Diameter 模型。**
6. **几何/QC 组成不能单独解释主要直径排序。D500 denominator 低值若在所有限制分析中持续存在，应把它视为需要进一步做 source-ellipse 内部信号分布审计的信号层特征，而不是通过重新设计 source ROI 去消除。**

## 9. 输出文件

- `common_grid_volume_metrics.csv`: 4×5 主网格的 source / 100–500 um tail / RI。
- `diameter_contrasts_by_flow.csv`: 每个固定流速下的直径差值和百分比差。
- `flow_variation_by_diameter.csv`: 每个固定直径在五个共同流速中的变化范围。
- `source_tail_ratio_decomposition.csv`: source、absolute 500 um tail、RI 并排。
- `depth_profile_summary.csv`, `depth_selected_common_grid.csv`, `depth_diameter_contrasts.csv`, `depth_band_summary.csv`: RI(r) 深度结果。
- `slow_axis_segment_summary.csv`, `slow_axis_robustness_details.csv`: volume 内空间稳健性。
- `geometry_qc_summary.csv`: X1/X4/z_top/source area 与 source/tail/RI 的描述性关系。
- `restricted_subset_summary.csv`, `restricted_subset_diameter_summary.csv`, `restricted_subset_contrasts.csv`: 预定义限制子集。
- `d500_denominator_audit.csv`: D500 denominator 的 flow/segment/geometry/QC 审计。
- `d500_extra_flow_extension.csv`: D500 的 2/9/12 mm/s 次级扩展。
- PNG figures、`validation.json`、`provenance.json`。

Validation status: **passed**. No t-test, ANOVA, p-value, mixed model, regression surface, nonlinear fit or detection-depth estimate was computed.
