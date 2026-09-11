# SV 固定真实直径 ROI v2：结果解释与后续定义迁移分析说明

日期：2026-09-11

本文件记录对 `analysis/sv_fixed_diameter_selective_v2/` 冻结结果的科学解释，并明确下一步与师兄拖尾指标做“定义迁移”比较的目的。它不修改原始结果、正式 SV 定义或既有主分析合同。

## 一、当前结果的核心解释

### 1. 直径效应是目前最稳定的主发现

在 D128、D235、D285 三种真实直径和 1/3/5/7/10 mm/s 五个匹配流速下，`tail100`、`tail500`、`RI100`、`RI500` 均满足：

`D128 > D235 > D285`。

这说明小直径血管不仅具有更高的 tail/source 相对比值，绝对 tail raw-SV 本身也更高，因此不能把 RI 的下降仅解释为 source 分母增大。

Source 本身不满足单一的直径排序：D128 通常低于 D235/D285，D285 相对 D235 在 1/3/5/7 mm/s 较高、10 mm/s 较低。因此当前数据不支持“直径越大，所有 SV 信号都越低”的概括。

四个固定 ±1 pixel 位置扰动下，RI500 的 `D128 > D235 > D285` 排序仍在全部五个匹配流速中保持，说明该主排序对小幅系统定位偏移较稳健。

### 2. Source 与 tail 存在有限但稳定的局部空间对应

在同一 scan volume 的 slow-axis 上，source 与其正下方 tail 的局部空间起伏存在正向对应。固定候选区域 `S3–S6 × T1–T4` 的 C16 在 15 卷中均为正，51-frame 去趋势后的范围约为 0.084–0.221。

三个代表量 `C100`、`W100`、`W500` 的 51-frame lag specificity 均为正。D128 的三个代表量在五个流速中峰值均位于 lag=0；D235/D285 的部分峰位于 -2、-1 或 +1 帧，整体集中在零位附近而非远端错位。

这一结果支持“同位置附近存在局部空间对应”，不代表 source 决定全部 tail，也不代表每个相邻 B-scan 都同步增减。相邻差分结果总体更弱，并出现接近零和负值，说明该对应主要体现为若干帧尺度的局部空间模式。

完整 6×20 地图中仍存在大量负相关 cell，因此不能把固定候选区域的正方向推广为全深度普遍正相关。

### 3. 匹配背景支持真实区域存在额外 coupling，但该特征并非全深度独有

15 卷的 C16 区域中位 excess 均为正，说明真实血管区域的固定近端 source-tail 关联整体强于左右匹配 pseudo control。

同时，完整地图中仍有 540 个 cell 的 real-minus-control excess ≤ 0，固定 C16 区域内也有少量非正 cell；pseudo ROI 自身也可出现正相关和正 specificity。因此“相关为正”本身不能作为真实拖尾的充分证据。

### 4. 深部 raw SV 非零不能解释为拖尾可检测到该深度

局部背景标准化对比随深度下降。475–500 μm 的最后一个 band 中，D235 和 D285 的局部 Z 在五个流速中均为负，D128 接近零。

因此 `tail500`/`RI500` 是固定 0–500 μm 窗口的描述性定量指标，不是校准 detection depth。非零 raw SV 不等于该深度存在可与局部背景区分的真实拖尾。

### 5. Raw SV 与结构 OCT 的共变是解释 coupling 时的重要限制

补充结构诊断中，real ROI 的 structural–SV 51-frame rho 约为 0.543–0.803，说明 raw SV 与局部结构信号存在明显空间共变。

因此当前观测到的 source-tail coupling 可包含 vessel-related、structure-related 及 local-background/spatial components。现有数据不能把这些成分定量分解，也不能据 Spearman rho 推导机制贡献比例。

## 二、对当前证据强度的排序

1. **最强、最稳健：** 128–285 μm 范围内 absolute tail 与 RI 的直径排序 `D128 > D235 > D285`。
2. **次一级但有明确空间证据：** source 与 tail 的零位/近零位局部空间对应。
3. **需要背景和结构约束解释：** coupling 的精确幅度、深度分布及直径间排序。
4. **当前不能定义：** 独立下壁真值、PSF、仪器噪声底、校准 detection depth、独立散斑数和因果机制贡献。

## 三、与师兄方法比较时已经确认的定义差异

师兄交付包的正式量化流程使用 OMAG `p_bld_ed`；当前项目主分析使用 `SV_raw = var(abs(IMG),1,3)`。

师兄冻结拖尾指标的关键定义为：

- 横向 central ROI：真实直径 D 对应横向宽度的 40%；
- 每个深度对 central ROI 取横向 median，得到 profile；
- 左右背景 ROI 各宽约 D/3，与血管边缘间隔约 D/4；左右背景像素合并后取 depth-wise median；
- excess：`max(profile - background, 0)`；
- 操作性血管下边界：`z_upper + round(D/6.7)`；
- 主 guard：2 px；
- 主 absolute tail：从 tail start 起 0–200 μm 的 excess AUC；
- 主 denominator：血管上边界以下 0.10D–0.45D 核心区 excess 的 P95；
- secondary normalized endpoint：`raw AUC / P95 denominator`；
- scan-level 主值：符合师兄 assessability 定义且 AUC 有效帧的 pooled median。

因此师兄的 raw AUC、P95 denominator 和 normalized AUC 与当前 `tail mean`、`whole-source mean`、`RI=tail/source` 不是同一指标。

## 四、下一步比较的科学定位

下一轮计划保持：

- 同一批 SV 数据；
- 同一正式 SV 算法；
- 同一 D128/D235/D285 × 1/3/5/7/10 mm/s 15 卷；
- 不引入 D500；

仅把 **ROI/背景/guard/窗口/denominator/scan-level aggregation 等统计定义迁移为师兄冻结定义**，建立“SV signal + senior metric definitions”的交叉分析。

该分析的目标不是让结果复现师兄曲线，也不作为当前 v2 主分析的替代，而是判断：

1. 采用师兄指标后，SV 是否也出现类似的 raw AUC 直径形态；
2. P95 vessel-core denominator 是否随直径明显增大；
3. normalized AUC 是否呈随直径下降的趋势；
4. 如果结果与当前 `tail mean/RI` 不同，差异主要来自 metric definition 还是 SV signal 本身。

如果未来再在同一数据上补充 OMAG + 当前指标，即可形成完整的 `signal algorithm × metric definition` 交叉比较；当前下一轮只执行 `SV + senior metric definitions`。
