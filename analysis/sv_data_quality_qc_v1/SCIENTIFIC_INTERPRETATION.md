# SV-OCTA 数据质量 QC：科学解读与当前结论边界

本文件记录对 `sv_data_quality_qc_v1` 结果的科学解释，目的在于区分：

1. **数据事实**：计算直接得到的结果；
2. **统计/测量解释**：这些结果能支持什么；
3. **当前不能下的结论**：避免把相关性或 QC 指标误写成因果机制。

本文件不修改任何正式 raw SV、冻结 geometry 或既有 source–tail coupling 结果。

---

## 1. 本轮最核心的问题

前一轮二维 source-depth × tail-depth 分析发现：

- D128–D285 存在较稳定的局部 source–tail 空间共变；
- D500 的 detrended mean–mean coupling 明显减弱、跨 flow map 稳定性也下降。

因此本轮 QC 的关键不是简单问“哪个数据更好”，而是问：

> **D500 的弱 coupling，究竟可能有多少来自 measurement detectability 下降、背景波动、geometry jitter 或系统级共同变化？**

当前答案是：

> **D500 确实存在明显更低的 source 与 proximal-tail SV detectability，这很可能会削弱可观测 coupling；但 D500 并不是整卷数据全面失效，背景稳定性、有效帧率和 matched pseudo-control 并未出现同等级恶化，因此不能用“D500 数据质量差”简单解释全部 coupling 差异。**

---

## 2. robust Z 的正确含义

本轮主 QC 指标为：

\[
Z_{robust}=
\frac{Median(Signal)-Median(Background)}{1.4826\,MAD(Background)}
\]

它表示：

> **目标区域的 SV 信号相对于匹配局部背景，高出了多少个 robust background fluctuation scale。**

因此它更适合称为：

- local-background standardized contrast；
- SV detectability；
- robust signal-background separation。

不能直接称为 calibrated instrument SNR，因为本轮没有 provenance 明确的 dark/air/noise-floor acquisition。

同样，`Z>3` 只能作为描述性 QC 阈值，不能被解释成经过验证的 OCTA detection threshold。

---

## 3. Whole-source detectability：D500 明显下降

五个 common-flow volumes 的 whole-source robust Z 中位数：

| Diameter | Whole-source Z |
| --- | ---: |
| D128 | 3.716 |
| D235 | 4.841 |
| D285 | 4.504 |
| D500 | 1.516 |

关键事实：

- D235/D285 的 source detectability 最好；
- D500 在 5/5 matched flows 中均低于 D235/D285；
- 因此并不是“直径越大，质量逐渐线性变差”，而是到了 D500 出现明显下降。

数据层解释：

> **D500 的血管 SV 本体相对于 matched local background 明显更难区分。**

但这还不能单独证明 D500 coupling 弱是由 detectability 下降导致。

---

## 4. Proximal-tail detectability：随直径增大明显下降

Pooled 0–100 μm tail robust Z：

| Diameter | Tail 0–100 μm Z |
| --- | ---: |
| D128 | 1.699 |
| D235 | 1.119 |
| D285 | 0.997 |
| D500 | 0.258 |

分层结果：

| Diameter | 0–25 μm | 25–50 μm | 50–75 μm | 75–100 μm |
| --- | ---: | ---: | ---: | ---: |
| D128 | 2.333 | 1.833 | 1.554 | 1.301 |
| D235 | 1.596 | 1.268 | 0.992 | 0.766 |
| D285 | 1.436 | 1.126 | 0.872 | 0.683 |
| D500 | 0.397 | 0.293 | 0.206 | 0.152 |

两个非常稳定的模式：

1. **同一直径内，tail 越深，detectability 越低；**
2. **同一 proximal depth 下，D500 明显最弱。**

到 475–500 μm，四种直径的 Z 均已接近 0 或略低于 0。

因此必须继续坚持：

> **deep raw SV 非零，不等于该深度仍具有可靠的 local-background detectability。**

这也支持以后把“tail effective depth”与单纯的非零 RI 或非零 raw SV 区分开。

---

## 5. 近端 0–100 μm 为什么既是 coupling 热区，也是 QC 较好的区域

前一轮 source-depth × tail-depth 分析发现 D128–D285 的主要 coupling hot region 位于 proximal tail，尤其约 0–100 μm。

本轮 QC 又发现：

> **0–100 μm 同时也是 tail 最容易区别于 local background 的区域。**

因此这一区域是目前最适合研究 local source–tail covariance 的 tail 深度范围。

但不能据此说：

> “coupling 高只是因为 Z 高。”

因为 source-depth 层面的结果并不支持这种简单解释，见下一节。

---

## 6. Coupling hotspot ≠ detectability hotspot

各直径 source detectability 最高的 bin：

- D128：S6；
- D235：S2；
- D285：S2；
- D500：S1。

而上一轮最大 coupling source bin 与最大 detectability bin 的一致卷数，四种直径均为：

\[
\boxed{0/5}
\]

这非常重要。

如果 coupling topology 只是由“哪个层信号最清楚”决定，那么最高 detectability bin 应该经常与最高 coupling bin 重合。

实际并没有。

因此：

> **source-depth coupling pattern 不能简单归结为 source-depth SNR/detectability pattern。**

这增强了前一轮“source axial coupling topology 具有独立空间信息”的可信度。

---

## 7. D500 source 内部存在明显的轴向 detectability 梯度

D500 的 S1–S6 robust Z：

| Bin | Z |
| --- | ---: |
| S1 | 6.428 |
| S2 | 3.712 |
| S3 | 1.418 |
| S4 | 0.945 |
| S5 | 0.763 |
| S6 | 0.612 |

这说明 D500 不是“整个 source 都差”。

更准确的图景是：

> **上部 source 仍然很清楚，但越往中部/下部，SV detectability 快速衰减。**

这和当前科学问题有直接关系，因为前一轮我们最关心的正是中部至中下部 source 与 proximal tail 的局部 coupling。

因此：

> **D500 中下部较差的 measurement detectability，很可能会削弱我们观测到的 local source–tail coupling。**

但这只能称为 plausible contributor，不能称为完全解释。

---

## 8. D500 并不是“整卷数据全面变坏”

如果 D500 coupling 变弱只是因为 acquisition / pipeline 整体失效，我们应该看到多个 QC 指标同步恶化。

实际并不是这样。

### 8.1 Background RCV

| Diameter | Background RCV |
| --- | ---: |
| D128 | 0.084 |
| D235 | 0.050 |
| D285 | 0.059 |
| D500 | 0.060 |

D500 与 D285 基本相同。

因此：

> **D500 的主要问题不是 background fluctuation amplitude 特别高。**

### 8.2 Valid-frame fraction

| Diameter | Valid-frame fraction |
| --- | ---: |
| D128 | 0.972 |
| D235 | 0.896 |
| D285 | 0.948 |
| D500 | 0.980 |

D500 反而最高。

因此：

> **D500 弱 coupling 不能解释成“大量帧无效”。**

---

## 9. D500 的 slow-axis drift 更明显

Background relative drift：

- D128：−0.125；
- D235：−0.064；
- D285：−0.096；
- D500：−0.151。

D500 的 slow-axis background decline 更明显。

这与上一轮结果一致：

- raw source–tail correlation 看起来较强；
- 51-frame detrending 后 mean–mean coupling 大幅下降。

因此更合理的统计解释是：

> **D500 含有更明显的大尺度 slow-axis trend，所以 raw correlation 更容易被共同空间趋势抬高。**

这也再次证明 detrending 是必要的。

---

## 10. Background–tail 负对照接近 0

51-frame detrended 后，BG–real-tail 的 volume-level median rho：

| Diameter | BG–real-tail rho |
| --- | ---: |
| D128 | 0.027 |
| D235 | −0.014 |
| D285 | −0.006 |
| D500 | 0.007 |

基本都在 0 附近。

这说明：

> **real tail 的局部 fluctuation 并不是简单跟附近 matched background 一起同步变化。**

因此前一轮 real source–tail coupling 不太可能只是“整幅图同一帧一起亮/一起暗”的系统背景效应。

---

## 11. Matched pseudo-source / pseudo-tail 负对照支持 real spatial specificity

左右 pseudo-source / pseudo-tail 完全复制真实 ROI geometry，只做 lateral translation。

Pseudo map 的典型 detrended correlation 约为：

\[
\rho\sim0.02-0.03
\]

而 real high-region rho 中位数约为：

| Diameter | Real high-region rho |
| --- | ---: |
| D128 | 0.237 |
| D235 | 0.210 |
| D285 | 0.167 |
| D500 | 0.117 |

Real minus matched-control excess：

| Diameter | Control-excess rho |
| --- | ---: |
| D128 | 0.204 |
| D235 | 0.156 |
| D285 | 0.148 |
| D500 | 0.098 |

而 real high-region cells 中 rho 高于左右两个 matched controls 的比例，五卷中位数四个直径均为 1.000。

这提供了很重要的负对照证据：

> **真实 vessel 位置的 source–tail local covariance 明显强于几何完全匹配的 lateral pseudo controls。**

因此前一轮 coupling 不是简单由 ROI 形状、全局空间背景或算法本身自动制造出来的。

---

## 12. D500 即使较弱，仍有正的 real-minus-control excess

D500 real high-region rho 约为 0.117，低于 D128–D285。

但 D500 的 control-excess 仍约为：

\[
0.098
\]

并保持正值。

因此当前不能把 D500 描述成：

> “不存在任何 real source–tail spatial relationship。”

更准确的是：

> **D500 的真实 spatial correspondence 仍高于 matched pseudo controls，但其强度明显减弱。**

---

## 13. Geometry QC：D500 主要是 X4 lateral centre 更抖

Median |ΔX4|：

- D128：0.537 px；
- D235：0.360 px；
- D285：0.412 px；
- D500：0.770 px。

所以 D500 lateral centre 的 slow-axis variation 较大。

但其他 geometry 指标没有同等级恶化：

- X1 median |Δ|：所有 common volumes 均为 2 px；
- z_top median |Δ|：所有直径均为 0；
- z_top p95：均约 1 px = 6.7 μm。

因此：

> **D500 存在较大的 lateral-centre jitter，但不能说 frozen geometry 整体失控。**

---

## 14. Structural OCT 也提示 D500 local contrast 较低

Structural source/local-background SBR 中位数：

- D128：0.968；
- D235：0.869；
- D285：0.828；
- D500：0.445。

D500 同样最低。

但由于缺少 calibrated dark/air/noise-floor reference：

> **这只能称为 structural source-to-local-background contrast，不能称为 true instrument-noise-referenced OCT SNR。**

---

## 15. 当前如何重新理解 D500

目前最不应该写成：

\[
D500\rightarrow low\ SNR\rightarrow weak\ coupling
\]

因为我们没有证据确定单一因果链及其解释比例。

更准确的描述是：

- D500 source detectability 下降；
- D500 proximal-tail detectability 大幅下降；
- D500 中下部 source detectability 特别低；
- D500 slow-axis background drift 更明显；
- D500 X4 lateral jitter 更高；
- D500 detrended local mean–mean coupling 更弱；
- 但 D500 background RCV 并未明显更差；
- valid-frame fraction 反而较高；
- pseudo-control coupling 没有明显升高；
- real-minus-control coupling 仍为正。

因此当前最严谨的结论是：

> **较低 measurement detectability 很可能是 D500 较弱 observed coupling 的贡献因素之一，但它不能单独、充分解释全部 diameter-dependent coupling difference。**

---

## 16. 对 D128–D285 前一轮 coupling 结论的影响

本轮 QC 反而增强了 D128–D285 局部 source–tail covariance 的可信度，因为：

1. real coupling 明显高于 matched pseudo coupling；
2. BG–real-tail detrended correlation 近 0；
3. coupling hotspot 不等同于 detectability hotspot；
4. 前一轮 lag=0 spatial localization 稳定；
5. proximal-tail 区域仍具有相对更好的 local-background detectability。

因此：

> **D128–D285 的 local source–tail spatial covariance 很难仅用系统背景、ROI geometry 或简单 detectability bias 解释。**

仍然不能将其升级为物理因果结论。

---

## 17. 这一轮最重要的三个结论

### 结论 1

\[
\boxed{\text{D500 的 source 与 proximal-tail SV detectability 明显降低，尤其是 source 中下部。}}
\]

### 结论 2

\[
\boxed{\text{D500 不是整卷数据全面失效：background stability、valid-frame rate 与 matched negative controls 总体正常。}}
\]

### 结论 3

\[
\boxed{\text{Real source–tail coupling 明显强于 matched pseudo controls，因此已有 local spatial covariance 不是简单背景伪相关。}}
\]

---

## 18. 当前推荐的论文级安全表述

可以写：

> Across D128–D285, the observed local source–tail covariance remained spatially specific relative to matched lateral pseudo-source/pseudo-tail controls and was not reproduced by local background fluctuations. D500 exhibited substantially lower source and proximal-tail SV detectability, particularly in deeper source strata, which may contribute to the weaker observed local coupling. However, background variability, valid-frame fraction, and matched pseudo-control coupling did not show a corresponding global deterioration, indicating that reduced detectability alone does not fully account for the diameter-dependent coupling differences.

不能写：

- “D500 低 SNR 导致 coupling 消失”；
- “D500 发生了确定的物理机制转换”；
- “Z<3 表示不可检测”；
- “pseudo controls 证明了 causality”；
- “structural SBR 就是真正 OCT SNR”。

---

## 19. 当前阶段科学图景

综合 source-depth × tail-depth coupling 与本轮 QC，目前最合理的整体图景是：

> **D128–D285 中，血管中部至中下部的局部 SV 波动与 proximal tail 存在跨 flow、近零 lag 的空间共变；这一关系明显强于 matched lateral pseudo controls。随着直径增大到 D500，局部 mean–mean coupling 明显减弱，同时 source 中下部和 proximal tail 的 SV detectability 大幅下降、slow-axis drift 与 X4 jitter 增加。现有数据支持 measurement quality 参与 D500 弱 coupling，但不能证明其为唯一或充分解释。**

这应作为当前阶段的解释边界，后续如继续研究机制，应保持“测量可检测性变化”与“真实 source–tail spatial topology 变化”两个假设并行，而不是预先认定其中一个。
