# SV 四直径 source 信号与轴向结构本轮分析总结 v1

日期：2026-09-08

## 1. 本轮范围

本轮只分析已有 SV 数据，不讨论新实验，也不扩展到新的物理机制模型。

比较对象为共同流速 1 / 3 / 5 / 7 / 10 mm/s 下的四个物理直径：

- D128
- D235
- D285
- D500

冻结定义保持不变：

- `sv_raw = var(abs(E),1,3)`；
- source ROI 为冻结 X1 × 物理直径的 fractional ellipse；
- X4、X1、continuity-first v2.1 z_top 均不改变；
- dx = 12.7 μm，dz = 6.7 μm；
- ellipse supersample = 16；
- 无 log、normalization、background subtraction、gain、clipping 或 positive truncation。

source ellipse 内部继续使用完全相同的三等分：upper / middle / lower，各占真实物理直径的 1/3。

## 2. 数据身份与验证

四直径三段分析全部来自真实二维 `sv_raw` 数组：

- D128：冻结 release 回放，2422 个有效 B-scan；
- D235：retained MAT，2225 个有效 B-scan；
- D285：retained MAT，2371 个有效 B-scan；
- D500：retained MAT，3913 个有效 B-scan（含 2/9/12 mm/s secondary extension）。

合计 10931 个有效空间 B-scan。

D128 release ZIP/NPZ、D235/D285/D500 retained MAT 均通过已有 SHA gate。完整 frozen source 回放和 upper+middle+lower Q 重组误差均在 `1e-10` 门限以内，实际最大误差为机器精度量级（约 1e-16）。

B-scan 和 slow-axis segment 仍是 volume 内空间位置，不作为独立实验重复。

## 3. 先明确“血管总信号”的两个不同含义

### 3.1 `source_mean_raw`

这是冻结 source ellipse 内的加权平均 raw SV，也是主指标 `RI_tail = tail_mean_raw / source_mean_raw` 的分母。

它回答的是：

> 血管 ROI 内平均每单位面积有多强的 SV 信号？

### 3.2 `source_q_raw`

这是整个 source ellipse 的积分 raw-SV 总量。

它回答的是：

> 把整个血管 ROI 内的 raw SV 全部积分起来，一共是多少？

由于直径增大时 source area 同时大幅增加，因此 `source_mean_raw` 和 `source_q_raw` 不能混用。

本轮下面的正式四直径趋势以 `source_mean_raw` 和三段内部结构为主。此前用“median area × median source mean”得到的积分量只是帮助理解面积效应的粗略示意，不作为正式 `source_q_raw` 汇总结果录入本轮结论。若后续需要严格比较四直径 `source_q_raw`，应直接从已有 framewise `source_q_raw` 做 volume-level 汇总。

## 4. 四直径 source 平均强度关系

五个共同 flow 的 volume-level `source_mean_raw` 中位数，再取跨 flow 中位数：

| Diameter | source_mean_raw | 相邻直径变化 |
|---:|---:|---:|
| 128 μm | 1.54205×10^8 | — |
| 235 μm | 1.83653×10^8 | +19.10% |
| 285 μm | 1.87972×10^8 | +2.35% |
| 500 μm | 1.02259×10^8 | −45.60% |

因此 source 平均强度与直径不是单调关系，而是：

`D128 ↑ D235 ≈ D285 ↓↓ D500`

也就是说：

1. 128→235：平均 vessel/source signal 明显增强；
2. 235→285：基本进入平台；
3. 285→500：平均 vessel/source signal 大幅下降，D500 成为四者中最低。

## 5. 四直径 source 内部轴向结构

为保证四直径完全同口径，下面的 `middle/upper` 与 `lower/upper` 使用每个 volume 的 band raw-SV 中位数之比，然后再取五个共同 flow 的中位数。

| Diameter | Upper | Middle | Lower | Middle/Upper | Lower/Upper | Upper Q frac | Middle Q frac | Lower Q frac |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 128 | 1.50860×10^8 | 1.57701×10^8 | 1.34899×10^8 | 1.0763 | 0.9337 | 0.2863 | 0.4361 | 0.2623 |
| 235 | 3.03001×10^8 | 1.62454×10^8 | 0.89014×10^8 | 0.5327 | 0.2938 | 0.4852 | 0.3685 | 0.1422 |
| 285 | 3.37897×10^8 | 1.56517×10^8 | 0.79349×10^8 | 0.4655 | 0.2385 | 0.5232 | 0.3474 | 0.1239 |
| 500 | 2.10310×10^8 | 0.68352×10^8 | 0.46149×10^8 | 0.3311 | 0.2212 | 0.5899 | 0.2782 | 0.1301 |

Q fraction 为逐帧 Q fraction 的 volume 中位数再做跨-flow中位数，因此三项中位数不要求严格相加为 1；逐帧 upper+middle+lower Q 已通过精确重组门禁。

## 6. `upper > middle > lower` 从什么时候开始稳定出现？

共同 5 个 flow：

| Diameter | volume-level strict order | slow-axis segments strict order |
|---:|---:|---:|
| 128 | 0/5 | 8/25 |
| 235 | 5/5 | 25/25 |
| 285 | 5/5 | 25/25 |
| 500 | 5/5 | 25/25 |

D500 若计入 2/9/12 mm/s extension，则为 8/8 volumes、40/40 segments。

因此：

- D128 没有稳定的 `upper > middle > lower` 轴向梯度；多数 flow 中 middle 略高于 upper，flow 7 甚至出现 lower > middle > upper。
- D235 开始突然建立非常稳定的 `upper > middle > lower`。
- D285 延续并进一步强化该结构。
- D500 继续保持 upper-dominant，但 middle 相对不足更加明显。

## 7. 三个直径跃迁分别发生了什么？

以下百分比使用四直径跨-flow band 中位数做描述性比较。

### 7.1 128→235 μm：第一次“大重排”

- Upper：+100.85%
- Middle：+3.01%
- Lower：−34.01%
- 整体 `source_mean_raw`：+19.10%

最突出的不是整个 vessel 同比例变亮，而是：

`Upper 大幅增强 + Middle 基本不变 + Lower 明显下降`

因此 source 从 D128 的近似均匀 / middle-dominant 状态，转为 D235 的强 upper-dominant 状态。

对应形状比值：

- Middle/Upper：1.076 → 0.533
- Lower/Upper：0.934 → 0.294

这是四直径中最大的内部结构转折。

### 7.2 235→285 μm：已有梯度进一步增强

- Upper：+11.52%
- Middle：−3.65%
- Lower：−10.86%
- 整体 `source_mean_raw`：+2.35%

此阶段 source 平均强度已经接近平台，但 upper-dominant 结构继续小幅强化：

- Middle/Upper：0.533 → 0.466
- Lower/Upper：0.294 → 0.239

### 7.3 285→500 μm：整体 source 下降，并出现更突出的 middle deficit

- Upper：−37.76%
- Middle：−56.33%
- Lower：−41.84%
- 整体 `source_mean_raw`：−45.60%

D500 三段都降低，但 middle 的相对下降最大。

形状变化：

- Middle/Upper：0.466 → 0.331，继续明显下降；
- Lower/Upper：0.239 → 0.221，只小幅下降；
- 因为 middle 降得比 lower 更明显，Lower/Middle 在 D500 反而回升。

因此 D500 的特殊之处不能简单描述为“越深越弱得更多”，更准确的是：

> upper-dominant 结构仍然存在，但 middle 相对于 upper 出现额外不足，而 lower/upper 已接近平台。

## 8. 和主 RI_tail 变化的对应关系

此前主结果的跨-flow RI_tail 中位数约为：

- D128：0.25947
- D235：0.17529
- D285：0.16095
- D500：0.25764

即：

- 128→235：RI −32.44%
- 235→285：RI −8.18%
- 285→500：RI +60.07%

与本轮 source 结果对应：

### 128→235

RI 大幅下降时，恰好也是 source 内部第一次发生巨大重排、并且 `source_mean_raw` 明显升高的阶段。数学上，分母增加会推动 RI 下降；与此同时 absolute tail 本身也在下降，因此两者共同作用。

### 235→285

source 平均强度只增加约 2.35%，内部结构也只是继续小幅强化，因此 RI 只进一步小幅下降。

### 285→500

D500 的 absolute 500 μm tail 并没有增强，反而较 D285 更低；但 `source_mean_raw` 大幅下降约 45.6%，所以低 denominator 是 RI rebound 的主要来源。

这些是比值分解和同步变化关系，不等同于物理因果机制。

## 9. 本轮最简洁的数据图景

四直径 source 可以概括成两个阶段：

### 阶段 A：D128→D235

从“内部近似均匀 / middle 较强”突然转成“upper 明显主导”。这是第一次、也是最大的轴向结构重排。

### 阶段 B：D235→D285→D500

upper-dominant 结构持续存在并强化；到 D500 时，整个 source 平均强度明显下降，同时 middle 相对 upper 的不足进一步突出，而 lower/upper 变化已经较小。

因此现有数据不支持简单的“直径越大，整个血管信号越强”或“直径越大，三层同比例衰减”。更符合数据的描述是：

> **直径变化同时改变了 source 的平均强度和 source 内部轴向信号分布，而且最明显的结构转折发生在 128→235 μm。**

## 10. 解释边界

本轮只登记数据事实和数学关系：

- 不把 upper / middle / lower 的变化归因于具体光学机制；
- 不把 z_top 共变解释成因果；
- 不做新回归、p-value 或机制拟合；
- 不把 spatial segment 当成独立实验重复；
- 不用粗略的 `median area × median mean` 替代正式 `source_q_raw` 汇总。

后续如继续分析已有数据，最自然的两个方向是：

1. 直接从已有 framewise `source_q_raw` 严格汇总四直径积分总量；
2. 把 source 的三段重排与对应 tail absolute signal 的变化逐 flow 对齐，进一步拆解 RI 变化来自 source 与 tail 的各自贡献。
