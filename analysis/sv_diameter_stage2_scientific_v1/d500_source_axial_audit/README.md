# D500 冻结 source ellipse 轴向局部数组审计

本轮结果显示，D500 的低 source denominator 分布在上、中、下三个轴向区段，中部的相对下降最突出。它并非只出现在某一小段 slow axis。与 z_top 的体积内空间共变则以上部最强；这两个观察分别描述跨直径差异和体积内共变，不应混作同一个效应。

## D500 与 D285：共同流速比较

下表为 all-valid 子集各 band 原始 SV 均值的体积中位数差异，百分比为 100×(D500−D285)/D285。每个条件只有一个 scan volume，以下为描述性对比。

| 流速 (mm/s) | 上部 | 中部 | 下部 |
|---:|---:|---:|---:|
| 1 | -33.21% | -59.44% | -45.46% |
| 3 | -35.33% | -59.91% | -46.18% |
| 5 | -39.39% | -57.50% | -43.39% |
| 7 | -41.96% | -53.91% | -41.36% |
| 10 | -36.68% | -53.55% | -39.80% |

五个共同流速中，三段均低于 D285；中部相对下降 53.55%–59.91%，上部下降 33.21%–41.96%，下部下降 39.80%–46.18%。相对下降以中部最大，但 band 均值的绝对差以上部最大（约 1.08×10^8–1.52×10^8 raw-SV units）。不能把“相对降幅最大”直接等同于对完整分母下降的最大绝对贡献。

## D500 内部轴向结构与 spatial robustness

全部 8 个 D500 volume 均呈上部 > 中部 > 下部。逐帧 lower/upper 比值的体积中位数为 0.197–0.225。上部 Q 占比的逐帧中位数约为 0.583–0.621，中部为 0.254–0.284，下部为 0.123–0.131；各占比中位数不要求相加恰为 1，逐帧 Q 重组由数值门禁验证。

8 卷 × 5 个 slow-axis segments 的全部 40 个分段均保留该轴向排序。在五个共同流速的 25 个同位置分段对比中，D500 三段均低于 D285，中部的相对下降均最大。因此，这一描述模式贯穿所分析 volume 的 slow axis，幅度仍随空间位置变化。

## Restricted subsets

all valid、direct only、central geometry、direct + central geometry 四种预定义子集均保留三段低于 D285且中部相对下降最大的结果；D500 的上 > 中 > 下排序也在每个子集的全部 8 卷中保留。central geometry 按每卷 X1 apparent width 与 z_top 的第 10–90 百分位交集定义，包含边界。direct only 沿用 frozen direct_candidate 标签。这使几何极端值或 short-gap composition 难以单独解释所见模式。

## z_top 的体积内空间共变

Spearman 由平均秩的 Pearson 相关直接计算，未计算 p-value。all-valid D500 的上部 ρ=0.888–0.915，中部 ρ=0.525–0.659，下部 ρ=0.283–0.495；每卷均以上部最强。direct + central 子集内，上部仍为 0.853–0.913，并在全部 8 卷中强于其余两段。

这说明冻结 source 与 z_top 的空间关联主要体现在上部信号；它是共变强度，不能解释为 z_top 对信号的因果作用或物理灵敏度。跨直径相对降低以中部最明显，体积内 z_top 关联以上部最明显，两者可以同时成立。

## 输入、定义与验证

使用原 retained MAT：5 个 D285 volume 与 8 个 D500 volume，正式有效帧共 6284（D285 2371、D500 3913）。D500 的 2/9/12 mm/s 仅作 within-D500 secondary extension。无需重新 export；每个输入 MAT 的 SHA 均与 frozen framewise_primary.csv 完全一致，shape=351×500 且全部 finite。
信号为 sv_raw = var(abs(E),1,3)，分母 N；没有 log、normalization、gain、clipping、positive truncation 或 background subtraction。
冻结几何保持 X4 centre、X1 apparent width、v2.1 z_top 和真实物理直径，dx=12.7 μm、dz=6.7 μm、ellipse supersample=16。仅在同一个 ellipse 的 subpixel samples 内按相对物理轴向位置划分三等分；没有新矩形或新 source ROI。
完整 source area/Q/mean 回放最大相对误差 3.307e-16，三段 Q 重组最大相对误差 5.279e-16，均严格小于 1e-10。7 项合成分区/门禁/相关诊断测试通过。525 个既有结果及配置文件保持哈希不变，唯一允许修改的既有文件为本次审计脚本。

原脚本的 geometry_from_row、band_weights、weighted_stats、qtls 经 AST 比较完全一致。补充了原先缺失的 band-Q 失败停止、所要求的 z_top 诊断和 provenance；详见 implementation_check.json。网页端状态文件中的数组缺失描述对应当时环境，本记录确认本地真实数组已恢复并完成执行。

## 解释边界

本轮支持的数据描述是：D500 的 source reduction 涉及全部轴向三段，并伴随中部相对不足和稳定的上高下低分布。没有据此识别光学衰减、多重散射、血流物理或系统灵敏度机制。scan volume 是实验单位；B-scan 与 slow-axis segment 均非独立重复。未重新计算 Stage 2 主 RI_tail、depth profile 或主对比，未进行回归、显著性检验或机制拟合。

## 文件

- source_axial_band_framewise.csv：6284 帧的三段面积、Q、raw mean、Q fraction、比值、几何标记、z_top 和 MAT SHA。
- source_axial_band_volume_summary.csv：13 卷 × 4 子集，52 行。
- d500_vs_d285_axial_band_contrasts.csv：5 flows × 4 子集 × 3 bands，60 行。
- d500_axial_band_slow_axis_summary.csv：D500 8 卷 × 5 segments，40 行。
- d500_axial_band_z_top_spearman.csv：D500 8 卷 × 4 子集 × 3 bands，96 行。
- d500_vs_d285_axial_band_segment_contrasts.csv：5 flows × 5 segments × 3 bands，75 行；由同一 framewise band 表机械汇总。
- validation.json、provenance.json、implementation_check.json、protected_sha256.json、input_sha256.csv、output_sha256.csv：门禁与身份记录。

执行入口：`python analysis/sv_diameter_stage2_scientific_v1/d500_source_axial_local_array_audit.py --mat-root <retained-MAT-root>`。主机绝对定位信息仅保存在忽略的本地执行回执。
