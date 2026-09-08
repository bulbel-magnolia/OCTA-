# D128 / D235 frozen source ellipse 轴向三等分数据总结

本次新增 D128 的 2,422 帧和 D235 的 2,225 帧，共 10 个 volume、4,647 帧及 50 个 volume × slow-axis segment 汇总。D285/D500 未重新计算；原有 946 个已跟踪文件的 SHA 保持不变。

## 固定定义与汇总口径

直接导入并调用 `../d500_source_axial_local_array_audit.py` 中的 `band_weights`、`weighted_stats`、`require_band_reconstruction`。每帧使用已有 frozen localization geometry，同一个 source ellipse 内按相对物理轴向坐标 `u=(z-z_top)*dz/diameter` 分为 upper `[0,1/3)`、middle `[1/3,2/3)`、lower `[2/3,1]`。使用原有 16×16 subpixel samples，三部分权重逐帧精确重构原 ellipse。

读取冻结导出的 `sv_raw = var(abs(E),1,3)`（方差分母 N），不重新导出或变换信号。`dx=12.7 μm`，`dz=6.7 μm`，`supersample=16`；无 log、normalization、background subtraction 或新 ROI。均值为 band 内 fractional-area-weighted raw SV，Q 为 raw SV 按实际面积积分，Q fraction 为该 band Q / 完整 source Q。

只分析冻结主表中已有有效帧，不追加筛帧。帧编号从 0 开始；slow-axis segment=`frame_index_0based // 100`，五段依次覆盖 000–099、100–199、200–299、300–399、400–499。D128 的 legacy `flowXX` 标识同时保留，并映射为 `D128_FXX_V01`。

Volume 和 segment 的各数值列分别取帧间中位数；ratio median 是逐帧比值的中位数，不是两个 band 中位数的比值。`strict_upper_middle_lower` 在逐帧表中判断该帧均值，在汇总表中判断三个 band 中位数；另列逐帧排序计数、比例及是否全部帧满足。分别取中位数后的 Q fraction 列不要求加和恰为 1，逐帧 Q fractions 的和已校验。所有结果均为描述性汇总，不计算推断统计量。

## 1. D128 五个 flow 的排序

五个 flow 中 **0/5** 的 band 中位数满足 `upper > middle > lower`。Flow 1、3、5、10 为 `middle > upper > lower`；flow 7 为 `lower > middle > upper`。逐帧共有 567/2,422 帧（23.41%）满足 upper > middle > lower。

## 2. D235 五个 flow 的排序

五个 flow 中 **5/5** 的 band 中位数满足 `upper > middle > lower`。逐帧为 2,213/2,225 帧（99.46%）；flow 1、3、5、7、10 分别有 8、2、0、1、1 帧不满足。因此 volume 排序一致不代表每一帧都满足。

下表 upper/middle/lower 为 band raw SV 帧间中位数，数值以 `×10^8` 展示；CSV 保留未缩放数值。M/U、L/U、L/M 为逐帧比值的中位数。

| Diameter (μm) | Flow (mm/s) | 帧数 | upper ×10^8 | middle ×10^8 | lower ×10^8 | M/U | L/U | L/M | band 中位数 U>M>L |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---:|
| 128 | 1 | 486 | 1.479069 | 1.542477 | 1.395322 | 1.082228 | 0.966705 | 0.911398 | 否 |
| 128 | 3 | 468 | 1.512742 | 1.630230 | 1.338595 | 1.067881 | 0.880654 | 0.855568 | 否 |
| 128 | 5 | 491 | 1.549211 | 1.577008 | 1.348991 | 1.038320 | 0.861737 | 0.852965 | 否 |
| 128 | 7 | 493 | 1.033971 | 1.269548 | 1.334704 | 1.229880 | 1.280572 | 1.040641 | 否 |
| 128 | 10 | 484 | 1.508602 | 1.623678 | 1.408570 | 1.082803 | 0.946494 | 0.863740 | 否 |
| 235 | 1 | 448 | 2.873780 | 1.624536 | 0.877416 | 0.571444 | 0.305890 | 0.532832 | 是 |
| 235 | 3 | 448 | 3.022723 | 1.628667 | 0.901628 | 0.542127 | 0.296040 | 0.554458 | 是 |
| 235 | 5 | 424 | 3.123428 | 1.602537 | 0.874019 | 0.520029 | 0.281605 | 0.550446 | 是 |
| 235 | 7 | 444 | 3.030010 | 1.614107 | 0.890144 | 0.531319 | 0.295884 | 0.551253 | 是 |
| 235 | 10 | 461 | 3.155958 | 1.657916 | 0.898539 | 0.522685 | 0.286586 | 0.545140 | 是 |

## 3. D128 的比值范围

五个 volume 的 middle/upper ratio median 为 **1.038320–1.229880**，lower/upper 为 **0.861737–1.280572**。这两个范围取自五个 volume 的比值中位数，不是所有帧的极值。

## 4. D235 的比值范围

五个 volume 的 middle/upper ratio median 为 **0.520029–0.571444**，lower/upper 为 **0.281605–0.305890**。

## 5. D128 → D235 的轴向结构变化

在相同的五个 flow 下，D128 的 middle 接近或高于 upper，flow 7 的 lower 也高于 upper；D235 的 upper 均高于 middle 和 lower。逐一匹配 flow 后，D235 的 upper band 中位数为 D128 的 1.943–2.930 倍，middle 为 0.999–1.271 倍，lower 为 0.629–0.674 倍。对应的 middle/upper 和 lower/upper 比值中位数在五个 flow 下均降低。这描述了当前两组冻结数据中更强的上部信号占优，不作直径变化的因果推断。

| Diameter (μm) | Flow (mm/s) | upper Q fraction median | middle Q fraction median | lower Q fraction median | 逐帧 U>M>L |
|---:|---:|---:|---:|---:|---:|
| 128 | 1 | 0.283352 | 0.428457 | 0.272419 | 120/486 |
| 128 | 3 | 0.289908 | 0.439321 | 0.259888 | 127/468 |
| 128 | 5 | 0.294359 | 0.436071 | 0.254661 | 163/491 |
| 128 | 7 | 0.249060 | 0.426481 | 0.312803 | 31/493 |
| 128 | 10 | 0.286338 | 0.437453 | 0.262303 | 126/484 |
| 235 | 1 | 0.471187 | 0.382256 | 0.142839 | 440/448 |
| 235 | 3 | 0.483234 | 0.370806 | 0.143806 | 446/448 |
| 235 | 5 | 0.493492 | 0.364736 | 0.138600 | 424/424 |
| 235 | 7 | 0.485151 | 0.368486 | 0.142248 | 443/444 |
| 235 | 10 | 0.491829 | 0.365767 | 0.140948 | 460/461 |

## 6. Slow-axis 五段的稳定性

D235 的 **25/25** 个 segment 都保持 band 中位数 `upper > middle > lower`，每个 flow 均为 5/5。排序稳定，但比值随段变化：全部段的 M/U 中位数为 0.459928–0.651765，L/U 为 0.247469–0.357005；五个 flow 的末段 M/U 和 L/U 均高于首段。

D128 仅 **8/25** 个 segment 满足该排序，各 flow 依次为 2/5、1/5、4/5、0/5、1/5，五段未保持统一排序。全部段的 M/U 中位数为 0.905246–1.328673，L/U 为 0.651638–1.473224。

下表按 segment 0 → 4 列出 band 中位数的实际排序；U=upper，M=middle，L=lower。

| Diameter (μm) | Flow (mm/s) | seg 0 | seg 1 | seg 2 | seg 3 | seg 4 |
|---:|---:|:---:|:---:|:---:|:---:|:---:|
| 128 | 1 | M>U>L | U>M>L | U>M>L | M>U>L | L>M>U |
| 128 | 3 | M>U>L | U>M>L | M>L>U | M>U>L | M>U>L |
| 128 | 5 | U>M>L | U>M>L | U>M>L | U>M>L | M>L>U |
| 128 | 7 | L>M>U | L>M>U | L>M>U | M>L>U | L>M>U |
| 128 | 10 | M>U>L | U>M>L | M>L>U | M>U>L | L>M>U |
| 235 | 1 | U>M>L | U>M>L | U>M>L | U>M>L | U>M>L |
| 235 | 3 | U>M>L | U>M>L | U>M>L | U>M>L | U>M>L |
| 235 | 5 | U>M>L | U>M>L | U>M>L | U>M>L | U>M>L |
| 235 | 7 | U>M>L | U>M>L | U>M>L | U>M>L | U>M>L |
| 235 | 10 | U>M>L | U>M>L | U>M>L | U>M>L | U>M>L |

## 输入身份与验证

D128 来自冻结 Release `formal-sv-d128-v21-run001`：25/25 ZIP 的字节数和 SHA-256 均与冻结清单一致，所分析 2,422/2,422 NPZ 的字节数和 SHA-256 均匹配。使用冻结 `localization.csv`，完整 source mean 与 D128 no-background 主表比较，最大相对误差 **3.239476120439153e-16**；同时对照 observed raw source area 和 Q，两者最大相对误差均为 0。

D235 来自五个 volume 的 retained MAT；2,225/2,225 MAT 的 SHA-256 与 `analysis/formal_sv_diameter_v1/framewise_primary.csv::input_mat_sha256` 完全一致。完整 source area、Q、mean 的最大回放相对误差均为 **0**。

全部 4,647 个 `sv_raw` 数组均为 351×500 且 finite；band 权重重构最大绝对误差为 0，band Q 重构最大相对误差为 **4.71412883780283e-16**。所有 replay 和 Q 重构误差均严格小于 `1e-10`。CSV 回读后重新核对冻结帧集合、input SHA、比值、Q fraction、volume/segment 汇总及排序。

代码和输入固定于提交 `eb92d7a5f8ba93d39e5cb68a1cd66ee65cf80588`，具体 SHA 见 `provenance.json`。该 SHA 表示计算基线；结果提交是包含本目录的 Git commit。

## 输出文件

- `d128_d235_source_axial_framewise.csv`：4,647 帧，含 band area/Q/mean、比值、Q fractions、冻结 source mean、几何、输入身份和逐帧误差。
- `d128_d235_source_axial_volume_summary.csv`：10 个 volume。
- `d128_d235_source_axial_segment_summary.csv`：50 个 volume × slow-axis segment。
- `validation.json`：SHA、shape/finite、数值回放、权重/Q 重构和派生表验证记录。
- `provenance.json`：冻结来源、函数、参数、软件版本、基线提交及输入文件 SHA。
- `output_sha256.csv`：上述输出和本 README 的字节数及 SHA-256；不包含自身。

本目录仅发布派生 CSV / JSON / README，供后续四直径统一比较使用。原始 MAT、NPZ、ZIP 保留本地。
