# 数据字段与分析入口

分析前按 `scan_id + frame_index_0based` 一对一连接表格。`frame_index_0based` 为 0–499；MATLAB 帧号为 1–500。每个 scan 是一个完整扫描卷，卷内帧为相邻空间位置。实验独立单位字段未知时为 NA。

| 字段 | 含义与单位 |
|---|---|
| valid | 正式指标纳入开关；仅 True 纳入 |
| raw_input_available / input_mapping_valid | 原始重建图存在 / 文件、帧号、尺寸、重复数与裁剪元数据核对通过 |
| z_upper_px / seed_z_upper_px | 冻结最终上边界 / 原候选位置，0 基像素中心坐标 |
| geometry_qc_valid | 冻结定位几何 QC 是否通过 |
| localization_source | direct_candidate_supported、short_gap_assisted 或 missing |
| new_tracking_class | 原追踪器分类，不等价于最终 valid |
| assessability | 原 vessel_presence_prediction 的明确别名 |
| q_vessel | 扣背景后在完整源椭圆内积分；原生幅值平方 × μm² |
| source_area_um2 | 分数权重计算的椭圆面积；μm² |
| source_mean | q_vessel / source_area；原生幅值平方 |
| q_tail | 扣背景后在 0/500 μm 矩形窗口内积分；原生幅值平方 × μm² |
| ratio_tail_to_vessel | q_tail / q_vessel，无量纲；不截断负值 |
| diagnostic_* | 被排除帧仍可计算时的诊断积分；不能代替正式纳入开关 |
| invalid_reason / localization_invalid_reason | 定量失效原因 / 更具体的原定位失效原因 |
| background_qc_valid / window_qc_valid | 双侧背景完整 / source 与 500 μm tail 窗口完整 |
| detection_status | 本批均 not_evaluated_no_matched_blank |
| detected / detectable_length_um | 本批 NA，不能用 0 或 false 替换 |

完整 profile 的 z_index_0based 是裁剪 B-scan 内原始深度；z_um=z×6.7 μm，r_um 相对于血管最深下缘。V 是固定血管宽度的加权行均值；B_left/B_right 是跳 3 取 5 的两侧行均值，B 为双侧合并均值；T=V−B；P=有效横向宽度×T，单位为原生幅值平方×μm。tail_z_fraction 记录轴向尾窗覆盖权重。`validity` 是该深度行有限性；正式纳入仍须 `frame_valid=True`。

二维数组保存源区覆盖权重 source_weights；尾区二维权重等于 tail_z_weights[:,None]×tail_x_weights[None,:]。原始 sv_raw 与背景列索引足以重建原深度背景和所有积分。负残差必须保留；PNG 中的亮度变换不用于计算。信号单位未标定到光学功率或物理血流量，不应自行解释为这些单位。

敏感性表的 variant 分别为 x_minus_1px、x_plus_1px、z_minus_1px、z_plus_1px、background_skip_plus_2px。先保留 primary_valid 与每个扰动自身 valid，避免在失效扰动上直接计算百分比。既有规则未按流速趋势或显著性挑选参数。
