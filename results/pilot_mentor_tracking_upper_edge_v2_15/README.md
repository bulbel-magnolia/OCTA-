# 15 帧上边界 v2 候选试运行

日期：2026-09-06

本目录用于比较 `legacy_connected_component_v1` 与候选
`persistent_core_paired_edge_v2`。它不覆盖
`results/pilot_mentor_tracking_15` 中已经冻结的 v1 结果。

## 本次实际改变

- 横向 X4/X1 仍沿用师兄的完整 500 帧自动定位链；
- SV 原始数据、128 μm 椭圆、双侧背景、500 μm 矩形拖尾和积分方法不变；
- 只升级 z_upper：加入横向覆盖、连续多行信号、顶缘内外对比和物理直径
  平衡；
- 峰值阈值由 15% 提高到 20%，背景噪声阈值由 3 倍提高到 4 倍；
- 原有 QC 代码和阈值未修改，没有人工点击或逐帧手工调整。

## 15 帧结果

15/15 帧通过原有定位、几何、背景和窗口 QC，并得到有效 Q_vessel、Q_tail
与 R。此前指出的五帧变化如下，正值表示顶缘向图像下方/深部移动：

| 面板 | 扫描与帧 | v1 顶缘 | v2 顶缘 | 变化 |
|---|---|---:|---:|---:|
| 1行3列 | flow01, B-scan 499 | 188 | 192 | +4 px / +26.8 μm |
| 2行2列 | flow03, B-scan 249 | 184 | 189 | +5 px / +33.5 μm |
| 3行1列 | flow05, B-scan 0 | 199 | 211 | +12 px / +80.4 μm |
| 4行1列 | flow07, B-scan 0 | 202 | 214 | +12 px / +80.4 μm |
| 4行3列 | flow07, B-scan 499 | 191 | 192 | +1 px / +6.7 μm |

[同图旧/新边界对照](comparison_vs_legacy_upper_edge_15.png)中，黄色虚线为
v1 顶缘，红线为 v2 顶缘，橙线为 v2 的 128 μm 下边界，青色为 v2 横向
范围。可用
[v2 单独总览](mentor_tracking_upper_edge_v2_overview_15.png)查看不叠加旧线
的结果。

## 不能忽略的全卷异常

[5 卷轨迹图](mentor_tracking_upper_edge_v2_trajectories_5x500.png)显示 flow03
约第 425–465 帧的 v2 上边界错误向深部偏离，最大相对 v1 为 +46 pixel。
该段合格候选稀疏，并混入少量更深的假候选，原轨迹补全器把它们连接起来。
这段不包含本次抽取的 0、249、499 帧，所以不影响上表 15 帧是否成功运行；
它说明 v2 目前只能作为候选供目测复核，不能宣称 2500 帧均已精确定位。

## 参数选择与复现

`legacy_threshold_matrix_15.csv` 与 `persistent_threshold_matrix_15.csv` 保存
峰值比例 0.15/0.20/0.25 × 噪声倍数 3/4/5 的完整 15 帧矩阵。噪声 5 倍
使 2/15 帧缺少局部种子，因此候选选择 0.20/4。配置分别见：

- `../../config/tracking_config.upper_edge_v2.a020_n4.json`
- `../../config/run_config.mentor_tracking.upper_edge_v2.pilot.json`

矩阵可由 `../../scripts/compare_upper_edge_thresholds.py` 重跑。原始 OCT、Flow
DICOM 与大体积数组仍保存在仓库外；本目录只保存可审计的派生表格和图片。
