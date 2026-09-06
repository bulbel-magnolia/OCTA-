# 师兄自动定位与 SV 矩形拖尾指标接入说明

日期：2026-09-06

## 1. 方案边界

当前方案把“寻找血管”和“计算拖尾”分成两条清楚的链。

定位链读取每个扫描对应的完整 500 帧 Flow DICOM，执行师兄的慢轴连续
跟踪。定量链读取我们从原始 OCT 重建的线性 sv_raw 单帧，执行物理椭圆
源区、逐行双侧背景、500 μm 矩形拖尾积分和敏感性分析。Flow DICOM 不参与
Q_vessel、Q_tail 或 R 的数值积分，师兄原有的 AUC 指标也不参与本项目结果。

这一路径没有表面检测、固定表面 z、折射率或表面到血管距离参数。

## 2. 自动定位流程

每卷输入数组顺序为 frame × z × x，尺寸为 500 × 351 × 500。

1. 师兄 tracking_core 在 500 帧上计算横向与轴向候选，并用完整慢轴序列
   做 Viterbi 连续轨迹约束。
2. 顶缘同时计算 alpha 0.10、0.15、0.20，正式定位使用主 alpha 0.15 的
   z_upper。
3. 每帧在 Viterbi 局部窗内运行 X1。X1 的 75% 分位数连续超阈值段提供
   血管横向宽度。
4. 同一局部窗计算 X2 稳健质心；权重经过背景扣除和 90% 分位截断，避免
   单个极亮像素控制中心。
5. X4 仅修正三帧序列中的孤立跳点，保留连续的真实空间变化。X4 是最终
   中央 A-line。
6. CNR、连续性、宽度一致性、轴向完整性和相邻帧支持共同生成可评估性。
   分类为 assessable、uncertain 或 not_assessable。

以上是已冻结的 `legacy_connected_component_v1`。2026-09-06 增加的候选
`persistent_core_paired_edge_v2` 保留第 1、3–6 步，只升级轴向顶缘：

1. 用 20% 峰值且至少 4 倍背景噪声的强阈值确认血管核心；
2. 用较低滞回阈值从核心向上寻找可能顶缘，但要求信号横向覆盖一个按
   `dx` 换算的血管宽度，且沿 z 连续存在；
3. 候选顶缘按上方到内部的亮度跃迁、内部相对外部的横向支持、固定
   128 μm 窗内的持续比例，以及上下半径平衡共同排序；
4. 不使用表面 z 或表面到血管距离，不改变原有 assessable/uncertain/
   not_assessable QC 规则。

v2 的完整参数冻结在
`config/tracking_config.upper_edge_v2.a020_n4.json`。默认配置仍指向 v1，
避免旧命令在未声明版本时产生不同结果。

## 3. 转换为我们的几何

正式几何使用以下固定映射：

| 几何量 | 来源 |
|---|---|
| 中央 A-line | 师兄 X4 |
| 左右边界间宽度 | 师兄 X1 连续段像素数 |
| 左右边界位置 | 保持 X1 宽度，以 X4 为中心重新放置 |
| 上边界 | 主 alpha=0.15 的 z_upper 像素中心上方 0.5 pixel |
| 下边界 | 上边界加 128 μm / 6.7 μm·pixel⁻¹，不取整 |
| 椭圆纵轴 | 128 μm |
| 拖尾矩形长度 | 500 μm |

横向和纵向始终分别使用 12.7 μm/A-line 与 6.7 μm/pixel。椭圆、矩形和
边界像素采用分数面积权重。

正式 source QC 要求跟踪未失败、X1 局部体有效、X1 未使用备用段，且帧分类
为 assessable。uncertain 与 not_assessable 仍保存坐标和诊断，R 不报告。

## 4. 真实数据接口验证

师兄原始 tracking_core.py 与仓库副本的 SHA-256 均为：

    214CD12F5AF958EB5AD3875E20226D85D5ABE5F7F5A3735EAC490DBDFF9C3292

5 个 Flow DICOM 均解码为 500 × 351 × 500，2500 帧全部形成连续轨迹：

| 扫描 | 高置信 | 模型辅助 | 跟踪失败 | assessable | uncertain | not_assessable | X4 修正 |
|---|---:|---:|---:|---:|---:|---:|---:|
| flow01 | 338 | 162 | 0 | 490 | 10 | 0 | 0 |
| flow03 | 312 | 188 | 0 | 490 | 7 | 3 | 4 |
| flow05 | 316 | 184 | 0 | 490 | 10 | 0 | 6 |
| flow07 | 328 | 172 | 0 | 480 | 11 | 9 | 7 |
| flow10 | 315 | 185 | 0 | 488 | 12 | 0 | 3 |

前、中、后所选的 15 帧全部为 assessable，X1 均未使用备用段。15 帧全部
通过定位 QC、双侧背景 QC、源区/拖尾窗完整性 QC 和最终指标有效性检查。

相对固定表面历史试运行，新中心的横向变化范围为 -1.94 至 +0.73 pixel；
新 z_upper 比旧单帧顶缘浅 1 至 19 pixel。该差异来自师兄整卷轴向轨迹，
不是表面距离换算。

## 5. 数值与评价边界

X4 会产生小数像素中心。sv_raw 的原始基线可达约 10⁹，直接以“血管均值减
背景均值”计算 T(z) 会放大浮点相消误差。实现改为先逐像素扣除同一背景，
再做加权均值；两式数学等价，且稳定满足 P(z)=有效宽度×T(z)。对应回归
测试覆盖大基线与小数边界。

当前 Flow DICOM 是定位图。正式 SV 指标始终来自原始 OCT 重建的线性
sv_raw。实验身份与部分采集时间元数据仍为空，运行表会如实记录完整性。
当前没有匹配空白剖面，探索性检出长度记录为未运行；Q_vessel、Q_tail 与
R 不受此项影响。

### 5.1 v2 候选试运行边界

v2 对 15 个前/中/后样本全部产生有效几何，原 QC 下仍为 15/15 有效。相对
v1，顶缘变化范围为 -2 至 +12 pixel，中位数为 0 pixel。五个先前目测偏上
的样本分别变化 +4、+5、+12、+12、+1 pixel（正值表示向深部移动）。

阈值矩阵同时测试了峰值比例 0.15/0.20/0.25 与背景噪声倍数 3/4/5。最终
候选为 0.20 与 4；噪声 5 倍在 15 帧中出现 2 帧无局部种子，所以没有采用。
这不是以 QC 通过率挑选最有利结果，旧 QC 代码未改。

全 500 帧检查发现 flow03 第 425–465 帧附近出现向深部偏离的轨迹，最大
相对 v1 为 +46 pixel。该段的合格候选稀疏，并混入少量更深的假候选，随后
被原轨迹补全器连接。它没有覆盖本次选取的 0/249/499 帧，但说明 v2 尚不能
直接替代整卷正式定位。冻结 v1 结果保留不变，v2 以候选结果发布供看图复核。

## 6. 可复现命令

生成一卷定位表：

    svrecttail mentor-track --flow-dicom scan-Flow_ed.dcm --scan-id flow01 --diameter-um 128 --output outputs/mentor_tracking/flow01

生成 v2 候选定位表时显式传入：

    svrecttail mentor-track --flow-dicom scan-Flow_ed.dcm --scan-id flow01 --diameter-um 128 --tracking-config config/tracking_config.upper_edge_v2.a020_n4.json --output outputs/mentor_tracking_upper_edge_v2_a020_n4/flow01

运行 15 帧混合先导：

    svrecttail run --config config/run_config.mentor_tracking.pilot.json --manifest data/manifest.csv --output outputs/pilot_mentor_tracking_15

v2 候选定量使用：

    svrecttail run --config config/run_config.mentor_tracking.upper_edge_v2.pilot.json --manifest outputs/local_run_inputs/manifest_mentor_tracking_upper_edge_v2.csv --output outputs/pilot_mentor_tracking_upper_edge_v2_15

定位表保存 500 帧主表、三个 alpha 轨迹表、输入哈希、师兄代码哈希和有效
配置。原始 OCT 与 DICOM 继续放在仓库外，GitHub 只保存代码、清单、哈希和
可审计的派生结果。
