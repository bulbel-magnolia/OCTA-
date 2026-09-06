# 128 μm v2.1：五卷完整 SV 数据

已完成 5×500=2500 帧的真实 OCT 重建、固定区域定量和二维归档。定量有效 2422 帧，排除 78 帧；全部计划帧均在主表和数组包中。检出深度未运行，原因是没有匹配空白。

| 扫描卷 | 计划 | 有 z | 几何通过 | 定量有效 | 排除 | 纳入直接支持 | 纳入辅助定位 |
|---|---:|---:|---:|---:|---:|---:|---:|
| flow01 | 500 | 498 | 486 | 486 | 14 | 320 | 166 |
| flow03 | 500 | 476 | 468 | 468 | 32 | 300 | 168 |
| flow05 | 500 | 494 | 491 | 491 | 9 | 325 | 166 |
| flow07 | 500 | 499 | 493 | 493 | 7 | 327 | 166 |
| flow10 | 500 | 486 | 484 | 484 | 16 | 302 | 182 |

直接候选支持与短缺口补全的总数、纳入数分别见 `scan_summary.csv`。定位总覆盖率与定量有效率分别计算，分母始终为计划的 2500 帧。`exclusions.csv` 保留全部排除帧；`exclusion_detail_counts.csv` 将 source QC 失败展开到原定位原因。

定位、几何、背景及正式 SV 定义均冻结于基线 `ee3012c04be0d02289520a34409d0855a3e5ca58`。配置复制为 `tracking_config.json` 与 `run_config.json`，原文件分别为 `config/tracking_config.continuity_first_v2_1.a020_n4.json` 和 `config/run_config.mentor_tracking.continuity_first_v2_1.pilot.json`。新增代码版本见 `software_versions.json`；冻结代码、输入和定位表哈希见 `collection_run.json`、`input_sha256.csv`、`frozen_inputs_sha256.csv`。

正式信号是三次位置内重复复数 OCT 的 `var(abs(E),1,3)`，分母 N；无归一化、对数变换或跨帧信号平滑。dx=12.7、dz=6.7 μm。源区采用 X4 中心、X1 宽度和 128 μm 高度的完整椭圆。尾区为源区最深下缘起始的等宽 500 μm 矩形，guard=0；左右各跳过 3 列、取 5 列、均值合并。每个原始深度扣同深度背景，保留负残差与分数像素权重。跟踪配置内旧 2 px/300 μm 窗口不用于正式积分。

`frame_results.csv` 保存全部位置的状态与指标；只纳入 `valid=True`。`diagnostic_*` 是失效帧仍可计算时的诊断值，不能作为正式分析值。`localization.csv` 和 `tracking/` 保留冻结几何及来源；`assessability` 对应原表 `vessel_presence_prediction`。`profiles/` 为完整 351 深度的 V、左右 B、B、T、P、分数尾权重与有效性；无几何帧的这些量为 NA。`sensitivity_results.csv` 保留 x±1、z±1、背景间隔+2 px 的五项扰动。`exclusions.csv` 与 `scan_summary.csv` 分别记录排除和描述性汇总。

`arrays/<scan_id>/frame_NNN.npz` 保存完整 351×500 的 `sv_raw`、`stru_amp`、`omag_raw`、定位 Flow DICOM 像素；有效几何还包括二维 `source_weights`、可外积恢复二维尾权重的 `tail_z_weights`/`tail_x_weights`、背景列索引及 B。数组无损保存，使用 `numpy.load(..., allow_pickle=False)`。`metadata_json` 记录同一帧的来源、几何、重建参数、0 基坐标及 `z,x` 轴顺序。完整 B-scan 的裁剪原点为 [0,0]，对应重建 FFT 的 MATLAB 50:400 行；数组 z=0 对应 FFT 第 50 行。原始 OCT、DICOM 和临时 MAT 不在包内，只有原文件名与 SHA-256。

原 15 帧接入判据已通过，见同级 bridge15 结果目录。全卷与首包的原 15 帧指标逐项一致，2500 行冻结坐标与原定位表一对一核对通过，见 `frozen_mapping_validation.json`。`archive_validation.json` 保存独立从 NPZ 二维矩阵重算积分和背景的检查结果。`preselected_qc_manifest.csv` 列出全部 91 个固定 QC 位置，包括原 15 帧、均匀位置、flow03 原错误分支与长缺口/重锁定邻域。

打开 `qc_gallery.html` 查看固定图及每帧纳入状态。有几何但未通过 QC 的帧仍画模板用于诊断，画出椭圆并不表示纳入分析。已复看 flow03 第 420、423、432、434、436、447、469–471 帧的同帧叠加；没有新增坐标映射错误。第 432、447 帧因冻结 assessability 不通过被排除。此轮没有重估定位准确率或调整核心追踪。

所有输入信号均来自该帧自身的三次重复，辅助定位只影响坐标。无匹配空白，`detected` 和长度保留 NA；不存在“长度为零”的未运行结果。实验身份、独立重复、时间及慢轴物理标定未知处保留 NA；每卷 500 个空间位置不作为 500 次独立实验。该批支持空间描述，未进行显著性检验或参数筛选。

大数组按每卷每 100 帧一个 ZIP，共 25 包。每包含该段的完整二维数据、相应 profile/QC 和共同表格/冻结配置；全批表格为全 2500 行，数组只含包名对应段。Release 入口为 [formal-sv-d128-v21-run001](https://github.com/bulbel-magnolia/OCTA-/releases/tag/formal-sv-d128-v21-run001)。逐包名称、大小和 SHA-256 在 `download_packages.csv`，逐帧在 `arrays_sha256.csv`。首包 bridge15 是另一个自包含包，不替换全卷数据。

`output_sha256.csv` 和 `frozen_inputs_sha256.csv` 对应 ZIP 内实际字节；Git 可能规范化文本换行，Git 浏览文件的字节哈希因此可能不同。冻结算法同时记录实际字节和 Git LF 哈希。`collection_run.json` 是准备阶段的冻结快照，最终执行状态以 `execution_status.json` 和验证报告为准。原 OCT/DICOM/MAT 留本地，不上传原始采集文件。
