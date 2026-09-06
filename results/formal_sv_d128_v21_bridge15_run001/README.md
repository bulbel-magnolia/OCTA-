# 128 μm v2.1：原 15 帧 SV 接入包

本包已完成真实 OCT 重建和确定性 SV 定量。原 15 个位置全部保留：13 帧有效，flow01/499 与 flow07/499 无可靠几何，正式指标为 NA，二维原图仍保留。11 帧为直接候选支持、2 帧为短缺口辅助定位。此处帧号为 0 基；MATLAB 帧号加 1。

定位、几何、背景及正式 SV 定义均冻结于基线 `ee3012c04be0d02289520a34409d0855a3e5ca58`。配置复制为 `tracking_config.json` 与 `run_config.json`，原文件分别为 `config/tracking_config.continuity_first_v2_1.a020_n4.json` 和 `config/run_config.mentor_tracking.continuity_first_v2_1.pilot.json`。新增代码版本见 `software_versions.json`；冻结代码、输入和定位表哈希见 `collection_run.json`、`input_sha256.csv`、`frozen_inputs_sha256.csv`。

正式信号是三次位置内重复复数 OCT 的 `var(abs(E),1,3)`，分母 N；无归一化、对数变换或跨帧信号平滑。dx=12.7、dz=6.7 μm。源区采用 X4 中心、X1 宽度和 128 μm 高度的完整椭圆。尾区为源区最深下缘起始的等宽 500 μm 矩形，guard=0；左右各跳过 3 列、取 5 列、均值合并。每个原始深度扣同深度背景，保留负残差与分数像素权重。跟踪配置内旧 2 px/300 μm 窗口不用于正式积分。

`frame_results.csv` 保存全部位置的状态与指标；只纳入 `valid=True`。`diagnostic_*` 是失效帧仍可计算时的诊断值，不能作为正式分析值。`localization.csv` 和 `tracking/` 保留冻结几何及来源；`assessability` 对应原表 `vessel_presence_prediction`。`profiles/` 为完整 351 深度的 V、左右 B、B、T、P、分数尾权重与有效性；无几何帧的这些量为 NA。`sensitivity_results.csv` 保留 x±1、z±1、背景间隔+2 px 的五项扰动。`exclusions.csv` 与 `scan_summary.csv` 分别记录排除和描述性汇总。

`arrays/<scan_id>/frame_NNN.npz` 保存完整 351×500 的 `sv_raw`、`stru_amp`、`omag_raw`、定位 Flow DICOM 像素；有效几何还包括二维 `source_weights`、可外积恢复二维尾权重的 `tail_z_weights`/`tail_x_weights`、背景列索引及 B。数组无损保存，使用 `numpy.load(..., allow_pickle=False)`。`metadata_json` 记录同一帧的来源、几何、重建参数、0 基坐标及 `z,x` 轴顺序。完整 B-scan 的裁剪原点为 [0,0]，对应重建 FFT 的 MATLAB 50:400 行；数组 z=0 对应 FFT 第 50 行。原始 OCT、DICOM 和临时 MAT 不在包内，只有原文件名与 SHA-256。

所有 15 张 `qc/*_mapping.png` 均已检查同帧四图叠加：无需镜像或平移，未见新增坐标偏移及明显背景材料错误。OMAG raw 与显示 DICOM 的数值不相等，因为后者经过显示滤波；原方向的秩相关为 0.647–0.672，均高于翻转诊断。`bridge_validation.json` 记录接入判据，`archive_validation.json` 记录从归档数组独立复算 13 帧积分、15 个哈希与 5265 行曲线通过。这是技术接入验证，不是新的人工定位准确率测量。

无匹配空白，检出深度状态为 `not_evaluated_no_matched_blank`，长度和 detected 为 NA。血管、phantom、session、独立重复身份及时间/慢轴标定未知，清单保留 NA；每卷 500 个空间位置不视为独立血管，3 次位置内重复不视为独立完整采集。辅助定位帧使用该帧自己的 SV，分析端可另做“仅直接支持”稳健性描述，不能按结果选择有利版本。

二维数据采用 GitHub Release 附件，入口为 [formal-sv-d128-v21-run001](https://github.com/bulbel-magnolia/OCTA-/releases/tag/formal-sv-d128-v21-run001)。实际包名、字节数和 SHA-256 见 `download_packages.csv`；逐帧哈希见 `arrays_sha256.csv`。下载后解压即可得到本目录表格、配置、图片和 arrays，路径全部相对。
