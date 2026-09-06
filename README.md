# OCTA SV rectangular-tail quantification

## 当前正式先导路径（2026-09-06）

当前先导分析采用混合方案：师兄的完整 500 帧 Flow 体数据定位链负责寻找
血管，我们的线性 sv_raw、物理椭圆源区、双侧背景和 500 μm 矩形拖尾指标
负责定量。定位依次使用慢轴 Viterbi 轨迹、X1 局部连续宽度、X2 稳健质心、
X4 孤立跳点修正、主 alpha=0.15 的 z_upper 和帧可评估性。此路径不读取
表面位置、折射率或表面到血管的距离。

固定表面 z=176 与 200 μm 距离模式只保留为历史对照。新方案的接口、
坐标转换、真实数据试运行结果和运行命令见
[师兄定位与 SV 指标接入说明](docs/MENTOR_TRACKING_INTEGRATION.md)。
冻结的 15 帧真实数据结果、逐帧 QC 图和 5 卷完整定位轨迹见
[师兄自动定位 × SV 矩形拖尾先导结果](results/pilot_mentor_tracking_15/README.md)。

### 上边界 v2 候选试运行（尚未替代冻结结果）

针对顶缘偏上的问题，仓库增加了版本化的
`persistent_core_paired_edge_v2`：强信号核心只负责确认血管存在，较低的
滞回阈值向上寻找边界，同时要求横向覆盖、连续多行信号、边界内外对比和
按 `dx`/`dz` 换算后的轴向平衡。候选配置把峰值阈值由 15% 提高到 20%，
背景噪声阈值由 3 倍提高到 4 倍；原 QC 判据没有修改，旧 v1 仍是无配置时
的默认方法。

同一批 15 帧试运行全部得到有效指标。先前指出的 1行3列、2行2列、3行1列、
4行1列、4行3列顶缘分别向下移动 4、5、12、12、1 pixel。全 500 帧轨迹
审计同时发现 flow03 约第 425–465 帧存在错误的向下分支，因此 v2 当前是
供人工看图确认的候选结果，不能据此宣称整卷定位已经解决。图、逐帧差异和
复现说明见
[上边界 v2 候选试运行](results/pilot_mentor_tracking_upper_edge_v2_15/README.md)。

本仓库实现 OCTA 散斑方差（SV）矩形拖尾定量流程。MATLAB 保留现有
`.oct` 重建与 OMAG 生成路径，Python 负责局部几何定位、固定双侧背景、
物理面积积分、检出长度、敏感性分析和 QC 输出。

## 方法状态

- 协议版本：`SV_Rectangle_v1_pilot`
- 正式信号：`sv_raw = var(abs(E), 1, 3)`，方差分母为 `N`
- 管径：128 μm
- 标定：横向 12.7 μm/A-line，轴向 6.7 μm/pixel
- 正式拖尾窗：间隔 0 μm，长度 500 μm
- 先导样本框架：5 个流速 × 前/中/后 3 个位置，共 15 帧

当前代码已通过 44 项自动测试，并完成 5 卷 Flow 体数据的完整定位和 15 帧
真实 SV 定量试运行。正式先导配置使用师兄定位表中的 X4、X1 和 z_upper；
清单锚点与固定表面深度模式仅保留为历史对照，不用于当前正式结果。未知的
实验身份及采集时间信息保持为空，结果会标记其完整性。

## 快速开始

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -e ".[test]"

svrecttail validate `
  --config config/run_config.pilot.json `
  --manifest data/manifest_template.csv `
  --allow-template-gaps

svrecttail run `
  --config config/run_config.pilot.json `
  --manifest data/manifest.csv `
  --output outputs/pilot_run_001
```

如有匹配的空白噪声剖面，可增加
`--blank-profiles data/interim/blank_profiles.npz`。NPZ 中每个键使用对应
`scan_id` 或通用键 `all`，数组为经相同背景扣除的 `T(r)`，形状是
“空白样本数 × 拖尾深度行数”。

## MATLAB 中间图导出

```matlab
addpath('matlab');
export_sv_omag_frame('scan.oct', 1, 'data/interim/scan_b000.mat', 2);
```

`bscan_index` 在 MATLAB 中为 1 基；清单中的 `x_anchor_center_px` 与
`z_anchor_center_px` 均为 0 基像素中心坐标。导出文件只包含线性
`sv_raw`、`sv_cv2`、`omag_raw`、`stru_amp` 和重建元数据。正式积分始终只读
`sv_raw`。

`config/run_config.surface_z176.pilot.json` 固定表面参考为 0 基
`z=176`，按 `200 μm × 1.12 / 6.7 μm/px` 得到血管顶缘先验
`z≈209.43`。程序只在该深度带内用 OMAG 自动寻找横向粗定位，再执行局部
X1 边界和顶缘重定位；此模式不读取清单中的旧横向/轴向锚点。

## 运行产物

每次成功运行生成：

- `frame_results.csv`：逐帧 `Q_vessel`、`Q_tail`、`R` 与完整 QC 状态；
- `scan_summary.csv`：以真实 `scan_id` 汇总帧数、失败计数和中位数/IQR；
- `profiles.csv`：逐深度 `V`、`B_left`、`B_right`、`B`、`T`、`P`；
- `profiles/`：逐帧完整 raw profile，避免仅依赖合并表；
- `localization.csv`：X1 边界、中央线顶缘及定位诊断；
- `sensitivity_results.csv`：横向 ±1 px、顶缘 ±1 px、背景间隔 +2 px；
- `detection_results.csv`、`detection_bins.csv`：5 行分箱、连续检出和右删失；
- `arrays/`：逐帧 MAT，含输入图、背景像素、背景校正图、掩膜和全部权重；
  启用检测时另存实际 `T(r)`、全部匹配空白、阈值参数及输入哈希；
- `qc/`：每帧 QC01–QC03 六面板图，检测可用时另存阈值图；
- `logs/`：运行元数据、单帧处理错误和清单预声明的人工调整记录；
- `run_config.json`、`manifest.csv`、`run_complete.json`：冻结输入、哈希和软件版本。

输出先写入被 Git 忽略的 `outputs/`。完成 QC 后，将冻结运行复制到
`data/processed/<run_id>/` 并提交。大型 MAT、HDF5 和 NumPy 二进制文件由
Git LFS 管理。

## 原始数据管理

原始 `.oct`/`.dcm` 保存在仓库外的本地数据目录，不进入 Git 历史；仓库用
`data/raw_inventory.csv` 保存文件名、大小、实验条件和 SHA-256，从而精确确认
每次分析使用的原始文件。机器相关的绝对路径只写入被忽略的
`outputs/local_run_inputs/manifest.csv`，不提交到公开仓库。

如果需要通过 GitHub 分发完整原始数据，应把冻结数据集作为 GitHub Release
附件上传，不要作为普通 Git 文件提交。详细规则见 `data/README.md`。

## 验证

```powershell
python -m pytest -q
matlab -batch "addpath(fullfile(pwd,'matlab')); addpath(fullfile(pwd,'tests','matlab')); test_compute_sv_maps; test_reconstruction_dependencies"
```

完整协议见 `docs/METHOD_SPECIFICATION.md`，实施审计见
`docs/IMPLEMENTATION_AUDIT.md`。`examples/synthetic_demo/` 提供纯合成的数值
CSV 与代表 QC，用于核对输出接口，不作为实验结果。
