# OCTA SV rectangular-tail quantification

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

当前代码已通过合成测试和 5 个历史重建帧的端到端冒烟测试。15 帧正式
先导运行需要在 `data/manifest_template.csv` 中补齐原始文件、B-scan、定位
锚点和重建版本；未知的实验身份及采集时间信息保持为空，结果会标记其完整性。

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
`data/processed/<run_id>/` 并提交。`.oct`、`.mat`、HDF5 和 NumPy 二进制
文件由 Git LFS 管理。

## 验证

```powershell
python -m pytest -q
matlab -batch "addpath(fullfile(pwd,'matlab')); addpath(fullfile(pwd,'tests','matlab')); test_compute_sv_maps; test_reconstruction_dependencies"
```

完整协议见 `docs/METHOD_SPECIFICATION.md`，实施审计见
`docs/IMPLEMENTATION_AUDIT.md`。`examples/synthetic_demo/` 提供纯合成的数值
CSV 与代表 QC，用于核对输出接口，不作为实验结果。
