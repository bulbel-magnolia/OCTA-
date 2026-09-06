# 数据清单填写规范

`data/manifest_template.csv` 已展开为 15 个先导帧。每一行对应一个实际
B-scan 帧，不以文件名推断管径、流速、血管身份或采集顺序。

## 计算必填字段

- `scan_id`：原始扫描标识，同一扫描的三个 B-scan 共用该值；
- `source_file`：MAT/HDF5/NPZ 中间图路径，可相对清单所在目录；
- `diameter_um`、`flow_speed_mm_s`、`dx_um`、`dz_um`：物理条件与标定；
- `position_label`：`front`、`middle` 或 `rear`；
- `bscan_index`：建议使用 0 基索引；
- `reconstruction_version`：生成中间图的代码或提交版本；
- `x_anchor_center_px`、`z_anchor_center_px`：使用 `manifest_anchor` 模式时，
  填写导师轨迹给出的 0 基像素中心；固定表面自动模式允许这两列留空；
- `geometry_source`：定位来源，例如 `mentor_tracking` 或
  `fixed_surface_global_x`。

正式验证同时确认 `source_file` 实际存在。

## 信息字段与未知值

`vessel_id`、`session_id` 应在来源可确认时填写，`phantom_id` 在已知时填写。
`slow_axis_position_um`、`temporal_repeat_id`、`temporal_repeat_count`、
`scan_time_interval_s` 和 `acquisition_order` 按真实采集记录填写。未知信息留空，
不从文件名推断，也不为每个 B-scan 人工创建不同 `vessel_id`。批处理在逐帧表
和运行元数据中报告身份及采集元数据完整性。

若一侧背景存在明确空间或结构问题，同时填写：

- `background_excluded_side`：`left` 或 `right`；
- `background_exclusion_reason`：可审计的事实原因。

两个字段必须同时出现。留空时执行固定双侧背景，程序不依据亮暗自动选边。
`notes` 记录重扫、运动、遮挡等其他帧级事实。

## 索引转换

MATLAB 导出函数接收 1 基 `bscan_index`。进入清单时，横向和轴向定位锚点
统一转换为 0 基像素中心：

```text
x_manifest = x_matlab - 1
z_manifest = z_matlab - 1
```

转换只执行一次。批处理不自动猜测索引基数。

一条测量记录由 `scan_id + bscan_index + temporal_repeat_id` 唯一确定。没有
时间重复时 `temporal_repeat_id` 留空。

## 运行前检查

```powershell
svrecttail validate `
  --config config/run_config.pilot.json `
  --manifest data/manifest.csv
```

正式验证会拒绝不存在的源文件，以及缺失物理参数、流速、B-scan、重建版本
或几何来源的行，并检查每行的管径和像素标定是否与冻结配置一致。
`manifest_anchor` 模式还会拒绝缺失定位锚点的行。
