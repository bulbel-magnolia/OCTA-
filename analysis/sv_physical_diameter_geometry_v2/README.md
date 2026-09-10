# 固定真实直径 ROI v2：已批准的几何与积分输入更新

2026-09-10 用户批准将 source 改为固定真实直径，并明确选择“D500 不再进入后续科学分析，但保留数据及排除记录”。本目录已落实该决定；分析入口政策为 `analysis/ACTIVE_ANALYSIS_POLICY.json`。

## 新版范围

D128、D235、D285，各共同流速1/3/5/7/10 mm/s，共15卷、7,018个原有效帧。7,500个名义位置全部保留，482个原无效geometry位置仍为缺失。没有依据新结果恢复或新增筛除任何帧。

D500全部8卷（含2/9/12 mm/s扩展）不进入后续科学分析；原始数据、3,913个有效帧的已有结果、QC、截图及排除记录全部保留。排除决定发生在查看已有图像、QC及coupling之后，属于事后分析范围调整，不是预注册排除，也不表示仪器采集失败已被证实。后续结论的直径范围限于128–285 μm。

## 已执行的几何更新

- 保留每个原有效帧的冻结X4和z_top。
- Source物理宽度=真实D、物理高度=真实D，圆形截面先验；像素坐标中的椭圆宽D/12.7、高D/6.7，不取整。
- 圆心z=z_top+D/(2×6.7)，左右边界X4±D/(2×12.7)。下界仍为z_top+D/6.7。
- Tail横向宽度同步为真实D，从同一物理下界开始，guard=0，无扩散锥；保留0–100/0–500 μm整体窗口和20个25 μm bands。
- Source保留6个归一化轴向bins和16×16 fractional weighting。
- X1原值改以original_apparent_width_um保存，作为原算法的表观宽度记录；不将其冒充真实直径。
- 正式SV仍为var(abs(IMG),1,3)，分母N，linear raw SV，无gain/log/normalization/background subtraction。

## 完成范围与后续边界

已从同身份raw arrays完成15卷新版source/tail mean/Q/area与6×20深度输入的重新积分，另保存逐卷新旧差值和有效像素支持。旧版两个coupling目录、旧QC及原几何结果均未改写。

本次没有重跑新的coupling或新的matched-background QC。旧coupling与旧QC属于X1几何，不能贴上v2标签继续使用；后续科学分析应读取本目录输入并按v2重算。旧版本脚本为复现旧结果而保留，不应直接作为新版科学分析入口。固定D也不自动修正上边界或中心的潜在误差。

## 为什么记录D500为范围调整

此前旧几何QC中，D500五个共同流速的S1 robust Z中位数为6.428，S6为0.612，whole-source为1.516，proximal tail为0.258；上部与深部可检测性存在明显差异。该事实与上部亮信号集中的观察相容，不等同于只有上部有血流，不能证明整组数据无效。用户决定将其排除于后续分析；排除原因和时点在cohort_manifest及d500_exclusion_record中公开记录。

## 覆盖与宽度

| diameter_um | volumes | nominal_frames | valid_frames |
| --- | --- | --- | --- |
| 128 | 5 | 2500 | 2422 |
| 235 | 5 | 2500 | 2225 |
| 285 | 5 | 2500 | 2371 |

| diameter_um | median | min | max | new_width_um |
| --- | --- | --- | --- | --- |
| 128 | 177.79999999999998 | 76.19999999999999 | 254.0 | 128 |
| 235 | 342.9 | 228.6 | 431.79999999999995 | 235 |
| 285 | 406.4 | 241.29999999999998 | 508.0 | 285 |

## 验证与复现

所有7,018个原有效帧的数组SHA256与正式输入逐一一致，D128所用25个release ZIP验证通过。先用同一数组重放原source/tail mean/Q核对输入，再积分新版ROI。新旧X4、z_top、物理bottom一致，source/tail均完整落在FOV内，新增无效帧=0。45个预先固定帧（每卷首、中、末有效帧）的全图source权重独立对照通过；各项数值误差见validation.json。

```powershell
python analysis/sv_physical_diameter_geometry_v2/build_fixed_diameter_inputs.py
python analysis/sv_physical_diameter_geometry_v2/finalize_geometry_migration.py
```

主要文件：analysis_policy.json与全局ACTIVE_ANALYSIS_POLICY.json；cohort_manifest.csv；d500_exclusion_record.csv；geometry_manifest.csv.gz；三个framewise_fixed_diameter_D*.csv.gz；volume_geometry_change_summary.csv；array_identity_audit.csv.gz；input_manifest.csv；validation.json；provenance.json；output_sha256.csv。
