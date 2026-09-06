# 15 帧师兄自动定位 × SV 矩形拖尾先导结果

本目录是 2026-09-06 冻结的真实数据试运行。定位使用师兄完整 500 帧慢轴
跟踪，指标使用本项目的线性 sv_raw、物理椭圆源区、双侧背景与 500 μm
矩形拖尾。运行对应代码提交 ac2d2791622464f7fc2c57dfea45bbb731246c01。

## 结论

- 5 卷各 500 帧均成功生成连续定位轨迹；
- 所选 15 个前/中/后 B-scan 全部为 assessable；
- 15/15 帧定位 QC、几何 QC、背景 QC 和 500 μm 窗口 QC 全部通过；
- 15/15 帧得到有效 Q_vessel、Q_tail 与 R；
- 0 帧使用 X1 备用位置，0 次人工调整；
- 定位过程不使用表面位置、折射率或表面到血管距离；
- 当前未提供匹配空白剖面，探索性检出长度明确记录为未运行。

## 看图

[15 帧边界总览](mentor_tracking_localization_overview_15.png)中，绿色竖线为
X4 中央 A-line，青色竖线为按 X1 宽度放置的左右边界，红色横线为 z_upper
上边界，橙色横线为按 128 μm 物理内径计算的下边界。

[5 卷完整轨迹](mentor_tracking_trajectories_5x500.png)显示每卷全部 500 帧的
Viterbi、X4 与 z_upper。黄色圆点是本次选择的 15 帧；橙色小点是整卷中的
uncertain 或 not_assessable 帧。

comparison_vs_surface_guided.csv 保存新旧定位逐帧差异。新 X4 相对旧单帧
中心变化 -1.94 至 +0.73 pixel；新 z_upper 比旧顶缘浅 1 至 19 pixel。

## 文件说明

- frame_results.csv：15 帧正式指标及全部 QC；
- localization.csv：X4、X1 宽度、z_upper、可评估性与最终几何；
- scan_summary.csv：按流速扫描汇总；
- profiles.csv 与 profiles/：完整 V、B、T、P 深度剖面；
- sensitivity_results.csv：横向和轴向 ±1 pixel、背景间隔变化结果；
- qc/：15 张六面板定位与定量 QC 图；
- tracking/：5 卷主定位表、三个 alpha 轨迹表及去除机器绝对路径后的元数据；
- arrays_sha256.csv：本地逐帧 MAT 归档的大小与 SHA-256。MAT 总计约 96 MB，
  未纳入普通 Git 提交；
- run_config.json、manifest.csv、run_complete.json：冻结配置、运行清单、
  代码提交、输入哈希和软件环境。

原始 OCT、Flow DICOM 和线性单帧 MAT 保存在仓库外。data/raw_inventory.csv
记录 10 个原始/定位输入的文件名、大小和 SHA-256。tracking 中的派生 CSV
足以审计 2500 帧定位；正式 SV 数值复算仍需要原始 OCT 或对应线性 MAT。

完整方法边界见 ../../docs/MENTOR_TRACKING_INTEGRATION.md。
