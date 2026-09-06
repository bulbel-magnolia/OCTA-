# 128 μm SV v2.1 数据交接（run001）

本轮按已批准的采集计划完成原 15 帧接入和五卷完整 2500 帧提取。冻结科学基线为 `ee3012c04be0d02289520a34409d0855a3e5ca58`；完整新增采集/校验代码快照为 `c967409`；全卷数据提交为 `e1eb0d71536df7c93c821eb0da928eee2576753d`。逐帧实际执行代码版本另存于 frame_results.csv 和 software_versions.json；不同提交之间只有批量调度、校验、输出和打包的增补，冻结科学代码哈希均一致。

- [全卷 README 与文件目录](https://github.com/bulbel-magnolia/OCTA-/tree/e1eb0d71536df7c93c821eb0da928eee2576753d/results/formal_sv_d128_v21_full2500_run001)
- [2500 行主结果表](https://github.com/bulbel-magnolia/OCTA-/blob/e1eb0d71536df7c93c821eb0da928eee2576753d/results/formal_sv_d128_v21_full2500_run001/frame_results.csv)
- [字段说明](https://github.com/bulbel-magnolia/OCTA-/blob/e1eb0d71536df7c93c821eb0da928eee2576753d/results/formal_sv_d128_v21_full2500_run001/DATA_DICTIONARY.md)
- [原 15 帧结果](https://github.com/bulbel-magnolia/OCTA-/tree/3d226d8/results/formal_sv_d128_v21_bridge15_run001)
- [真实二维数组 Release](https://github.com/bulbel-magnolia/OCTA-/releases/tag/formal-sv-d128-v21-run001)
- [25 个全卷包的名称、大小和 SHA-256](https://github.com/bulbel-magnolia/OCTA-/blob/e1eb0d71536df7c93c821eb0da928eee2576753d/results/formal_sv_d128_v21_full2500_run001/download_packages.csv)

| 卷 | 计划 | 有 z | 几何通过/定量有效 | 排除 | 纳入直接支持 | 纳入短缺口辅助 |
|---|---:|---:|---:|---:|---:|---:|
| flow01 | 500 | 498 | 486 | 14 | 320 | 166 |
| flow03 | 500 | 476 | 468 | 32 | 300 | 168 |
| flow05 | 500 | 494 | 491 | 9 | 325 | 166 |
| flow07 | 500 | 499 | 493 | 7 | 327 | 166 |
| flow10 | 500 | 486 | 484 | 16 | 302 | 182 |
| 合计 | 2500 | 2453 | 2422 | 78 | 1574 | 848 |

定位覆盖率 98.12%，定量纳入比例 96.88%；两者均以全部 2500 计划帧为分母。这些比例不是人工测得的定位准确率。47 帧无 z，另 31 帧未通过几何 QC；其中 flow03/431 同时存在源信号非正，不重复计数。所有排除帧保留原始二维矩阵及逐帧记录。原 15 帧有 13 帧定量有效；flow01/499 与 flow07/499 保持 NA，没有替换抽样位置。

已完成 59 项 Python 测试、2 项 MATLAB 重建/信号测试；2500 个 NPZ 哈希和帧身份核对、2453 帧独立二维积分复算、877500 行完整深度曲线检查全部通过。首包与全卷的原 15 帧指标一致，原始定位文件及全部冻结坐标核对通过。没有修改追踪器、SV 定义、区域、背景估计或 QC 标准。

首包为 75,258,930 字节，已从公开 GitHub 下载链接完整下载并校验：`efe23ab8ffb215995796eaed9c64fe3ad4c7f1d05458c40a1a84c1486d21d707`。全卷 25 包共 10,900,978,108 字节。每包核对 GitHub 服务器完整 SHA-256，并从公开链接读取首尾字节与本地 ZIP 比较；没有将全部 10.9 GB 再完整下载一遍。最终远端验证记录见全卷目录的 release_assets_verified.csv 和 release_validation.json。GitHub Release 的初始 tag 保留在首包发布提交，完整数据提交由本页明确链接，未改写该 tag 或合并 main。

分析端先读主表、字段说明与 README，以 valid=True 为主结果；保留“包含合格辅助定位”与“仅直接候选支持”的预定比较。每帧使用自己的 SV，未跨帧平滑信号或指标。无匹配空白，检出深度为 not_evaluated，长度和 detected 保持 NA。实验独立单位和时间/慢轴标定未知处仍为 NA。未进行显著性检验、参数筛选或流速总体推断；后续统计解释交给分析端。
