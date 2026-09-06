# 数据管理

raw_inventory.csv 以 asset_type 区分 raw_oct 与
localization_flow_dicom。当前 5 个原始 OCT 和 5 个定位 DICOM 均记录精确
字节数与 SHA-256；DICOM 只用于师兄自动定位，不参与 SV 指标积分。

## 默认方案

原始 `.oct` 和配套 `.dcm` 保存在仓库外，不进入 Git 历史。仓库中的
`raw_inventory.csv` 记录相对文件名、字节数、实验条件和 SHA-256；运行前可用
这些字段确认本机文件与项目输入完全一致。

- `outputs/local_interim/`：由原始 `.oct` 重建的线性图，仅本机使用；
- `outputs/local_run_inputs/`：包含本机绝对路径的运行清单，仅本机使用；
- `data/processed/<run_id>/`：通过 QC 后可提交的冻结配置、表格、剖面和图；
- 大型 MAT/HDF5/NumPy 处理产物：确需版本化时使用 Git LFS。

每个已提交运行都必须包含实际使用的 manifest、运行配置和原始文件 SHA-256。
`.gitignore` 已阻止误提交 `.oct`、`.dcm` 及本机中间数据。

## 在 GitHub 上发布完整原始数据

当前 5 个 `.oct` 各约 1.43 GiB，总计约 7.15 GiB。普通 GitHub 仓库会拒绝
超过 100 MiB 的单文件。Git LFS 虽可容纳这些单文件，整批数据会占用大部分
免费 LFS 存储与下载流量。

需要在 GitHub 分发时，优先创建单独的 GitHub Release，把每个 `.oct` 作为
Release asset 上传。每个文件必须小于 2 GiB；不要把这些文件加入提交或分支。
发布前确认仓库公开范围、数据许可和脱敏状态。

GitHub 官方限制：

- https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-large-files-on-github
- https://docs.github.com/en/repositories/releasing-projects-on-github/about-releases
- https://docs.github.com/en/billing/concepts/product-billing/git-lfs
