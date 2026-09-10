# SV-OCTA 固定真实直径 ROI v2：选择性分析执行任务

本文件是下一轮的执行计划草案。收到用户明确的执行指令后，在指定本地工作区执行；生成本计划的当前对话不运行科学分析、不修改分析代码、不提交或推送。

## 1. 工作区、输入与只读边界

仓库：bulbel-magnolia/OCTA-
分支：analysis/sv-diameter-stage2
本地工作区：C:/Users/zby/Desktop/2026血管三维重建项目/算法差异分析/_worktrees/OCTA-/sv-diameter-stage2

正式输入目录：analysis/sv_physical_diameter_geometry_v2/
有效政策：analysis/ACTIVE_ANALYSIS_POLICY.json
新输出目录：analysis/sv_fixed_diameter_selective_v2/

输入包括：
- framewise_fixed_diameter_D128.csv.gz、D235、D285；
- geometry_manifest.csv.gz、coverage.csv、cohort_manifest.csv；
- volume_geometry_change_summary.csv、d500_exclusion_record.csv；
- analysis_policy.json、input_manifest.csv、array_identity_audit.csv.gz；
- validation.json、provenance.json、output_sha256.csv。

只读历史目录：
- analysis/sv_diameter_stage2_scientific_v1/
- analysis/sv_source_tail_framewise_coupling_v1/
- analysis/sv_source_tail_depth_coupling_v1/
- analysis/sv_data_quality_qc_v1/

先检查实际工作区、HEAD、工作区差异、政策和文件哈希。允许既有未提交的新版本地成果存在，完整记录其身份。不要自动 checkout、pull、reset、clean；不要假定远端已经包含新版成果。现有新旧结果目录只读，后续代码和新结果写入新目录。新目录若已有冲突成果，停止覆盖并报告。

从 manifest 解析 retained SV 数组的实际位置；不得猜测路径，不下载替代版本，不从 .oct 默认重建。需要读取的数组在读取时与已有身份清单核对；没有身份变化，不安排另一个全量完整性审计前置步骤。

## 2. 固定队列与信号

仅 D128/D235/D285 × 1、3、5、7、10 mm/s，共 15 卷。
每卷保留原始 frame_index=0…499，共 7500 个名义位置。
原有效帧分别为 2422、2225、2371，共 7018；482 个原无效 geometry 位置保持缺失。

D500 全部 8 卷只保留排除记录和历史档案，不读取其数组进行任何后续科学计算，不进入主分析、次分析、敏感性、对照、图或统计汇总。记录该决定发生在查看旧图像/QC之后，属于事后范围调整；不称为预注册排除，不断言无血流或采集失败。

正式 SV_raw=var(abs(IMG),1,3)，方差分母 N，linear raw SV。
Q=sum(SV_raw*w)*dx*dz；area=sum(w)*dx*dz；mean=Q/area。
dx=12.7 μm/pixel，dz=6.7 μm/pixel；不能据此推导 PSF 或 slow-axis 步距。
不以 log、CV²、gain、normalization、clipping、positive truncation、背景相减替代正式信号。

真实 D×D 的 source、D 宽 tail、guard=0、冻结 X4/z_top、由 z_top+D/dz 确定的下界、source supersample=16、6 source bins、20 个 25 μm tail bands 均不改写。下界来自几何先验，不称为独立检测到的下壁。

## 3. 在运行新版科学统计之前固定的分析合同

创建 analysis_plan.json、selection_log.md、column_mapping.json、run_status.json。
标记 informed_by_historical_results=true、independent_validation=false。
记录所有参数、输出定义和条件扩展；不得根据新版排序、相关系数或图形修改合同。

三个核心问题：
1. 三直径在匹配流速条件下的 absolute source/tail 与相对 tail/source 关系。
2. 同卷空间关联的深度分布与零位对应。
3. 匹配局部背景、空间负对照和定位敏感性对上述解释的约束。

主标量：RI500；同时固定报告 RI100、source mean、tail100 mean、tail500 mean，20-band tail mean 与 RI profile。
逐帧 RI_H=T_H/S；S=0 时 NaN，不加 epsilon。

固定候选区域 C={(S3,S4,S5,S6)×(T1,T2,T3,T4)}，共16 cells。
该区域来自历史同数据观察，不能称为未见数据前确定或独立验证。
C16=16 个 lag0 rho 的中位数，保留全部正、零、负值；另存每 cell 值、正值数和可定义数。主 C16 要求16项均可定义，否则输出 NaN 和原因，不以剩余 cells 替代完整区域。
同时报告每个 source bin 对 T1–T4 的4项 rho 中位数，展示全部六层，不自动寻找或替换热点。

固定三个代表配对：
W100：whole-source mean 对 tail0–100 mean；
W500：whole-source mean 对 tail0–500 mean；
C100：pooled source S3–S6 mean 对 tail0–100 mean。
Pooled source mean=sum(Q_S3…S6)/sum(area_S3…S6)，不能直接平均四个 bin mean。
C100 的标量相关不等于 C16 的 cell-rho 中位数，两者分列。

主深度图：每卷 real raw 与 real 51-frame detrended 的完整6×20 lag0 mean–mean图。
输出明确保留source的u边界/中点，以及tail band的上下界/中点（12.5、37.5、…、487.5 μm）；最后一band覆盖475–500 μm，不冒充在500 μm处的点测量。
保留全部15卷图和逐 cell 数值，包括候选区域以外、深部、零值与负值；不得只展示共识图。
跨流速汇总仅展示五个原值、median、min–max、正/负/零及可定义卷数；不得称为同条件五次重复。
直径共识图按cell对该直径五个flow卷的系数取中位数，同时保存全部flow值、min/max、可定义数和符号计数；不是先合并帧或把五个条件视作同条件重复。

## 4. 通用统计与缺失合同

以 geometry_manifest 的原始7500位置为连接骨架。按 scan_id+frame_index 连接，不按行序连接，不把缺失补零、插值或压缩距离。
原 geometry valid flag 永不改变；背景 FOV、扰动 FOV、零分母或统计不可定义均使用独立诊断状态列。

去趋势：在原0–499坐标上，centered moving median，主窗口51，敏感性31/101；端点截短、忽略邻域 NaN、min_periods=1，无效中心仍 NaN。
保存每个窗口内的实际支持数量。残差只用于统计诊断。
Spearman：当前有效配对集合内 average ranks，再 Pearson；至少3对且双方非常量；否则 NaN+原因。不调用输出 p-value 的检验接口。3对只是可计算条件，不是可靠性门槛。

Lag=-50…50，Source(i) 对 Tail(i+lag)，不循环，不把 dx/dz 用来换算 lag。
每个 lag 按原坐标重新配对并保存 n_pairs。
Specificity=rho(0)-median[rho(lag),30≤|lag|≤50]，远端固定42项。rho(0)或任何远端项不可定义，主 specificity 为 NaN，记录原因和可定义 lag 数。
Peak 是最大有符号 rho；完全并列先选最小 |lag|，再选较小有符号 lag。
zero_lag_rank=1+count(rho(lag)>rho(0))。
主 peak/rank 要求全部101项可定义，否则 NaN；曲线中已有值仍保存。
±50边界峰原样报告，不自动扩大 lag 范围，不推断范围外峰位。

51-frame 的三个代表配对另作配对坐标匹配核查：
对每个远端 lag，取 I_lag={i: S(i)、T(i)、T(i+lag)均可用}；
delta_lag=corr[S(i),T(i)|I_lag]-corr[S(i),T(i+lag)|I_lag]；
matched_specificity=42个 delta_lag 的中位数。任一项不可定义则主摘要 NaN。
它是配对组成诊断，不替换历史 specificity，不是显著性检验。

相邻差分只对原坐标 i、i+1 均有效的相邻帧计算：dS=S(i+1)-S(i)，dT同理；仅三个代表配对，使用 raw 差分，不叠加去趋势。

## 5. 阶段一：读取既有 v2 表，完成基础描述

复用已有 real mean/Q/area 和6×20积分，不重新生产这些正式输入。
核对六层 source 的 Q/area 重构 whole-source，前4个 tail bands 重构100 μm，20 bands 重构500 μm。
每卷先计算逐帧 RI，再取卷内中位数；不得用两个卷中位数的比替代 median(T/S)。

输出 source、tail100、tail500、RI100、RI500 的逐卷摘要，以及全20-band raw tail/RI profile。
按相同 flow 比较128–235、128–285、235–285三对直径，展示具体值、差值及分母有效时的百分比变化，不拟合共同流速单调趋势。
保留每卷固定空间块0–99、100–199、200–299、300–399、400–499的基础指标摘要；仅作为空间覆盖/异质性诊断，不当独立重复或跨卷解剖同位配对。

面积理论：source=πD²/4；tail_H=D*H；每个 source bin 面积是固定的几何常数。
保留实际 fractional rasterized area 及其相对解析面积的偏差；不能要求16×16离散圆面积与πD²/4机器精度相等。
一次性核查 raw/51 的三个代表配对中 actual Q–Q 与 mean–mean 的差异；另以理论常量面积乘 mean 验证正比例缩放的相关等价性。
此为数值检查，不输出科学 Q–Q 地图/lag/面积相关，不把面积抖动作为协变量。

## 6. 阶段二：新增 matched-background 像素统计与最小定位 QC

尽量在一次 retained-array 遍历中完成本阶段及第8节所需扰动积分。
已验证 real mean/Q/area 继续从表读取。新读像素仅用于新分布统计、pseudo ROI、新扰动和指定 spot-check。

匹配背景：
x_L=X4-1.5*D/dx；x_R=X4+1.5*D/dx。
左右 source/tail 与 real 连续几何严格同形同深度平移，source仍为D×D，tail仍D宽、0–500 μm、相同 bins；边到边间距0.5D。
旧 ±1.5X1 的位置、像素分布、Z 和 pseudo coupling 数值不得继承为v2。
每侧完整 source+0–500 μm tail 在 FOV 内才可用，不挪回、不缩小、不按亮度删除或换背景。
连续几何平移正确不等于离散像素权重逐元素相同。

新增 real/L/R 的分布统计区域：whole-source、6个source bins、pooled S3–S6、20个tail bands、pooled tail100、pooled tail500。
按 fractional weights 计算 weighted median/MAD；pseudo 另计算 mean/Q/area。
weighted median 为累计权重达到一半处；恰为一半时取相邻取值均值；weighted MAD 围绕对应 weighted median 计算。

主 combined BG 要求该帧双侧均完整可用，拼接左右像素及权重，不平均左右中位数或左右Z。
Z_local=(weighted_median_real-weighted_median_BG)/(1.4826*weighted_MAD_BG)。
MAD=0、空ROI等输出NaN和原因，不加epsilon，不设Z>1/2/3阈值，不称仪器噪声底、检验Z或校准检测深度。
Pooled100/500及pooled source的中位数/MAD从相应完整像素分布计算，不平均分区统计量。
保存 real-minus-BG mean 作为有单位的诊断差，不改写正式SV/RI。

左右各自可用值保留。主双侧对照在共同支持 M=real/L/R均完整有效 的原0–499位置上执行；先将三者按M置缺失，再分别去趋势，保证对照的窗口支持一致。
正式 real-all 结果保留；real-matched 单独保存，不能与不同帧集合的 pseudo rho 直接相减。
单侧可用不能静默代替主双侧背景；输出覆盖情况与独立状态。

背景漂移仅对 whole-source、tail100、tail500 的 L/R/combined 像素中位数空间曲线，输出：
- RCV=1.4826*MAD_across_frames(b)/median_across_frames(b)；
- median(b[449:499])-median(b[0:50]) 的绝对及相对差；
- 两端实际覆盖数和 raw/51趋势曲线。
零分母或空端块为NaN，不移动端块寻找可用帧。
同时保留左右 signed asymmetry=2*(mean_L-mean_R)/(mean_L+mean_R)，分母0为NaN。

像素支持：优先读取已经保存且定义一致的 Neff；缺少 sum(w²) 时由几何权重计算，不用 area 推算。
Neff=(sum w)²/sum(w²)，另存非零权重像素数、面积误差和bin轴向像素跨度。
Neff只是加权像素支持，不是独立散斑数、光学分辨率或统计独立n。
对6×20权重计算 source-bin/tail-band 的共享像素数与 overlap_index=sum(w_s*w_t)/sqrt(sum(w_s²)*sum(w_t²))，作为几何共享指标，不当作测得的噪声相关或校正系数。

图像核查：优先恢复既有45个固定帧ID，每卷首/中/末；无法恢复名单时，按原有效帧排序选择首、floor((n-1)/2)、末帧并记录这是新固定名单。
并列展示 retained stru_amp 与 linear SV、冻结X4/z_top、几何先验下界、D圆和左右pseudo ROI，按真实物理长宽比显示。
同卷统一显示范围，不按结果选亮图。标记“与可见轮廓相容/明显错位/图像不足以判断”，不把先验下界或信号峰当真值。
复用同身份三直径每卷X4/z_top相邻抖动事实；必要时从冻结几何按真正相邻帧恢复，不再做X1筛选。

## 7. 阶段三：主空间关联与负对照

Real：完整6×20 raw/51 lag0地图；固定C16、六层proximal行摘要；三个代表配对raw/51的lag曲线、specificity、peak/rank。
相邻差分和配对坐标匹配specificity按第4节执行。

Pseudo：L/R各一套完整6×20的51-frame lag0地图；三个代表配对各自的51-frame lag与matched_specificity。不跑每个cell的lag。

同一M上逐cell定义 control_rho=(rho_L+rho_R)/2；excess_cell=rho_real_matched-control_rho。
输出全部cell的L/R/real-matched/excess；C16_excess=16个excess_cell的中位数，不用两个区域中位数相减替代。
另存C16的real-minus-L、real-minus-R，防止左右平均掩盖相反判断。
不在real或pseudo地图重新选top10%区域，不以任何新峰点确定对照位置。

增加两个小型共变诊断，均无新积分：
1. 20个同深度real-tail mean与combined BG-tail mean的raw/51 lag0 rho；combined mean由左右Q/area合并。
2. 三个代表配对的real-source对L/R pseudo-tail的51-frame lag0 rho，作为crossed control。

所有相减仅是系数诊断，不是background-corrected正式coupling，不是无偏因果效应。

## 8. 阶段四：限定敏感性

去趋势尺度：仅三个代表配对重算31/101的完整lag；仅固定16 cells重算31/101的lag0。不重跑31/101全深度地图或pseudo全套尺度组合。

位置扰动：全部15卷、全部原有效帧，固定四个单轴刚性平移：
(dx_pixel,dz_pixel)=(+1,0)、(-1,0)、(0,+1)、(0,-1)。
物理幅度分别12.7 μm横向、6.7 μm轴向；每种偏移在整卷一致，不采样随机抖动，不搜最优位置。
横向移动source和tail；轴向移动z_top、圆心、所有source bins、先验bottom和tail起点。D、guard、window、supersample不变。
仅新增real whole-source/6 bins、tail T1–T4、tail100/500所需积分；不生成扰动版20-band全套分析，不再积分其pseudo背景。

输出每个偏移的source/tail/RI卷摘要、三个代表配对51-frame rho与specificity、16-cell候选摘要及逐cell差值。
以baseline与四偏移共同可用的原坐标掩码作配对比较，重新在该共同掩码去趋势；正式baseline-all不变。
扰动越FOV仅标记该诊断缺失，不能把原有效帧变为无效。
四个方向全部报告，不选择最有利方向；它们检验小幅系统偏移敏感性，不给出定位误差真值，也不代表逐帧随机抖动已验证。

## 9. 阶段五：限定新旧几何比较

只比较相同三直径、15卷、7018同身份有效帧。
直接复用既有 volume_geometry_change_summary 的迁移摘要，不重新计算已存在数值。
仅补充其中缺少的RI100/500及指定coupling对照：旧版已保存的whole-source到tail100/500的raw/51 rho，以及从旧6×20已保存rho提取的固定16-cell摘要。
不重跑旧积分、不恢复旧high region、不比较旧/new背景Z或excess来声称单一real ROI改动的效应。

逐帧新旧比值仅在所需分母非零时计算；可核对 RI_new/RI_old=(T_new/T_old)/(S_new/S_old)。
不得把该逐帧恒等式直接套到分别求得的卷中位数，也不量化“分母解释百分比”。
主科学结果一律来自v2；新旧差异标记为同数据几何敏感性，不是独立验证或优劣评分。

## 10. 只允许以下触发式扩展

每项先记录 trigger、scan_id、相关表行/帧ID、固定扩展内容和状态；原主分析不改动。

A. 新增权重检测到source–T1共享像素（overlap_index>0）：从现有表补S3–S6×T2–T4的12-cell摘要和pooled S3–S6对tail25–100的51-frame诊断。基线T1保留，不修改guard，不将该检查声称为排除PSF串扰。

B. 图像核查记录明确错位或可见界面穿越，且现有45帧不足以判读：对受影响卷补最大|相邻ΔX4|及最大|相邻Δz_top|对应的两帧对；并列取最早帧。不修改定位或新增帧筛除。

C. 图像记录背景结构不匹配，或同卷C16 real-minus-L与real-minus-R异号、提示结论依赖背景侧：先检查FOV/覆盖/几何。若两侧完整可用，再对受影响卷一次性增加±2.5D对照，固定同形同深度规则，仅做固定区域/三个代表量和全深度背景对比；不搜索其他距离，不替换±1.5D主结果，不按亮度挑选。

D. 图像显示结构界面/梯度随ROI位置变化并涉及source或tail解释：对受影响卷补固定D的real/L/R structural mean，限whole-source、tail100、tail500；计算同ROI structural–SV的raw/51相关及结构source–tail相关。只作描述，不算仪器SNR，不做因果解释比例或复杂partial模型。

E. 只有在报告将依赖精确单个source bin的位置，而规定扰动中该位置发生变化时，补一套6→3合并bin诊断：相邻两层Q和area先求和再取mean，51-frame lag0，完整20个tail bands；该分箱直接从表恢复。它只检验更粗轴向描述，不提高分辨率。不启动4/8 bins。

出现窗口敏感、负相关、excess≤0、旧排序消失或边界lag峰，按观察报告，不以此为理由继续搜索参数或扩大lag。
本轮不启动normalized-tail、4/8 bins、D500、逐cell全套lag、高区重选、map相似度竞赛、留一流速“独立验证”、QC回归或机制模型。

## 11. 验证与完成标准

复用没有变更的身份/完整性/覆盖证据；只对新模块补核查。
输入及旧结果前后哈希一致，队列、原valid身份和正式几何不变。
新增积分以同16×16定义在固定45帧独立核查real参考、pseudo形状/FOV、分层重构；同定义非零积分相对误差容限1e-10，零参考单列检查。数值容限不是生物可靠性阈值，解析圆面积近似误差另列。
weighted median/MAD独立测试包括非均匀权重、正好半权重、零MAD、空ROI；不引入epsilon。
统计独立核查至少覆盖各直径一卷和所有新统计模块：原坐标lag配对、缺失内部/端点、ties、constant、已知shift、positive scaling、common-mask对照、差分。相同统计定义两种实现的rho差容限1e-12；失败停止受影响模块并记录。

主地图每模式15×120个cell均占位；NaN保留原因。所有15卷、全部flow、深部/负值/对照结果均交付。图与表的数值/方向/单位核对，相关图固定[-1,1]、excess图固定[-2,2]；不把空间离散范围画成独立样本置信区间。

输出至少包括以下文件组，允许合理长表拆分但不漏定义：
- analysis_plan.json / selection_log.md / column_mapping.json / run_status.json；
- input_manifest.csv / provenance.json / validation.json / output_sha256.csv；
- volume_metrics / diameter_contrasts_by_flow / raw_depth_profiles / spatial_block_summary；
- depth_maps_raw51 / fixed_region_summary / whole_coupling / representative_lag / adjacent_difference；
- framewise_local_background_qc / local_contrast_profiles / pseudo_maps / fixed_region_control_excess / crossed_controls / background_drift；
- geometry_support / shared_pixel_support / localization_review及固定抽查图；
- position_sensitivity / detrend_sensitivity / q_mean_numerical_check；
- geometry_version_comparison / conditional_extension_log / README。

README逐项回答三个核心问题，分别列出保持、减弱、反向或不可定义的历史候选观察；不写机制Discussion。
未取得retained数组时，继续完成表驱动阶段；像素背景/扰动标记BLOCKED，不虚构已完成，也不改用旧QC数值填补。
没有跨scan采集设置记录时，absolute SV按记录仪器单位报告，标注设置一致性未核实；不自动归一化。
slow-axis物理步距未知时使用帧单位；无独立真值/PSF/噪声参考时不补造。

验收依据是定义、覆盖、身份、数值正确性和完整呈现，不是系数高低、符号或旧排序是否保留。
独立实验单位是scan volume；B-scan、cell、空间段、10对flow maps均不得提升独立n。禁止pooled-frame correlation、frame-level p-value、bootstrap帧当独立样本、复杂模型、因果归因。
完成后停在用户审核，不执行git commit/push、远端发布或论文发布。
