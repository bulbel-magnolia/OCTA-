"""Generate data-only report from verified tables; no mechanism or exclusion rules."""
import json
import numpy as np
import pandas as pd
from derive_qc import OUT,ROOT,DEP,ID,read,csv,js,digest
DS=[128,235,285,500]
def fmt(v,n=3):return f'{v:.{n}f}' if np.isfinite(v) else 'NaN'
def mm(x):
    x=pd.Series(x).dropna();return f'{x.median():.3f} [{x.min():.3f}, {x.max():.3f}]'
def table(df):
    lines=['| '+' | '.join(map(str,df.columns))+' |','| '+' | '.join(['---']*len(df.columns))+' |']
    lines+=['| '+' | '.join(str(v).replace('|','/') for v in row)+' |' for row in df.itertuples(index=False,name=None)]
    return '\n'.join(lines)
def main():
    val=json.loads((OUT/'validation.json').read_text());assert val['status']=='passed'
    v=read(OUT/'volume_qc_summary.csv');v=v[v.grid_scope.eq('common_grid')]
    r=read(OUT/'volume_region_detectability_summary.csv');r=r[r.grid_scope.eq('common_grid')]
    n=read(OUT/'negative_control_volume_summary.csv');n=n[n.grid_scope.eq('common_grid')]
    bc=read(OUT/'background_real_tail_correlation.csv');bc=bc[bc.grid_scope.eq('common_grid')&bc['mode'].eq('detrended_51')&bc.side.eq('bg')]
    bgvol=bc.groupby(ID).rho.median().reset_index();csv('background_tail_volume_median_correlation.csv',bgvol)
    j=read(OUT/'geometry_jitter_summary.csv');j=j[j.grid_scope.eq('common_grid')]
    a=read(OUT/'qc_vs_coupling_summary.csv');raw=read(OUT/'raw_sv_array_audit.csv.gz');extra=read(OUT/'additional_invalid_geometry_raw_audit.csv.gz')
    allraw=pd.concat([raw,extra]);csv('raw_sv_volume_audit_summary.csv',allraw.groupby(ID).agg(arrays=('frame_index','size'),pixels=('n_pixels','sum'),finite_min=('finite_fraction','min'),nan_count=('nan_count','sum'),inf_count=('inf_count','sum'),negative_count=('negative_count','sum'),min=('min','min'),max=('max','max'),zero_fraction_mean=('zero_fraction','mean')).reset_index())
    def by(col):return '；'.join(f'D{d} {mm(v[v.diameter_um.eq(d)][col])}' for d in DS)
    def region(name,col='z_median'):return '；'.join(f'D{d} {mm(r[r.diameter_um.eq(d)&r.region.eq(name)][col])}' for d in DS)
    def neg(col):return '；'.join(f'D{d} {mm(n[n.diameter_um.eq(d)][col])}' for d in DS)
    sourcebins=[];prox=[];main=[];negs=[];jrows=[]
    for d in DS:
        g=v[v.diameter_um.eq(d)];rg=r[r.diameter_um.eq(d)];ng=n[n.diameter_um.eq(d)]
        sourcebins.append({'D':d}|{f'S{i}':mm(rg[rg.region.eq(f's{i}')].z_median) for i in range(1,7)})
        prox.append({'D':d}|{label:mm(rg[rg.region.eq(name)].z_median) for name,label in [('t01','0–25'),('t02','25–50'),('t03','50–75'),('t04','75–100'),('tail100','pooled 0–100')]})
        main.append({'D':d,'source Z':mm(g.source_z),'source SBR':mm(g.source_sbr),'BG RCV':mm(g.background_rcv),'BG relative drift':mm(g.background_drift),'valid fraction':mm(g.valid_fraction)})
        negs.append({'D':d,'real high':mm(ng.real_high_median),'L at real high':mm(ng.left_high_at_real_region),'R at real high':mm(ng.right_high_at_real_region),'cellwise excess high':mm(ng.excess_high_median),'real > both fraction':mm(ng.high_fraction_real_gt_both),'BG–tail rho':mm(bgvol[bgvol.diameter_um.eq(d)].rho)})
        for par in ['X4','X1','z_top']:
            z=j[j.diameter_um.eq(d)&j.parameter.eq(par)&j.unit.eq('pixel')]
            jrows.append({'D':d,'geometry':par,'median |delta| pixel':mm(z['median']),'p95 |delta| pixel':mm(z.p95),'median |delta| um':mm(j[j.diameter_um.eq(d)&j.parameter.eq(par)&j.unit.eq('um')]['median']),'flag fraction':mm(z.outlier_flag_fraction)})
    qs={}
    qs['Q1']='Whole-source robust Z 的五卷 median [min,max]：'+by('source_z')+'。'
    qs['Q2']='D500 与 D128/D235/D285 的五个匹配 flow 比较，source Z 均更低（各 5/5）；具体值见 Q1。'
    qs['Q3']='直接从完整 0–100 μm ROI 计算的 pooled tail robust Z：'+by('proximal_z')+'。四个 25 μm band 见下表。'
    hot=[]
    for d in DS:
        g=v[v.diameter_um.eq(d)];hot.append(f'D{d} {mm(g.high_tail_z_median)}')
    qs['Q4']='将既有每卷 real high-region cells 对应的 tail-band volume Z 取中位数（重复 tail bin 按入选 cell 保留），结果为 '+'；'.join(hot)+'。这是相对局部背景的描述性分离度，Z>3 未经校准，不能据此称为“可靠检测”或统计显著。'
    qs['Q5']='T1（0–25 μm）为 '+region('t01')+'；T20（475–500 μm）为 '+region('t20')+'。完整 20 层见 Figure 2；不假定逐层严格单调。'
    qs['Q6']='S1–S6 的五卷 Z 见下表；各直径最高 median source-bin 为 '+ '；'.join(f'D{d} S{int(r[r.diameter_um.eq(d)&r.region.str.fullmatch("s[1-6]")].groupby("region").z_median.median().idxmax()[1:])}' for d in DS)+'。'
    qs['Q7']='每卷最大 coupling cell 所在 source bin 与最大 source detectability bin 完全一致的卷数：'+'；'.join(f'D{d} {int((v[v.diameter_um.eq(d)].source_peak_coupling_bin==v[v.diameter_um.eq(d)].source_peak_detectability_bin).sum())}/5' for d in DS)+'。两种位置不作等同或因果解释。'
    qs['Q8']='Source matched BG median curve 的 RCV：'+by('background_rcv')+'；首尾 51 原始帧块 median 差 / 全卷 median：'+by('background_drift')+'。'
    qs['Q9']='Valid geometry fraction：'+by('valid_fraction')+'。三项 jitter 的 median 与 p95 见表，仅使用真正相邻原始 indices；没有以 QC 删除帧或卷。'
    qs['Q10']='每卷 20 个 real-tail / combined-BG-tail detrended rho 的中位数，再汇总五卷：'+'；'.join(f'D{d} {mm(bgvol[bgvol.diameter_um.eq(d)].rho)}' for d in DS)+'。原始和 detrended 各 band/侧结果均保留。'
    qs['Q11']='Pseudo whole-map median rho：L '+neg('left_map_median')+'；R '+neg('right_map_median')+'。Real-vs-L map Spearman similarity：'+neg('left_similarity')+'；Real-vs-R：'+neg('right_similarity')+'。'
    qs['Q12']='Real high region 的 median rho：'+neg('real_high_median')+'；同一组 cells 的逐 cell real-minus-control，再取 median：'+neg('excess_high_median')+'。Real 大于双侧 controls 的 cell fraction：'+neg('high_fraction_real_gt_both')+'。这些不是 background-corrected formal coupling。'
    ass=[]
    for p in ['source_z','proximal_z','background_rcv','z_top_jitter','X1_jitter','valid_fraction','negative_control_rho']:
        aa=a[a.predictor.eq(p)&a.outcome.eq('real_high_rho')&a.diameter.astype(str).isin([str(d) for d in DS])]
        vals=[float(aa[aa.diameter.astype(str).eq(str(d))].spearman.iloc[0]) for d in DS]
        ass.append(dict(predictor=p,**{f'D{d}':fmt(x) for d,x in zip(DS,vals)}))
    qs['Q13']='七项 QC 对 real high-region rho 的 within-diameter descriptive Spearman 见表，每项 n=5；常量指标的 rho 为 NaN。其余三个 coupling outputs 和 pooled-20 exploratory coefficients 保存在 qc_vs_coupling_summary.csv；pooled 结果可能受 diameter 混杂。'
    qs['Q14']='Mixed evidence。D500 source 与 proximal-tail detectability 在五个 matched common flows 中均低于 D235/D285；geometry valid fraction、jitter、background stability 和 negative controls 需分别按表读取，不能合并成一个 SNR。较低局部可检测性是较弱观测 coupling 的可能贡献因素；本轮描述性数据不能量化其解释比例，也不能断言完全解释或无关。'
    qs['Q15']='本轮不建立 post-hoc exclusion threshold。Z>0/1/2/3 与 jitter robust flag 仅为预先固定的描述性标记，全部 23 卷保留。'
    js('data_answers.json',qs)
    facts=[qs['Q1'],qs['Q3'],qs['Q8'],qs['Q9'],qs['Q10'],qs['Q11'],qs['Q12'],qs['Q7'],
       '20 common-grid volumes / 9,456 valid frames；3 extension volumes / 1,475 valid frames；共 11,500 原始位置中 569 个 geometry 缺失。',
       f'全部 {len(allraw):,} retained raw SV arrays 为 351×500，{int(allraw.n_pixels.sum()):,} pixels；NaN={int(allraw.nan_count.sum())}，Inf={int(allraw.inf_count.sum())}，负值={int(allraw.negative_count.sum())}。',
       '10,931 个正式有效帧逐一 SHA256 通过；左右 matched ROI 在全部有效帧均完整覆盖，无 source/tail 重叠。',
       '23 卷 10,931 个有效帧均有同形状、有限的 retained stru_amp；只计算 structural local contrast/SBR，没有 calibrated noise-floor reference。']
    js('key_data_facts.json',facts)
    errors=pd.DataFrame([dict(metric=k,**vv) for k,vv in val['formal_reintegration'].items()])
    text='# SV-OCTA 数据质量 / 背景可检测性 / 负对照 QC 审计 v1\n\n'
    text+='本目录为独立诊断层。正式 raw SV、continuity-first v2.1 geometry 及两个冻结 coupling 目录保持原样；没有背景相减、重算 formal coupling、自动排除或物理机制解释。独立实验单位为 diameter × flow × scan volume，B-scans 是卷内空间位置。\n\n'
    text+=f'Starting HEAD：`{val["starting_head"]}`；branch：`analysis/sv-diameter-stage2`。已执行 fetch、checkout、pull --ff-only；未 merge main。最终提交 SHA 由 Git commit 与交付回复标识，避免自引用。\n\n'
    text+='## 固定方法与统计约定\n\n'
    text+='SV = var(abs(IMG),1,3)，分母 N，全部正式 mean/Q 使用 linear raw SV。Source 是冻结 X4 中心的 X1 × physical-D ellipse，dx=12.7 μm、dz=6.7 μm、16×16 supersample；tail 从 true physical bottom 起，guard=0，X1 宽，无 cone/spreading。左右中心只移动 ±1.5 X1_px，边到边间距 0.5 X1。连续几何严格平移，仍在固定图像子像素格点积分；分数像素平移会使 rasterized 有效面积略有变化，不能要求左右离散权重逐元素相同。每侧需完整 source 与 0–500 μm tail 在 FOV 内才有效，不挪回、不缩小、不因高 SV 删除。\n\n'
    text+='ROI mean、median、MAD 与 population SD 均采用 fractional weights；weighted median 是累计质量过半处，恰好一半时取相邻值均值。Combined BG 拼接有效侧像素和权重，按总像素权重计算，非左右 Z 的平均。Z=(weighted real median−weighted BG median)/(1.4826 weighted BG MAD)；CNR-like=(mean difference)/weighted BG SD。SBR=real mean/BG mean，SV SBR dB=10log10；contrast=(real−BG)/(real+BG)。Z/CNR 零分母记 NaN，有计数；不加 epsilon。Asymmetry 保留 signed log ratio、absolute log ratio、absolute difference 及 2|L−R|/(L+R)，log epsilon 固定 float64 tiny。Z>0/1/2/3 的 fraction 同时给出可评价帧分母和全部有效 geometry 帧分母。Kish Neff=(sum w)^2/sum(w²)，所有 region/side 都保留。\n\n'
    text+='6 source bins 严格沿用 frozen normalized bin boundaries；20 tail bands 为 25 μm 等宽。Pooled 0–100 μm 直接整块 ROI 重算 median/MAD/Z，不平均四个 Z。背景稳定性分别保存每个 ROI 的 left/right/combined、pixel mean/median 的 slow-axis curve。CV 用 population SD/mean；RCV=1.4826 MAD/median；dynamic range=(p95−p5)/median。Drift 为固定原始 0–50 与 449–499 两块的 median 差，另有相对值、每原始 frame slope、51-frame trend max−min 及 residual MAD，无 p-value。\n\n'
    text+='所有 detrend 在原始 0…499 网格上做 centered 51-frame moving median，端点截短，min_periods=1，忽略邻域 NaN，缺失中心保持缺失。Pseudo mean–mean 6×20 的 Spearman、lag ±50、signed peak/tie、specificity 直接调用冻结数值函数，不运行其正式分析入口。Lag=Δ 为 Source(i) vs Tail(i+Δ)，不循环。Specificity=rho(0)−median(rho at 30≤|lag|≤50)。High region 与上一轮相同：正 rho cells 中最高 ceil(10% positive cells)，截止 ties 全纳入。Pseudo 自身 high region 与在 frozen real high region 上的值分别保存；real/control excess 比较使用同一 real region，避免不同 cells 的混淆。左右 rho 的 cellwise median 形成 control rho，Δrho 仅 diagnostic。\n\n'
    text+='Geometry jitter 只在原始 indices 相差 1 时计算；|Δg| > median+5×1.4826 MAD 仅 flag。每直径五个 common flows 1/3/5/7/10 独立卷；D500 2/9/12 单独 extension，不进入 common 汇总。无 frame-level tests、p-values、mixed models 或复杂预测模型。\n\n'
    text+='## 输入与验证\n\n'+table(errors)+'\n\n'
    text+=f'10,931 valid arrays（D128 NPZ=2,422，MAT=8,509）与上一轮 identity 完全一致；25 release ZIP SHA 全部核对。另审计 569 个 invalid-geometry raw arrays：D128 使用冻结 release NPZ SHA，其他 MAT 保存新 SHA receipt（这些无效帧不在既有 coupling identity 表中，未伪称其 SHA 为旧表验证）。全部 11,500 arrays 可用，共 {int(allraw.n_pixels.sum()):,} pixels，finite=100%，NaN/Inf/negative=0。逐数组 shape、zero fraction、min、median、p95/p99/p99.9/max 保存在 raw audit CSV。\n\n'
    text+=f'背景左右 coverage 均 10,931/10,931；最大 overlap weight={val["max_overlap_weight"]}；连续 translation 参数最大误差={val["translation_parameter_max_error"]:.3g}；full-image ellipse spot-check 权重误差=0。独立 scalar lag spot-check={val["independent_lag_spotcheck_count"]} 条，最大 rho 差={val["independent_lag_max_absolute_error"]:.3g}，原始 index 配对计数通过。完整冻结文件 SHA 检查 {val["n_frozen_files_checked"]} 项通过。\n\n'
    text+='## Common-grid 数字汇总\n\n所有单元格均为五卷值的 median [min,max]；五个原值另保存在 CSV。\n\n'+table(pd.DataFrame(main))+'\n\nProximal-tail robust Z：\n\n'+table(pd.DataFrame(prox))+'\n\nSource-bin robust Z：\n\n'+table(pd.DataFrame(sourcebins))+'\n\nGeometry：\n\n'+table(pd.DataFrame(jrows))+'\n\nNegative controls（51-frame detrended）：\n\n'+table(pd.DataFrame(negs))+'\n\nQC → real high-region rho 的 descriptive Spearman：\n\n'+table(pd.DataFrame(ass))+'\n\n'
    text+='## Q1–Q15：纯数据答案\n\n'+'\n\n'.join(f'**{k}**　{x}' for k,x in qs.items())+'\n\n'
    text+='## Structural OCT 与噪声参考\n\n23 卷 retained arrays 均含共配准、同形状、全有限 stru_amp。该量按代码为 repeat amplitude mean，只报告 structural source-to-local-background contrast、local SBR 和 BG variability；structural ratio 未转换为 dB。没有在 retained cohort 的 manifest、README、signal implementation 或 metadata 中识别到有明确 provenance 的独立 dark/air/noise-floor reference。No calibrated instrument-noise-referenced OCT SNR could be computed from the available retained data。结构 QC 与 SV/coupling 关系仅为 same-cohort volume-level descriptive Spearman。\n\n'
    text+='## 异常、缺失与范围边界\n\n569 个 geometry 缺失在名义 11,500 行中保留；没有 QC 新增 exclusion。无 FOV coverage 缺失，无 raw nonfinite/negative。左右不对称按事实保留，没有按值删除“疑似其他血管”；asymmetry 本身不标定污染来源。零分母计数、constant-input NaN 和 Neff 全保留 CSV。仪器 noise reference 不可用，未计算 true instrument SNR；本轮不能从描述性 n=5 关系确定 measurement detectability 对 coupling 降低的因果解释比例。全部任务验证以 validation.json 为准。\n\n'
    text+='## 文件与复现\n\n逐帧文件按 D128/D235/D285/D500 与 D500_extension 分组 gzip；含名义帧、valid flag、全 source/bin/tail 的 real/L/R/BG mean/median/MAD/SD/Q/area/Neff、contrast、阈值、structural 指标。volume_region_detectability_summary 为所有 ROI 的汇总；用户指定的 source/tail/proximal/geometry/background/negative-control/relationship CSV 均单独输出。完整 lag 曲线按 scan×side×mode gzip；slow-axis raw/trend/residual 另有 gzip。Figure 1–10 同时提供 PNG/PDF，mean/median 统计不将 B-scans 当误差重复；所有 flow 全展示，coupling 色轴 [-1,1]，excess [-2,2]，同 metric 共享含零尺度。\n\n'
    text+='在仓库根目录执行（需 input_manifest 内相同 SHA retained arrays；文件路径可按已有 identities 解析）：\n\n```powershell\npython analysis/sv_data_quality_qc_v1/derive_qc.py\npython analysis/sv_data_quality_qc_v1/audit_additional_arrays.py\npython analysis/sv_data_quality_qc_v1/analyze_qc.py\npython analysis/sv_data_quality_qc_v1/make_figures.py\npython analysis/sv_data_quality_qc_v1/verify_qc.py\npython analysis/sv_data_quality_qc_v1/build_report.py\n```\n\n'
    text+='## 最重要的 12 条数据事实\n\n'+'\n\n'.join(f'{i}. {s}' for i,s in enumerate(facts,1))+'\n'
    supplemental=[]
    for d in DS:
        g=v[v.diameter_um.eq(d)];b=r[r.diameter_um.eq(d)&r.region.eq('source')]
        supplemental.append({'D':d,'source BG pixel median':mm(b.bg_median_median),'source BG pixel MAD':mm(b.bg_mad_median),
           'BG absolute log asymmetry':mm(b.bg_asym_log_median),'structural source local SBR':mm(g.structural_source_sbr)})
    text+='\n## 背景不对称与结构对照补充\n\n'+table(pd.DataFrame(supplemental))+'\n\n'
    text+='左右背景不对称完整保留，包括每卷 median/p90/p95；上述数值没有定位其来源，也没有据此认定或移除其他血管。D500 structural source/local-background SBR 低于其他三直径，且与更低 SV detectability 同时出现；这些是 depth-matched local contrast，不是噪声参考的仪器 SNR。noise-reference 的 metadata 检索范围是预先固定每卷最接近 index 249 的有效帧，共 23 个记录，加上 cohort manifests、README 和 source implementation；完整检索记录另存 CSV，不声称穷尽每个原始 acquisition metadata。\n\n'
    text+='Common-grid 的 X1 median absolute jitter 全为 2 pixels，z_top median absolute jitter 全为 0；两项作为候选 predictor 在每直径五卷内恒定，相关系数不可定义（NaN），不是漏算。X1/z_top 的 p90/p95/p99/max/MAD 和 flag fraction 仍完整保留；z_top p95 为 1 pixel。零 MAD 时的 jitter flag 按原先固定公式执行，没有额外放宽。\n'
    text+='\n最后执行 `python analysis/sv_data_quality_qc_v1/finalize_delivery.py` 核对交付范围及输出身份，生成 output_sha256.csv（排除其自身）；intermediate checkpoint 和缓存不提交。所有文本固定 UTF-8/LF，确保 Git checkout 后 SHA 一致。\n'
    (OUT/'README.md').write_text(text,encoding='utf-8',newline='\n')
    p=json.loads((OUT/'provenance.json').read_text());p['delivered_script_sha256']={x.name:digest(x) for x in OUT.glob('*.py')};p['input_manifest_sha256']=digest(OUT/'input_manifest.csv')
    p['read_policy_note']='CSV parser low_memory=False fixes mixed boolean inference; no numerical definition changed';p['calibrated_noise_reference']='none documented; unavailable, no true instrument SNR'
    js('provenance.json',p)
    print(table(pd.DataFrame(main)));print(table(pd.DataFrame(prox)));print(table(pd.DataFrame(negs)))

if __name__=='__main__':main()
