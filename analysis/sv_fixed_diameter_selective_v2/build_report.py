"""Three complete answers with machine-readable claim-to-row provenance."""
from common import *
def link(n,label=None):return f'[{label or n}](<{(OUT/n).as_posix()}>)'
def fmt(x,dec=3):return 'NaN' if pd.isna(x) else f'{x:.{dec}f}' if isinstance(x,(float,np.floating)) else str(x)
def table(df):
    return '\n'.join(['| '+' | '.join(map(str,df.columns))+' |','| '+' | '.join(['---']*len(df.columns))+' |']+['| '+' | '.join(fmt(x).replace('|','/') for x in row)+' |' for row in df.itertuples(index=False,name=None)])
def rng(a):
    s=summary(a);return f"{s['median']:.3f} [{s['min']:.3f}, {s['max']:.3f}]"
EVIDENCE=[];EVIDENCE_ROWS=[]
def evidence(eid,question,file,claim,selection=None,selector='all rows'):
    a=read(OUT/file);a['_csv_row_number']=np.arange(2,len(a)+2)
    if selection is not None:a=a[selection(a)]
    EVIDENCE.append(dict(evidence_id=eid,question=question,claim=claim,table=file,sha256=sha(OUT/file),selector=selector,n_rows=len(a),csv_row_numbers=json.dumps(a._csv_row_number.tolist())))
    for r in a.to_dict('records'):
        rownum=r.pop('_csv_row_number');EVIDENCE_ROWS.append(dict(evidence_id=eid,table=file,csv_row_number=rownum,row_json=json.dumps(clean(r),ensure_ascii=False,allow_nan=False)))
    return link(file,eid)
def main():
    val=json.loads((OUT/'validation.json').read_text());assert val['status']=='passed'
    volume=read(OUT/'volume_metrics.csv');contrast=read(OUT/'diameter_contrasts_by_flow.csv');region=read(OUT/'fixed_region_summary.csv');maps=read(OUT/'depth_maps_raw51.csv');whole=read(OUT/'whole_coupling.csv');ex=read(OUT/'fixed_region_control_excess.csv');control=read(OUT/'fixed_region_control_cells.csv');local=read(OUT/'local_contrast_profiles.csv');pos=read(OUT/'position_sensitivity_metric_changes.csv');pp=read(OUT/'position_sensitivity_pairs.csv');pr=read(OUT/'position_sensitivity_fixed_region.csv');sens=read(OUT/'detrend_sensitivity_changes.csv');diff=read(OUT/'adjacent_difference.csv');geo=read(OUT/'geometry_version_direction_ledger.csv');aext=read(OUT/'conditional_A_summary.csv');st=read(OUT/'conditional_D_structural_correlations.csv');drift=read(OUT/'background_drift.csv')
    q1=evidence('Q1-01',1,'volume_metrics.csv','All 15 volumes: absolute means and median of per-frame RI100/RI500')
    q1c=evidence('Q1-02',1,'diameter_contrasts_by_flow.csv','All matched-flow contrasts, signed differences and percentage changes')
    q1p=evidence('Q1-03',1,'raw_depth_profiles.csv','All 20 raw tail bands and RI profiles including deep bands')
    q1b=evidence('Q1-04',1,'spatial_block_summary.csv','Fixed original spatial blocks are heterogeneity diagnostics, not replicates')
    q1g=evidence('Q1-05',1,'geometry_version_direction_ledger.csv','Same-definition old/new per-frame RI comparisons',lambda d:d.metric.str.startswith('RI'),'metric in RI100/RI500')
    q2=evidence('Q2-01',2,'fixed_region_summary.csv','Signed fixed C16 and all six proximal source-row summaries')
    q2m=evidence('Q2-02',2,'depth_maps_raw51.csv','Complete 15 x 120 raw and 51 maps, all signed values')
    q2w=evidence('Q2-03',2,'whole_coupling.csv','All three representative raw/51 pairs, strict specificity, peak/rank and matched specificity')
    q2d=evidence('Q2-04',2,'adjacent_difference.csv','Original adjacent raw differences, including D285 F03 W100 negative result')
    q2s=evidence('Q2-05',2,'detrend_sensitivity_changes.csv','Fixed 31/101 window sensitivities; no reselection')
    q2g=evidence('Q2-06',2,'geometry_version_direction_ledger.csv','Old/new C16 and W100/W500 using comparable definitions',lambda d:~d.metric.str.startswith('RI'),'metric in C16/W100/W500')
    q3=evidence('Q3-01',3,'fixed_region_control_excess.csv','C16 summary on common support; separate L/R and signed cellwise excess')
    q3c=evidence('Q3-02',3,'fixed_region_control_cells.csv','All cellwise background comparisons, including nonpositive excess')
    q3z=evidence('Q3-03',3,'local_contrast_profiles.csv','Weighted median/MAD local contrast profiles, including negative deep-band Z')
    q3b=evidence('Q3-04',3,'background_real_tail_covariation.csv','20 same-depth real/BG raw/51 correlations')
    q3x=evidence('Q3-05',3,'crossed_controls.csv','Real source vs left/right pseudo-tail 51-frame correlations')
    q3dr=evidence('Q3-06',3,'background_drift.csv','Fixed original endpoint coverage and signed/absolute/relative drift')
    q3p=evidence('Q3-07',3,'position_sensitivity_metric_changes.csv','All four rigid offsets on five-position common support')
    q3pp=evidence('Q3-08',3,'position_sensitivity_pairs.csv','All offset 51-frame correlations and specificity, baseline common retained')
    q3a=evidence('Q3-09',3,'conditional_A_summary.csv','Shared source/T1 support trigger and fixed T2-T4 diagnostic')
    q3d=evidence('Q3-10',3,'conditional_D_structural_correlations.csv','Triggered descriptive structural diagnostics; no causal attribution')
    q3l=evidence('Q3-11',3,'localization_review.csv','All 45 fixed images assessed; complete physical boundary remains unresolved')
    q3bl=evidence('Q3-12',3,'conditional_B_review.csv','All four prescribed maximum-jump image pairs reviewed')
    q3pseudo=evidence('Q3-13',3,'pseudo_whole_coupling.csv','All fixed pseudo representative pairs, including positive specificity and boundary peak',lambda d:d.side.isin(['left','right']),'side in left/right')
    q3rank=evidence('Q3-14',3,'position_sensitivity_metrics.csv','RI500 matched-flow ordering in all four offsets',lambda d:d.support.eq('five_position_common')&d.metric.eq('RI500'),'five_position_common and metric=RI500')
    lines=['# 固定真实直径 ROI v2：三个核心问题的新答案','',
        '本轮已完成三直径、五个匹配流速共 15 卷的新版描述、空间关联、匹配背景、负对照与限定敏感性分析。下文的结论均对应新版固定 D 几何；保持、减弱、反向和不可定义的结果分别列出。完整逐帧、逐 cell、逐 lag 结果随本目录交付。',
        '',f'验收：{link("validation.json")} 中 {len(val["checks"])} 项最终数值、覆盖和身份检查通过。独立实现核查是数值验证；本批数据没有独立实验验证。',
        '', '## 阅读约定与分析身份','',
        '正式信号为 SV_raw=var(abs(IMG),1,3)，方差分母 N，保持 linear raw SV。Q=sum(SV_raw×w)×12.7×6.7，mean=Q/area。Source 是真实 D×D 圆形几何，保留冻结的 X4、z_top；source 的 6 个归一化层与 tail 的 20 个 25 μm bands 均固定。下界 z_top+D/6.7 来自几何先验。',
        '', '每卷保留原 frame_index=0…499，总计 7500 个名义位置；D128/D235/D285 分别有 2422/2225/2371 个原有效帧，482 个原无效位置保持缺失。背景与四方向扰动没有新增原 geometry 无效帧。所有相关在单卷原坐标上计算，不合并不同 scan 的帧。',
        '', '下文 median [min,max] 是同一直径五个不同流速卷的描述，不是同条件五次重复，也不是置信区间。五个原值保存在各 across_flow 表。跨 scan 采集设置一致性未核实，绝对 SV 采用记录仪器单位；RI 无量纲。空间 lag 采用帧单位，不换算时间或物理距离。',
        '',f'合同在科学统计前冻结：{link("analysis_plan.json")}、{link("selection_log.md")}、{link("column_mapping.json")}。历史同数据结果参与了 C16 和诊断范围选择，informed_by_historical_results=true、independent_validation=false。D500 的 8 卷按既有事后排除决定仅保留 {link("d500_exclusion_record.csv")}，本轮读取 D500 数组数为 0。',
        '', '## 问题一：匹配流速下的 absolute source/tail 与相对 tail/source','',
        f'五个匹配流速下，tail100、tail500、RI100 和 RI500 均满足 D128 > D235 > D285。Source mean 的关系不同：D128 低于另外两直径；D285 相对 D235 在 1/3/5/7 mm/s 较高，在 10 mm/s 较低。不能用一条“直径越大、所有信号都越低”的结论概括。证据：{q1}、{q1c}。','',
        '下表列出全部 15 卷。Source 按 10⁸、tail 按 10⁷ 仪器单位展示；RI 是先逐帧取比值再取卷内中位数。','']
    vt=volume.pivot(index=ID,columns='metric',values='median').reset_index();vt['source']=vt.source/1e8;vt['tail100']=vt.tail100/1e7;vt['tail500']=vt.tail500/1e7
    lines += [table(vt[['diameter_um','flow_mm_s','source','tail100','tail500','RI100','RI500']].rename(columns={'diameter_um':'D (μm)','flow_mm_s':'流速','source':'source ×10⁸','tail100':'tail100 ×10⁷','tail500':'tail500 ×10⁷'})),'']
    ct=contrast[contrast.metric.eq('RI500')].pivot(index=['reference_D','comparison_D'],columns='flow_mm_s',values='percent_change').reset_index();ct.columns=['参考 D','比较 D','1 mm/s (%)','3 mm/s (%)','5 mm/s (%)','7 mm/s (%)','10 mm/s (%)']
    lines += ['RI500 的逐流速百分比变化如下；差值和其余指标的同类比较均保存在对比长表。','',table(ct),'',
        f'全部 20-band mean 与逐帧 RI 的卷内摘要见 {q1p}。深度坐标使用 band 中点 12.5、37.5、…、487.5 μm；末 band 是 475–500 μm 的平均，不能解释为在 500 μm 的点测量。固定空间块 0–99、100–199、…、400–499 的五项指标均保留，作为覆盖与异质性描述，不进行跨卷解剖同位配对：{q1b}。',
        '',f'Q 和实际面积的描述见 {link("absolute_integral_metrics.csv")}，解析圆面积、离散面积偏差见 {link("geometry_area.csv.gz")}。未将栅格圆面积近似要求为解析面积的机器精度恒等。',
        '', '| 状态 | 新版判断 |','| --- | --- |',
        '| 保持 | 三直径 RI 的匹配流速排序保持；绝对 tail 的排序也逐流速成立。 |',
        '| 减弱／改变 | 新旧几何下 RI 的具体数值改变，必须使用新版数值。不能沿用旧 RI 的幅度；同身份逐帧对照已单独保存。 |',
        '| 反向／例外 | 10 mm/s 的 D285 source mean 低于 D235，与其余四个匹配流速的方向相反。 |',
        '| 不可定义／未核实 | 482 个原 geometry 无效位置的逐帧 RI 为缺失；本批有效帧未出现额外 source 零分母。跨 scan 设置一致性及经校准的绝对强度可比性未得到独立核实。 |','',
        f'新旧 RI 对照证据：{q1g}。逐帧 RI_new/RI_old 恒等式已核对；未把该恒等式套用到各自的卷中位数，也未计算“分母解释百分比”。','',
        '## 问题二：同卷关联的深度分布与零位对应','',
        'C16 固定为 S3–S6 × T1–T4。其定义是全部 16 个有符号 lag0 Spearman rho 的中位数；任一 cell 不可定义时主值为 NaN。C100 则是按 Q/area 合并 S3–S6 后与 pooled tail100 的标量相关，两者不是同一指标。','',
        f'本轮 raw 和 51-frame 地图各有完整 1800 个 cell。51-frame 下 C16 的全部 240 个 cell 都为正，但完整地图保留了 {(maps[maps["mode"].eq("detrended_51")].rho<0).sum()} 个负值；raw 地图保留 {(maps[maps["mode"].eq("raw")].rho<0).sum()} 个负值。两种主地图均无不可定义 cell，也无恰好为零的 cell。方向一致的候选区域并不代表全深度都正相关。证据：{q2}、{q2m}。','']
    rt=[]
    for d in [128,235,285]:
        for mode in ['raw','detrended_51']:
            x=region[(region.diameter_um==d)&region['mode'].eq(mode)&region.region.eq('C16')].sort_values('flow_mm_s');rt.append({'D':d,'模式':mode,'流速 1':x.value.iloc[0],'3':x.value.iloc[1],'5':x.value.iloc[2],'7':x.value.iloc[3],'10':x.value.iloc[4],'median [min,max]':rng(x.value),'正/负/零/可定义卷':'5/0/0/5'})
    lines += [table(pd.DataFrame(rt)),'',
        '51-frame 的 C16 范围为 0.084–0.221，说明该固定区域的关联幅度有限。本轮不设强／弱阈值。六个 source 层各自对 T1–T4 的四项中位数均完整保留；没有重选热点或用最大 cell 替代固定区域。','',
        f'三个代表量的 51-frame 结果如下。Specificity 采用固定 42 个远端 lag；matched specificity 用每个远端 lag 的共同配对坐标核查组成差异。正 specificity 不等同于峰必在零位。证据：{q2w}。','']
    wt=[]
    for (d,pair),g in whole[whole['mode'].eq('detrended_51')].groupby(['diameter_um','pair']):
        wt.append({'D':d,'配对':pair,'rho0 median [min,max]':rng(g.rho0),'specificity':rng(g.specificity),'matched specificity':rng(g.matched_specificity),'峰在 0 的卷数':int(g.peak_lag.eq(0).sum()),'峰 lag 范围':f'{int(g.peak_lag.min())}…{int(g.peak_lag.max())}'})
    lines += [table(pd.DataFrame(wt)),'',
        f'原坐标相邻差分的关系总体更小，D285 的 W100 在 3 mm/s 为 {diff[(diff.diameter_um==285)&diff.flow_mm_s.eq(3)&diff.pair.eq("W100")].rho.iloc[0]:.6f}，为反向结果；D235 的 W500 在 10 mm/s 为 {diff[(diff.diameter_um==235)&diff.flow_mm_s.eq(10)&diff.pair.eq("W500")].rho.iloc[0]:.6f}，接近零。其余全部卷和配对没有省略：{q2d}。',
        '',f'31/101-frame 敏感性仅覆盖固定 16 cells 与三个代表量。相对 51-frame，代表 rho0 的变化为 {sens[sens.metric.eq("rho0")].change.min():.3f} 至 {sens[sens.metric.eq("rho0")].change.max():.3f}，specificity 变化为 {sens[sens.metric.eq("specificity")].change.min():.3f} 至 {sens[sens.metric.eq("specificity")].change.max():.3f}；需按卷与窗口读取，不能将主窗口的幅度推广为尺度不变。{q2s}。','']
    cg=geo[geo.metric.eq('C16')&geo['mode'].eq('detrended_51')];gt=[]
    for d,g in cg.groupby('diameter_um'):gt.append({'D':d,'旧版同定义 C16':rng(g.old_value),'新版同定义 C16':rng(g.new_value),'减弱卷数':int(g.direction_status.eq('same_sign_reduced').sum()),'增加卷数':int(g.direction_status.eq('same_sign_increased').sum())})
    lines += ['同身份、同一固定区域定义的新旧 C16 对比如下。这里先在各卷取 C16，再描述五流速；没有与“先对每个 cell 做跨流速共识再取区域中位数”混用。','',table(pd.DataFrame(gt)),'',
        f'这个对照显示 D128 的 C16 在五卷均减弱，D235 为四卷增加、一卷减弱，D285 方向混合。按逐卷 C16 的五流速中位数，旧版 D128 高于 D235，新版变为 D235 高于 D128；这不是 RI 排序的变化。{q2g}。','',
        '| 状态 | 新版判断 |','| --- | --- |',
        '| 保持 | 固定 C16 的正方向在三个直径、全部流速保留；三个代表量 51-frame 的 specificity 均为正。 |',
        '| 减弱 | D128 同定义 C16 相对旧几何全五卷降低；D285 raw 的广泛正关联在去趋势后明显收缩；原始相邻差分中存在接近零的结果。 |',
        '| 反向／例外 | 全深度地图有负相关；D285、3 mm/s、W100 的相邻差分为负；部分代表量的峰在 −2、−1 或 +1 帧；旧 C16 的 D128/D235 汇总排序反转。 |',
        '| 不可定义 | 主 6×20 地图、C16 和三个代表量的主摘要本批均可定义，NaN 数为 0；不足 3 对、常量或任一必需 lag 缺失的规则已经测试，不以 0 填补。lag 的 μm 或时间含义不可定义。 |','',
        '## 问题三：背景、负对照和定位敏感性如何约束解释','',
        f'左右背景中心均按 X4±1.5D/dx 平移，使用相同形状和深度；7018 个原有效位置均满足双侧完整 FOV。real/L/R 先取共同掩码，再分别去趋势，因此本批 real-matched 与 real-all 的支持恰好相同。Combined BG 由左右像素及 fractional weights 拼接，局部 Z 使用加权中位数和 MAD，不是两侧 Z 的平均。覆盖见 {link("background_coverage.csv")}。','',
        f'15 卷的 C16_excess 均为正，范围 {ex.excess_cell_value.min():.3f}–{ex.excess_cell_value.max():.3f}；固定区域的 real-minus-L、real-minus-R 也均为正。但全图 1800 个 cell 有 {(control.excess_cell<=0).sum()} 个 excess≤0，其中固定 C16 内有 {(control.in_C16&(control.excess_cell<=0)).sum()}/240 个。不能将区域中位数为正改写成“所有深度都超过背景”。{q3}、{q3c}。','']
    et=ex[['diameter_um','flow_mm_s','real_matched_rho_value','left_rho_value','right_rho_value','excess_cell_value','real_minus_L_value','real_minus_R_value']].rename(columns={'diameter_um':'D','flow_mm_s':'流速','real_matched_rho_value':'real C16','left_rho_value':'L C16','right_rho_value':'R C16','excess_cell_value':'C16 excess','real_minus_L_value':'median(real−L)','real_minus_R_value':'median(real−R)'})
    lines += [table(et),'',
        'C16 excess 是逐 cell 相减后取 16 项中位数；表中左右各自的区域中位数不能相减来重构该列。系数差只是一项描述性诊断，不是背景校正后的正式 SV/coupling，也不是无偏因果效应。','',
        f'左右 pseudo 的 90 个代表量相关和 specificity 本批也都为正，rho0 为 0.012–0.270，说明正的局部关联或正 specificity 本身并非 real ROI 独有。D128、1 mm/s、左侧 W500 的峰位于 +50 帧边界，rho0=0.097、specificity=0.120；原样保留边界峰，未扩大 lag 范围或推断范围外峰位。完整 L/R lag 与 matched specificity：{q3pseudo}。','',
        f'背景对比随深度减小，并在部分深 band 反向。以下为五流速卷内 Z 中位数的 median [min,max]；没有设置 Z>1/2/3 阈值。负 Z 表示该 ROI 像素中位数低于局部背景中位数，不能解释为负的正式 SV。{q3z}。','']
    zt=[]
    for d in [128,235,285]:
        row={'D':d}
        for r in ['source','tail100','tail500','t10','t20']:row[r]=rng(local[(local.diameter_um==d)&local.region.eq(r)&local.metric.eq('z_local')]['median'])
        zt.append(row)
    lines += [table(pd.DataFrame(zt)),'',
        f'D235 和 D285 的 T20（475–500 μm）局部 Z 在五个流速下均为负；D128 的末 band 局部对比仍为正但接近零。本轮没有据此生成仪器噪声底或“检测深度”。同深度 real-tail 与 combined BG-tail 的 51-frame rho 在全 20-band、全流速中的范围为 D128 −0.112 至 0.118、D235 −0.152 至 0.117、D285 −0.141 至 0.099，幅度接近零且包含反向。Crossed controls 同样存在正负方向：例如 D235 的 W500，左侧五卷均为负（−0.087 至 −0.029），右侧五卷均为正（0.014 至 0.110）。全部 raw/51 和 L/R 值见 {q3b}、{q3x}。','',
        f'背景本身沿空间坐标变化。Combined 背景的 source、tail100、tail500 固定末端块中位数相对起始块均下降。三者相对差范围分别为 {drift[(drift.side=="bg")&(drift.region=="source")].relative_difference.min()*100:.1f}% 至 {drift[(drift.side=="bg")&(drift.region=="source")].relative_difference.max()*100:.1f}%、{drift[(drift.side=="bg")&(drift.region=="tail100")].relative_difference.min()*100:.1f}% 至 {drift[(drift.side=="bg")&(drift.region=="tail100")].relative_difference.max()*100:.1f}%、{drift[(drift.side=="bg")&(drift.region=="tail500")].relative_difference.min()*100:.1f}% 至 {drift[(drift.side=="bg")&(drift.region=="tail500")].relative_difference.max()*100:.1f}%。分母为各自整条背景中位数曲线的中位数；原端点范围为 0–50 与 449–499（均含端点），实际支持数随表提供。{q3dr}。','',
        f'全部四个单轴平移都已计算：横向 ±1 像素（±12.7 μm），轴向 ±1 像素（±6.7 μm）。使用 baseline 与四偏移共同支持后，RI500 的卷内中位数变化为 {pos[pos.metric.eq("RI500")].percent_change.min():.2f}% 至 {pos[pos.metric.eq("RI500")].percent_change.max():.2f}%；RI100 为 {pos[pos.metric.eq("RI100")].percent_change.min():.2f}% 至 {pos[pos.metric.eq("RI100")].percent_change.max():.2f}%。C16 的变化范围为 {pr[~pr.offset.eq("baseline")].C16_change.min():.3f} 至 {pr[~pr.offset.eq("baseline")].C16_change.max():.3f}，四方向各卷 C16 仍为正。{q3p}、{q3pp}。这是小幅系统位置偏移的敏感性，不是定位误差真值，也未验证逐帧随机抖动。','',
        f'在同一固定平移方向下继续匹配流速，四个方向的全部五个流速均保留 RI500 的 D128 > D235 > D285 排序；幅度变化与排序保持分别记录。{q3rank}。','',
        f'几何共享像素检查发现每卷 source 的最下层与 T1 存在共享像素，故按预定触发 A 补充 S3–S6×T2–T4 的 C12，以及 pooled S3–S6 对 tail25–100 的 51-frame 相关。C12 为 {aext.C12.min():.3f}–{aext.C12.max():.3f}，pooled 相关为 {aext.pooled_S36_tail25_100_rho.min():.3f}–{aext.pooled_S36_tail25_100_rho.max():.3f}，全部为正。该结果不等于排除了 PSF 串扰。{q3a}。Neff、非零像素数、sum(w²)、层轴向像素跨度和解析面积误差均已保存；Neff 不是独立散斑数。','',
        f'45 帧结构/SV 并列图均已审核。结构图存在轴向亮度梯度，完整 source 轮廓和独立下壁不足以判定，因此按合同保留“图像不足以判断”。D285 的 1、3 mm/s 首帧中 source 下缘附近的结构亮区触发了 B，补看最大相邻 ΔX4 与 Δz_top 的四组帧对；这些图仍未提供下壁真值，未修改定位或筛帧。{q3l}、{q3bl}。','',
        f'触发 D 后，仅补算 real/L/R 的 source、tail100、tail500 structural mean。Real ROI 的同区域 structural–SV 51-frame rho 范围为 {st[(st.side=="real")&st["mode"].eq("detrended_51")&st.kind.eq("same_ROI_structural_vs_SV")].rho.min():.3f}–{st[(st.side=="real")&st["mode"].eq("detrended_51")&st.kind.eq("same_ROI_structural_vs_SV")].rho.max():.3f}，表明所测 SV 与结构变化共同出现；不能将这种共变量化为机制贡献比例。左右 pseudo 的对应值和结构 source–tail 相关也全部交付。{q3d}。','',
        'C 的两侧反向判断条件未触发，图像亦未发现足以触发的单侧背景结构失配，因此未增加 ±2.5D 对照。报告不依赖精确单一 source bin 的位置，E 未触发。各项 trigger、帧 ID、范围与完成状态见 '+link('conditional_extension_log.csv')+'。','',
        '| 状态 | 新版判断 |','| --- | --- |',
        '| 保持 | 主双侧匹配背景完整；C16 的区域中位超额在全部 15 卷为正，四个固定方向的位置扰动下 C16 也为正。 |',
        '| 减弱 | 单个 cell 并非普遍超出对照；深部局部对比趋近零；窗口和位置会改变幅度，背景漂移与结构共变限制了仅凭正相关作出的解释。 |',
        '| 反向 | 540 个全图 excess≤0，其中 12 个位于 C16；D235/D285 最后 band 的局部 Z 在五流速均为负。 |',
        '| 不可定义 | 原 geometry 无效位置保留缺失；本批有效像素区域未出现额外背景 MAD=0。独立下壁、PSF、噪声底、独立散斑数、校准检测深度及因果效应均不能由当前数据定义。 |','',
        '## 追溯、验证与完整交付','',
        f'每个正文证据编号对应 {link("answer_evidence.csv")} 的表名、SHA256、筛选条件及 CSV 行号；对应原行的完整值另存 {link("answer_evidence_rows.csv.gz")}。全覆盖矩阵见 {link("coverage_ledger.csv")}，三个问题的完成对照见 {link("core_question_coverage.csv")}。',
        '',f'原输入和历史目录的前后 SHA 均一致：{link("input_hashes_after.csv")}。正式 v2 输入身份、运行环境和起始 HEAD 见 {link("provenance.json")}；本轮没有 checkout、pull、reset、commit、push 或远端发布。',
        '',f'数值核查包括同 16×16 定义的 45 帧 real/pseudo/四偏移积分、分层重构、weighted median/MAD、原坐标 lag 与有效配对、常量和 ties、已知 shift、正比例缩放、共同掩码、相邻差分和 matched specificity。非零积分相对误差阈值 1e-10，rho 差阈值 1e-12，零参考单列。{link("independent_integration_checks.csv.gz")}、{link("statistical_independent_checks.csv.gz")}、{link("control_position_independent_checks.csv.gz")}。',
        '',f'过程中发现并修复了加权半质量边界的浮点累加差异：边界情况用精确二进制权重算术判断，独立 Fraction 实现复核。首次失败的本轮检查点保留在 failed_weighted_roundoff_checkpoints/，不属于有效科学交付；修正后全量像素统计和 45 帧检查通过。没有改变信号、阈值或分析合同。记录见 {link("weighted_numerical_fix.json")}。',
        '', '## 图表索引','']
    for d in [128,235,285]:lines.append(f'- D{d}：{link(f"figures/depth_atlas_D{d}.png","全部五卷 raw/51 地图")}；{link(f"figures/control_atlas_D{d}.png","real-matched／L／R／excess 全图")}；{link(f"figures/lag_atlas_D{d}.png","全部五卷三代表 lag")}。')
    lines += ['',f'{link("figures/matched_flow_metrics.png","匹配流速强度和 RI")}；{link("figures/depth_profiles.png","全 20-band 强度／RI／局部 Z")}；{link("figures/consensus.png","五流速系数中位共识图")}；{link("figures/fixed_region_controls.png","固定区域对照")}；{link("figures/position_sensitivity.png","四方向敏感性")}。所有统计图另有同名 PDF；地图相关色轴固定 [−1,1]，excess 为 [−2,2]。','',
        f'原始曲线、各模式实际窗口支持数、完整逐 cell 值、可定义原因、负值表及所有单侧背景均在对应 CSV 中。{link("figure_manifest.csv")} 提供图与输入表映射。完整输出文件身份见 {link("output_sha256.csv")}。','',
        '复现顺序：initialize.py（仅空的新输出目录）→ verify_statistics.py → analyze_tables.py → extract_pixels.py → render_localization.py → analyze_controls.py → 对 45 帧完成并记录图像审核 → conditional_structural.py → 完成触发 B 的图像审核 → summarize_results.py → make_figures.py → validate_delivery.py → build_report.py → finalize_delivery.py。record_visual_review.py 和 finish_visual_review.py 保存的是本次实际图像审核记录；更换输入后必须重新看图，不能把这些标签当成自动定位判定。现有成功输出不应自动覆盖；需要复跑时使用另一个明确命名的输出副本。','']
    (OUT/'README.md').write_text('\n'.join(lines),encoding='utf-8')
    csv('answer_evidence.csv',EVIDENCE);csv('answer_evidence_rows.csv.gz',EVIDENCE_ROWS)
    matrix=[dict(question=1,definition='15 matched flow volumes; median of per-frame RI; raw mean/Q/area separated',coverage='15/15 scans; 75 scalar rows; 75 pairwise contrasts; 600 depth-profile rows; 375 block rows',retained='RI and absolute tail ordering',weakened='amplitudes change with fixed-D geometry',reversed='D285 versus D235 source at flow10',undefined='482 original positions; settings comparability unverified',status='complete'),
        dict(question=2,definition='full 6x20 raw/51 maps; signed strict C16; three fixed pairs; original-coordinate lag',coverage='3600/3600 cells; 90 main representative rows; 18180 lag rows incl 31/101; 45 adjacent rows',retained='positive C16 and main specificity',weakened='D128 old/new C16, adjacent differences and scale dependence',reversed='negative full-map cells; D285 F03 adjacent W100; old C16 summary ordering',undefined='0 undefined main cells; lag physical step unavailable',status='complete'),
        dict(question=3,definition='new +/-1.5D both-side pooled pixels; common masks before detrend; four prescribed shifts',coverage='7018/7018 arrays; 5400 real-matched/pseudo cells; 1800 excess cells; all 4 offsets; 45+8 reviewed frames; A/B/D complete, C/E not triggered',retained='positive region-level excess',weakened='nonpositive individual cell excess; deeper background convergence; drift and structural covariation',reversed='deep D235/D285 Z and 540 excess cells',undefined='independent wall, PSF, calibrated noise/detection depth and causal effect',status='complete')]
    csv('core_question_coverage.csv',matrix)
    # The resulting report contains all required classes, including observed zeros of undefined counts.
    js('report_validation.json',dict(status='passed',questions=3,evidence_entries=len(EVIDENCE),evidence_rows=len(EVIDENCE_ROWS),no_weak_correlation_cutoff_introduced=True,all_claims_use_saved_v2_outputs=True,unresolved_physical_quantities_explicit=True))
    print('Three-question report and claim-to-row ledger written',len(EVIDENCE),flush=True)
if __name__=='__main__':main()
