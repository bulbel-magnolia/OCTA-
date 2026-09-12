"""Endpoint-specific intersections, independent medians and scientific figures."""
from run_omag import *
import math,statistics
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
M=s.METRICS;D=s.DIMS;F=s.FLOWS
def eligible(t,m):
 q=t.senior_metric_assessable.astype(bool)&t.auc_defined.astype(bool)
 if 'senior_full_assessable' in t:q &= t.senior_full_assessable.astype(bool)
 if m!='raw_AUC_200':q &= t.denominator_defined.astype(bool)
 if m=='normalized_AUC_um':q &= t.normalized_defined.astype(bool)
 return q
def loadstage(p):return pd.concat([read(f) for f in sorted(p.glob('framewise_metrics_*.csv.gz'))],ignore_index=True)
def ranking(a):return ' < '.join('D'+str(D[i]) for i in np.argsort(a)) if len(set(a))==3 and np.isfinite(a).all() else 'ties/undefined: '+s.direction(a)
def median(values):return statistics.median(list(values)) if len(values) else np.nan
def comparisons(tables):
 summary=[];checks=[];ledger=[];pct=[];oatr=[]
 for pair,left,right in [('SV-A_vs_OMAG-A','SV-A','OMAG-A'),('SV-B_vs_OMAG-B','SV-B','OMAG-B'),('OMAG-A_vs_OMAG-B','OMAG-A','OMAG-B')]:
  rows=[]
  for sid,a in tables[left].groupby('scan_id'):
   a=a.set_index('frame_index').sort_index();b=tables[right].query('scan_id==@sid').set_index('frame_index').sort_index();assert a.index.equals(b.index)
   for m in M:
    ma=eligible(a,m);mb=eligible(b,m);common=ma&mb
    # Independently reconstruct set intersection and medians using keyed records.
    def keys(g):
     out={}
     for fi,r in g.iterrows():
      valid=bool(r.senior_metric_assessable and r.auc_defined and r.get('senior_full_assessable',True))
      if m!='raw_AUC_200':valid=valid and bool(r.denominator_defined)
      if m=='normalized_AUC_um':valid=valid and bool(r.normalized_defined)
      if valid:out[int(fi)]=float(r[m])
     return out
    ka=keys(a);kb=keys(b);kk=set(ka)&set(kb);assert kk==set(a.index[common])
    ac=float(a.loc[common,m].median());bc=float(b.loc[common,m].median());ia=median([ka[k] for k in kk]);ib=median([kb[k] for k in kk]);passed=np.allclose([ac,bc],[ia,ib],rtol=1e-14,atol=0,equal_nan=True);assert passed
    ch=dict(pair=pair,scan_id=sid,metric=m,n_common=len(kk),frame_support_sha256=hashlib.sha256(','.join(map(str,sorted(kk))).encode()).hexdigest(),left_independent_median=ia,right_independent_median=ib,passed=bool(passed));checks.append(ch)
    rows.append(dict(pair=pair,left_algorithm=left,right_algorithm=right,scan_id=sid,diameter_um=int(a.diameter_um.iloc[0]),flow_mm_s=int(a.flow_mm_s.iloc[0]),metric=m,n_left_native=int(ma.sum()),n_right_native=int(mb.sum()),n_common=int(common.sum()),n_left_only=int((ma&~mb).sum()),n_right_only=int((mb&~ma).sum()),left_native_median=float(a.loc[ma,m].median()),right_native_median=float(b.loc[mb,m].median()),left_common_median=ac,right_common_median=bc,frame_support_sha256=ch['frame_support_sha256']))
  tab=pd.DataFrame(rows);summary+=rows
  csv('comparison/'+{'SV-A_vs_OMAG-A':'svA_omagA_comparison.csv','SV-B_vs_OMAG-B':'svB_omagB_comparison.csv','OMAG-A_vs_OMAG-B':'stageOA_stageOB_common_support.csv'}[pair],tab)
  for (f,m),g in tab.groupby(['flow_mm_s','metric']):
   g=g.set_index('diameter_um').loc[D]
   for support in ['native','common']:
    a=g['left_'+support+'_median'].to_numpy();b=g['right_'+support+'_median'].to_numpy()
    lr=dict(pair=pair,support=support,flow_mm_s=int(f),metric=m,left_algorithm=left,right_algorithm=right,left_direction=s.direction(a),right_direction=s.direction(b),left_full_rank=ranking(a),right_full_rank=ranking(b),same_full_rank=ranking(a)==ranking(b),same_adjacent_directions=s.direction(a)==s.direction(b))
    for i,j in [(0,1),(1,2),(0,2)]:
     pa=100*(a[j]-a[i])/a[i] if a[i]!=0 else np.nan;pb=100*(b[j]-b[i])/b[i] if b[i]!=0 else np.nan
     lr[f'left_pct_{D[i]}_{D[j]}']=pa;lr[f'right_pct_{D[i]}_{D[j]}']=pb;lr[f'same_direction_{D[i]}_{D[j]}']=bool(np.sign(a[j]-a[i])==np.sign(b[j]-b[i]))
     for alg,x,p in [(left,a,pa),(right,b,pb)]:pct.append(dict(pair=pair,support=support,algorithm=alg,flow_mm_s=int(f),metric=m,from_diameter_um=D[i],to_diameter_um=D[j],percent_change=p,reason='ok' if np.isfinite(p) else 'zero_or_undefined_reference'))
    ledger.append(lr)
    if pair=='OMAG-A_vs_OMAG-B':oatr.append(lr)
 csv('comparison/common_support_summary.csv',summary);csv('validation/common_support_checks.csv',checks);csv('comparison/algorithm_direction_ledger.csv',ledger);csv('comparison/diameter_percent_change_by_algorithm.csv',pct);csv('comparison/stageOA_stageOB_support_trends.csv',oatr)
 return pd.DataFrame(ledger),pd.DataFrame(summary)
def validate(tables):
 rows=[];geom=read(R/'analysis/sv_physical_diameter_geometry_v2/geometry_manifest.csv.gz')
 for label,stage in [('OMAG-A','stageOA'),('OMAG-B','stageOB')]:
  t=tables[label];assert len(t)==7500 and not t.duplicated(['scan_id','frame_index']).any()
  z=t.merge(geom,on=['scan_id','frame_index'],validate='one_to_one');assert z.original_geometry_valid.equals(z.valid_geometry)
  assert np.allclose(z.x4_frozen,z.X4,rtol=0,atol=0,equal_nan=True) and np.allclose(z.z_top_frozen,z.z_top,rtol=0,atol=0,equal_nan=True)
  if stage=='stageOA':assert t.loc[~t.original_geometry_valid,M].isna().all().all()
  v=read(O/stage/'volume_metrics.csv').set_index('scan_id')
  for sid,g in t.groupby('scan_id'):
   assert set(g.frame_index)==set(range(500))
   for m in M:
    vals=g.loc[eligible(g,m),m].tolist();calc=median(vals);stored=v.loc[sid,'median_'+m];passed=np.isclose(calc,stored,rtol=1e-14,atol=0,equal_nan=True);assert passed and len(vals)==v.loc[sid,'n_pooled_'+m]
    rows.append(dict(stage=stage,scan_id=sid,metric=m,n_frames=len(vals),independent_sorted_median=calc,stored_median=stored,passed=bool(passed)))
 csv('validation/scan_median_reconstruction.csv',rows)
 for col in ['x4_frozen','z_top_frozen','original_geometry_valid']:
  a=tables['OMAG-A'].sort_values(['scan_id','frame_index'])[col].to_numpy();b=tables['SV-A'].sort_values(['scan_id','frame_index'])[col].to_numpy();assert np.allclose(a,b,rtol=0,atol=0,equal_nan=True)
 identity=read(O/'validation/retained_omag_identity.csv.gz');assert len(identity)==7500 and identity.status.eq('PASS').all()
 fixed=read(O/'validation/metric_formula_checks.csv');assert fixed.passed.all() and fixed.omag_stage.eq('OA').sum()==45
 after=[]
 for row in read(O/'validation/frozen_before.csv.gz').itertuples():
  h=sha(R/row.path);assert h==row.sha256;after.append(dict(path=row.path,sha256=h,unchanged=True))
 csv('validation/frozen_after.csv.gz',after)
 rawafter=[]
 for row in read(O/'raw_oct_input_manifest.csv').itertuples():
  h=sha(row.raw_path);assert h==row.raw_sha256;rawafter.append(dict(scan_id=row.scan_id,raw_sha256_after=h,unchanged=True))
 csv('validation/raw_oct_after.csv',rawafter)
 return dict(status='PASS',raw_identity_volumes=15,raw_sha_before_after_unchanged=15,common_preprocessing_frames=45,independent_omag_frames=45,retained_omag_files_sha_metadata_verified=7500,nominal_positions_per_stage=7500,original_geometry_valid=7018,original_invalid_preserved=482,independent_metric_OA=int(fixed.omag_stage.eq('OA').sum()),independent_metric_OB=int(fixed.omag_stage.eq('OB').sum()),independent_scan_medians=len(rows),frozen_files_unchanged=len(after))
def figures(summary):
 (O/'figures').mkdir(exist_ok=True)
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'savefig.dpi':220})
 colors=['#1670a9','#dc7c22','#23836b','#ba4565','#7655a2'];titles=['Raw tail AUC','Vessel-core P95','Normalized AUC']
 def save(fig,name):
  fig.savefig(O/'figures'/f'{name}.png');fig.savefig(O/'figures'/f'{name}.pdf');plt.close(fig)
 v=read(O/'stageOA/volume_metrics.csv');fig,axes=plt.subplots(1,3,figsize=(14,4.8))
 for ax,m,title,unit in zip(axes,M,['A  Raw tail AUC (nominal 0–200 μm)','B  Vessel-core P95','C  Normalized AUC'],['OMAG instrument units × μm','OMAG instrument units','AUC / P95 (μm)']):
  for f,c in zip(F,colors):
   g=v[v.flow_mm_s.eq(f)].set_index('diameter_um').loc[D];ax.plot(D,g['median_'+m],'-o',color=c,label=f'{f} mm/s')
  ax.set(xlabel='Diameter (μm)',ylabel=unit,title=title,xticks=D);ax.grid(axis='y',alpha=.2)
 axes[0].legend(frameon=False);fig.suptitle('OMAG O-A | Frozen SV geometry + senior metrics');fig.text(.5,.015,'Pooled frame medians; one scan per diameter × flow. Nominal 0–200 μm window: 29 pixels = 194.3 μm.',ha='center',fontsize=9);fig.tight_layout(rect=[0,.06,1,.94]);save(fig,'omag_senior_metrics_stageOA')
 for pair,name,title in [('SV-A_vs_OMAG-A','sv_vs_omag_same_geometry','Primary comparison | same geometry and endpoint-specific common frames'),('SV-B_vs_OMAG-B','sv_vs_omag_full_pipeline','Full pipeline | separate localization and assessability; common eligible frames')]:
  fig,axes=plt.subplots(2,3,figsize=(14,7),sharex=True,sharey='col')
  for j,(m,mt) in enumerate(zip(M,titles)):
   for i,alg in enumerate(['left','right']):
    ax=axes[i,j]
    for f,c in zip(F,colors):
     g=summary[summary.pair.eq(pair)&summary.metric.eq(m)&summary.flow_mm_s.eq(f)].set_index('diameter_um').loc[D];x=g[alg+'_common_median'].to_numpy();y=100*x/x[0] if x[0]!=0 else x*np.nan;ax.plot(D,y,'-o',color=c,label=f'{f} mm/s')
    ax.axhline(100,color='gray',lw=.7,alpha=.5);ax.set(title=('SV' if i==0 else 'OMAG')+' · '+mt,ylabel='Within algorithm: D128 = 100%',xticks=D);ax.grid(axis='y',alpha=.2)
    if i==1:ax.set_xlabel('Diameter (μm)')
  axes[0,0].legend(frameon=False,fontsize=8);fig.suptitle(title);fig.text(.5,.012,'Each line uses its own algorithm and flow D128 reference. No cross-algorithm amplitude ratio.',ha='center',fontsize=9);fig.tight_layout(rect=[0,.035,1,.94]);save(fig,name)
 t=read(O/'stageOB/tracking_summary.csv');va=read(O/'stageOA/coverage.csv');vb=read(O/'stageOB/coverage.csv');fig,ax=plt.subplots(1,3,figsize=(15,4.7));idx=np.arange(15)
 ax[0].bar(idx-.18,va.n_auc_defined/500*100,.36,label='O-A');ax[0].bar(idx+.18,vb.n_auc_defined/500*100,.36,label='O-B');ax[0].set(ylabel='Eligible AUC frames (%)',ylim=(0,105),title='Coverage');ax[0].legend(loc='lower left',frameon=True,facecolor='white',framealpha=.95)
 ax[1].bar(idx,t.median_abs_dx);ax[1].set(ylabel='Median |Δx| (pixels)',title='O-B vs frozen X4');ax[2].bar(idx,t.median_abs_dz);ax[2].set(ylabel='Median |Δz| (pixels)',title='O-B vs frozen z top')
 for a in ax:a.set_xticks(idx,[x.replace('_V01','') for x in t.scan_id],rotation=90,fontsize=7)
 fig.suptitle('QC | Localization and support differences');fig.tight_layout();save(fig,'omag_localization_coverage_qc')
def report(ledger,summary,val):
 ta=read(O/'stageOA/trend_ledger.csv');tb=read(O/'stageOB/trend_ledger.csv');primary=ledger.query("pair=='SV-A_vs_OMAG-A' and support=='common'");au=primary[primary.metric.eq(M[0])];diff=int((~au.same_adjacent_directions).sum())
 counts={stage:{m:int(t.loc[t.metric.eq(m),'target_shape_met'].sum()) for m in M} for stage,t in [('OA',ta),('OB',tb)]}
 lines=['# 同数据 OMAG × 师兄指标','',f'共同 preprocessing identity gate：**PASS**。15 卷 raw OCT 的前后 SHA 与正式来源一致；45 帧重建 SV 与 retained SV 逐像素完全相同（最大绝对误差 0）。同一 complex IMG 的独立读取检查及 raw OMAG 固定帧检查通过。','',f'O-A absolute AUC 与 SV-A 的相邻直径方向在 **{diff}/5** 个 flow 不同，以下为完全相同 geometry、endpoint-specific common support 的主比较。','', '| Flow (mm/s) | SV-A raw AUC | OMAG O-A raw AUC | SV-A P95 | O-A P95 | SV-A normalized | O-A normalized |','|---|---|---|---|---|---|---|']
 for f in F:
  p=primary[primary.flow_mm_s.eq(f)].set_index('metric');lines.append('| '+str(f)+' | '+' | '.join(str(p.loc[m,c]) for m in M for c in ['left_direction','right_direction'])+' |')
 lines+=['','预先冻结趋势的满足流速数：','', '| 预定义形态 | O-A | O-B |','|---|---|---|']
 for m,target in [(M[0],'raw AUC: D128<D235>D285'),(M[1],'P95: D128<D235<D285'),(M[2],'normalized: D128>D235>D285')]:lines.append(f"| {target} | {counts['OA'][m]}/5 | {counts['OB'][m]}/5 |")
 lines+=['','以上统计采用各阶段 native support；完整排序（包括非相邻直径的次序）、native/common support 和逐算法直径百分比变化保存在 comparison/。','', '## 判读','']
 if counts['OA'][M[0]]:
  lines.append('在同一原始 OCT、共同预处理、固定坐标、相同 senior metric 下，OMAG 与 SV 呈现不同的 absolute-tail diameter response；出现中间峰的具体 flow 见表。这支持当前数据中 absolute-tail 直径关系具有 algorithm dependence，不能据此证明机制不同。')
 elif ta[ta.metric.eq(M[0])].direction.eq('D128>D235>D285').all():
  lines.append('OMAG O-A 的 absolute AUC 也在全部五个流速呈 D128>D235>D285。师兄历史中间峰不能由 OMAG algorithm alone 在当前数据中重现。D128 与历史 D185 的直径范围、原始采集及其他 historical dataset differences 仍需区分；本轮未检验这些因素。')
 else:lines.append('OMAG O-A 的实际直径方向见表。当前数据并未稳定重现五个流速的中间峰；解释按逐流速观察限定，不以目标曲线筛选或调整结果。')
 commonob=ledger.query("pair=='OMAG-A_vs_OMAG-B' and support=='common'")
 for m,label in zip(M,['absolute AUC','P95 denominator','normalized AUC']):
  q=commonob[commonob.metric.eq(m)];lines.append(f"O-A → O-B 在 common support 下改变 {int((~q.same_adjacent_directions).sum())}/5 个 flow 的 {label} 相邻直径方向。")
 lines+=['','O-B 与 SV-B 的比较包含信号算法、算法依赖的定位和纳入变化，属于 full-pipeline comparison。O-A/O-B 的 common-support 结果用于检查定位对 denominator 排序的影响，不把覆盖变化和算法差异混成纯 signal algorithm effect。','', 'Normalized 方向相同不能证明算法等效；应同时读取 numerator 与 P95 的方向及各自的直径百分比变化。本轮不计算 OMAG AUC/SV AUC，也不以跨算法绝对单位的比值描述“高多少倍”。','', '## 四条 P95 轨迹与定位','', '| Flow | SV-A P95 | O-A P95 | SV-B P95 | O-B P95 | O-A common P95 | O-B common P95 |','|---|---|---|---|---|---|---|']
 for f in F:
  pa=primary[primary.metric.eq(M[1])&primary.flow_mm_s.eq(f)].iloc[0];pb=ledger.query("pair=='SV-B_vs_OMAG-B' and support=='native' and metric=='vessel_core_P95' and flow_mm_s==@f").iloc[0];po=commonob[commonob.metric.eq(M[1])&commonob.flow_mm_s.eq(f)].iloc[0]
  lines.append(f'| {f} | {pa.left_direction} | {pa.right_direction} | {pb.left_direction} | {pb.right_direction} | {po.left_direction} | {po.right_direction} |')
 va=read(O/'stageOA/volume_metrics.csv');vb=read(O/'stageOB/volume_metrics.csv')
 lines+=['','## 数据、计算与核验','',f'O-A 保留 7500 nominal positions、原 geometry-valid 7018、原 invalid 482。有效 AUC/P95/normalized pooled 帧数为 {int(va.n_pooled_raw_AUC_200.sum())}/{int(va.n_pooled_vessel_core_P95.sum())}/{int(va.n_pooled_normalized_AUC_um.sum())}。O-B 对应为 {int(vb.n_pooled_raw_AUC_200.sum())}/{int(vb.n_pooled_vessel_core_P95.sum())}/{int(vb.n_pooled_normalized_AUC_um.sum())}。0 AUC 保留，P95 必须 finite 且 >1e-12；不因低信号或趋势新增排除。','', '15 卷 OMAG 采用已验证 retained export 的 raw 第二输出无损复用：每个来源文件 SHA、frame index、OMAG NumEV、共同预处理元数据均核验，并由每卷三帧原始 OCT 独立重建校验。大型中间数组为 float64 NPY，逻辑顺序 z,x,frame，frame=0…499；SHA 和路径见 omag_reconstruction/。它们未进入 Git。没有声称本轮从 OCT 重新导出了全部 7500 帧。','', 'O-A 直接调用上一轮冻结 metric；O-B 使用原 senior tracking/X4/body assessability（alpha=.15、阈值≥.60、uncertain exclusion、原 numerical floors），定位工作数组按原 loader 转 float32，无 SV-specific scale mapping。指标保留原 raw float64，与已有 SV-A/SV-B adapter 的精度一致。原始 senior exporter 的 single 存储转换没有施加到无损 retained raw 上；这一点在合同中明确，避免额外量化成为主比较的变化项。','', '中央/单侧背景/gap 宽度分别为 4/8/9、3/6/7、2/5/6 px。先 round(D/dx)，z_top 按 Python ties-to-even；D/dz 为19/35/43，core offsets [2,9)/[4,16)/[5,20)，P95 为 linear quantile。P、合并双侧 B 按深度取 median，E=max(P−B,0)；AUC=[tail_start+1,tail_start+30)，29×6.7=194.3 μm。Normalized=AUC/P95，单位 μm。','',f'核验通过：45 帧共同预处理、45 帧 OMAG 独立重建、O-A/O-B 指标独立复算 {val["independent_metric_OA"]}/{val["independent_metric_OB"]} 帧、90 个 volume endpoint medians、135 个 endpoint-specific common-support 组合；15 卷原始 OCT 前后 SHA 不变，{val["frozen_files_unchanged"]} 个冻结文件前后 SHA 不变。所有 7500 个 OMAG frame finite。','', '这是同数据、同共同 preprocessing、同 senior metric 的 algorithm comparison within this dataset。每个 diameter×flow 只有一卷，空间帧不作为独立重复；不做帧级显著性检验或因果机制归因。历史最小直径为 D185，本轮为 D128，仅可比较 qualitative middle-diameter-peak shape，不能称 exact replication。','', '## 复现','', '在本目录依次运行 `python -B prepare.py`；MATLAB R2023a `rebuild_fixed`；`python -B validate_gate.py`；`python -B run_omag.py`；`python -B compare_report.py`；`python -B finalize.py`。原始路径和 retained 路径记录于 manifests；需在新的空输出副本复现，脚本拒绝覆盖完成的中间数组/汇总。依赖为 Python/NumPy/pandas/SciPy/Matplotlib 与 MATLAB 原函数；版本及源文件 SHA 见 provenance/contract。','', 'BLOCKED/MISSING：无。']
 (O/'README.md').write_text('\n'.join(lines)+'\n',encoding='utf-8');return counts
def main():
 tables={'SV-A':loadstage(OLD/'stageA'),'SV-B':loadstage(OLD/'stageB'),'OMAG-A':loadstage(O/'stageOA'),'OMAG-B':loadstage(O/'stageOB')}
 print('Independent medians, source hashes and frozen geometry',flush=True);val=validate(tables)
 print('Endpoint-specific common support',flush=True);ledger,summary=comparisons(tables);val['independent_common_support_checks']=len(read(O/'validation/common_support_checks.csv'));assert val['independent_common_support_checks']==135
 for stage in ['stageOA','stageOB']:
  tr=read(O/stage/'trend_ledger.csv');tr['full_rank_ascending']=[ranking(np.array([r.D128,r.D235,r.D285])) for r in tr.itertuples()];csv(stage+'/trend_ledger.csv',tr)
 figures(summary);counts=report(ledger,summary,val);val['trend_counts']=counts;js('validation.json',val);js('run_status.json',dict(status='COMPLETE_VALIDATED',common_preprocessing_gate='PASS',omag_fixed_checks='PASS',completed_volumes=15,blocked=[],missing=[]));print(json.dumps(val),flush=True)
if __name__=='__main__':main()
