"""Independent table/endpoint validation, direction comparisons and standalone figures."""
from run_transfer import *
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

def eligible(g,m):
 mask=g.senior_metric_assessable.astype(bool)&g.auc_defined.astype(bool)
 if 'senior_full_assessable' in g:mask &= g.senior_full_assessable.astype(bool)
 if m!='raw_AUC_200':mask &= g.denominator_defined.astype(bool)
 if m=='normalized_AUC_um':mask &= g.normalized_defined.astype(bool)
 return mask

def endpoint_checks():
 rows=[]
 for d,vw,cw,bw,gap,dp,cs,ce in [(128,10,4,3,2,19,2,9),(235,19,8,6,5,35,4,16),(285,22,9,7,6,43,5,20)]:
  for x,z in [(250.,180.),(250.5,180.5),(251.5,181.5)]:
   r=metric(np.arange(351*500,dtype=float).reshape(351,500),x,z,d);cx=round(x);zz=round(z);vl=cx-vw//2
   expected={'central_x_start':cx-cw//2,'central_x_end':cx-cw//2+cw,'background_L_start':vl-gap-bw,'background_L_end':vl-gap,'background_R_start':vl+vw+gap,'background_R_end':vl+vw+gap+bw,'z_lower':zz+dp,'tail_start':zz+dp+2,'auc_z_start':zz+dp+3,'tail_end':zz+dp+32,'core_start':zz+cs,'core_end':zz+ce}
   passed=all(r[k]==v for k,v in expected.items())
   rows.append(dict(test='literal_diameter_and_half_integer_endpoints',diameter_um=d,x=x,z=z,D_over_dx=d/DX,D_over_dz=d/DZ,vessel_width=vw,central_width=cw,background_width=bw,gap=gap,core_pixels=ce-cs,window_pixels=29,actual_integral_um=194.3,passed=passed,**expected));assert passed
 # Boundary tests explicitly preserve Python ties-to-even, including 0.4D and fraction-derived widths.
 for value,expected in [(10.5,10),(11.5,12),(np.nextafter(10.5,-np.inf),10),(np.nextafter(10.5,np.inf),11),(200/DZ,30),(0/DZ,0)]:
  passed=round(float(value))==expected;rows.append(dict(test='round_half_even_or_window_ratio',input_value=value,expected=expected,actual=round(float(value)),passed=passed));assert passed
 frame=np.ones((351,500));frame[:,248:252]=10
 cases=[('valid',frame,250,180,True,True,'ok'),('left_out_of_bounds',frame,0,180,False,False,'lateral_roi_incomplete_or_bg_less_than_6'),('right_out_of_bounds',frame,499,180,False,False,'lateral_roi_incomplete_or_bg_less_than_6'),('tail_partial',frame,250,330,False,True,'tail_out_of_bounds'),('zero_denominator',np.ones((351,500)),250,180,True,False,'denominator_nonpositive')]
 for name,im,x,z,av,dv,reason in cases:
  r=metric(im,x,z,128);passed=r['auc_defined']==av and r['denominator_defined']==dv and reason in [r['auc_reason'],r['denominator_reason']];assert passed,(name,r)
  if not av:assert np.isnan(r['raw_AUC_200'])
  if not dv:assert np.isnan(r['vessel_core_P95']) and np.isnan(r['normalized_AUC_um'])
  rows.append(dict(test=name,passed=passed,reason=reason))
 # Nonfinite samples propagate; positive truncation includes genuine zero AUC.
 for which in ['tail','core']:
  im=frame.copy();im[203 if which=='tail' else 182,250]=np.nan;r=metric(im,250,180,128)
  key='auc_defined' if which=='tail' else 'denominator_defined';assert not r[key];rows.append(dict(test=which+'_nan_propagates',passed=True,reason=r['auc_reason' if which=='tail' else 'denominator_reason']))
 im=np.zeros((351,500));im[:,248:252]=EPS/2;r=metric(im,250,180,128);assert not r['denominator_defined'] and r['denominator_reason']=='denominator_at_or_below_epsilon';rows.append(dict(test='positive_epsilon_guard',passed=True,reason=r['denominator_reason']))
 # Below-background tail is a real zero, not invalid/missing.
 im=frame.copy();im[202:231,248:252]=0;r=metric(im,250,180,128);assert r['auc_defined'] and r['raw_AUC_200']==0 and r['normalized_AUC_um']==0;rows.append(dict(test='negative_excess_truncates_to_valid_zero',passed=True,reason='ok'))
 csv('validation/rounding_endpoint_checks.csv',rows);return len(rows)

def load_stage(stage):return pd.concat([read(p) for p in sorted((O/stage).glob('framewise_metrics_*.csv.gz'))],ignore_index=True)

def all_validation(A,B):
 originals=read(G/'geometry_manifest.csv.gz');median_checks=[];preserved=[];nan_summary=[]
 for stage,t in [('stageA',A),('stageB',B)]:
  assert len(t)==7500 and not t.duplicated(['scan_id','frame_index']).any();assert t.groupby('scan_id').frame_index.apply(lambda x:set(x)==set(range(500))).all()
  merged=t.merge(originals,on=['scan_id','frame_index'],suffixes=('','_orig'),validate='one_to_one')
  assert (merged.original_geometry_valid==merged.valid_geometry).all()
  assert np.allclose(merged.x4_frozen,merged.X4,equal_nan=True,rtol=0,atol=0) and np.allclose(merged.z_top_frozen,merged.z_top,equal_nan=True,rtol=0,atol=0)
  assert t.original_geometry_valid.sum()==7018
  if stage=='stageA':assert t.loc[~t.original_geometry_valid,METRICS].isna().all().all()
  v=read(O/stage/'volume_metrics.csv').set_index('scan_id')
  for sid,g in t.groupby('scan_id'):
   assert v.loc[sid,'n_metric_assessable']==int(g.senior_metric_assessable.sum())
   assert v.loc[sid,'n_nominal']==len(g) and v.loc[sid,'n_geometry_valid']==int(g.original_geometry_valid.sum())
   for flag,key in [('auc_defined','n_auc_defined_all_positions'),('denominator_defined','n_denominator_defined_all_positions'),('normalized_defined','n_normalized_defined_all_positions')]:assert v.loc[sid,key]==int(g[flag].sum())
   for m in METRICS:
    arr=np.asarray(g.loc[eligible(g,m),m],float);arr=arr[np.isfinite(arr)];s=np.sort(arr);n=len(s);med=float((s[(n-1)//2]+s[n//2])/2) if n else np.nan
    stored=v.loc[sid,'median_'+m];passed=bool(np.isclose(stored,med,rtol=1e-14,atol=0,equal_nan=True));assert passed
    assert n==v.loc[sid,'n_pooled_'+m]
    median_checks.append(dict(stage=stage,scan_id=sid,metric=m,n_independent=n,independent_sorted_median=med,stored_median=stored,passed=passed))
  for m,flag,reason in [('raw_AUC_200','auc_defined','auc_reason'),('vessel_core_P95','denominator_defined','denominator_reason'),('normalized_AUC_um','normalized_defined','normalized_reason')]:
   assert np.array_equal(np.isfinite(t[m]),t[flag]);assert t.loc[~t[flag],reason].notna().all() and t.loc[~t[flag],reason].ne('ok').all()
   for rs,gg in t.groupby(reason):nan_summary.append(dict(stage=stage,metric=m,reason=rs,n=len(gg),n_missing=int(gg[m].isna().sum())))
  preserved.append(dict(stage=stage,n_nominal=7500,n_original_valid=7018,n_original_invalid=482,frozen_values_exact=True,passed=True))
 for row in read(O/'input_manifest.csv').itertuples():assert sha(R/row.input_path)==row.sha256,row.input_path
 identity=read(O/'validation/array_identity_on_read.csv.gz');assert len(identity)==7500 and identity.matched_original_export_manifest.all() and identity.matched_geometry_v2.sum()==7018
 fixed=read(O/'validation/fixed_frame_checks.csv');assert fixed.passed.all();assert (fixed.stage=='A').sum()==45
 scale=read(O/'validation/stageB_scale_guards.csv');assert len(scale)==15 and scale.status.eq('passed').all() and scale.n_positive_body_sigma_floor.eq(0).all()
 # Export manifests also checked against frozen intake result's hash ledger.
 exp_hashes=read(R/'results/diameter_intake_v1/output_sha256.csv').set_index('relative_path')
 manifest_receipts=[]
 for sid in sorted(A.scan_id.unique()):
  if sid.startswith('D128'):continue
  rel=f'results/diameter_intake_v1/exports/{sid}_frames.csv';h=sha(R/rel);assert h==exp_hashes.loc[rel,'sha256'];manifest_receipts.append(dict(input_path=rel,sha256=h,matched_frozen_export_ledger=True))
 csv('validation/export_manifest_checks.csv',manifest_receipts);csv('validation/scan_median_reconstruction.csv',median_checks);csv('validation/frozen_geometry_checks.csv',preserved);csv('validation/nan_reason_counts.csv',nan_summary)
 after=read(O/'validation/array_identity_after_dtype_correction.csv.gz') if (O/'validation/array_identity_after_dtype_correction.csv.gz').exists() else identity
 cross=identity.merge(after,on=['scan_id','frame_index'],suffixes=('_initial','_final'),validate='one_to_one');assert len(cross)==7500 and (cross.sha256_initial==cross.sha256_final).all()
 nend=endpoint_checks()
 return dict(status='passed',nominal_positions_per_stage=7500,original_geometry_valid=7018,original_invalid_preserved=482,raw_arrays_matched_original_exports=7500,raw_arrays_matched_fixed_v2=7018,fixed_frame_checks=int(len(fixed)),stageA_fixed_frame_checks=45,stageB_fixed_frame_checks=int((fixed.stage=='B').sum()),endpoint_invalid_checks=nend,independent_scan_medians=len(median_checks),input_and_history_files_unchanged=len(read(O/'input_manifest.csv')),export_manifest_checks=len(manifest_receipts),scale_guard_scans_passed=15,no_reconstruction=True,no_inferential_tests=True,no_d500=True,validation_utc=now())

def comparisons(A,B):
 common=[]
 for sid,a in A.groupby('scan_id'):
  b=B[B.scan_id.eq(sid)].set_index('frame_index');a=a.set_index('frame_index');assert a.index.equals(b.index)
  for m in METRICS:
   ma=eligible(a,m);mb=eligible(b,m);both=ma&mb
   am=float(a.loc[ma,m].median());bm=float(b.loc[mb,m].median());ac=float(a.loc[both,m].median());bc=float(b.loc[both,m].median())
   row=dict(scan_id=sid,diameter_um=int(a.diameter_um.iloc[0]),flow_mm_s=int(a.flow_mm_s.iloc[0]),metric=m,n_A=int(ma.sum()),n_B=int(mb.sum()),n_common=int(both.sum()),n_A_only=int((ma&~mb).sum()),n_B_only=int((mb&~ma).sum()),A_native_median=am,B_native_median=bm,A_common_median=ac,B_common_median=bc,common_difference_B_minus_A=bc-ac,common_percent_difference=100*(bc-ac)/ac if ac!=0 else np.nan,common_difference_reason='ok' if ac!=0 and np.isfinite([ac,bc]).all() else 'zero_or_undefined_A_common',delta_native_B_minus_A=bm-am,delta_A_native_to_common=ac-am,delta_B_common_to_native=bm-bc)
   # Exact additive median-difference decomposition (description, no causal attribution).
   assert np.isclose(row['delta_native_B_minus_A'],row['delta_A_native_to_common']+row['common_difference_B_minus_A']+row['delta_B_common_to_native'],rtol=1e-12,atol=1e-8)
   common.append(row)
 common=pd.DataFrame(common);csv('comparison/stageA_stageB_common_support.csv',common)
 ct=[]
 for (m,f),g in common.groupby(['metric','flow_mm_s']):
  g=g.set_index('diameter_um')
  for scope,col in [('A_native','A_native_median'),('B_native','B_native_median'),('A_common','A_common_median'),('B_common','B_common_median')]:
   x=g.loc[DIMS,col].to_numpy();ct.append(dict(metric=m,flow_mm_s=f,support=scope,direction=direction(x),D128=x[0],D235=x[1],D285=x[2]))
 csv('comparison/stageA_stageB_support_trends.csv',ct)
 current=read(V/'volume_metrics.csv');records=[];ledger=[]
 for f in FLOWS:
  row={'flow_mm_s':f}
  for name,m in [('current_tail100','tail100'),('current_tail500','tail500'),('current_RI100','RI100'),('current_RI500','RI500'),('current_source','source')]:
   g=current[current.flow_mm_s.eq(f)&current.metric.eq(m)].set_index('diameter_um');x=g.loc[DIMS,'median'].to_numpy();row[name+'_direction']=direction(x);row[name+'_rank_ascending']=' < '.join('D'+str(DIMS[i]) for i in np.argsort(x));records.append(dict(flow_mm_s=f,method='current_SV_v2',metric=m,unit=g.unit.iloc[0],direction=direction(x),D128=x[0],D235=x[1],D285=x[2]))
  for stage,prefix in [('stageA','stageA'),('stageB','stageB')]:
   v=read(O/stage/'volume_metrics.csv');v=v[v.flow_mm_s.eq(f)].set_index('diameter_um')
   for m,key,unit in zip(METRICS,['AUC','denominator','normalized'],['SV instrument units * um','SV instrument units','um']):
    x=v.loc[DIMS,'median_'+m].to_numpy();row[prefix+'_'+key+'_direction']=direction(x);row[prefix+'_'+key+'_rank_ascending']=' < '.join('D'+str(DIMS[i]) for i in np.argsort(x));records.append(dict(flow_mm_s=f,method=stage,metric=m,unit=unit,direction=direction(x),D128=x[0],D235=x[1],D285=x[2]))
  ledger.append(row)
 csv('comparison/current_vs_senior_definition.csv',records);csv('comparison/metric_definition_direction_ledger.csv',ledger)
 return pd.DataFrame(ledger),pd.DataFrame(ct)

def figures(ledger):
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'savefig.dpi':220})
 colors=['#1670a9','#dc7c22','#23836b','#ba4565','#7655a2']
 for stage,suffix in [('stageA',''),('stageB','_stageB')]:
  v=read(O/stage/'volume_metrics.csv');fig,axes=plt.subplots(1,3,figsize=(14,4.8))
  titles=['A  Raw tail AUC, nominal 0–200 μm','B  Vessel-core P95 denominator','C  Normalized AUC']
  ys=['SV instrument units × μm','SV instrument units','AUC / core P95 (μm)']
  for ax,m,title,ylabel in zip(axes,METRICS,titles,ys):
   for f,col in zip(FLOWS,colors):
    g=v[v.flow_mm_s.eq(f)].set_index('diameter_um');ax.plot(DIMS,g.loc[DIMS,'median_'+m],marker='o',lw=1.8,ms=5,color=col,label=f'{f} mm/s')
   ax.set(xlabel='Vessel diameter (μm)',ylabel=ylabel,title=title,xticks=DIMS);ax.grid(axis='y',alpha=.18);ax.ticklabel_format(axis='y',style='sci',scilimits=(-3,4),useMathText=True)
  axes[0].legend(frameon=False,fontsize=9)
  fig.suptitle('SV signal + senior metric definitions  |  '+('Stage A: frozen localization' if stage=='stageA' else 'Stage B: senior localization and assessability'),fontsize=14,y=.98)
  fig.text(.5,.015,'One scan per diameter × flow; pooled frame medians. Senior code integrates 29 pixels (194.3 μm); no independent frame replicates.',ha='center',fontsize=9)
  fig.tight_layout(rect=[0,.045,1,.92]);fig.savefig(O/f'figures/senior_style_sv_metrics{suffix}.png');fig.savefig(O/f'figures/senior_style_sv_metrics{suffix}.pdf');plt.close(fig)
 fields=['current_tail100','current_tail500','current_RI100','current_RI500','current_source','stageA_AUC','stageA_denominator','stageA_normalized','stageB_AUC','stageB_denominator','stageB_normalized'];labels=['Current SV v2 · tail100','Current SV v2 · tail500','Current SV v2 · RI100','Current SV v2 · RI500','Current SV v2 · whole source','Stage A · AUC','Stage A · core P95','Stage A · normalized AUC','Stage B · AUC','Stage B · core P95','Stage B · normalized AUC']
 categories={'D128>D235>D285':0,'D128<D235>D285':1,'D128<D235<D285':2,'D128>D235<D285':3,'UNDEFINED':4};mat=np.array([[categories[ledger.loc[j,k+'_direction']] for j in range(5)] for k in fields]);fig,ax=plt.subplots(figsize=(11,6.6));ax.imshow(mat,cmap=ListedColormap(['#cfe4f1','#fae2bf','#d5eade','#e7ddf0','#eeeeee']),vmin=0,vmax=4,aspect='auto')
 for i,k in enumerate(fields):
  for j in range(5):ax.text(j,i,ledger.loc[j,k+'_direction'].replace('D128','128 ').replace('D235',' 235 ').replace('D285',' 285'),ha='center',va='center',fontsize=10)
 ax.set_xticks(range(5),[f'{f} mm/s' for f in FLOWS]);ax.set_yticks(range(len(labels)),labels);ax.set_title('Matched-flow diameter directions | SV v2 vs definition/process transfer',pad=15);ax.tick_params(length=0)
 for y in [4.5,7.5]:ax.axhline(y,color='white',lw=4)
 fig.text(.5,.018,'Directions only. Different metric definitions and units are not compared by absolute magnitude.',ha='center',fontsize=9);fig.tight_layout(rect=[0,.04,1,1]);fig.savefig(O/'figures/definition_direction_comparison.png');fig.savefig(O/'figures/definition_direction_comparison.pdf');plt.close(fig)

def report(ledger,ct,val):
 ta=read(O/'stageA/trend_ledger.csv');tb=read(O/'stageB/trend_ledger.csv');va=read(O/'stageA/volume_metrics.csv');vb=read(O/'stageB/volume_metrics.csv');common=read(O/'comparison/stageA_stageB_common_support.csv')
 lines=['# SV + 师兄拖尾指标定义迁移分析','', 'Stage A 与 Stage B 均已完成，使用与正式 SV v2 完全相同的 15 卷 retained SV_raw。正式信号仍为 `var(abs(IMG),1,3)`，方差分母 N；dx=12.7 μm/px，dz=6.7 μm/px。旧结果和 geometry 均保持只读。','', '## 主要结果','', '| 指标与执行前冻结形态 | Stage A 满足流速 | Stage B 满足流速 |','|---|---|---|']
 for m,target in [('raw_AUC_200','D128 < D235 > D285'),('vessel_core_P95','D128 < D235 < D285'),('normalized_AUC_um','D128 > D235 > D285')]:
  a=ta[ta.metric.eq(m)];b=tb[tb.metric.eq(m)];af=a.loc[a.target_shape_met,'flow_mm_s'].tolist();bf=b.loc[b.target_shape_met,'flow_mm_s'].tolist()
  lines.append(f'| {m}: {target} | {len(af)}/5；{af} mm/s | {len(bf)}/5；{bf} mm/s |')
  lines += []
 lines+=['','Stage A 未满足各形态的流速：']
 for m in METRICS:lines.append(f"- {m}：{ta.loc[ta.metric.eq(m)&~ta.target_shape_met,'flow_mm_s'].tolist()} mm/s。")
 lines+=['','## 与正式 SV v2 的方向比较','', '正式 tail100、tail500、RI100、RI500 在五个匹配流速均为 D128 > D235 > D285。下表保留新方法每一个流速的方向；不同单位之间不比较数值大小。','', '| Flow (mm/s) | A: AUC | A: P95 | A: normalized | B: AUC | B: P95 | B: normalized |','|---|---|---|---|---|---|---|']
 for r in ledger.itertuples():lines.append(f'| {r.flow_mm_s} | {r.stageA_AUC_direction} | {r.stageA_denominator_direction} | {r.stageA_normalized_direction} | {r.stageB_AUC_direction} | {r.stageB_denominator_direction} | {r.stageB_normalized_direction} |')
 lines+=['','Whole-source mean 与 vessel-core P95 的范围、横向汇总和背景处理不同。正式 source 在 1/3/5/7 mm/s 递增，10 mm/s 为 D128 < D235 > D285。Stage A 的 P95 相比 source 仅在 7 mm/s 改变这一直径方向；Stage B 又使 7/10 mm/s 的 P95 均变为递增。完整三直径排序见 `comparison/metric_definition_direction_ledger.csv`。', '', '## 对四个问题的判断', '', '迁移指标定义后，absolute AUC 的中间直径峰值在五个流速均未出现；进一步迁移定位与 assessability 后仍未出现。AUC 与 normalized AUC 的直径排序均保持正式 tail/RI 的递减方向。因此本轮不能把既有 absolute-tail 形态差异解释为单纯的指标定义差异。', '', 'Denominator 对定义与定位流程更敏感：Stage A 在 7 mm/s 相比当前 whole-source 排序发生改变；Stage B 在 7/10 mm/s 相比 Stage A 改变排序，并且相同 common support 下仍然改变。这支持定位对 denominator 方向有影响，不能仅归于帧覆盖变化；它不证明某个因素就是历史差异的原因。', '', 'Normalized endpoint 虽与师兄同向，却也与当前正式 RI 同向，不能据此称迁移后完整三联图得到复现。对于未出现的 raw-AUC 中间峰，在同一数据和 SV 信号下迁移两阶段流程仍不能得到师兄形态，使 signal algorithm / dataset 对差异的贡献更值得考虑；本轮不能单独区分 OMAG 算法与两套原始数据，也不进行机制归因。', '', '## 定位与覆盖比较','',f'Stage A 保留全部 7500 个 nominal positions，其中原 geometry valid 7018，原 invalid 482；原 invalid 的指标仍为 NaN。Stage A 新增 ROI 可评估数 {int(va.n_metric_assessable.sum())}，AUC 有效数 {int(va.n_auc_defined.sum())}，denominator 有效数 {int(va.n_denominator_defined.sum())}，normalized 有效数 {int(va.n_normalized_defined.sum())}。', '',f'Stage B 在原始 0–499 空间序列上运行师兄 alpha=0.15 定位、X4 与 body assessability，未覆盖冻结字段。完整流程 assessable 数 {int(vb.n_full_assessable.sum())}，AUC 有效纳入数 {int(vb.n_auc_defined.sum())}，normalized 有效纳入数 {int(vb.n_normalized_defined.sum())}。', '', 'Common support 分别针对 AUC、P95、normalized endpoint 取 A/B 同时纳入的位置。`comparison/stageA_stageB_common_support.csv` 保留各自 median、共同帧 median、差值和纳入变化；`stageA_stageB_support_trends.csv` 对两类 support 分别判断趋势。']
 for m in METRICS:
  c=common[common.metric.eq(m)];same=[];different=[]
  for f in FLOWS:
   p=ct[ct.metric.eq(m)&ct.flow_mm_s.eq(f)].set_index('support').direction
   (same if p['A_common']==p['B_common'] else different).append(f)
  lines+=['',f'{m}：共同帧集合总计 {int(c.n_common.sum())}，每卷 {int(c.n_common.min())}–{int(c.n_common.max())} 个位置。共同 support 下 A/B 方向保持的 flow 为 {same}，改变的 flow 为 {different} mm/s。']
 lines+=['', '各自 support 的总差值在 CSV 中分解为 A-native→A-common、A-common→B-common、B-common→B-native 三项。这是描述性差值分解：共同帧上的变化对应定位改变后的指标，support 项对应帧集合改变，不能据此进行因果归因。','', '## 代码审计中必须保留的定义','', '- 师兄正式入口调用 `formal_9scan_freeze._interval_result`，其 `_window_pixels(0)` 返回 1。因此名义 0–200 μm 实际为 `[tail_start+1, tail_start+30)`，共 29 像素、194.3 μm。操作性 `tail_start=z_lower+2` 保留，实际首个积分样本相对 z_lower 为 +3 px。未修正该历史实现。','- 旧 generic/development helper 使用 0 起点，得到 30 像素；本次使用正式 freeze 入口。图中明确为 nominal 0–200 μm。','- 先 round(D/dx)，再计算 central/background 宽度及 gap。D128/D235/D285 central 宽度分别为 4/8/9 px；单侧 background 为 3/6/7 px，gap 为 2/5/6 px。','- 左右背景像素合并后逐深度取 median；central 同样取 median，E=max(P−B,0)。AUC 为 sum×dz。','- D/dz 取整为 19/35/43 px；core 相对整数化 z_top 为 [2,9)、[4,16)、[5,20)，P95 使用 NumPy 默认 linear quantile。','- 原 frozen z_top 以边界坐标保存，可能为半整数。新方法严格使用 Python ties-to-even round，另存 z_upper_metric；x4_frozen/z_top_frozen 原值未改。','- denominator 必须 finite 且 >1e-12。无效主字段为 NaN，候选值另存 denominator_candidate。normalized=AUC/P95，单位 μm，正式 runner 不加 epsilon。','- Stage A 只使用几何/有限值/ROI/窗口可定义性，不提前施加 Stage B 的 body assessability score。扫描值遵循 pooled median；denominator 主汇总限于 AUC 纳入帧且 P95 有效。','- Stage B 使用相对 SNR、原几何参数及 >=0.60 assessability 阈值，uncertain 不纳入主值。numeric floor 原值保留并核查；未做 SV 尺度映射或调参。MATLAB Thr=85 dB 属于 OMAG 重建，本次不执行。','- `senior_source/` 保存实际读取代码的字节相同快照与配置，完整 source SHA 和文档差异见 `senior_method_contract.json`。', '- Stage B 定位与 assessability 使用师兄加载器原有的 float32 内存副本；指标仍使用原 retained SV_raw 精度。补齐加载精度后重算了全部 Stage B，修正前表仅在 validation/ 的 DIAGNOSTIC_ONLY 目录保留，不进入最终分析。初始合同原文、字段映射和修正记录见 `senior_code_audit.md`。','', '## 分析身份与解释边界','', '这是 method-definition transfer / sensitivity analysis，不是新的独立实验验证。历史师兄结果在分析前已知，因此不是盲验证或预注册。D128 与师兄 D185 不同，相似方向只能称 qualitative shape similarity，不能称 exact replication。','', 'B-scan/frame 是 slow-axis 空间位置，不是时间重复。15 volumes 是 15 个 diameter×flow 条件扫描；每个条件仅一卷，不同 flow 不是同条件重复，也不是每个直径五次独立重复。本次不做 frame-level p-value、复杂模型或机制归因。','', '归一化结果必须结合 raw AUC 和 P95 同时解释。定义改变后趋势是否改变只能说明对定义/流程的敏感性；不能单独证明哪种定义是差异原因，也不能区分 OMAG 算法与两套原始数据的贡献。','', '## 验证与文件','',f"验证状态：{val['status']}。7500 个数组均匹配原 export SHA，其中 7018 个逐一匹配固定 v2 的数组 SHA；45 个 Stage A 固定帧和 {val['stageB_fixed_frame_checks']} 个 Stage B 固定帧独立重算公式；{val['independent_scan_medians']} 个 scan median 从 framewise 表独立重建；{val['endpoint_invalid_checks']} 项 rounding/NaN/FOV 检查通过。冻结输入与历史文件 {val['input_and_history_files_unchanged']} 个 SHA 保持不变。",'', '- `stageA/`、`stageB/`：7500 行 framewise、15 卷 summary、coverage、逐 flow contrasts/trends。','- `comparison/`：正式 v2 与迁移定义的方向 ledger，以及 A/B common-support 数值与趋势。','- `figures/`：SV 三联图与方向比较图，PNG/PDF。','- `validation/`、`validation.json`：数组身份、固定帧、rounding、NaN reason、独立 median 和尺度保护核验。','- `output_sha256.csv`：交付文件 SHA256（不含自身）。','', '复跑命令：在新 checkout 中保留相同输入布局及合同/来源快照，将本轮结果另行归档后，在固定 analysis/sv_senior_metric_transfer_v1 路径运行 `python -B run_transfer.py`、`python -B validate_compare_report.py`。最终 run_transfer 已包含 float32 定位步骤，无需再运行一次性 correct_stageB_loader_fidelity.py；脚本拒绝覆盖已有 volume_metrics。输入路径由 manifest 给出，不运行原始 OCT/OMAG 导出。','', '![SV Stage A](figures/senior_style_sv_metrics.png)','', '![SV Stage B](figures/senior_style_sv_metrics_stageB.png)','', '![Direction comparison](figures/definition_direction_comparison.png)','']
 (O/'README.md').write_text('\n'.join(lines),encoding='utf-8')

def main():
 A=load_stage('stageA');B=load_stage('stageB');validation=all_validation(A,B);ledger,ct=comparisons(A,B);figures(ledger);report(ledger,ct,validation);js('validation.json',validation);js('run_status.json',dict(status='analysis_validated_figures_pending_visual_review',stageA='success',stageB='success',missing=[],blocked=[],validation='passed'))
 print(json.dumps(validation,ensure_ascii=False),flush=True)
 for st in ['stageA','stageB']:print(st,read(O/st/'trend_ledger.csv').groupby('metric').target_shape_met.sum().to_dict())
if __name__=='__main__':main()
