"""SV + audited senior metric definitions. No source reconstruction or outcome tuning."""
from pathlib import Path
import sys, json, hashlib, io, zipfile, importlib.util, platform, datetime, subprocess
sys.dont_write_bytecode=True
import numpy as np
import pandas as pd
import scipy
from scipy.io import loadmat
O=Path(__file__).resolve().parent;R=O.parents[1];G=R/'analysis/sv_physical_diameter_geometry_v2';V=R/'analysis/sv_fixed_diameter_selective_v2'
S=O/'senior_source/python';sys.path.insert(0,str(S))
from tailq import metrics as sm
DIMS=[128,235,285];FLOWS=[1,3,5,7,10];DX=12.7;DZ=6.7;EPS=1e-12
METRICS=['raw_AUC_200','vessel_core_P95','normalized_AUC_um']
CFG=json.loads((O/'senior_source/config/example_config.json').read_text(encoding='utf-8'))
def loadmodule(name,path):
 spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m
freeze=loadmodule('senior_freeze',S/'13_formal_9scan/formal_9scan_freeze.py')
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def js(n,x):(O/n).write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def csv(n,x):
 p=O/n;p.parent.mkdir(parents=True,exist_ok=True)
 pd.DataFrame(x).to_csv(p,index=False,float_format='%.17g',lineterminator='\n',compression={'method':'gzip','mtime':0} if str(n).endswith('.gz') else None)
def read(p):return pd.read_csv(p,float_precision='round_trip',low_memory=False)
def metric(frame,x,z,d):
 """Execute only the primary audited interval and denominator; all bounds half open."""
 row={k:np.nan for k in METRICS+['normalized_AUC','x_metric','z_upper_metric','central_x_start','central_x_end','background_L_start','background_L_end','background_R_start','background_R_end','z_lower','tail_start','tail_end','auc_z_start','core_start','core_end','denominator_candidate','n_effective_depth_pixels','actual_integration_depth_um']}
 row.update(senior_metric_assessable=False,assessability_reason='missing_localization',auc_defined=False,denominator_defined=False,normalized_defined=False,auc_reason='missing_localization',denominator_reason='missing_localization',normalized_reason='missing_localization',lateral_fov_complete=False,tail_fov_complete=False,core_fov_complete=False,tail_n_finite_pixels=0,nominal_window_um=200)
 if not np.isfinite([x,z]).all():return row
 p=sm.extract_profiles(frame,x_center=x,diameter_um=d,lateral_um_per_px=DX,metrics_cfg=CFG['metrics'])
 zint=int(round(float(z)));dp=max(1,int(round(d/DZ)));lo=zint+dp;ts=lo+2;a=ts+sm._window_pixels(0,DZ);b=ts+sm._window_pixels(200,DZ)
 c=zint+int(np.ceil(.10*dp));e=zint+int(np.floor(.45*dp))+1
 row.update(x_metric=int(round(float(x))),z_upper_metric=zint,central_x_start=p.central_bounds[0],central_x_end=p.central_bounds[1],background_L_start=p.left_bounds[0],background_L_end=p.left_bounds[1],background_R_start=p.right_bounds[0],background_R_end=p.right_bounds[1],z_lower=lo,tail_start=ts,auc_z_start=a,tail_end=b,core_start=c,core_end=e,lateral_fov_complete=p.roi_complete,tail_fov_complete=bool(0<=a<b<=len(frame)),core_fov_complete=bool(0<=c<e<=len(frame)))
 if not p.roi_complete:
  row.update(assessability_reason='lateral_roi_incomplete_or_bg_less_than_6',auc_reason='lateral_roi_incomplete_or_bg_less_than_6',denominator_reason='lateral_roi_incomplete_or_bg_less_than_6',normalized_reason='lateral_roi_incomplete_or_bg_less_than_6');return row
 row.update(senior_metric_assessable=True,assessability_reason='ok')
 item=freeze._interval_result(p,ts,0,200)
 den,dvalid,cc,ee=sm._normalizer(p.excess,zint,dp,CFG['metrics']['core_fraction'],'p95',EPS)
 assert (cc,ee)==(c,e)
 row['denominator_candidate']=den
 if not row['core_fov_complete']:dr='core_out_of_bounds'
 elif not np.isfinite(p.excess[c:e]).all():dr='core_nonfinite'
 elif not np.isfinite(den):dr='denominator_nonfinite'
 elif den<=0:dr='denominator_nonpositive'
 elif den<=EPS:dr='denominator_at_or_below_epsilon'
 else:dr='ok'
 av=bool(item['complete'] and np.isfinite(item['raw_auc']))
 ar='ok' if av else ('tail_out_of_bounds' if not row['tail_fov_complete'] else 'tail_nonfinite')
 nv=av and bool(dvalid)
 row.update(raw_AUC_200=item['raw_auc'] if av else np.nan,vessel_core_P95=den if dvalid else np.nan,normalized_AUC_um=item['raw_auc']/den if nv else np.nan,normalized_AUC=item['raw_auc']/den if nv else np.nan,auc_defined=av,denominator_defined=bool(dvalid),normalized_defined=nv,auc_reason=ar,denominator_reason=dr,normalized_reason='ok' if nv else ';'.join(s for s in [ar,dr] if s!='ok'),tail_n_finite_pixels=int(np.isfinite(p.excess[max(0,a):min(len(frame),b)]).sum()),n_effective_depth_pixels=b-a if av else 0,actual_integration_depth_um=(b-a)*DZ if av else np.nan)
 return row

def independent_check(frame,row,stage):
 """Separate direct implementation using sorted medians, explicit linear P95 and fsum."""
 import math
 d=int(row['diameter_um']);x=row['x4_frozen'] if stage=='A' else row['x_senior'];z=row['z_top_frozen'] if stage=='A' else row['z_top_senior']
 vw=round(d/DX);cw=max(3,round(vw*.4));bw=max(1,round(vw/3));gap=max(1,round(vw/4));center=round(float(x));vl=center-vw//2;vr=vl+vw;c0=center-cw//2;c1=c0+cw;l0=vl-gap-bw;l1=vl-gap;r0=vr+gap;r1=r0+bw
 def med(vals):
  a=sorted(float(v) for v in vals);n=len(a);return (a[(n-1)//2]+a[n//2])/2
 pp=np.array([med(a[c0:c1]) for a in frame]);bb=np.array([med(list(a[l0:l1])+list(a[r0:r1])) for a in frame]);ex=np.array([max(float(p-b),0.) for p,b in zip(pp,bb)])
 zz=round(float(z));dp=round(d/DZ);ts=zz+dp+2;aa=ts+1;ee=ts+30;cs=zz+math.ceil(.1*dp);ce=zz+math.floor(.45*dp)+1
 raw=DZ*math.fsum(ex[aa:ee]);core=sorted(ex[cs:ce]);rank=(len(core)-1)*.95;il=math.floor(rank);ih=math.ceil(rank);den=core[il]+(core[ih]-core[il])*(rank-il);norm=raw/den if den>EPS else np.nan
 orig=sm.extract_profiles(frame,x_center=x,diameter_um=d,lateral_um_per_px=DX,metrics_cfg=CFG['metrics'])
 rec=dict(stage=stage,scan_id=row['scan_id'],frame_index=row['frame_index'],max_profile_error=float(np.max(abs(pp-orig.profile))),max_background_error=float(np.max(abs(bb-orig.background))),max_excess_error=float(np.max(abs(ex-orig.excess))))
 rec['bounds_passed']=tuple(row[k] for k in ['central_x_start','central_x_end','background_L_start','background_L_end','background_R_start','background_R_end','auc_z_start','tail_end','core_start','core_end'])==(c0,c1,l0,l1,r0,r1,aa,ee,cs,ce)
 for m,v in zip(METRICS,[raw,den if den>EPS else np.nan,norm]):
  rec['independent_'+m]=v;rec['output_'+m]=row[m];rec[m+'_passed']=bool(np.isclose(row[m],v,rtol=1e-12,atol=0,equal_nan=True))
 rec['passed']=rec['bounds_passed'] and all(rec[m+'_passed'] for m in METRICS) and rec['max_excess_error']==0
 assert rec['passed'],rec
 return rec

def volume_summary(table):
 rows=[]
 for sid,g in table.groupby('scan_id',sort=True):
  eligible=g.senior_metric_assessable.astype(bool)
  if 'senior_full_assessable' in g:eligible &= g.senior_full_assessable.astype(bool)
  rawmask=eligible & g.auc_defined;denmask=rawmask & g.denominator_defined;normmask=rawmask & g.normalized_defined
  rec=dict(scan_id=sid,diameter_um=int(g.diameter_um.iloc[0]),flow_mm_s=int(g.flow_mm_s.iloc[0]),n_nominal=len(g),n_geometry_valid=int(g.original_geometry_valid.sum()),n_metric_assessable=int(g.senior_metric_assessable.sum()),n_primary_eligible=int(eligible.sum()),n_auc_defined=int((eligible&g.auc_defined).sum()),n_denominator_defined=int((eligible&g.denominator_defined).sum()),n_normalized_defined=int(normmask.sum()),n_denominator_pooled=int(denmask.sum()),n_auc_defined_all_positions=int(g.auc_defined.sum()),n_denominator_defined_all_positions=int(g.denominator_defined.sum()),n_normalized_defined_all_positions=int(g.normalized_defined.sum()))
  if 'senior_full_assessable' in g:rec['n_full_assessable']=int(g.senior_full_assessable.sum())
  for m,mask in zip(METRICS,[rawmask,denmask,normmask]):
   vals=g.loc[mask,m];rec['median_'+m]=vals.median();rec['min_'+m]=vals.min();rec['max_'+m]=vals.max();rec['n_pooled_'+m]=int(vals.count())
  rows.append(rec)
 return pd.DataFrame(rows)

def trend_tables(v):
 rows=[];contrasts=[]
 for f in FLOWS:
  sub=v[v.flow_mm_s.eq(f)].set_index('diameter_um')
  for m in METRICS:
   x=np.array([sub.loc[d,'median_'+m] for d in DIMS]);ok=bool(np.isfinite(x).all());a,b,c=x
   shape={'raw_AUC_200':a<b and b>c,'vessel_core_P95':a<b<c,'normalized_AUC_um':a>b>c}[m] if ok else False
   rows.append(dict(flow_mm_s=f,metric=m,D128=a,D235=b,D285=c,evaluable=ok,target_shape_met=bool(shape),direction=direction(x),reason='ok' if ok else 'one_or_more_volume_medians_undefined'))
   for i,j in [(0,1),(1,2),(0,2)]:
    finite=np.isfinite(x[[i,j]]).all();pctok=finite and x[i]!=0
    contrasts.append(dict(flow_mm_s=f,metric=m,from_diameter_um=DIMS[i],to_diameter_um=DIMS[j],from_value=x[i],to_value=x[j],absolute_difference=x[j]-x[i] if finite else np.nan,percent_difference=100*(x[j]-x[i])/x[i] if pctok else np.nan,reason='ok' if pctok else ('zero_reference' if finite else 'undefined_volume_median')))
 return pd.DataFrame(rows),pd.DataFrame(contrasts)
def direction(x):
 if not np.isfinite(x).all():return 'UNDEFINED'
 return 'D128'+('<' if x[0]<x[1] else '>' if x[0]>x[1] else '=')+'D235'+('<' if x[1]<x[2] else '>' if x[1]>x[2] else '=')+'D285'

def input_volume(sid,geom,ids,release,receipts):
 g=geom[geom.scan_id.eq(sid)];d=int(g.diameter_um.iloc[0]);f=int(g.flow_mm_s.iloc[0]);valid_ids=ids[ids.scan_id.eq(sid)].set_index('frame_index');first=valid_ids.iloc[0];volume=np.empty((500,351,500),dtype=np.float64);cache={}
 if d==128:
  expected=release[release.scan_id.eq(f'flow{f:02d}')].set_index('frame_index_0based')
  packages={Path(p).name:Path(p) for p in valid_ids.input_path.unique()}
 else:
  mp=R/f'results/diameter_intake_v1/exports/{sid}_frames.csv';expected=read(mp).set_index('frame_index');assert len(expected)==500
  parent=Path(first.input_path).parent
 for fi in range(500):
  if d==128:
   member=f'arrays/flow{f:02d}/frame_{fi:03d}.npz';start=fi//100*100;name=f'formal_sv_d128_v21_full2500_run001_flow{f:02d}_{start:03d}_{start+99:03d}.zip';p=packages[name]
   if name not in cache:cache[name]=zipfile.ZipFile(p)
   blob=cache[name].read(member)
   with np.load(io.BytesIO(blob),allow_pickle=False) as z:arr=np.asarray(z['sv_raw'],float)
  else:
   member='';p=parent/f'frame_{fi:03d}.mat';blob=p.read_bytes();mat=loadmat(io.BytesIO(blob),variable_names=['sv_raw','metadata'],simplify_cells=True);arr=np.asarray(mat['sv_raw'],float);meta=mat['metadata'];assert meta['bscan_index_matlab_1based']==fi+1 and meta['formal_signal_definition']=='var(abs(E), 1, 3)' and meta['variance_denominator']=='N'
  h=hashlib.sha256(blob).hexdigest();assert h==expected.loc[fi,'sha256'],(sid,fi,'release identity')
  original=fi in valid_ids.index
  if original:assert h==valid_ids.loc[fi,'sha256'],(sid,fi,'v2 identity')
  assert arr.shape==(351,500) and np.isfinite(arr).all() and (arr>=0).all()
  volume[fi]=arr
  receipts.append(dict(scan_id=sid,frame_index=fi,diameter_um=d,flow_mm_s=f,input_path=str(p),member=member,sha256=h,matched_original_export_manifest=True,matched_geometry_v2=original,identity_scope='v2+original_export' if original else 'original_export; original_geometry_invalid',array_dtype='float64_loaded',status='passed'))
 for z in cache.values():z.close()
 return volume

class ScaleGuardBlocked(RuntimeError):pass

def senior_tracking(volume,sid,d,track,dev):
 """Original primary tracking + X4 + body assessability, with observational scale checks."""
 # Both senior HDF5 loaders cast Flow to float32 before tracking/features.
 # This is an in-memory localization working copy; metric() still sees original SV_raw.
 volume=volume.astype(np.float32,copy=False)
 cfg=track.merge_tracking_config(CFG['tracking']);audit={'positive_sigma_min':float('inf'),'n_numerical_sigma_zero':0}
 orig=track.robust_sigma
 def observed(values,axis=None):
  sig=orig(values,axis=axis);a=np.asarray(sig);pos=a[np.isfinite(a)&(a>0)];audit['n_numerical_sigma_zero']+=int((a==0).sum())
  if len(pos):
   audit['positive_sigma_min']=min(audit['positive_sigma_min'],float(pos.min()))
   if (pos<=np.finfo(np.float32).eps).any():raise ScaleGuardBlocked('positive robust sigma reaches absolute float32 epsilon; no scale mapping')
  return sig
 track.robust_sigma=observed
 try:
  xt=track.locate_lateral_track(volume,float(d),cfg);profiles=track.extract_profiles(volume,xt['x_center'],float(d),cfg);peak=track.locate_peak_path(profiles,cfg)
  tt=track.track_one_alpha(scan_id=sid,diameter_um=float(d),alpha=.15,xtrack=xt,peak_path=peak,profiles=profiles,cfg=cfg)
 finally:track.robust_sigma=orig
 feats=[]
 for i in range(500):feats.append(dev._local_body_features(volume[i],tt.iloc[i]))
 ff=pd.DataFrame(feats)
 if 'local_body_peak_cnr' not in ff:raise track.TrackingError('No local body features')
 ff['neighbor_peak_cnr']=ff.local_body_peak_cnr.shift(1).combine(ff.local_body_peak_cnr.shift(-1),lambda l,r:np.nanmedian([l,r])).fillna(ff.local_body_peak_cnr)
 ff=dev._add_assessability_score(ff);corrected,changed=dev._isolated_jump_correct(ff.x2_robust_centroid_px.to_numpy(float));ff['x_senior']=corrected;ff['x4_jump_corrected']=changed;ff['frame_index']=range(500);ff['scan_id']=sid
 floor=ff.local_body_sigma.le(1e-6)&ff.local_body_peak_excess.gt(0)
 if floor.any():raise ScaleGuardBlocked('local body positive signal activates 1e-6 sigma floor; no objective mapping')
 audit.update(scan_id=sid,n_positive_body_sigma_floor=int(floor.sum()),min_body_sigma=float(ff.local_body_sigma.min()),status='passed',parameters_unchanged=True)
 return tt,ff,audit

def main():
 plan=json.loads((O/'analysis_plan.json').read_text(encoding='utf-8'));prov=json.loads((O/'provenance.json').read_text(encoding='utf-8'));assert sha(O/'analysis_plan.json')==prov['analysis_plan_sha256'];assert sha(O/'senior_method_contract.json')==prov['contract_sha256']
 for s in json.loads((O/'senior_method_contract.json').read_text(encoding='utf-8'))['sources']:assert sha(O/s['snapshot_path'])==s['sha256']
 assert not (O/'stageA/volume_metrics.csv').exists(),'Do not overwrite completed analysis'
 geom=read(G/'geometry_manifest.csv.gz');ids=read(G/'array_identity_audit.csv.gz').fillna('');vids=read(V/'array_identity_on_read.csv.gz');assert len(geom)==7500 and geom.valid_geometry.sum()==7018 and len(ids)==7018
 chk=ids.merge(vids,on=['scan_id','frame_index'],suffixes=('_g','_v'),validate='one_to_one');assert (chk.sha256_g==chk.sha256_v).all()
 release=read(R/'results/formal_sv_d128_v21_full2500_run001/arrays_sha256.csv');fixed=set(map(tuple,read(G/'independent_integration_spotchecks.csv')[['scan_id','frame_index']].values))
 track=loadmodule('senior_tracking',S/'02_track_vessel/tracking_core.py');dev=loadmodule('senior_development',S/'10_xroi_refinement/develop_xroi_assessability.py')
 AA=[];BB=[];receipts=[];checks=[];scale=[];errors=[]
 for sid,g in geom.groupby('scan_id',sort=True):
  print(now(),sid,'loading and hashing 500 retained arrays',flush=True);volume=input_volume(sid,geom,ids,release,receipts);d=int(g.diameter_um.iloc[0]);flow=int(g.flow_mm_s.iloc[0]);arows=[]
  for gr in g.sort_values('frame_index').itertuples():
   fi=int(gr.frame_index);base=dict(frame_index=fi,scan_id=sid,diameter_um=d,flow_mm_s=flow,original_geometry_valid=bool(gr.valid_geometry),x4_frozen=gr.X4,z_top_frozen=gr.z_top)
   m=metric(volume[fi],gr.X4,gr.z_top,d) if gr.valid_geometry else metric(volume[fi],np.nan,np.nan,d)
   if not gr.valid_geometry:
    for k in ['assessability_reason','auc_reason','denominator_reason','normalized_reason']:m[k]='original_geometry_invalid'
   row={**base,**m};arows.append(row)
   if (sid,fi) in fixed:checks.append(independent_check(volume[fi],row,'A'))
  at=pd.DataFrame(arows);csv(f'stageA/framewise_metrics_{sid}.csv.gz',at);AA.append(at)
  print(now(),sid,'Stage A done; original senior tracking alpha=.15',flush=True)
  try:
   tt,ff,audit=senior_tracking(volume,sid,d,track,dev);scale.append(audit);csv(f'stageB/tracking_{sid}.csv.gz',tt);csv(f'stageB/features_{sid}.csv.gz',ff);brows=[]
   for fi in range(500):
    base={k:arows[fi][k] for k in ['frame_index','scan_id','diameter_um','flow_mm_s','original_geometry_valid','x4_frozen','z_top_frozen']};x=ff.iloc[fi].x_senior;z=tt.iloc[fi].z_upper_px
    m=metric(volume[fi],x,z,d);row={**base,**m,'x_senior':x,'z_top_senior':z,'senior_tracking_status':tt.iloc[fi].tracking_class,'senior_tracking_qc_valid':bool(tt.iloc[fi].qc_valid),'senior_assessability_class':ff.iloc[fi].vessel_presence_prediction,'senior_assessability_score':ff.iloc[fi].assessability_score,'senior_full_assessable':bool(ff.iloc[fi].vessel_presence_prediction=='assessable'),'senior_full_assessability_reason':'ok' if ff.iloc[fi].vessel_presence_prediction=='assessable' else str(ff.iloc[fi].vessel_presence_prediction)}
    brows.append(row)
    if (sid,fi) in fixed and row['auc_defined'] and row['denominator_defined']:checks.append(independent_check(volume[fi],row,'B'))
   bt=pd.DataFrame(brows);csv(f'stageB/framewise_metrics_{sid}.csv.gz',bt);BB.append(bt)
  except (ScaleGuardBlocked,track.TrackingError) as e:
   errors.append(dict(scan_id=sid,type=type(e).__name__,reason=str(e)));print('Stage B blocked',sid,str(e),flush=True)
  csv('validation/fixed_frame_checks.csv',checks);csv('validation/array_identity_on_read.csv.gz',receipts);csv('validation/stageB_scale_guards.csv',scale)
  js('run_status.json',dict(status='running',stageA_scans=len(AA),stageB_scans=len(BB),stageB_errors=errors));del volume
 for stage,tables in [('stageA',AA),('stageB',BB)]:
  if not tables:continue
  table=pd.concat(tables,ignore_index=True);v=volume_summary(table);csv(stage+'/volume_metrics.csv',v);csv(stage+'/coverage.csv',v[[k for k in v if k.startswith('n_') or k in ['scan_id','diameter_um','flow_mm_s']]])
  if len(v)==15:
   tr,con=trend_tables(v);csv(stage+'/trend_ledger.csv',tr);csv(stage+'/diameter_contrasts_by_flow.csv',con)
 if errors:(O/'stageB/BLOCKED.md').write_text('# Stage B incomplete\n\n'+json.dumps(errors,ensure_ascii=False,indent=2),encoding='utf-8')
 prov.update(software=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__),analysis_script_sha256=sha(__file__),completed_metrics_utc=now());js('provenance.json',prov)
 js('run_status.json',dict(status='metrics_complete_validation_pending',stageA='success',stageB='success' if len(BB)==15 else 'BLOCKED_PARTIAL',stageA_scans=15,stageB_scans=len(BB),stageB_errors=errors))
 print('All executable metrics completed; independent validation/report next',flush=True)
if __name__=='__main__':main()
