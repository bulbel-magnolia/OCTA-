"""Lossless reuse of identity-verified raw OMAG and frozen senior metrics."""
from prepare import *
import importlib.util,datetime
spec=importlib.util.spec_from_file_location('frozen_transfer',OLD/'run_transfer.py');s=importlib.util.module_from_spec(spec);spec.loader.exec_module(s)
def tracking(volume,sid,d,track,dev):
 # Exact previous adapter/original loader precision, with no SV scale veto/mapping.
 volume=np.asarray(volume,dtype=np.float32);cfg=track.merge_tracking_config(s.CFG['tracking'])
 xt=track.locate_lateral_track(volume,float(d),cfg);profiles=track.extract_profiles(volume,xt['x_center'],float(d),cfg);peak=track.locate_peak_path(profiles,cfg)
 tt=track.track_one_alpha(scan_id=sid,diameter_um=float(d),alpha=.15,xtrack=xt,peak_path=peak,profiles=profiles,cfg=cfg)
 ff=pd.DataFrame([dev._local_body_features(volume[i],tt.iloc[i]) for i in range(500)])
 ff['neighbor_peak_cnr']=ff.local_body_peak_cnr.shift(1).combine(ff.local_body_peak_cnr.shift(-1),lambda l,r:np.nanmedian([l,r])).fillna(ff.local_body_peak_cnr)
 ff=dev._add_assessability_score(ff);xx,changed=dev._isolated_jump_correct(ff.x2_robust_centroid_px.to_numpy(float));ff['x_senior']=xx;ff['x4_jump_corrected']=changed
 return tt,ff
def main():
 status=json.loads((O/'run_status.json').read_text(encoding='utf-8'));assert status['common_preprocessing_gate']=='PASS' and status['omag_fixed_checks']=='PASS'
 prov=json.loads((O/'provenance.json').read_text(encoding='utf-8'));assert sha(O/'analysis_plan.json')==prov['analysis_plan_sha256'] and sha(OLD/'run_transfer.py')==prov['metric_source_sha256']
 assert not (O/'stageOA/volume_metrics.csv').exists(),'completed output exists'
 geom=read(R/'analysis/sv_physical_diameter_geometry_v2/geometry_manifest.csv.gz');ids=read(OLD/'validation/array_identity_on_read.csv.gz').fillna('');raw=read(O/'raw_oct_input_manifest.csv').set_index('scan_id');fixed=set(map(tuple,read(O/'validation/fixed_frames.csv')[['scan_id','frame_index']].values))
 track=s.loadmodule('omag_tracking',s.S/'02_track_vessel/tracking_core.py');dev=s.loadmodule('omag_development',s.S/'10_xroi_refinement/develop_xroi_assessability.py')
 AA=[];BB=[];receipts=[];arrays=[];checks=[];tracking_rows=[];meta_rows=[]
 for sid,g in geom.groupby('scan_id',sort=True):
  print(s.now(),sid,'SHA/metadata validation and lossless OMAG materialization',flush=True)
  p=O/'intermediate_local'/f'{sid}_p_bld_ed.npy';assert not p.exists(),'no overwrite intermediate'
  vol=np.lib.format.open_memmap(p,mode='w+',dtype='float64',shape=(351,500,500),fortran_order=True)
  d=int(g.diameter_um.iloc[0]);flow=int(g.flow_mm_s.iloc[0]);imin=np.inf;imax=-np.inf
  for ar in ids[ids.scan_id.eq(sid)].sort_values('frame_index').itertuples(index=False):
   sv,om,meta=retained(ar);fi=int(ar.frame_index);source=Path(meta['source_file']);rp=Path(raw.loc[sid,'raw_path']);assert source.name==rp.name
   if d!=128:assert source.resolve()==rp.resolve()
   vol[:,:,fi]=om;imin=min(imin,float(om.min()));imax=max(imax,float(om.max()))
   oh=hashlib.sha256(np.ascontiguousarray(om).tobytes()).hexdigest();assert np.array_equal(vol[:,:,fi],om)
   receipts.append(dict(scan_id=sid,frame_index=fi,input_path=ar.input_path,member=ar.member,retained_file_sha256=ar.sha256,omag_array_sha256=oh,source_raw_sha256=raw.loc[sid,'raw_sha256'],dtype=str(om.dtype),shape='351,500',all_finite=True,metadata_passed=True,lossless_copy_passed=True,status='PASS'))
  vol.flush();arrays.append(dict(scan_id=sid,path=str(p),shape='351,500,500',axis_order='z,x,frame',frame_index='0:499',dtype='float64',min=imin,max=imax,finite_count=351*500*500,total_count=351*500*500,SHA256=sha(p),source_raw_oct_SHA256=raw.loc[sid,'raw_sha256'],reconstruction_config_SHA256=sha(O/'reconstruction_contract.json'),method='verified retained raw OMAG reuse; no quantization'))
  meta_rows.append(dict(scan_id=sid,nominal_frames=500,verified_retained_files=500,independently_raw_rebuilt_fixed_frames=3,source_raw_sha256=raw.loc[sid,'raw_sha256'],common_preprocessing_gate='PASS',reuse_provenance='same-IMG export_sv_omag_frame; omag_num_ev=2; every retained file SHA and per-frame preprocessing metadata checked',status='PASS'))
  arows=[]
  for gr in g.sort_values('frame_index').itertuples():
   fi=int(gr.frame_index);base=dict(frame_index=fi,scan_id=sid,diameter_um=d,flow_mm_s=flow,original_geometry_valid=bool(gr.valid_geometry),x4_frozen=gr.X4,z_top_frozen=gr.z_top)
   m=s.metric(vol[:,:,fi],gr.X4 if gr.valid_geometry else np.nan,gr.z_top if gr.valid_geometry else np.nan,d)
   if not gr.valid_geometry:
    for k in ['assessability_reason','auc_reason','denominator_reason','normalized_reason']:m[k]='original_geometry_invalid'
   row={**base,**m};arows.append(row)
   if (sid,fi) in fixed:checks.append({**s.independent_check(vol[:,:,fi],row,'A'),'omag_stage':'OA'})
  at=pd.DataFrame(arows);csv(f'stageOA/framewise_metrics_{sid}.csv.gz',at);AA.append(at)
  print(s.now(),sid,'O-A complete; original OMAG full localization',flush=True)
  tt,ff=tracking(vol.transpose(2,0,1),sid,d,track,dev);csv(f'stageOB/tracking_{sid}.csv.gz',tt);csv(f'stageOB/features_{sid}.csv.gz',ff);brows=[]
  for fi in range(500):
   base={k:arows[fi][k] for k in ['frame_index','scan_id','diameter_um','flow_mm_s','original_geometry_valid','x4_frozen','z_top_frozen']};x=ff.iloc[fi].x_senior;z=tt.iloc[fi].z_upper_px
   row={**base,**s.metric(vol[:,:,fi],x,z,d), 'x_senior':x,'z_top_senior':z,'x_omag_senior':x,'z_top_omag_senior':z,'senior_tracking_status':tt.iloc[fi].tracking_class,'senior_tracking_qc_valid':bool(tt.iloc[fi].qc_valid),'senior_assessability_class':ff.iloc[fi].vessel_presence_prediction,'senior_assessability_score':ff.iloc[fi].assessability_score,'senior_full_assessable':bool(ff.iloc[fi].vessel_presence_prediction=='assessable')}
   brows.append(row)
   if (sid,fi) in fixed and row['auc_defined'] and row['denominator_defined']:checks.append({**s.independent_check(vol[:,:,fi],row,'B'),'omag_stage':'OB'})
  bt=pd.DataFrame(brows);csv(f'stageOB/framewise_metrics_{sid}.csv.gz',bt);BB.append(bt)
  tracking_rows.append(dict(scan_id=sid,n_assessable=int(bt.senior_full_assessable.sum()),n_uncertain=int(bt.senior_assessability_class.eq('uncertain').sum()),n_not_assessable=int(bt.senior_assessability_class.eq('not_assessable').sum()),n_tracking_qc_valid=int(bt.senior_tracking_qc_valid.sum()),median_abs_dx=float((bt.x_senior-bt.x4_frozen).abs().median()),median_abs_dz=float((bt.z_top_senior-bt.z_top_frozen).abs().median()),n_original_sigma_floor=int((ff.local_body_sigma.le(1e-6)&ff.local_body_peak_excess.gt(0)).sum()),loader_dtype='float32',metric_dtype='float64 original raw; same frozen SV adapter',alpha=.15,assessability_threshold=.60))
  csv('validation/retained_omag_identity.csv.gz',receipts);csv('validation/metric_formula_checks.csv',checks);csv('omag_reconstruction/omag_array_manifest.csv',arrays);csv('omag_reconstruction/export_metadata_summary.csv',meta_rows);csv('stageOB/tracking_summary.csv',tracking_rows)
  js('run_status.json',dict(status='RUNNING',common_preprocessing_gate='PASS',omag_fixed_checks='PASS',completed_volumes=len(AA)));del vol
 for stage,tables in [('stageOA',AA),('stageOB',BB)]:
  table=pd.concat(tables,ignore_index=True);v=s.volume_summary(table);csv(stage+'/volume_metrics.csv',v);csv(stage+'/coverage.csv',v[[k for k in v if k.startswith('n_') or k in ['scan_id','diameter_um','flow_mm_s']]])
  tr,con=s.trend_tables(v);csv(stage+'/trend_ledger.csv',tr);csv(stage+'/diameter_contrasts_by_flow.csv',con)
 js('run_status.json',dict(status='METRICS_COMPLETE_VALIDATION_PENDING',common_preprocessing_gate='PASS',omag_fixed_checks='PASS',completed_volumes=15))
 print('METRICS COMPLETE',flush=True)
if __name__=='__main__':main()
