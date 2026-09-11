"""One-time fidelity correction; reuse Stage A, preserve all provisional Stage B outputs."""
from run_transfer import *
import shutil

def main():
 arch=O/'validation/stageB_before_loader_dtype_correction_DIAGNOSTIC_ONLY'
 assert not arch.exists();shutil.copytree(O/'stageB',arch)
 shutil.copyfile(O/'validation/stageB_scale_guards.csv',arch/'stageB_scale_guards.csv')
 shutil.copyfile(O/'validation/fixed_frame_checks.csv',arch/'fixed_frame_checks_A_B.csv')
 correction=dict(recorded_utc=now(),reason='Late source-loader audit: both senior tracking_core.load_flow_volume and tailq.io.load_flow cast incoming Flow to float32 before localization/features. Add this unchanged original dtype behavior.',original_plan_and_trend_definitions_unchanged=True,stageA_reused=True,metric_input='unchanged retained SV_raw float64; no dtype change for endpoint measurement',stageB_localization_input='float32 in-memory working copy as original loader; no saved replacement raw arrays',parameter_tuning=False,previous_outputs=str(arch.relative_to(O)),additional_audit=['run_single_scan._make_features emits vessel_presence_prediction; _make_metrics expects assessability_class. Adapter uses original prediction field directly; class threshold unchanged.','Formal freeze does not add tracking qc_valid as a separate exclusion; it uses frozen assessability class plus primary metric validity.'])
 js('validation/loader_fidelity_correction.json',correction)
 track=loadmodule('senior_tracking',S/'02_track_vessel/tracking_core.py');dev=loadmodule('senior_development',S/'10_xroi_refinement/develop_xroi_assessability.py')
 geom=read(G/'geometry_manifest.csv.gz');ids=read(G/'array_identity_audit.csv.gz').fillna('');release=read(R/'results/formal_sv_d128_v21_full2500_run001/arrays_sha256.csv');fixed=set(map(tuple,read(G/'independent_integration_spotchecks.csv')[['scan_id','frame_index']].values));checks=read(O/'validation/fixed_frame_checks.csv').query("stage=='A'").to_dict('records');receipts=[];scale=[];BB=[]
 for sid,g in geom.groupby('scan_id',sort=True):
  print(now(),sid,'Stage B original float32 loader fidelity',flush=True);volume=input_volume(sid,geom,ids,release,receipts);d=int(g.diameter_um.iloc[0]);a=read(O/f'stageA/framewise_metrics_{sid}.csv.gz').set_index('frame_index');tt,ff,audit=senior_tracking(volume,sid,d,track,dev);scale.append(audit);csv(f'stageB/tracking_{sid}.csv.gz',tt);csv(f'stageB/features_{sid}.csv.gz',ff);rows=[]
  for fi in range(500):
   base={k:a.loc[fi,k] for k in ['scan_id','diameter_um','flow_mm_s','original_geometry_valid','x4_frozen','z_top_frozen']};base['frame_index']=fi;x=ff.iloc[fi].x_senior;z=tt.iloc[fi].z_upper_px;m=metric(volume[fi],x,z,d)
   row={**base,**m,'x_senior':x,'z_top_senior':z,'senior_tracking_status':tt.iloc[fi].tracking_class,'senior_tracking_qc_valid':bool(tt.iloc[fi].qc_valid),'senior_assessability_class':ff.iloc[fi].vessel_presence_prediction,'senior_assessability_score':ff.iloc[fi].assessability_score,'senior_full_assessable':bool(ff.iloc[fi].vessel_presence_prediction=='assessable'),'senior_full_assessability_reason':'ok' if ff.iloc[fi].vessel_presence_prediction=='assessable' else str(ff.iloc[fi].vessel_presence_prediction)};rows.append(row)
   if (sid,fi) in fixed and row['auc_defined'] and row['denominator_defined']:checks.append(independent_check(volume[fi],row,'B'))
  bt=pd.DataFrame(rows);BB.append(bt);csv(f'stageB/framewise_metrics_{sid}.csv.gz',bt);del volume
 v=volume_summary(pd.concat(BB,ignore_index=True));csv('stageB/volume_metrics.csv',v);csv('stageB/coverage.csv',v[[k for k in v if k.startswith('n_') or k in ['scan_id','diameter_um','flow_mm_s']]]);tr,con=trend_tables(v);csv('stageB/trend_ledger.csv',tr);csv('stageB/diameter_contrasts_by_flow.csv',con)
 csv('validation/stageB_scale_guards.csv',scale);csv('validation/fixed_frame_checks.csv',checks);csv('validation/array_identity_after_dtype_correction.csv.gz',receipts)
 before=read(arch/'volume_metrics.csv');comp=before.merge(v,on='scan_id',suffixes=('_before','_final'),validate='one_to_one');csv('validation/stageB_loader_correction_volume_comparison.csv',comp)
 prov=json.loads((O/'provenance.json').read_text(encoding='utf-8'));prov['analysis_script_sha256']=sha(O/'run_transfer.py');prov['loader_fidelity_correction']='validation/loader_fidelity_correction.json';prov['completed_metrics_utc']=now();js('provenance.json',prov)
 print('Stage B fidelity correction complete',tr.groupby('metric').target_shape_met.sum().to_dict(),flush=True)
if __name__=='__main__':main()
