"""Freeze design, inventory raw inputs and prepare fixed retained references."""
from pathlib import Path
import sys,json,hashlib,zipfile,io,shutil,platform
sys.dont_write_bytecode=True
import numpy as np
import pandas as pd
from scipy.io import loadmat,savemat
O=Path(__file__).resolve().parent;R=O.parents[1];OLD=R/'analysis/sv_senior_metric_transfer_v1'
def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def csv(n,x):
 p=O/n;p.parent.mkdir(parents=True,exist_ok=True);pd.DataFrame(x).to_csv(p,index=False,compression={'method':'gzip','mtime':0} if str(n).endswith('.gz') else None)
def js(n,x):(O/n).write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
def read(p):return pd.read_csv(p,float_precision='round_trip',low_memory=False)
def retained(row):
 p=Path(row.input_path)
 if row.member:
  with zipfile.ZipFile(p) as z:blob=z.read(row.member)
  with np.load(io.BytesIO(blob),allow_pickle=False) as z:
   sv=z['sv_raw'];om=z['omag_raw'];meta=json.loads(z['metadata_json'].item())['reconstruction']['export_metadata']
 else:
  blob=p.read_bytes();m=loadmat(io.BytesIO(blob),simplify_cells=True);sv=m['sv_raw'];om=m['omag_raw'];meta=m['metadata']
 assert hashlib.sha256(blob).hexdigest()==row.sha256
 assert meta['bscan_index_matlab_1based']==int(row.frame_index)+1 and meta['omag_num_ev']==2
 a=meta['reconstruction'];assert list(a['IMcropRg'])==list(range(50,401)) and a['registration_Nsub']==10 and not a['registration_Colshift'] and a['phase_threshold_db']==85 and a['do_median_bulk_phase_shift'] and a['do_registration'] and a['do_phase_compensation']
 assert all(not a[k] for k in ['subtract_reference','useDVV','useKES','use_dispersion_coefficient'])
 assert sv.shape==om.shape==(351,500) and np.isfinite(sv).all() and np.isfinite(om).all()
 return sv,om,meta
def main():
 (O/'intermediate_local').mkdir(exist_ok=True);(O/'validation').mkdir(exist_ok=True)
 (O/'.gitignore').write_text('intermediate_local/\n__pycache__/\n',encoding='utf-8');(O/'.gitattributes').write_text('* -text\n',encoding='utf-8')
 geom=read(R/'analysis/sv_physical_diameter_geometry_v2/geometry_manifest.csv.gz');ids=read(OLD/'validation/array_identity_on_read.csv.gz').fillna('')
 fixed=[]
 for sid,g in geom.groupby('scan_id'):
  ff=sorted(g.loc[g.valid_geometry,'frame_index'].astype(int));fixed += [dict(scan_id=sid,frame_index=f,role=role) for f,role in zip([ff[0],ff[(len(ff)-1)//2],ff[-1]],['first_valid','lower_median_valid','last_valid'])]
 plan=dict(baseline_commit='b365c276c20a16d22ac51afd891b8d924819db86',diameters=[128,235,285],flows=[1,3,5,7,10],nominal_frames=7500,original_geometry_valid=7018,fixed_frames=fixed,frame_selection='first/lower-median/last ordered original geometry-valid frame; frozen before reconstruction',identity_tolerance=dict(rtol=1e-10,atol=1e-6,explanation='Float64 reconstruction on same MATLAB release; tight relative FFT roundoff allowance; no post-result tuning. Dtype casts reported separately.'),metrics_source='../sv_senior_metric_transfer_v1/run_transfer.py',primary='SV-A vs O-A endpoint-specific intersection, with native support retained',secondary='SV-B vs O-B full algorithm-dependent localization and support',trends=dict(raw_AUC_200='D128<D235>D285',vessel_core_P95='D128<D235<D285',normalized_AUC_um='D128>D235>D285'),raw_storage='float64 retained raw, lossless; same precision as SV comparator; original float32 working loader for full-pipeline localization',prohibited=['D500','parameter tuning','new exclusions','cross-algorithm amplitude ratios','frame p values'])
 js('analysis_plan.json',plan);csv('validation/fixed_frames.csv',fixed)
 before=[]
 for name in ['sv_fixed_diameter_selective_v2','sv_senior_metric_transfer_v1','sv_physical_diameter_geometry_v2']:
  for p in sorted((R/'analysis'/name).rglob('*')):
   if p.is_file():before.append(dict(path=str(p.relative_to(R)),sha256=sha(p),bytes=p.stat().st_size))
 csv('validation/frozen_before.csv.gz',before)
 src=R.parents[2]/'师兄方法/Bscan_tail_auto_quantification_for_junior_20260903(1)/matlab'
 (O/'reconstruction_source').mkdir(exist_ok=True);deps=[]
 for n in ['iniset.m','ini2struct.m','OCTA_F_SubPixReg.m','SSOCT_F_PhCompV3.m','OCTA_F_ED_Clutter_EigFeed.m','config.ini']:
  target=O/'reconstruction_source'/n;shutil.copyfile(src/n,target)
  deps.append(dict(name=n,source_path=str(src/n),actual_path=str(target),sha256=sha(target),retained_export_dependency_path=str(R/'matlab'/n),retained_export_dependency_sha256=sha(R/'matlab'/n),text_equal_ignoring_line_endings=(target.read_bytes().replace(b'\r\n',b'\n')==(R/'matlab'/n).read_bytes().replace(b'\r\n',b'\n'))))
 contract=dict(dependencies=deps,Thr_dB=85,Thr_linear=10**(85/20),crop_matlab=[50,400],Nsub=10,Colshift=False,NumEV=2,doMedianShift=True,repeats=3,phase_second_arg=1,signal='second output OCTA_F_ED_Clutter_EigFeed(IMG,2), raw linear unfiltered',axis_order='z,x,frame',sv_formula='var(abs(IMG),1,3)',plan_sha256=sha(O/'analysis_plan.json'),metric_contract_sha256=sha(OLD/'senior_method_contract.json'))
 js('reconstruction_contract.json',contract)
 manifest=read(R/'data/diameter_stage_v1/volume_manifest.csv');records=[];jobs=[]
 for sid,g in geom.groupby('scan_id'):
  rr=manifest[manifest.scan_id.eq(sid)].iloc[0];raw=Path('D:/OCTA原始数据')/rr.oct_relative_path
  print(sid,'hashing original OCT',flush=True);h=sha(raw);assert h==rr.oct_sha256 and raw.stat().st_size==rr.oct_size_bytes
  selected=[f for f in fixed if f['scan_id']==sid]
  for f in selected:
   ar=ids[ids.scan_id.eq(sid)&ids.frame_index.eq(f['frame_index'])].iloc[0];sv,om,meta=retained(ar)
   assert Path(meta['source_file']).name==raw.name
   if int(rr.diameter_um)!=128:assert Path(meta['source_file']).resolve()==raw.resolve()
   ref=O/'intermediate_local'/f"{sid}_{f['frame_index']:03d}_reference.mat";savemat(ref,dict(sv_retained=sv,omag_retained=om))
   jobs.append(dict(scan_id=sid,frame_index=f['frame_index'],raw_path=str(raw),reference_path=str(ref),output_path=str(O/'intermediate_local'/f"{sid}_{f['frame_index']:03d}_rebuilt.mat")))
  records.append(dict(scan_id=sid,diameter_um=int(rr.diameter_um),flow_mm_s=int(rr.flow_speed_mm_s),raw_path=str(raw),raw_bytes=raw.stat().st_size,raw_sha256=h,expected_raw_sha256=rr.oct_sha256,source_provenance='data/diameter_stage_v1/volume_manifest.csv; retained export_metadata; frozen array_identity_on_read',retained_sv_mapping=f"../sv_senior_metric_transfer_v1/validation/array_identity_on_read.csv.gz#{sid}",n_retained_frames=int(ids.scan_id.eq(sid).sum()),status='PASS'))
  csv('raw_oct_input_manifest.csv',records)
 js('intermediate_local/jobs.json',jobs)
 js('provenance.json',dict(baseline_commit=plan['baseline_commit'],analysis_plan_sha256=sha(O/'analysis_plan.json'),reconstruction_contract_sha256=sha(O/'reconstruction_contract.json'),raw_manifest_sha256=sha(O/'raw_oct_input_manifest.csv'),python=platform.python_version(),numpy=np.__version__,metric_source_sha256=sha(OLD/'run_transfer.py'),frozen_input_manifest_sha256=sha(O/'validation/frozen_before.csv.gz')))
 js('run_status.json',dict(status='PREPARED_IDENTITY_PENDING',raw_identity='15/15 PASS'))
if __name__=='__main__':main()
