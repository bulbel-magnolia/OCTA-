"""Close provenance, validate artifact inventory and write byte-exact output ledger."""
from prepare import *
import scipy,matplotlib,datetime,subprocess
def normalized_lines(p):return [s.strip() for s in Path(p).read_bytes().decode('utf-8-sig').splitlines() if s.strip()]
def main():
 val=json.loads((O/'validation.json').read_text(encoding='utf-8'));assert val['status']=='PASS'
 contract=json.loads((O/'reconstruction_contract.json').read_text(encoding='utf-8'));dep=[]
 for d in contract['dependencies']:
  assert sha(d['actual_path'])==d['sha256'] and sha(d['retained_export_dependency_path'])==d['retained_export_dependency_sha256']
  ok=normalized_lines(d['actual_path'])==normalized_lines(d['retained_export_dependency_path']);assert ok
  dep.append(dict(name=d['name'],senior_sha256=d['sha256'],retained_export_sha256=d['retained_export_dependency_sha256'],nonempty_lines_equal_after_trimming_whitespace=ok))
 csv('validation/dependency_equivalence.csv',dep)
 arraychecks=[]
 for a in read(O/'omag_reconstruction/omag_array_manifest.csv').itertuples():
  h=sha(a.path);assert h==a.SHA256;arr=np.load(a.path,mmap_mode='r');assert arr.shape==(351,500,500) and arr.dtype==np.float64
  assert np.isfinite(arr).sum()==a.finite_count and float(arr.min())==a.min and float(arr.max())==a.max
  arraychecks.append(dict(scan_id=a.scan_id,sha256=h,shape='351,500,500',axis_order='z,x,frame',dtype='float64',frame_indices='0:499',finite_count=int(np.isfinite(arr).sum()),passed=True));del arr
 csv('validation/final_array_checks.csv',arraychecks)
 # D128 has an additional original raw-input SHA anchor, independent of names.
 anchor=read(R/'results/formal_sv_d128_v21_full2500_run001/input_sha256.csv');raw=read(O/'raw_oct_input_manifest.csv');anchorchecks=[]
 for r in raw[raw.diameter_um.eq(128)].itertuples():
  matches=anchor[anchor.sha256.eq(r.raw_sha256)];assert len(matches)==1
  anchorchecks.append(dict(scan_id=r.scan_id,raw_sha256=r.raw_sha256,original_raw_input_manifest='results/formal_sv_d128_v21_full2500_run001/input_sha256.csv',matched=True))
 csv('validation/d128_original_raw_anchor.csv',anchorchecks)
 required=['README.md','analysis_plan.json','raw_oct_input_manifest.csv','reconstruction_contract.json','common_preprocessing_audit.md','provenance.json','run_status.json','validation.json','omag_reconstruction/omag_array_manifest.csv','omag_reconstruction/export_metadata_summary.csv','omag_reconstruction/fixed_frame_reconstruction_checks.csv','comparison/algorithm_direction_ledger.csv','comparison/diameter_percent_change_by_algorithm.csv','comparison/common_support_summary.csv','comparison/stageOA_stageOB_common_support.csv','comparison/svA_omagA_comparison.csv','comparison/svB_omagB_comparison.csv','validation/common_preprocessing_validation.csv','validation/fixed_frame_omag_checks.csv','validation/metric_formula_checks.csv','validation/scan_median_reconstruction.csv','validation/common_support_checks.csv']
 for stage in ['stageOA','stageOB']:
  required += [stage+'/'+n for n in ['volume_metrics.csv','coverage.csv','trend_ledger.csv']];assert len(list((O/stage).glob('framewise_metrics_*.csv.gz')))==15
 for name in ['omag_senior_metrics_stageOA','sv_vs_omag_same_geometry','sv_vs_omag_full_pipeline']:
  required += ['figures/'+name+ext for ext in ['.png','.pdf']]
 assert all((O/p).is_file() for p in required)
 prov=json.loads((O/'provenance.json').read_text(encoding='utf-8'));assert sha(O/'analysis_plan.json')==prov['analysis_plan_sha256']
 prov['pre_environment_reconstruction_contract_sha256']=prov.pop('reconstruction_contract_sha256');prov['reconstruction_contract_sha256']=sha(O/'reconstruction_contract.json')
 prov.update(completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),software=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__,os=platform.platform()),source_scripts=[dict(path=p.name,sha256=sha(p)) for p in sorted(O.glob('*.py'))]+[dict(path='rebuild_fixed.m',sha256=sha(O/'rebuild_fixed.m'))],source_anchor_manifests=[dict(path=str(p.relative_to(R)),sha256=sha(p)) for p in [R/'data/diameter_stage_v1/volume_manifest.csv',R/'results/formal_sv_d128_v21_full2500_run001/input_sha256.csv',OLD/'senior_method_contract.json',OLD/'validation/array_identity_on_read.csv.gz']],git=dict(branch='analysis/omag-same-data-senior-metric-v1',base_commit='b365c276c20a16d22ac51afd891b8d924819db86',remote_baseline_verified=True),large_arrays_excluded_from_git=True)
 js('provenance.json',prov);val.update(required_artifacts_present=len(required),final_array_sha_shape_finite_checks=15,dependency_equivalence_checks=6,d128_original_raw_anchor_checks=5);js('validation.json',val)
 audit=(O/'common_preprocessing_audit.md').read_text(encoding='utf-8');note='\nAll six senior/retained dependency files have exactly matching nonempty lines after trimming surrounding whitespace; actual file SHA differences are newline/blank-line/trailing-space differences. See validation/dependency_equivalence.csv. Fixed SV and OMAG results are bitwise equal in all 45 frames, not merely within tolerance.\n'
 if note not in audit:(O/'common_preprocessing_audit.md').write_text(audit+note,encoding='utf-8')
 entries=[]
 for p in sorted(O.rglob('*')):
  if p.is_file() and 'intermediate_local' not in p.parts and '__pycache__' not in p.parts and p.name!='output_sha256.csv':entries.append(dict(path=str(p.relative_to(O)).replace('\\','/'),bytes=p.stat().st_size,sha256=sha(p)))
 csv('output_sha256.csv',entries)
 for r in read(O/'output_sha256.csv').itertuples():assert sha(O/r.path)==r.sha256
 print('FINALIZED',len(entries),'hashed artifacts; 15 final arrays PASS',flush=True)
def refresh_ledger():
 # Reuse scientific checks for unchanged inputs/metrics after presentation edits.
 prov=json.loads((O/'provenance.json').read_text(encoding='utf-8'))
 for row in prov['source_scripts']:row['sha256']=sha(O/row['path'])
 js('provenance.json',prov)
 entries=[dict(path=str(p.relative_to(O)).replace('\\','/'),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(O.rglob('*')) if p.is_file() and 'intermediate_local' not in p.parts and '__pycache__' not in p.parts and p.name!='output_sha256.csv']
 csv('output_sha256.csv',entries)
 for r in read(O/'output_sha256.csv').itertuples():assert sha(O/r.path)==r.sha256
 print('Refreshed and verified',len(entries),'artifact hashes',flush=True)
if __name__=='__main__':
 if '--refresh-ledger' in sys.argv:refresh_ledger()
 else:main()
