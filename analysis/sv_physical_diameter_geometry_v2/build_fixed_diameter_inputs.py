"""User-approved fixed physical diameter geometry, with explicit post-hoc D500 exclusion."""
from pathlib import Path
import sys,io,json,hashlib,zipfile,subprocess,platform
sys.dont_write_bytecode=True
import numpy as np
import pandas as pd
import scipy
from scipy.io import loadmat
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[1]
DEP=ROOT/'analysis/sv_source_tail_depth_coupling_v1';PREV=ROOT/'analysis/sv_source_tail_framewise_coupling_v1'
sys.path.insert(0,str(DEP))
from derive_depth_metrics import source_stats,tail_stats,source_local_weights,VesselGeometry,DX,DZ
from svrecttail.geometry import ellipse_weights,ellipse_is_complete,interval_is_complete
KEEP=[128,235,285];FLOWS=[1,3,5,7,10];EXCLUDE=[500]
def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def csv(n,x):x.to_csv(OUT/n,index=False,float_format='%.17g',lineterminator='\n',compression={'method':'gzip','mtime':0} if n.endswith('.gz') else None)
def js(n,x):(OUT/n).write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
def read(p):return pd.read_csv(p,float_precision='round_trip',low_memory=False)
def geometry(x4,z_top,diameter):
    return VesselGeometry(x4-diameter/(2*DX),x4+diameter/(2*DX),z_top,float(diameter),DX,DZ)
def main():
    start=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    refpaths=[PREV/'framewise_source_tail_metrics.csv',PREV/'validation.json',DEP/'derive_depth_metrics.py',ROOT/'src/svrecttail/geometry.py',DEP/'array_identity_audit.csv.gz',ROOT/'analysis/sv_data_quality_qc_v1/volume_qc_summary.csv']
    frozen={str(p.relative_to(ROOT)):sha(p) for p in refpaths}
    fw=read(refpaths[0]);assert sha(refpaths[0])==json.loads(refpaths[1].read_text())['framewise_sha256']
    source=fw[fw.diameter_um.isin(KEEP)&fw.flow_mm_s.isin(FLOWS)].copy()
    ids=read(DEP/'array_identity_audit.csv.gz').fillna('').set_index(['scan_id','frame_index'])
    assert source.scan_id.nunique()==15
    cohorts=fw.groupby(['scan_id','diameter_um','flow_mm_s','grid_scope']).size().reset_index(name='prior_valid_frames')
    cohorts['new_analysis_role']=np.where(cohorts.diameter_um.isin(KEEP),'primary_three_diameters','excluded_from_future_scientific_analysis')
    cohorts['nominal_frames']=500
    cohorts['decision_basis']=np.where(cohorts.diameter_um.eq(500),'user decision after image/QC review; upper signal concentration and lower deep detectability','retain common flows')
    cohorts['selection_timing']='post_hoc_after_previous_results_review'
    csv('cohort_manifest.csv',cohorts);csv('d500_exclusion_record.csv',cohorts[cohorts.diameter_um.eq(500)])
    js('analysis_policy.json',dict(version='physical_diameter_geometry_v2',decision_date='2026-09-10',
        user_geometry_approval='好的，我们就改进为这种方案吧',
        user_d500_approval='D500 不再进入后续科学分析，但保留数据及排除记录',
        primary_diameters_um=KEEP,common_flows_mm_s=FLOWS,excluded_diameters_um=EXCLUDE,excluded_volumes=8,
        exclusion_timing='post-hoc after image/QC inspection; not preregistered and not instrument failure certification',
        old_data_and_results='retain; do not overwrite or delete',
        source='physical D x D ellipse in pixel coordinates; same frozen X4 and z_top',
        tail='physical D lateral width; begins z_top+D/dz; guard 0; no cone; raw SV',
        original_X1='retain as original_apparent_width_um only; never relabel apparent width as physical D',
        invalid_geometry='preserve original missing frames; fixed D does not recover missing top/center',
        new_coupling='not computed in this geometry/input migration; old coupling not relabeled as new',
        old_qc='belongs to old X1 geometry; not relabeled as QC for new ROI'))
    manifest=[dict(input_path=p,sha256=h,role='frozen reference') for p,h in frozen.items()]
    # Every D128 package actually used is checked against its release manifest.
    prior=read(DEP/'input_manifest.csv');packages=set(ids.loc[list(zip(source.scan_id,source.frame_index))].query("member != ''").input_path)
    expected_zip={Path(r.input_path).name:r.sha256 for r in prior.itertuples() if str(r.input_path).endswith('.zip')}
    for s in sorted(packages):
        p=Path(s);p=p if p.is_absolute() else ROOT/p;h=sha(p);assert h==expected_zip[p.name]
        manifest.append(dict(input_path=str(p.resolve()),sha256=h,role='formal release ZIP'))
    errors={};out=[];arrayids=[];spots=[];summaries=[];geometry_rows=[]
    def checkerror(name,a,b,tol=1e-10):
        err=float(np.max(abs(np.asarray(a)-np.asarray(b))));rel=float(np.max(abs(np.asarray(a)-np.asarray(b))/np.maximum(abs(np.asarray(b)),np.finfo(float).tiny)))
        z=errors.setdefault(name,dict(max_absolute_error=0.,max_relative_error=0.));z['max_absolute_error']=max(err,z['max_absolute_error']);z['max_relative_error']=max(rel,z['max_relative_error']);assert rel<tol,name
    for scan,group in source.groupby('scan_id',sort=True):
        print('New geometry integration '+scan+f' ({len(group)} valid frames)',flush=True);cache={};volume=[]
        checkframes={int(group.frame_index.min()),int(min(group.frame_index,key=lambda i:(abs(i-249),i))),int(group.frame_index.max())}
        for r in group.itertuples():
            ir=ids.loc[(scan,r.frame_index)];p=Path(ir.input_path);p=p if p.is_absolute() else ROOT/p
            if ir.member:
                if str(p) not in cache:cache[str(p)]=zipfile.ZipFile(p)
                blob=cache[str(p)].read(ir.member)
                with np.load(io.BytesIO(blob),allow_pickle=False) as z:sv=np.asarray(z['sv_raw'],float)
            else:
                blob=p.read_bytes();sv=np.asarray(loadmat(io.BytesIO(blob),variable_names=['sv_raw'])['sv_raw'],float)
            h=hashlib.sha256(blob).hexdigest();assert h==ir.sha256==r.input_sha256
            assert sv.shape==(351,500) and np.isfinite(sv).all() and not (sv<0).any()
            arrayids.append(dict(scan_id=scan,frame_index=r.frame_index,input_path=str(p.resolve()),member=ir.member,sha256=h,identity_passed=True))
            if not ir.member:manifest.append(dict(input_path=str(p.resolve()),sha256=h,role='retained formal raw SV MAT'))
            old=VesselGeometry(r.x_left_edge_px,r.x_right_edge_px,r.z_top,float(r.diameter_um),DX,DZ)
            # Replay old scalar identity using exactly the same raw array, then integrate new geometry.
            oa,oq,om=source_stats(sv,old,6)
            checkerror('old_source_mean_replay',oq.sum()/oa.sum(),r.source_mean_raw)
            checkerror('old_source_q_replay',oq.sum(),r.source_q_raw)
            for dep in [100,500]:
                aa,qq,mm=tail_stats(sv,old,[0,dep]);checkerror(f'old_tail{dep}_mean_replay',mm[0],getattr(r,f'tail_mean_raw_{dep}um'));checkerror(f'old_tail{dep}_q_replay',qq[0],getattr(r,f'tail_q_raw_{dep}um'))
            g=geometry(r.X4,r.z_top,r.diameter_um)
            assert abs(g.lateral_width_um-r.diameter_um)<1e-9 and abs(g.x_center_px-r.X4)<1e-12 and g.z_top_edge_px==r.z_top
            assert abs(g.z_bottom_edge_px-r.z_bottom_edge_px)<1e-12
            assert ellipse_is_complete(sv.shape,g) and interval_is_complete(sv.shape[0],g.z_bottom_edge_px,g.z_bottom_edge_px+500/DZ)
            a,q,m=source_stats(sv,g,6);ta,tq,tm=tail_stats(sv,g,np.arange(21)*25)
            rec=dict(scan_id=scan,diameter_um=r.diameter_um,flow_mm_s=r.flow_mm_s,frame_index=r.frame_index,valid_geometry=True,
                analysis_role='primary_three_diameters',geometry_version='physical_diameter_geometry_v2',X4=r.X4,z_top=r.z_top,
                z_bottom_edge_px=g.z_bottom_edge_px,x_left_edge_px=g.x_left_edge_px,x_right_edge_px=g.x_right_edge_px,
                roi_width_um=r.diameter_um,roi_height_um=r.diameter_um,original_apparent_width_um=r.X1,
                source_area_um2=a.sum(),source_q_raw=q.sum(),source_mean_raw=q.sum()/a.sum(),input_sha256=h,
                old_source_mean_raw=r.source_mean_raw,old_source_q_raw=r.source_q_raw)
            for j in range(6):
                for key,val in [('area_um2',a[j]),('q_raw',q[j]),('mean_raw',m[j])]:rec[f'source_s{j+1}_{key}']=val
            for j in range(20):
                for key,val in [('area_um2',ta[j]),('q_raw',tq[j]),('mean_raw',tm[j])]:rec[f'tail_t{j+1:02}_{key}']=val
            for dep in [100,500]:
                aa,qq,mm=tail_stats(sv,g,[0,dep]);checkerror(f'new_tail{dep}_bin_sum',tq[:dep//25].sum(),qq[0]);checkerror(f'new_tail{dep}_physical_area',aa[0],r.diameter_um*dep)
                rec.update({f'tail_area_um2_{dep}um':aa[0],f'tail_q_raw_{dep}um':qq[0],f'tail_mean_raw_{dep}um':mm[0],
                    f'old_tail_mean_raw_{dep}um':getattr(r,f'tail_mean_raw_{dep}um'),f'old_tail_q_raw_{dep}um':getattr(r,f'tail_q_raw_{dep}um')})
            # Fractional support is preserved; do not round the diameter to integer pixels.
            zi,xi,ww,full=source_local_weights(sv.shape,g,6)
            rec['source_effective_pixels']=full.sum()**2/np.sum(full**2)
            for j in range(6):rec[f'source_s{j+1}_effective_pixels']=ww[j].sum()**2/np.sum(ww[j]**2)
            if r.frame_index in checkframes:
                fullref=ellipse_weights(sv.shape,g,supersample=16);refq=(fullref*sv).sum()*DX*DZ;refarea=fullref.sum()*DX*DZ
                checkerror('new_source_vs_independent_full_image_Q',q.sum(),refq);checkerror('new_source_vs_independent_full_image_area',a.sum(),refarea)
                spots.append(dict(scan_id=scan,frame_index=r.frame_index,source_full_image_q=refq,source_local_q=q.sum(),passed=True))
            out.append(rec);volume.append(rec)
        for z in cache.values():z.close()
        data=pd.DataFrame(volume)
        for metric in ['source_mean_raw','source_q_raw','tail_mean_raw_100um','tail_q_raw_100um','tail_mean_raw_500um','tail_q_raw_500um']:
            rel=(data[metric]/data['old_'+metric]-1)*100
            summaries.append(dict(scan_id=scan,diameter_um=int(group.diameter_um.iloc[0]),flow_mm_s=float(group.flow_mm_s.iloc[0]),metric=metric,
                n_valid_frames=len(data),old_volume_median=data['old_'+metric].median(),new_volume_median=data[metric].median(),
                paired_frame_change_percent_median=rel.median(),paired_frame_abs_change_percent_median=rel.abs().median(),paired_frame_abs_change_percent_p95=rel.abs().quantile(.95)))
        nominal=data.set_index('frame_index').reindex(range(500)).rename_axis('frame_index').reset_index()
        nominal['valid_geometry']=nominal.valid_geometry.eq(True)
        for col,v in [('scan_id',scan),('diameter_um',int(group.diameter_um.iloc[0])),('flow_mm_s',float(group.flow_mm_s.iloc[0])),('analysis_role','primary_three_diameters'),('geometry_version','physical_diameter_geometry_v2')]:nominal[col]=v
        geometry_rows.append(nominal)
    allrows=pd.concat(geometry_rows,ignore_index=True);valid=allrows[allrows.valid_geometry]
    assert len(allrows)==7500 and len(valid)==len(source) and set(zip(valid.scan_id,valid.frame_index))==set(zip(source.scan_id,source.frame_index))
    for d,g in allrows.groupby('diameter_um'):csv(f'framewise_fixed_diameter_D{int(d)}.csv.gz',g)
    columns=['scan_id','diameter_um','flow_mm_s','frame_index','valid_geometry','analysis_role','geometry_version','X4','z_top','z_bottom_edge_px',
             'x_left_edge_px','x_right_edge_px','roi_width_um','roi_height_um','original_apparent_width_um']
    csv('geometry_manifest.csv.gz',allrows[columns]);csv('volume_geometry_change_summary.csv',pd.DataFrame(summaries))
    csv('array_identity_audit.csv.gz',pd.DataFrame(arrayids));csv('input_manifest.csv',pd.DataFrame(manifest).drop_duplicates('input_path'))
    csv('independent_integration_spotchecks.csv',pd.DataFrame(spots))
    unchanged=all(sha(ROOT/p)==h for p,h in frozen.items());assert unchanged
    coverage=allrows.groupby('diameter_um').agg(volumes=('scan_id','nunique'),nominal_frames=('frame_index','size'),valid_frames=('valid_geometry','sum')).reset_index()
    csv('coverage.csv',coverage)
    js('validation.json',dict(status='passed',checks=errors,primary_volumes=15,primary_valid_frames=len(valid),primary_nominal_frames=7500,
       original_invalid_frames_preserved=7500-len(valid),excluded_d500_volumes=8,excluded_d500_valid_frames=int(cohorts[cohorts.diameter_um.eq(500)].prior_valid_frames.sum()),
       fixed_physical_dimensions_passed=True,frozen_centers_top_bottom_preserved=True,no_added_invalid_frames=True,
       raw_array_sha_all_passed=True,raw_finite_nonnegative_all=True,independent_full_image_checks=len(spots),frozen_references_unchanged=unchanged,
       formal_coupling_rerun=False,matched_background_qc_rerun=False,exclusion_is_posthoc=True))
    js('provenance.json',dict(starting_head=start,decision_date='2026-09-10',input_manifest_sha256=sha(OUT/'input_manifest.csv'),
       script_sha256=sha(__file__),source_reference_sha256=frozen,software=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__),
       signal='var(abs(IMG),1,3), denominator N; linear raw SV; no background subtraction',
       source_supersample=16,dx_um=DX,dz_um=DZ,source_normalized_bins=6,tail_band_edges_um=list(range(0,501,25)),
       scope='approved geometry migration and raw integration inputs only; new coupling and new matched-background QC are separate downstream analyses'))
    print(coverage.to_string(index=False),flush=True);print('Geometry migration and new raw integration inputs passed',flush=True)
if __name__=='__main__':main()
