#!/usr/bin/env python3
"""Observed/apparent tail intensity from frozen run001; no background subtraction."""
from __future__ import annotations
import argparse, hashlib, importlib.util, io, json, sys, tempfile, zipfile
from pathlib import Path
import numpy as np
import pandas as pd

EXPECTED_VALID=2422; EXPECTED_PACKAGES=25
RELEASE_TAG='formal-sv-d128-v21-run001'
RELEASE_BASE=f'https://github.com/bulbel-magnolia/OCTA-/releases/download/{RELEASE_TAG}'
ANCHORS=np.arange(0.,500.0001,10.)

def load_helper(path:Path):
    spec=importlib.util.spec_from_file_location('fixed128_helper',path); mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod); return mod

def wint(img,w,pix):
    m=w>0
    return float(np.sum(img[m]*w[m],dtype=np.float64)*pix) if np.any(m) and np.isfinite(img[m]).all() else np.nan

def rowmean(img,w):
    f=np.isfinite(img); den=(f*w[None,:]).sum(1,dtype=np.float64); out=np.full(img.shape[0],np.nan)
    np.divide((np.where(f,img,0.)*w[None,:]).sum(1,dtype=np.float64),den,out=out,where=den>0); return out

def qs(x):
    a=pd.to_numeric(pd.Series(x),errors='coerce').to_numpy(float); a=a[np.isfinite(a)]
    if not a.size:return dict(n=0,mean=np.nan,sd=np.nan,min=np.nan,q1=np.nan,median=np.nan,q3=np.nan,iqr=np.nan,max=np.nan)
    q1,med,q3=np.quantile(a,[.25,.5,.75]); return dict(n=int(a.size),mean=float(a.mean()),sd=float(a.std(ddof=1)) if a.size>1 else np.nan,min=float(a.min()),q1=float(q1),median=float(med),q3=float(q3),iqr=float(q3-q1),max=float(a.max()))

def frame_summary(df):
    metrics=['q_vessel_observed','source_mean_observed','q_tail_observed','tail_mean_observed','ri_tail_observed','ri_tail_corrected','source_background_fraction','tail_background_fraction']
    rows=[]
    for scan,g in df.groupby('scan_id',sort=False):
        base=dict(scan_id=scan,flow_mm_s=float(g.flow_mm_s.iloc[0]),n_valid_frames=len(g))
        for m in metrics: rows.append(base|dict(metric=m)|qs(g[m]))
    return pd.DataFrame(rows).sort_values(['flow_mm_s','metric'])

def depth_summary(df):
    rows=[]
    for (scan,t),g in df.groupby(['scan_id','target_r_um'],sort=False):
        o,c,v,b=qs(g.ri_r_observed),qs(g.ri_r_corrected),qs(g.V_observed),qs(g.B_frozen)
        rows.append(dict(scan_id=scan,flow_mm_s=float(g.flow_mm_s.iloc[0]),target_r_um=float(t),n_frames=len(g),selected_r_um_median=float(np.median(g.selected_r_um)),max_abs_sampling_error_um=float(g.abs_sampling_error_um.max()),ri_r_observed_q1=o['q1'],ri_r_observed_median=o['median'],ri_r_observed_q3=o['q3'],ri_r_corrected_q1=c['q1'],ri_r_corrected_median=c['median'],ri_r_corrected_q3=c['q3'],V_observed_median=v['median'],B_frozen_median=b['median']))
    return pd.DataFrame(rows).sort_values(['flow_mm_s','target_r_um'])

def write_readme(out,ss,val):
    med=ss[ss.metric.eq('ri_tail_observed')]
    lines=['# Observed SV tail intensity — no background subtraction','',
    'Question: **How bright does the tail appear in the SV image relative to the vessel itself?**','',
    '`RI_tail_observed = (Q_T_observed/A_T)/(Q_V_observed/A_V)` and `RI_observed(r)=V(r)/source_mean_observed`. Both source and tail use linear `sv_raw` without background subtraction. The previous background-corrected endpoint is retained separately.','',
    'Frozen X4+X1 geometry, v2.1 z geometry, 500 µm tail, guard=0, fractional weights and `valid == True` are unchanged. No p-values are computed; B-scans are spatial positions, not independent biological replicates.','',
    '## Validation','',f"- Valid frames: **{val['valid_frames']}**",f"- Release packages: **{val['release_packages_verified']}/{EXPECTED_PACKAGES}**",f"- NPZ SHA checks: **{val['valid_frame_npz_sha256_verified']}**",f"- Frozen corrected-metric replay max relative error: **{val['baseline_replay_max_relative_error']:.3e}**",f"- Area replay max relative error: **{val['area_replay_max_relative_error']:.3e}**",f"- Max depth-anchor sampling error: **{val['max_depth_anchor_sampling_error_um']:.3f} µm**",'',
    '## Scan-level observed RI_tail','', '| scan | flow (mm/s) | median | Q1 | Q3 |','|---|---:|---:|---:|---:|']
    for r in med.sort_values('flow_mm_s').itertuples(): lines.append(f'| {r.scan_id} | {r.flow_mm_s:g} | {r.median:.6g} | {r.q1:.6g} | {r.q3:.6g} |')
    lines += ['','`RI_tail_observed=0.20` means the average signal actually visible in the tail rectangle is 20% of the average signal actually visible in the vessel ROI. It intentionally includes local SV background.','',
    'Files: `observed_tail_intensity_framewise.csv`, `observed_depth_anchor_framewise.csv.gz`, `observed_tail_intensity_scan_summary.csv`, `observed_depth_summary.csv`, plus validation/provenance/SHA audits.']
    (out/'README.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--root',default='.'); ap.add_argument('--formal-dir',default='results/formal_sv_d128_v21_full2500_run001'); ap.add_argument('--task1-dir',default='analysis/formal_sv_d128_v21_run001/relative_tail_intensity_full2422'); ap.add_argument('--output-dir',default='analysis/formal_sv_d128_v21_run001/observed_tail_intensity_full2422'); ap.add_argument('--frozen-source-sha',default='acef5eb9d5ac356f1acf11aee885895da79a22e3'); ap.add_argument('--workflow-source-sha',default=''); a=ap.parse_args()
    root=Path(a.root).resolve(); formal=(root/a.formal_dir).resolve(); task1=(root/a.task1_dir).resolve(); out=(root/a.output_dir).resolve(); out.mkdir(parents=True,exist_ok=True)
    helper_path=root/'analysis/formal_sv_d128_v21_run001/fixed128_width_sensitivity_full2422/audit_fixed128_width.py'; H=load_helper(helper_path)
    sys.path.insert(0,str(root/'src')); from svrecttail.geometry import VesselGeometry,ellipse_weights,interval_overlap_weights; from svrecttail.quantification import quantify_frame
    base,_,cfg=H.prepare_inputs(root,formal,task1); packages=pd.read_csv(formal/'download_packages.csv'); arrays=pd.read_csv(formal/'arrays_sha256.csv')
    if len(base)!=EXPECTED_VALID or len(packages)!=EXPECTED_PACKAGES: raise AssertionError('Frozen input coverage mismatch')
    arrays['scan_id']=arrays.scan_id.astype(str); arrays['frame_index_0based']=pd.to_numeric(arrays.frame_index_0based).astype(int); ah={(r.scan_id,int(r.frame_index_0based)):str(r.sha256) for r in arrays.itertuples()}
    recs=[]; deps=[]; pkgs=[]; verified=0; replay_max=0.; area_max=0.; anchor_max=0.; pix=cfg['dx_um']*cfg['dz_um']
    with tempfile.TemporaryDirectory(prefix='sv-observed-') as td:
      td=Path(td)
      for pkg in packages.itertuples(index=False):
        name=str(pkg.file); scan,lo,hi=H.parse_package_name(name); path=td/name; H.download_file(f'{RELEASE_BASE}/{name}',path); zsha=H.sha256_file(path)
        if zsha!=str(pkg.sha256) or path.stat().st_size!=int(pkg.bytes): raise AssertionError(f'Package verification failed {name}')
        wanted=base[(base.scan_id.eq(scan))&base.frame_index.between(lo,hi)]; vin=0
        with zipfile.ZipFile(path) as zf:
          names=set(zf.namelist())
          for row in wanted.itertuples(index=False):
            fi=int(row.frame_index); member=f'arrays/{scan}/frame_{fi:03d}.npz'
            if member not in names: raise FileNotFoundError(member)
            data=zf.read(member)
            if H.sha256_bytes(data)!=ah[(scan,fi)]: raise AssertionError(f'NPZ SHA mismatch {scan}/{fi}')
            verified+=1; vin+=1
            with np.load(io.BytesIO(data),allow_pickle=False) as d: sv=np.asarray(d['sv_raw'],dtype=np.float64)
            geom=VesselGeometry(x_left_edge_px=float(row.baseline_x_left_edge_px),x_right_edge_px=float(row.baseline_x_right_edge_px),z_top_edge_px=float(row.z_top_edge_px),diameter_um=cfg['diameter_um'],dx_um=cfg['dx_um'],dz_um=cfg['dz_um'])
            replay=quantify_frame(sv,geom,tail_gap_um=cfg['tail_gap_um'],tail_length_um=cfg['tail_length_um'],background_skip_columns=cfg['skip_columns'],background_strip_width_columns=cfg['strip_width_columns'],background_excluded_side=None,ellipse_supersample=cfg['ellipse_supersample'],source_qc_valid=True)
            if not replay.valid: raise AssertionError(f'Baseline replay failed {scan}/{fi}: {replay.invalid_reason}')
            for n,o in [(replay.q_vessel,row.q_vessel),(replay.q_tail,row.q_tail),(replay.source_mean,row.source_mean),(replay.ratio_tail_to_vessel,row.ratio_tail_to_vessel)]: replay_max=max(replay_max,H.abs_rel_error(float(n),float(o)))
            if replay_max>1e-8: raise AssertionError(f'Baseline replay tolerance exceeded {replay_max}')
            sw=ellipse_weights(sv.shape,geom,supersample=cfg['ellipse_supersample']); xw=interval_overlap_weights(sv.shape[1],geom.x_left_edge_px,geom.x_right_edge_px); top=geom.z_bottom_edge_px+cfg['tail_gap_um']/cfg['dz_um']; zw=interval_overlap_weights(sv.shape[0],top,top+cfg['tail_length_um']/cfg['dz_um']); tw=np.multiply.outer(zw,xw)
            av=float(sw.sum()*pix); at=float(tw.sum()*pix); area_max=max(area_max,H.abs_rel_error(av,float(row.source_area_um2)),H.abs_rel_error(at,float(row.tail_area_um2)))
            if area_max>1e-10: raise AssertionError(f'Area replay tolerance exceeded {area_max}')
            qv=wint(sv,sw,pix); qt=wint(sv,tw,pix); svmean=H.safe_divide(qv,av); stmean=H.safe_divide(qt,at); ri=H.safe_divide(stmean,svmean); sb=qv-float(row.q_vessel); tb=qt-float(row.q_tail)
            recs.append(dict(scan_id=scan,frame_index=fi,frame_index_0based=fi,flow_mm_s=float(row.flow_mm_s),localization_source=row.localization_source,valid=True,q_vessel_observed=qv,q_tail_observed=qt,source_area_um2=av,tail_area_um2=at,source_mean_observed=svmean,tail_mean_observed=stmean,ri_tail_observed=ri,q_vessel_corrected=float(row.q_vessel),q_tail_corrected=float(row.q_tail),source_mean_corrected=float(row.source_mean),tail_mean_corrected=float(row.tail_mean),ri_tail_corrected=float(row.ri_tail),source_background_component=sb,tail_background_component=tb,source_background_fraction=H.safe_divide(sb,qv),tail_background_fraction=H.safe_divide(tb,qt)))
            vp=rowmean(sv,xw); rr=(np.arange(sv.shape[0])-geom.z_bottom_edge_px)*cfg['dz_um']; cand=np.flatnonzero((zw>0)&np.isfinite(vp))
            for target in ANCHORS:
              zi=int(cand[np.argmin(np.abs(rr[cand]-target))]); err=float(abs(rr[zi]-target)); anchor_max=max(anchor_max,err)
              if err>cfg['dz_um']/2+1e-9: raise AssertionError(f'Anchor tolerance exceeded {scan}/{fi}/{target}: {err}')
              deps.append(dict(scan_id=scan,frame_index=fi,flow_mm_s=float(row.flow_mm_s),target_r_um=float(target),selected_z_index_0based=zi,selected_r_um=float(rr[zi]),abs_sampling_error_um=err,tail_z_fraction=float(zw[zi]),V_observed=float(vp[zi]),B_frozen=float(replay.background.combined[zi]),T_corrected=float(replay.tail_contrast_profile[zi]),source_mean_observed=svmean,source_mean_corrected=float(row.source_mean),ri_r_observed=H.safe_divide(vp[zi],svmean),ri_r_corrected=H.safe_divide(replay.tail_contrast_profile[zi],float(row.source_mean))))
        pkgs.append(dict(package=name,scan_id=scan,frame_lo=lo,frame_hi=hi,expected_bytes=int(pkg.bytes),observed_bytes=path.stat().st_size,expected_sha256=str(pkg.sha256),observed_sha256=zsha,package_verified=True,valid_frames_in_package=len(wanted),valid_frame_npz_sha256_verified=vin)); path.unlink(); print(f'PROCESSED {name}: valid={len(wanted)}',flush=True)
    fw=pd.DataFrame(recs).sort_values(['flow_mm_s','frame_index']); dp=pd.DataFrame(deps).sort_values(['flow_mm_s','frame_index','target_r_um'])
    if len(fw)!=EXPECTED_VALID or verified!=EXPECTED_VALID or len(dp)!=EXPECTED_VALID*len(ANCHORS): raise AssertionError('Output coverage mismatch')
    if not np.isfinite(fw.ri_tail_observed).all() or not np.isfinite(dp.ri_r_observed).all(): raise AssertionError('Nonfinite observed relative intensity')
    ss=frame_summary(fw); ds=depth_summary(dp); fw.to_csv(out/'observed_tail_intensity_framewise.csv',index=False,float_format='%.17g'); dp.to_csv(out/'observed_depth_anchor_framewise.csv.gz',index=False,float_format='%.17g',compression='gzip'); ss.to_csv(out/'observed_tail_intensity_scan_summary.csv',index=False,float_format='%.17g'); ds.to_csv(out/'observed_depth_summary.csv',index=False,float_format='%.17g'); pd.DataFrame(pkgs).to_csv(out/'release_package_audit.csv',index=False)
    val=dict(valid_frames=len(fw),depth_anchor_rows=len(dp),depth_anchor_count_per_frame=len(ANCHORS),release_packages_verified=sum(bool(x['package_verified']) for x in pkgs),valid_frame_npz_sha256_verified=verified,baseline_replay_max_relative_error=float(replay_max),area_replay_max_relative_error=float(area_max),max_depth_anchor_sampling_error_um=float(anchor_max),observed_metrics_background_subtracted=False,source_and_tail_same_signal_scale=True,p_values_computed=False,bscans_treated_as_independent_biological_replicates=False,parameter_selection_by_flow_pattern=False); (out/'validation.json').write_text(json.dumps(val,indent=2)+'\n')
    paths=[formal/'frame_results.csv',formal/'download_packages.csv',formal/'arrays_sha256.csv',formal/'run_config.json',task1/'relative_tail_intensity_framewise.csv',helper_path,root/'src/svrecttail/geometry.py',root/'src/svrecttail/quantification.py']; pd.DataFrame([dict(input_path=p.relative_to(root).as_posix(),sha256=H.sha256_file(p),size_bytes=p.stat().st_size) for p in paths]).to_csv(out/'input_sha256.csv',index=False)
    prov=dict(analysis='Observed/apparent SV tail intensity without background subtraction',scientific_question='How bright does the tail appear in the SV image relative to the vessel itself?',frozen_input_directory=a.formal_dir,output_directory=a.output_dir,formal_release_tag=RELEASE_TAG,frozen_formal_source_sha=a.frozen_source_sha,workflow_source_sha=a.workflow_source_sha,entry_rule='valid == True only',geometry='frozen X4 centre + X1 apparent width + v2.1 z geometry',sv_signal='linear sv_raw = Var_t(|E|), N denominator',definitions={'RI_tail_observed':'(Q_T_observed/A_T)/(Q_V_observed/A_V)','RI_observed_r':'V(r)/source_mean_observed'},background_subtraction='none for observed source or tail',depth_summary='0..500 um every 10 um; nearest native axial sample; no interpolation or fit',statistical_scope='scan-volume descriptive distributions; no p-values'); (out/'provenance.json').write_text(json.dumps(prov,indent=2)+'\n'); write_readme(out,ss,val)
    key=ss[ss.metric.eq('ri_tail_observed')][['scan_id','flow_mm_s','median','q1','q3']]; print(json.dumps(dict(valid_frames=len(fw),depth_anchor_rows=len(dp),baseline_replay_max_relative_error=replay_max,max_depth_anchor_sampling_error_um=anchor_max,ri_tail_observed=key.to_dict('records')),indent=2))
if __name__=='__main__': main()
