#!/usr/bin/env python3
"""Promote no-background raw-SV relative tail intensity to the primary SV endpoint.

Primary scalar:
    RI_tail = mean(raw SV in tail ROI) / mean(raw SV in vessel ROI)
Primary depth curve:
    RI(r) = mean(raw SV across tail lateral span at native depth r) / mean(raw SV in vessel ROI)

Background is not subtracted anywhere in the primary endpoint. Historical background-corrected
outputs remain untouched. This script replays all 2422 valid frames from released 2-D sv_raw
arrays, audits fixed-128 lateral-width geometry, then recomputes flow summaries, spatial ACF,
and exhaustive systematic subsampling stability for the new primary endpoint.
"""
from __future__ import annotations

import argparse, hashlib, importlib.util, io, json, math, tempfile, zipfile
from pathlib import Path
import numpy as np
import pandas as pd

EXPECTED_COUNTS={"flow01":486,"flow03":468,"flow05":491,"flow07":493,"flow10":484}
EXPECTED_VALID=sum(EXPECTED_COUNTS.values()); EXPECTED_PACKAGES=25
RELEASE_TAG="formal-sv-d128-v21-run001"
RELEASE_BASE=f"https://github.com/bulbel-magnolia/OCTA-/releases/download/{RELEASE_TAG}"
STRIDES=(2,5,10,20,25,50); SELECTED_DEPTHS=(0,50,100,200,300,400,500)


def load_module(name,path):
    spec=importlib.util.spec_from_file_location(name,path); mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod); return mod

def sha256(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''): h.update(b)
    return h.hexdigest()

def safe_div(a,b):
    a=float(a); b=float(b)
    return a/b if np.isfinite(a) and np.isfinite(b) and b!=0 else np.nan

def weighted_integral(img,w,pixel_area):
    m=w>0
    if not np.any(m) or not np.isfinite(img[m]).all(): return np.nan
    return float(np.sum(img[m]*w[m],dtype=np.float64)*pixel_area)

def pct_delta(new,old): return 100.0*safe_div(float(new)-float(old),old)

def quantiles(x):
    a=pd.to_numeric(pd.Series(x),errors='coerce').to_numpy(float); a=a[np.isfinite(a)]
    if not len(a): return dict(n=0,q1=np.nan,median=np.nan,q3=np.nan,p95=np.nan,max=np.nan,mean=np.nan,sd=np.nan)
    q1,med,q3,p95=np.quantile(a,[.25,.5,.75,.95])
    return dict(n=len(a),q1=float(q1),median=float(med),q3=float(q3),p95=float(p95),max=float(np.max(a)),mean=float(np.mean(a)),sd=float(np.std(a,ddof=1)) if len(a)>1 else np.nan)

def threshold_cross(lags,acf,t):
    for i in range(1,len(lags)):
        a0,a1=acf[i-1],acf[i]
        if np.isfinite(a0) and np.isfinite(a1) and a0>t and a1<=t:
            return float(lags[i-1]+(a0-t)/(a0-a1)*(lags[i]-lags[i-1])) if a0!=a1 else float(lags[i])
    return np.nan

def acf_pairwise(series,max_lag=200):
    vals=series[np.isfinite(series)]; mu=float(vals.mean()); var=float(np.mean((vals-mu)**2))
    if len(vals)<3 or var<=0: raise RuntimeError('invalid ACF input')
    rows=[]; out=[]
    for lag in range(max_lag+1):
        l=series if lag==0 else series[:-lag]; r=series if lag==0 else series[lag:]
        m=np.isfinite(l)&np.isfinite(r); val=float(np.mean((l[m]-mu)*(r[m]-mu))/var) if m.any() else np.nan
        rows.append((lag,int(m.sum()),val)); out.append(val)
    return np.asarray(out,float),rows

def initial_positive_tau(acf):
    vals=[]; stop=np.nan
    for lag in range(1,len(acf)):
        if not np.isfinite(acf[lag]) or acf[lag]<=0: stop=float(lag); break
        vals.append(float(acf[lag]))
    return float(1+2*sum(vals)),stop

def write_readme(out,flow,geom,acf_cross,sub_cross,depth_selected,val):
    ri=flow[['scan_id','flow_mm_s','n_valid_frames','ri_tail_median','ri_tail_q1','ri_tail_q3','ri_delta_pct_vs_flow01','source_mean_raw_median','tail_mean_raw_median']]
    g=geom.set_index('metric')
    a=acf_cross.iloc[0]
    lines=[
      '# No-background Relative Tail Intensity — primary SV reanalysis v1','',
      'From this analysis forward, **RI_tail** means the no-background raw-SV ratio:','',
      '`RI_tail = mean(raw SV in tail ROI) / mean(raw SV in vessel ROI)`','',
      '`RI(r) = mean(raw SV across the vessel lateral span at depth r) / mean(raw SV in vessel ROI)`','',
      'No local background is subtracted. Historical background-corrected outputs remain in their old directories only and are not used by this primary analysis.','',
      '## Five-flow primary result','',
      '| scan | flow | valid | RI_tail median [Q1,Q3] | change vs flow01 | vessel raw mean | tail raw mean |','|---|---:|---:|---:|---:|---:|---:|']
    for r in ri.itertuples(index=False): lines.append(f'| {r.scan_id} | {r.flow_mm_s:g} | {r.n_valid_frames} | {r.ri_tail_median:.6f} [{r.ri_tail_q1:.6f},{r.ri_tail_q3:.6f}] | {r.ri_delta_pct_vs_flow01:+.2f}% | {r.source_mean_raw_median:.6g} | {r.tail_mean_raw_median:.6g} |')
    lines += ['', '## Fixed-128 geometry stress test','',
      f"- RI_tail framewise |Δ| median/P95: **{g.loc['ri_tail','abs_delta_pct_median']:.2f}% / {g.loc['ri_tail','abs_delta_pct_p95']:.2f}%**.",
      f"- Source raw mean |Δ| median/P95: **{g.loc['source_mean_raw','abs_delta_pct_median']:.2f}% / {g.loc['source_mean_raw','abs_delta_pct_p95']:.2f}%**.",'',
      '## Spatial autocorrelation of new RI_tail','',
      f"Across five volumes, median lag-1 ACF = **{a.acf_lag1_median:.3f}**; median 1/e crossing = **{a.lag_1e_median:.2f} frames**; median initial-positive correlation span = **{a.tau_median:.2f} frames**. This is spatial dependence, not biological n.",'',
      '## Systematic subsampling stability of new RI_tail','',
      '| stride | ~positions/500 | median |error| | P95 |error| | max |error| |','|---:|---:|---:|---:|---:|']
    for r in sub_cross.itertuples(index=False): lines.append(f'| {r.stride_frames} | {r.approx_positions_per_500:.0f} | {r.abs_error_pct_median:.2f}% | {r.abs_error_pct_p95:.2f}% | {r.abs_error_pct_max:.2f}% |')
    lines += ['', '## Selected RI(r) medians','', '| flow | 0 | 50 | 100 | 200 | 300 | 400 | 500 µm |','|---:|---:|---:|---:|---:|---:|---:|---:|']
    for flowv,gp in depth_selected.groupby('flow_mm_s'):
        d={int(r.target_r_um):r.ri_r_median for r in gp.itertuples(index=False)}
        lines.append('| '+str(int(flowv))+' | '+' | '.join(f'{d[x]:.4f}' for x in SELECTED_DEPTHS)+' |')
    lines += ['', '## Validation','',
      f"- valid frames: **{val['valid_frames']}**", f"- raw-array packages verified: **{val['release_packages_verified']}/{EXPECTED_PACKAGES}**", f"- per-frame NPZ SHA verified: **{val['valid_frame_npz_sha256_verified']}**",
      f"- baseline observed raw-SV replay max relative error: **{val['baseline_observed_replay_max_relative_error']:.3e}**", '- no interpolation, no zero-fill, no background subtraction, no p-values.', '- one scan volume per flow: all flow comparisons remain descriptive until independent-volume replication.','',
      '## Files','', '- `framewise_primary.csv`', '- `scan_flow_summary.csv`', '- `depth_selected.csv` / `depth_full_10um.csv`', '- `fixed128_observed_framewise.csv` / `fixed128_geometry_summary.csv`', '- `acf_ri_tail.csv` / `acf_summary.csv`', '- `subsampling_phase_results.csv` / `subsampling_summary.csv`', '- `validation.json`, `provenance.json`, `input_sha256.csv`']
    (out/'README.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--root',default='.'); ap.add_argument('--output-dir',default='analysis/formal_sv_d128_v21_run001/no_background_relative_tail_v1_full2422'); ap.add_argument('--source-sha',required=True); ap.add_argument('--workflow-source-sha',required=True); a=ap.parse_args()
    root=Path(a.root).resolve(); out=(root/a.output_dir).resolve(); out.mkdir(parents=True,exist_ok=True)
    formal=root/'results/formal_sv_d128_v21_full2500_run001'; task1=root/'analysis/formal_sv_d128_v21_run001/relative_tail_intensity_full2422'; observed=root/'analysis/formal_sv_d128_v21_run001/observed_tail_intensity_full2422'; helper_path=root/'analysis/formal_sv_d128_v21_run001/fixed128_width_sensitivity_full2422/audit_fixed128_width.py'
    H=load_module('fixed128_helper',helper_path)
    import sys; sys.path.insert(0,str(root/'src')); from svrecttail.geometry import VesselGeometry,ellipse_weights,interval_overlap_weights
    base,_,cfg=H.prepare_inputs(root,formal,task1)
    obs=pd.read_csv(observed/'observed_tail_intensity_framewise.csv'); depth=pd.read_csv(observed/'observed_depth_summary.csv'); packages=pd.read_csv(formal/'download_packages.csv'); arrays=pd.read_csv(formal/'arrays_sha256.csv')
    obs['scan_id']=obs.scan_id.astype(str); obs['frame_index']=pd.to_numeric(obs.frame_index).astype(int)
    if len(obs)!=EXPECTED_VALID or obs.groupby('scan_id').size().to_dict()!=EXPECTED_COUNTS: raise AssertionError('observed input coverage mismatch')
    arrays['scan_id']=arrays.scan_id.astype(str); arrays['frame_index_0based']=pd.to_numeric(arrays.frame_index_0based).astype(int); ah={(r.scan_id,int(r.frame_index_0based)):str(r.sha256) for r in arrays.itertuples(index=False)}
    merged=base.merge(obs[['scan_id','frame_index','q_vessel_observed','q_tail_observed','source_mean_observed','tail_mean_observed','ri_tail_observed']],on=['scan_id','frame_index'],validate='one_to_one')
    fixed=[]; pkg_rows=[]; verified=0; replay_max=0.; pix=cfg['dx_um']*cfg['dz_um']
    with tempfile.TemporaryDirectory(prefix='sv-nobg-primary-') as td:
      td=Path(td)
      for pkg in packages.itertuples(index=False):
        name=str(pkg.file); scan,lo,hi=H.parse_package_name(name); path=td/name; H.download_file(f'{RELEASE_BASE}/{name}',path)
        zsha=H.sha256_file(path)
        if zsha!=str(pkg.sha256) or path.stat().st_size!=int(pkg.bytes): raise AssertionError(f'package mismatch {name}')
        wanted=merged[(merged.scan_id.eq(scan))&merged.frame_index.between(lo,hi)]; vin=0
        with zipfile.ZipFile(path) as zf:
          for row in wanted.itertuples(index=False):
            fi=int(row.frame_index); member=f'arrays/{scan}/frame_{fi:03d}.npz'; data=zf.read(member)
            if H.sha256_bytes(data)!=ah[(scan,fi)]: raise AssertionError(f'NPZ SHA mismatch {scan}/{fi}')
            verified+=1; vin+=1
            with np.load(io.BytesIO(data),allow_pickle=False) as d: sv=np.asarray(d['sv_raw'],dtype=np.float64)
            def calc(width_um):
                x4=float(row.x4_px); wpx=width_um/cfg['dx_um']; geom=VesselGeometry(x_left_edge_px=x4-wpx/2,x_right_edge_px=x4+wpx/2,z_top_edge_px=float(row.z_top_edge_px),diameter_um=cfg['diameter_um'],dx_um=cfg['dx_um'],dz_um=cfg['dz_um'])
                sw=ellipse_weights(sv.shape,geom,supersample=cfg['ellipse_supersample']); xw=interval_overlap_weights(sv.shape[1],geom.x_left_edge_px,geom.x_right_edge_px); z0=geom.z_bottom_edge_px; zw=interval_overlap_weights(sv.shape[0],z0,z0+cfg['tail_length_um']/cfg['dz_um']); tw=np.multiply.outer(zw,xw)
                av=float(sw.sum()*pix); at=float(tw.sum()*pix); qv=weighted_integral(sv,sw,pix); qt=weighted_integral(sv,tw,pix); sm=safe_div(qv,av); tm=safe_div(qt,at); return qv,qt,sm,tm,safe_div(tm,sm),av,at
            b=calc(float(row.baseline_lateral_width_um)); f=calc(128.0)
            refs=[float(row.q_vessel_observed),float(row.q_tail_observed),float(row.source_mean_observed),float(row.tail_mean_observed),float(row.ri_tail_observed)]
            for x,y in zip(b[:5],refs): replay_max=max(replay_max,H.abs_rel_error(x,y))
            if replay_max>1e-10: raise AssertionError(f'raw observed replay tolerance exceeded {replay_max}')
            fixed.append(dict(scan_id=scan,frame_index=fi,flow_mm_s=float(row.flow_mm_s),baseline_width_um=float(row.baseline_lateral_width_um),fixed_width_um=128.0,q_vessel_raw_base=b[0],q_vessel_raw_fixed128=f[0],q_tail_raw_base=b[1],q_tail_raw_fixed128=f[1],source_mean_raw_base=b[2],source_mean_raw_fixed128=f[2],tail_mean_raw_base=b[3],tail_mean_raw_fixed128=f[3],ri_tail_base=b[4],ri_tail_fixed128=f[4],ri_tail_delta_pct=pct_delta(f[4],b[4]),source_mean_raw_delta_pct=pct_delta(f[2],b[2]),tail_mean_raw_delta_pct=pct_delta(f[3],b[3])))
        pkg_rows.append(dict(package=name,verified=True,valid_frames=len(wanted),npz_verified=vin)); path.unlink()
    fx=pd.DataFrame(fixed).sort_values(['flow_mm_s','frame_index']); fx.to_csv(out/'fixed128_observed_framewise.csv',index=False,float_format='%.17g')
    if len(fx)!=EXPECTED_VALID or verified!=EXPECTED_VALID: raise AssertionError('fixed128 coverage mismatch')
    # Primary framewise standardized output.
    primary=obs[['scan_id','frame_index','frame_index_0based','flow_mm_s','source_mean_observed','tail_mean_observed','ri_tail_observed']].rename(columns={'source_mean_observed':'source_mean_raw','tail_mean_observed':'tail_mean_raw','ri_tail_observed':'ri_tail'}).copy(); primary.to_csv(out/'framewise_primary.csv',index=False,float_format='%.17g')
    # Five-flow summaries and contrasts.
    flow_rows=[]
    for scan,g in primary.groupby('scan_id',sort=False):
        qr=quantiles(g.ri_tail); qs=quantiles(g.source_mean_raw); qt=quantiles(g.tail_mean_raw)
        flow_rows.append(dict(scan_id=scan,flow_mm_s=float(g.flow_mm_s.iloc[0]),n_valid_frames=len(g),ri_tail_q1=qr['q1'],ri_tail_median=qr['median'],ri_tail_q3=qr['q3'],source_mean_raw_median=qs['median'],tail_mean_raw_median=qt['median']))
    flow=pd.DataFrame(flow_rows).sort_values('flow_mm_s'); ref=float(flow.loc[flow.flow_mm_s.eq(1),'ri_tail_median'].iloc[0]); flow['ri_delta_pct_vs_flow01']=100*(flow.ri_tail_median-ref)/ref; flow.to_csv(out/'scan_flow_summary.csv',index=False,float_format='%.17g')
    # Depth new RI(r), no interpolation: reuse validated native-nearest 10um summary from observed raw-SV run.
    d=depth[['scan_id','flow_mm_s','target_r_um','n_frames','selected_r_um_median','max_abs_sampling_error_um','ri_r_observed_q1','ri_r_observed_median','ri_r_observed_q3']].rename(columns={'ri_r_observed_q1':'ri_r_q1','ri_r_observed_median':'ri_r_median','ri_r_observed_q3':'ri_r_q3'}); d.to_csv(out/'depth_full_10um.csv',index=False,float_format='%.17g'); ds=d[d.target_r_um.isin(SELECTED_DEPTHS)].copy(); ds.to_csv(out/'depth_selected.csv',index=False,float_format='%.17g')
    # Geometry sensitivity under the new no-background RI definition.
    geom_rows=[]
    for metric,col in [('ri_tail','ri_tail_delta_pct'),('source_mean_raw','source_mean_raw_delta_pct'),('tail_mean_raw','tail_mean_raw_delta_pct')]:
        q=quantiles(np.abs(fx[col])); geom_rows.append(dict(metric=metric,n=q['n'],abs_delta_pct_median=q['median'],abs_delta_pct_p95=q['p95'],abs_delta_pct_max=q['max']))
    geom=pd.DataFrame(geom_rows); geom.to_csv(out/'fixed128_geometry_summary.csv',index=False,float_format='%.17g')
    fx_scan=[]
    for scan,g in fx.groupby('scan_id'):
        fx_scan.append(dict(scan_id=scan,flow_mm_s=float(g.flow_mm_s.iloc[0]),ri_tail_base_median=float(g.ri_tail_base.median()),ri_tail_fixed128_median=float(g.ri_tail_fixed128.median()),delta_pct_of_scan_medians=pct_delta(g.ri_tail_fixed128.median(),g.ri_tail_base.median())))
    pd.DataFrame(fx_scan).sort_values('flow_mm_s').to_csv(out/'fixed128_scan_summary.csv',index=False,float_format='%.17g')
    # Spatial ACF for new RI_tail only.
    acf_rows=[]; acf_sum=[]
    for scan,g in primary.groupby('scan_id'):
        series=np.full(500,np.nan); series[g.frame_index_0based.astype(int).to_numpy()]=g.ri_tail.to_numpy(float); acf,rows=acf_pairwise(series,200)
        for lag,npair,val in rows: acf_rows.append(dict(scan_id=scan,flow_mm_s=float(g.flow_mm_s.iloc[0]),lag_frames=lag,pair_count=npair,acf=val))
        nonpos=next((i for i in range(1,len(acf)) if np.isfinite(acf[i]) and acf[i]<=0),np.nan); tau,stop=initial_positive_tau(acf)
        acf_sum.append(dict(scan_id=scan,flow_mm_s=float(g.flow_mm_s.iloc[0]),n_valid=len(g),acf_lag1=acf[1],acf_lag5=acf[5],acf_lag10=acf[10],acf_lag20=acf[20],lag_1e=threshold_cross(np.arange(len(acf)),acf,math.exp(-1)),lag_0p1=threshold_cross(np.arange(len(acf)),acf,.1),first_nonpositive=nonpos,initial_positive_tau=tau,descriptive_spatial_ess=len(g)/tau))
    acf_df=pd.DataFrame(acf_rows); acf_df.to_csv(out/'acf_ri_tail.csv',index=False,float_format='%.17g'); asum=pd.DataFrame(acf_sum).sort_values('flow_mm_s'); asum.to_csv(out/'acf_summary.csv',index=False,float_format='%.17g')
    acf_cross=pd.DataFrame([dict(acf_lag1_median=float(asum.acf_lag1.median()),lag_1e_median=float(asum.lag_1e.median()),lag_0p1_median=float(asum.lag_0p1.median()),first_nonpositive_median=float(asum.first_nonpositive.median()),tau_median=float(asum.initial_positive_tau.median()),spatial_ess_median=float(asum.descriptive_spatial_ess.median()))]); acf_cross.to_csv(out/'acf_cross_scan_summary.csv',index=False,float_format='%.17g')
    # Exhaustive systematic subsampling for the new RI_tail only.
    phase=[]
    for scan,g in primary.groupby('scan_id'):
        idx=g.frame_index_0based.astype(int).to_numpy(); full=float(g.ri_tail.median())
        for stride in STRIDES:
            for ph in range(stride):
                sub=g[(idx%stride)==ph]; est=float(sub.ri_tail.median()); err=pct_delta(est,full); phase.append(dict(scan_id=scan,flow_mm_s=float(g.flow_mm_s.iloc[0]),stride_frames=stride,phase=ph,n_selected=len(sub),full_median=full,subsample_median=est,error_pct=err,abs_error_pct=abs(err)))
    phase=pd.DataFrame(phase); phase.to_csv(out/'subsampling_phase_results.csv',index=False,float_format='%.17g')
    subrows=[]
    for stride,g in phase.groupby('stride_frames'):
        q=quantiles(g.abs_error_pct); subrows.append(dict(stride_frames=int(stride),approx_positions_per_500=500/int(stride),n_phase_records=len(g),selected_n_min=int(g.n_selected.min()),selected_n_median=float(g.n_selected.median()),selected_n_max=int(g.n_selected.max()),abs_error_pct_median=q['median'],abs_error_pct_p95=q['p95'],abs_error_pct_max=q['max'],fraction_within_5pct=float((g.abs_error_pct<=5).mean()),fraction_within_10pct=float((g.abs_error_pct<=10).mean())))
    sub=pd.DataFrame(subrows).sort_values('stride_frames'); sub.to_csv(out/'subsampling_summary.csv',index=False,float_format='%.17g')
    # Descriptive observed-effect context: new RI only. Background is absent by definition.
    geom_ri=geom.set_index('metric').loc['ri_tail']; sub100=sub[sub.stride_frames.eq(5)].iloc[0]; context=[]
    for r in flow.itertuples(index=False):
        if r.flow_mm_s==1: continue
        effect=abs(r.ri_delta_pct_vs_flow01); local_p95=float(sub100.abs_error_pct_p95); geom_p95=float(geom_ri.abs_delta_pct_p95)
        cls='above_sampling_and_geometry_stress' if effect>geom_p95 and effect>local_p95 else ('above_sampling_within_geometry_stress' if effect>local_p95 else 'within_sampling_sensitivity')
        context.append(dict(scan_id=r.scan_id,flow_mm_s=r.flow_mm_s,ri_delta_pct_vs_flow01=r.ri_delta_pct_vs_flow01,observed_abs_delta_pct=effect,subsampling_100_p95_pct=local_p95,fixed128_geometry_p95_pct=geom_p95,descriptive_context_class=cls))
    pd.DataFrame(context).to_csv(out/'flow_effect_context.csv',index=False,float_format='%.17g')
    # Validation/provenance.
    val=dict(valid_frames=len(primary),scan_valid_counts=primary.groupby('scan_id').size().to_dict(),release_packages_verified=len(pkg_rows),valid_frame_npz_sha256_verified=verified,baseline_observed_replay_max_relative_error=float(replay_max),primary_background_subtracted=False,primary_definition='mean raw SV tail / mean raw SV vessel',depth_definition='raw SV lateral mean / raw vessel mean',fixed128_background_used=False,invalid_frames_interpolated=False,invalid_frames_zero_filled=False,p_values_computed=False,bscans_treated_as_independent_biological_replicates=False)
    (out/'validation.json').write_text(json.dumps(val,indent=2)+'\n',encoding='utf-8')
    inputs=[observed/'observed_tail_intensity_framewise.csv',observed/'observed_depth_summary.csv',formal/'download_packages.csv',formal/'arrays_sha256.csv',formal/'run_config.json',helper_path,root/'src/svrecttail/geometry.py']
    pd.DataFrame([dict(input_path=p.relative_to(root).as_posix(),sha256=sha256(p),size_bytes=p.stat().st_size) for p in inputs]).to_csv(out/'input_sha256.csv',index=False)
    prov=dict(analysis='No-background Relative Tail Intensity primary reanalysis v1',source_sha=a.source_sha,workflow_source_sha=a.workflow_source_sha,formal_release=RELEASE_TAG,entry_rule='valid == True only; 2422/2500',geometry='frozen X4 center + X1 apparent width + v2.1 z; 500um tail; guard=0',primary_scalar='RI_tail = tail_mean_raw/source_mean_raw',primary_depth='RI(r)=raw tail lateral mean/source_mean_raw',background='not used in any primary calculation',legacy_background_corrected_outputs='retained historically only; not used',fixed128='geometry stress test only; no background involved',statistics='volume-level descriptive only; no p-values; B-scans are not independent biological replicates')
    (out/'provenance.json').write_text(json.dumps(prov,indent=2)+'\n',encoding='utf-8')
    write_readme(out,flow,geom,acf_cross,sub,ds,val)
    print(json.dumps({'valid_frames':len(primary),'flow_summary':flow[['scan_id','flow_mm_s','ri_tail_median','ri_delta_pct_vs_flow01']].to_dict('records'),'fixed128_ri_abs_delta_median_pct':float(geom_ri.abs_delta_pct_median),'fixed128_ri_abs_delta_p95_pct':float(geom_ri.abs_delta_pct_p95),'subsampling':sub.to_dict('records'),'output_dir':str(out.relative_to(root))},indent=2))

if __name__=='__main__': main()
