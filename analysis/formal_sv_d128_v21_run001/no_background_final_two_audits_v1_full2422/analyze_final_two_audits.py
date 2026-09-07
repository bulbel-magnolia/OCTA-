#!/usr/bin/env python3
"""Final two descriptive audits for no-background SV RI_tail.

1) Aggregation sensitivity: compare pre-specified scan-volume summaries.
2) Residual geometry/localization QC: audit remaining geometry and localization-source coupling.

No background subtraction. No p-values. B-scans are spatial samples within one volume,
not independent experimental replicates.
"""
from __future__ import annotations

import argparse, hashlib, importlib.util, json
from pathlib import Path
import numpy as np
import pandas as pd

EXPECTED_COUNTS={"flow01":486,"flow03":468,"flow05":491,"flow07":493,"flow10":484}
EXPECTED_VALID=sum(EXPECTED_COUNTS.values())
SEGMENTS=((0,99),(100,199),(200,299),(300,399),(400,499))


def sha256(path:Path)->str:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''): h.update(b)
    return h.hexdigest()

def safe_div(a,b):
    a=float(a); b=float(b)
    return a/b if np.isfinite(a) and np.isfinite(b) and b!=0 else np.nan

def pct(new,old): return 100.0*safe_div(float(new)-float(old),old)

def corr(x,y):
    x=np.asarray(x,float); y=np.asarray(y,float); m=np.isfinite(x)&np.isfinite(y)
    if m.sum()<3 or np.std(x[m])==0 or np.std(y[m])==0: return np.nan
    return float(np.corrcoef(x[m],y[m])[0,1])

def spearman(x,y):
    x=pd.Series(np.asarray(x,float)).rank(method='average').to_numpy(float)
    y=pd.Series(np.asarray(y,float)).rank(method='average').to_numpy(float)
    return corr(x,y)

def qstats(x):
    a=pd.to_numeric(pd.Series(x),errors='coerce').to_numpy(float); a=a[np.isfinite(a)]
    if not len(a): return dict(n=0,q1=np.nan,median=np.nan,q3=np.nan,mean=np.nan,min=np.nan,max=np.nan)
    q1,m,q3=np.quantile(a,[.25,.5,.75])
    return dict(n=int(len(a)),q1=float(q1),median=float(m),q3=float(q3),mean=float(a.mean()),min=float(a.min()),max=float(a.max()))

def load_module(name,path):
    spec=importlib.util.spec_from_file_location(name,path); mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod); return mod


def aggregation_audit(df:pd.DataFrame,out:Path):
    rows=[]
    for scan,g0 in df.groupby('scan_id',sort=False):
        g=g0.sort_values('frame_index_0based').copy(); flow=float(g.flow_mm_s.iloc[0])
        # Pre-specified methods.
        methods={}
        methods['primary_median_framewise_RI']=float(g.ri_tail.median())
        methods['ratio_of_median_raw_means']=safe_div(g.tail_mean_raw.median(),g.source_mean_raw.median())
        methods['ratio_of_mean_raw_means']=safe_div(g.tail_mean_raw.mean(),g.source_mean_raw.mean())
        methods['mean_framewise_RI']=float(g.ri_tail.mean())
        methods['pooled_area_weighted_raw_ratio']=safe_div(
            safe_div(g.q_tail_raw.sum(),g.tail_area_um2.sum()),
            safe_div(g.q_vessel_raw.sum(),g.source_area_um2.sum()))
        seg_medians=[]
        for lo,hi in SEGMENTS:
            sg=g[g.frame_index_0based.between(lo,hi)]
            if len(sg)==0: raise AssertionError(f'empty segment {scan} {lo}-{hi}')
            seg_medians.append(float(sg.ri_tail.median()))
        methods['median_of_5_spatial_segment_medians']=float(np.median(seg_medians))
        primary=methods['primary_median_framewise_RI']
        for method,value in methods.items():
            rows.append({'scan_id':scan,'flow_mm_s':flow,'n_valid':len(g),'aggregation_method':method,'ri_volume':value,
                         'delta_pct_vs_primary_within_scan':pct(value,primary)})
    ag=pd.DataFrame(rows).sort_values(['aggregation_method','flow_mm_s'])
    # Add flow01 reference and rank per method.
    outrows=[]
    for method,g in ag.groupby('aggregation_method',sort=False):
        ref=float(g[g.flow_mm_s.eq(1)].ri_volume.iloc[0])
        ranks=g.ri_volume.rank(method='min',ascending=False)
        for idx,r in g.iterrows():
            d=r.to_dict(); d['delta_pct_vs_flow01']=pct(r.ri_volume,ref); d['rank_desc']=int(ranks.loc[idx]); outrows.append(d)
    ag2=pd.DataFrame(outrows).sort_values(['aggregation_method','flow_mm_s'])
    ag2.to_csv(out/'aggregation_all_methods.csv',index=False,float_format='%.17g')

    summary=[]
    for method,g in ag2.groupby('aggregation_method',sort=False):
        f7=g[g.flow_mm_s.eq(7)].iloc[0]; others=g[~g.flow_mm_s.eq(7)]
        summary.append({'aggregation_method':method,'flow07_value':float(f7.ri_volume),'flow01_value':float(g[g.flow_mm_s.eq(1)].ri_volume.iloc[0]),
                        'other4_scan_median':float(others.ri_volume.median()),'flow07_vs_flow01_pct':float(f7.delta_pct_vs_flow01),
                        'flow07_vs_other4_median_pct':pct(f7.ri_volume,others.ri_volume.median()),'flow07_rank_desc':int(f7.rank_desc),
                        'max_abs_within_scan_delta_pct_vs_primary':float(np.max(np.abs(g.delta_pct_vs_primary_within_scan)))})
    sm=pd.DataFrame(summary); sm.to_csv(out/'aggregation_summary.csv',index=False,float_format='%.17g')
    # Per-scan dispersion across aggregation methods.
    disp=[]
    for scan,g in ag2.groupby('scan_id',sort=False):
        primary=float(g[g.aggregation_method.eq('primary_median_framewise_RI')].ri_volume.iloc[0])
        disp.append({'scan_id':scan,'flow_mm_s':float(g.flow_mm_s.iloc[0]),'primary':primary,'min_across_methods':float(g.ri_volume.min()),
                     'max_across_methods':float(g.ri_volume.max()),'range_pct_of_primary':100*(float(g.ri_volume.max())-float(g.ri_volume.min()))/primary,
                     'max_abs_delta_pct_vs_primary':float(np.max(np.abs(g.delta_pct_vs_primary_within_scan)))})
    pd.DataFrame(disp).sort_values('flow_mm_s').to_csv(out/'aggregation_dispersion_by_scan.csv',index=False,float_format='%.17g')
    return ag2,sm


def residual_qc(df:pd.DataFrame,out:Path):
    geom=['lateral_width_um','source_area_um2','z_top_edge_px','x4_px']
    metrics=['ri_tail','source_mean_raw','tail_mean_raw']
    corrrows=[]
    for scan,g in df.groupby('scan_id',sort=False):
        for m in metrics:
            for q in geom:
                corrrows.append({'scan_id':scan,'flow_mm_s':float(g.flow_mm_s.iloc[0]),'signal_metric':m,'geometry_metric':q,
                                 'spearman':spearman(g[q],g[m]),'pearson':corr(g[q],g[m]),'n':len(g)})
    cr=pd.DataFrame(corrrows).sort_values(['flow_mm_s','signal_metric','geometry_metric']); cr.to_csv(out/'geometry_signal_correlations.csv',index=False,float_format='%.17g')

    # Scan-level geometry context.
    gr=[]
    for scan,g in df.groupby('scan_id',sort=False):
        rec={'scan_id':scan,'flow_mm_s':float(g.flow_mm_s.iloc[0]),'n_valid':len(g)}
        for q in geom:
            s=qstats(g[q]); rec[f'{q}_median']=s['median']; rec[f'{q}_q1']=s['q1']; rec[f'{q}_q3']=s['q3']
        gr.append(rec)
    gs=pd.DataFrame(gr).sort_values('flow_mm_s'); gs.to_csv(out/'geometry_scan_summary.csv',index=False,float_format='%.17g')
    other=gs[~gs.flow_mm_s.eq(7)]; f7=gs[gs.flow_mm_s.eq(7)].iloc[0]
    ctx=[]
    for q in geom:
        v=float(f7[f'{q}_median']); o=float(other[f'{q}_median'].median()); f1=float(gs[gs.flow_mm_s.eq(1)][f'{q}_median'].iloc[0])
        ctx.append({'geometry_metric':q,'flow07_median':v,'flow01_median':f1,'other4_scan_median':o,
                    'flow07_minus_flow01':v-f1,'flow07_minus_other4_median':v-o,
                    'flow07_vs_flow01_pct':pct(v,f1) if f1!=0 else np.nan,'flow07_vs_other4_median_pct':pct(v,o) if o!=0 else np.nan})
    pd.DataFrame(ctx).to_csv(out/'geometry_flow07_context.csv',index=False,float_format='%.17g')

    # Localization source composition and within-category medians.
    counts=[]; cats=[]
    for scan,g in df.groupby('scan_id',sort=False):
        total=len(g)
        for loc,lg in g.groupby('localization_source',dropna=False,sort=True):
            cats.append(str(loc))
            counts.append({'scan_id':scan,'flow_mm_s':float(g.flow_mm_s.iloc[0]),'localization_source':str(loc),'n':len(lg),'fraction':len(lg)/total,
                           'ri_median':float(lg.ri_tail.median()),'source_mean_median':float(lg.source_mean_raw.median()),'tail_mean_median':float(lg.tail_mean_raw.median())})
    lc=pd.DataFrame(counts).sort_values(['flow_mm_s','localization_source']); lc.to_csv(out/'localization_source_summary.csv',index=False,float_format='%.17g')

    # Restricted-sample sensitivity. Restrictions are within each volume, so they do not make B-scans independent.
    subsetrows=[]
    direct_name='direct_candidate_supported'
    for scan,g0 in df.groupby('scan_id',sort=False):
        flow=float(g0.flow_mm_s.iloc[0]); g=g0.copy()
        zlo,zhi=g.z_top_edge_px.quantile([.1,.9]); wlo,whi=g.lateral_width_um.quantile([.1,.9])
        masks={
            'all_valid':np.ones(len(g),dtype=bool),
            'direct_only':g.localization_source.astype(str).eq(direct_name).to_numpy(),
            'central80_z':g.z_top_edge_px.between(zlo,zhi).to_numpy(),
            'central80_width':g.lateral_width_um.between(wlo,whi).to_numpy(),
            'central80_z_and_width':(g.z_top_edge_px.between(zlo,zhi)&g.lateral_width_um.between(wlo,whi)).to_numpy(),
            'direct_and_central80_z_width':(g.localization_source.astype(str).eq(direct_name)&g.z_top_edge_px.between(zlo,zhi)&g.lateral_width_um.between(wlo,whi)).to_numpy(),
        }
        for name,mask in masks.items():
            sg=g.loc[mask]
            if len(sg)<20:
                subsetrows.append({'scan_id':scan,'flow_mm_s':flow,'subset':name,'n':len(sg),'ri_median':np.nan,'source_mean_median':np.nan,'tail_mean_median':np.nan})
            else:
                subsetrows.append({'scan_id':scan,'flow_mm_s':flow,'subset':name,'n':len(sg),'ri_median':float(sg.ri_tail.median()),
                                   'source_mean_median':float(sg.source_mean_raw.median()),'tail_mean_median':float(sg.tail_mean_raw.median())})
    ss=pd.DataFrame(subsetrows).sort_values(['subset','flow_mm_s']); ss.to_csv(out/'restricted_subset_scan_summary.csv',index=False,float_format='%.17g')
    sc=[]
    for subset,g in ss.groupby('subset',sort=False):
        if not np.isfinite(g.ri_median).all():
            sc.append({'subset':subset,'all_5_scans_evaluable':False}); continue
        f7=g[g.flow_mm_s.eq(7)].iloc[0]; f1=g[g.flow_mm_s.eq(1)].iloc[0]; other=g[~g.flow_mm_s.eq(7)]
        ranks=g.ri_median.rank(method='min',ascending=False)
        sc.append({'subset':subset,'all_5_scans_evaluable':True,'flow07_ri_median':float(f7.ri_median),'flow01_ri_median':float(f1.ri_median),
                   'other4_scan_median_ri':float(other.ri_median.median()),'flow07_vs_flow01_pct':pct(f7.ri_median,f1.ri_median),
                   'flow07_vs_other4_median_pct':pct(f7.ri_median,other.ri_median.median()),'flow07_rank_desc':int(ranks.loc[f7.name]),
                   'flow07_source_vs_other4_pct':pct(f7.source_mean_median,other.source_mean_median.median()),
                   'flow07_tail_vs_other4_pct':pct(f7.tail_mean_median,other.tail_mean_median.median())})
    sctx=pd.DataFrame(sc); sctx.to_csv(out/'restricted_subset_flow07_context.csv',index=False,float_format='%.17g')
    return cr,gs,lc,ss,sctx


def write_readme(out,aggsum,cr,gs,lc,sctx,val):
    lines=['# Final two audits for no-background SV RI_tail v1','',
           'Primary endpoint remains `RI_tail = mean(raw SV in tail ROI) / mean(raw SV in vessel ROI)`. No background subtraction is used.','',
           '## 1. Aggregation sensitivity','',
           '| aggregation method | flow07 | flow07 vs flow01 | flow07 vs other4 median | rank | max within-scan change vs primary |',
           '|---|---:|---:|---:|---:|---:|']
    for r in aggsum.itertuples(index=False):
        lines.append(f'| {r.aggregation_method} | {r.flow07_value:.6f} | {r.flow07_vs_flow01_pct:+.2f}% | {r.flow07_vs_other4_median_pct:+.2f}% | {r.flow07_rank_desc}/5 | {r.max_abs_within_scan_delta_pct_vs_primary:.2f}% |')
    allrank=bool((aggsum.flow07_rank_desc==1).all())
    lines += ['',f'- Flow07 ranks first under **{int((aggsum.flow07_rank_desc==1).sum())}/{len(aggsum)}** pre-specified aggregation methods.' if allrank else f'- Flow07 ranks first under {int((aggsum.flow07_rank_desc==1).sum())}/{len(aggsum)} methods.','',
              '## 2. Residual geometry/localization QC','',
              '### Flow07 geometry medians','',
              '| metric | flow07 | flow01 | other4 scan median | flow07 vs other4 |','|---|---:|---:|---:|---:|']
    gctx=pd.read_csv(out/'geometry_flow07_context.csv')
    for r in gctx.itertuples(index=False): lines.append(f'| {r.geometry_metric} | {r.flow07_median:.6g} | {r.flow01_median:.6g} | {r.other4_scan_median:.6g} | {r.flow07_vs_other4_median_pct:+.2f}% |')
    lines += ['', '### RI_tail vs geometry: Spearman by scan','', '| scan | width | source area | z_top | X4 |','|---|---:|---:|---:|---:|']
    rcr=cr[cr.signal_metric.eq('ri_tail')]
    for scan,g in rcr.groupby('scan_id',sort=False):
        d={r.geometry_metric:r.spearman for r in g.itertuples(index=False)}
        lines.append(f"| {scan} | {d['lateral_width_um']:.3f} | {d['source_area_um2']:.3f} | {d['z_top_edge_px']:.3f} | {d['x4_px']:.3f} |")
    lines += ['', '### Restricted-subset sensitivity','', '| subset | flow07 RI | vs flow01 | vs other4 median | rank | vessel vs other4 | tail vs other4 |','|---|---:|---:|---:|---:|---:|---:|']
    for r in sctx.itertuples(index=False):
        if not bool(r.all_5_scans_evaluable): lines.append(f'| {r.subset} | NA | NA | NA | NA | NA | NA |'); continue
        lines.append(f'| {r.subset} | {r.flow07_ri_median:.6f} | {r.flow07_vs_flow01_pct:+.2f}% | {r.flow07_vs_other4_median_pct:+.2f}% | {int(r.flow07_rank_desc)}/5 | {r.flow07_source_vs_other4_pct:+.2f}% | {r.flow07_tail_vs_other4_pct:+.2f}% |')
    lines += ['', 'Localization-source composition and category-specific medians are in `localization_source_summary.csv`. Restrictions are QC sensitivity subsets only; they do not create independent replicates.','',
              '## Validation','',f"- valid frames: **{val['valid_frames']}**",f"- scan counts: `{val['scan_counts']}`",'- no background subtraction, interpolation, p-values, or B-scan pseudo-replication.','',
              '## Files','', '- `aggregation_all_methods.csv`, `aggregation_summary.csv`, `aggregation_dispersion_by_scan.csv`',
              '- `geometry_signal_correlations.csv`, `geometry_scan_summary.csv`, `geometry_flow07_context.csv`',
              '- `localization_source_summary.csv`, `restricted_subset_scan_summary.csv`, `restricted_subset_flow07_context.csv`',
              '- `validation.json`, `provenance.json`, `input_sha256.csv`']
    (out/'README.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--root',default='.'); ap.add_argument('--output-dir',default='analysis/formal_sv_d128_v21_run001/no_background_final_two_audits_v1_full2422'); ap.add_argument('--source-sha',required=True); ap.add_argument('--workflow-source-sha',required=True); a=ap.parse_args()
    root=Path(a.root).resolve(); out=(root/a.output_dir).resolve(); out.mkdir(parents=True,exist_ok=True)
    obs_path=root/'analysis/formal_sv_d128_v21_run001/observed_tail_intensity_full2422/observed_tail_intensity_framewise.csv'
    obs=pd.read_csv(obs_path)
    obs=obs.rename(columns={'q_vessel_observed':'q_vessel_raw','q_tail_observed':'q_tail_raw','source_mean_observed':'source_mean_raw','tail_mean_observed':'tail_mean_raw','ri_tail_observed':'ri_tail'})
    if len(obs)!=EXPECTED_VALID or obs.groupby('scan_id').size().to_dict()!=EXPECTED_COUNTS: raise AssertionError('observed coverage mismatch')
    helper_path=root/'analysis/formal_sv_d128_v21_run001/fixed128_width_sensitivity_full2422/audit_fixed128_width.py'
    H=load_module('fixed128_helper_final',helper_path)
    formal=root/'results/formal_sv_d128_v21_full2500_run001'; task1=root/'analysis/formal_sv_d128_v21_run001/relative_tail_intensity_full2422'
    base,_,_=H.prepare_inputs(root,formal,task1)
    geom=base[['scan_id','frame_index','x4_px','z_top_edge_px','baseline_lateral_width_um']].rename(columns={'frame_index':'frame_index_0based','baseline_lateral_width_um':'lateral_width_um'})
    df=obs.merge(geom,on=['scan_id','frame_index_0based'],how='inner',validate='one_to_one')
    if len(df)!=EXPECTED_VALID: raise AssertionError('geometry merge coverage mismatch')
    req=['source_mean_raw','tail_mean_raw','ri_tail','q_vessel_raw','q_tail_raw','source_area_um2','tail_area_um2','lateral_width_um','z_top_edge_px','x4_px']
    if not np.isfinite(df[req].to_numpy(float)).all(): raise AssertionError('nonfinite final audit input')
    ag,ags=aggregation_audit(df,out)
    cr,gs,lc,ss,sctx=residual_qc(df,out)
    val={'valid_frames':len(df),'scan_counts':df.groupby('scan_id').size().to_dict(),'aggregation_methods':int(ags.shape[0]),'background_subtracted':False,'p_values_computed':False,'bscans_independent_replicates':False}
    (out/'validation.json').write_text(json.dumps(val,indent=2)+'\n',encoding='utf-8')
    prov={'analysis':'aggregation sensitivity + residual geometry/localization QC','source_sha':a.source_sha,'workflow_source_sha':a.workflow_source_sha,
          'primary_definition':'mean raw SV tail / mean raw SV vessel','aggregation_methods':ags.aggregation_method.tolist(),
          'residual_qc':['lateral_width_um','source_area_um2','z_top_edge_px','x4_px','localization_source','pre-specified restricted subsets'],
          'scope':'descriptive within-volume QC; no flow inference'}
    (out/'provenance.json').write_text(json.dumps(prov,indent=2)+'\n',encoding='utf-8')
    inputs=[obs_path,helper_path,formal/'frame_results.csv',task1/'relative_tail_intensity_framewise.csv']
    pd.DataFrame([{'input_path':p.relative_to(root).as_posix(),'sha256':sha256(p),'size_bytes':p.stat().st_size} for p in inputs]).to_csv(out/'input_sha256.csv',index=False)
    write_readme(out,ags,cr,gs,lc,sctx,val)
    print(json.dumps({'aggregation':ags.to_dict('records'),'restricted':sctx.to_dict('records'),'localization_sources':sorted(df.localization_source.astype(str).unique().tolist())},indent=2))

if __name__=='__main__': main()
