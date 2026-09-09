"""Volume-level QC, exact frozen negative-control pipeline and descriptive joins."""
import sys, json
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from derive_qc import OUT,ROOT,DEP,PREV,ID,REGIONS,read,csv,js,digest
from analyze_depth_coupling import block_rho,lag_cube,curve_summary,cell_template,describe_map,high_region,similarity

def summarize(v):
    a=np.asarray(v,float);a=a[np.isfinite(a)]
    if len(a)==0:return dict.fromkeys(['n','median','min','max','p5','p10','p90','p95','p99','mad'],np.nan)
    med=np.median(a)
    return dict(n=len(a),median=med,min=a.min(),max=a.max(),p5=np.quantile(a,.05),p10=np.quantile(a,.1),p90=np.quantile(a,.9),p95=np.quantile(a,.95),p99=np.quantile(a,.99),mad=np.median(abs(a-med)))
def rho(a,b):return similarity(a,b)
def scalar(a,b):
    m=np.isfinite(a)&np.isfinite(b);a=rankdata(np.asarray(a)[m]);b=rankdata(np.asarray(b)[m])
    if len(a)<3 or np.std(a)==0 or np.std(b)==0:return np.nan
    return float(np.corrcoef(a,b)[0,1])
def series_stats(s):
    a=s.dropna();med=a.median();mad=(a-med).abs().median();trend=s.rolling(51,center=True,min_periods=1).median();res=s-trend
    # Endpoint blocks are fixed original frames 0..50 and 449..499, not compressed valid rows.
    first=s.iloc[:51].median();last=s.iloc[-51:].median()
    slope=np.polyfit(a.index.to_numpy(float),a.to_numpy(),1)[0] if len(a)>1 else np.nan
    return dict(n_valid=len(a),mean=a.mean(),sd=a.std(ddof=0),mad=mad,p5=a.quantile(.05),p50=med,p95=a.quantile(.95),
       cv=a.std(ddof=0)/a.mean() if a.mean()>0 else np.nan,rcv=1.4826*mad/med if med>0 else np.nan,
       dynamic_range=(a.quantile(.95)-a.quantile(.05))/med if med>0 else np.nan,
       first_block_median=first,last_block_median=last,first_to_last_difference=last-first,
       first_to_last_relative=(last-first)/med if med>0 else np.nan,linear_slope_per_frame=slope,
       trend_amplitude=trend.max()-trend.min(),trend_amplitude_relative=(trend.max()-trend.min())/med if med>0 else np.nan,
       residual_mad=(res-res.median()).abs().median()),trend,res
def main():
    ex=json.loads((OUT/'extraction_validation.json').read_text());assert ex['status']=='passed'
    data=pd.concat([read(p) for p in sorted(OUT.glob('framewise_qc_D*.csv.gz'))],ignore_index=True)
    realmap=read(DEP/'volume_absolute_depth_coupling_summary.csv')
    # Read required robustness and geometry-control context without modifying it.
    flow=read(DEP/'flow_robustness_summary.csv');gc=read(DEP/'geometry_control_summary.csv')
    allmetrics=[];geometry=[];jitters=[];coverage=[];stability=[];slow=[];support=[];bgcorr=[];pseudo=[];pseudo_summ=[];excess=[];negative=[];volumes=[];special=[];struct=[];spotchecks=[]
    for key,g in data.groupby(ID,sort=True):
        meta=dict(zip(ID,key));print('Analyze QC '+meta['scan_id'],flush=True)
        g=g.sort_values('frame_index').set_index('frame_index');assert np.array_equal(g.index,np.arange(500))
        valid=g.valid_geometry.eq(True);v=g[valid];missing=~valid.to_numpy();runs=[];run=0
        for b in missing:
            if b:run+=1
            elif run:runs.append(run);run=0
        if run:runs.append(run)
        geom=meta|dict(nominal_frames=500,valid_geometry_frames=len(v),invalid_frames=500-len(v),valid_fraction=len(v)/500,
                      longest_invalid_run=max(runs,default=0),invalid_gaps=len(runs))
        adjacent=np.diff(v.index)==1
        for name,scale in [('X4',12.7),('X1_px',12.7),('z_top',6.7)]:
            j=np.abs(np.diff(v[name]))[adjacent];sm=summarize(j);label='X1' if name=='X1_px' else name
            flag=(j>sm['median']+5*1.4826*sm['mad']).mean()
            for unit,factor in [('pixel',1),('um',scale)]:
                geometry_row=meta|dict(parameter=label,unit=unit,n_adjacent_pairs=len(j),outlier_flag_fraction=flag)
                geometry_row.update({k:(val if k=='n' else val*factor) for k,val in sm.items()});jitters.append(geometry_row)
            geom[label+'_jitter_median_px']=sm['median'];geom[label+'_jitter_p95_px']=sm['p95'];geom[label+'_jitter_flag_fraction']=flag
        geometry.append(geom)
        cov=meta|dict(valid_geometry_frames=len(v))
        for side in ['left','right','any','both']:
            cov[side+'_valid_fraction']=v[side+'_valid'].eq(True).mean();cov[side+'_valid_frames']=int(v[side+'_valid'].eq(True).sum())
        coverage.append(cov)
        metrics={}
        for name in REGIONS:
            row=meta|dict(region=name)
            for field in ['z','sbr','sbr_db','cnr','contrast','real_mean','real_median','real_mad','bg_mean','bg_median','bg_mad','bg_sd',
                          'bg_asym_abs','bg_asym_relative','bg_logratio','bg_asym_log']:
                vals=v.get(f'{name}_{field}',pd.Series(dtype=float));sm=summarize(vals)
                for k,val in sm.items():row[field+'_'+k]=val
            for t in [0,1,2,3]:
                z=v[name+'_z'];row[f'z_gt{t}_fraction_evaluable']=float((z.dropna()>t).mean());row[f'z_gt{t}_fraction_all_geometry']=float((z>t).mean())
            row['z_evaluable_frames']=int(v[name+'_z'].notna().sum());row['z_missing_frames']=int(v[name+'_z'].isna().sum())
            metrics[name]=row;allmetrics.append(row)
            for side in ['real','left','right','bg']:
                support.append(meta|dict(region=name,side=side)|summarize(v[f'{name}_{side}_neff']))
            for flag in ['robust_den_zero','sd_den_zero']:
                special.append(meta|dict(region=name,flag=flag,n_frames=int(v[name+'_'+flag].eq(True).sum())))
        vr=meta|dict(valid_frames=len(v),source_z=metrics['source']['z_median'],source_sbr=metrics['source']['sbr_median'],
                    source_z_p10=metrics['source']['z_p10'],source_z_p90=metrics['source']['z_p90'],source_z_gt3=metrics['source']['z_gt3_fraction_evaluable'],
                    proximal_z=metrics['tail100']['z_median'],valid_fraction=geom['valid_fraction'],
                    z_top_jitter=geom['z_top_jitter_median_px'],X1_jitter=geom['X1_jitter_median_px'],X4_jitter=geom['X4_jitter_median_px'])
        for name in REGIONS:
            for side in ['bg','left','right']:
                for stat in ['mean','median']:
                    col=f'{name}_{side}_{stat}';ss,tr,res=series_stats(g[col])
                    stability.append(meta|dict(region=name,side=side,pixel_statistic=stat)|ss)
                    slow.append(pd.DataFrame(meta|dict(region=name,side=side,pixel_statistic=stat,frame_index=np.arange(500),raw=g[col].values,trend51=tr.where(valid).values,residual51=res.values)))
                    if name=='source' and side=='bg' and stat=='median':
                        vr.update(background_median=ss['p50'],background_mad=ss['mad'],background_rcv=ss['rcv'],background_drift=ss['first_to_last_relative'],background_trend_amplitude=ss['trend_amplitude_relative'])
        for name in ['source','tail100']:
            for field in ['sbr','contrast','real_mean','bg_mean','bg_sd']:
                col=f'structural_{name}_{field}'
                if col in v:struct.append(meta|dict(region=name,metric=field)|summarize(v[col]))
        real=realmap[realmap.scan_id.eq(key[0])].sort_values(['source_bin','tail_bin']);realhigh,_=high_region(real)
        realkeys=set(zip(realhigh.source_bin,realhigh.tail_bin))
        curves=g[[f'{name}_{side}_mean' for name in [f's{i}' for i in range(1,7)]+[f't{i:02}' for i in range(1,21)] for side in ['real','bg','left','right']]]
        modes={'raw':curves,'detrended_51':curves-curves.rolling(51,center=True,min_periods=1).median()}
        control_maps={}
        for mode,cur in modes.items():
            for k in range(1,21):
                for side in ['bg','left','right']:
                    a=cur[f't{k:02}_{side}_mean'];b=cur[f't{k:02}_real_mean']
                    bgcorr.append(meta|dict(mode=mode,tail_bin=k,side=side,rho=rho(a,b),n_pairs=int((a.notna()&b.notna()).sum())))
            for side in ['left','right']:
                x=cur[[f's{i}_{side}_mean' for i in range(1,7)]].to_numpy();y=cur[[f't{i:02}_{side}_mean' for i in range(1,21)]].to_numpy()
                template=cell_template(meta,'absolute',6)
                if (np.isfinite(x).all(axis=1)&np.isfinite(y).all(axis=1)).sum()<3:
                    sm=template.copy();sm['rho']=np.nan;sm['mode']=mode;sm['side']=side;pseudo.append(sm);continue
                cube,ns=lag_cube(x,y);sm=curve_summary(cube,ns,template,mode);sm['rho']=sm.rho_lag0;sm['side']=side;pseudo.append(sm)
                lf=pd.DataFrame(dict(lag=np.repeat(np.arange(-50,51),120),source_bin=np.tile(np.repeat(np.arange(1,7),20),101),tail_bin=np.tile(np.arange(1,21),606),rho=cube.ravel(),n_pairs=np.repeat(ns,120)))
                csv(f'lag_{key[0]}_{side}_{mode}.csv.gz',lf)
                desc=describe_map(sm);mask=np.array([(a,b) in realkeys for a,b in zip(sm.source_bin,sm.tail_bin)])
                pseudo_summ.append(meta|dict(side=side,mode=mode,real_region_median_rho=sm.loc[mask,'rho'].median(),
                  real_region_median_specificity=sm.loc[mask,'zero_lag_specificity'].median(),map_median_specificity=sm.zero_lag_specificity.median(),
                  peak_pm2_count=int(sm.peak_within_pm2.sum()),peak_pm5_count=int(sm.peak_within_pm5.sum()),
                  similarity_to_real=rho(real.rho,sm.rho) if mode=='detrended_51' else np.nan)|desc)
                if mode=='detrended_51':control_maps[side]=sm
                # Independent scalar rank/lag verification, one cell per side/mode/volume, all 101 lags.
                for li,lag in enumerate(range(-50,51)):
                    ids=np.arange(500);target=ids+lag;ok=(target>=0)&(target<500);ids=ids[ok];target=target[ok]
                    aa=x[ids,2];bb=y[target,1];rr=scalar(aa,bb);nn=int((np.isfinite(aa)&np.isfinite(bb)).sum())
                    er=abs(rr-cube[li,2,1]);assert nn==ns[li] and (er<1e-12 or not np.isfinite(rr))
                    spotchecks.append(meta|dict(side=side,mode=mode,lag=lag,source_bin=3,tail_bin=2,n_pairs=nn,absolute_error=er))
        control=real.copy()
        for side in ['left','right']:control['pseudo_'+side+'_rho']=control_maps[side].rho.values if side in control_maps else np.nan
        control['control_rho']=control[['pseudo_left_rho','pseudo_right_rho']].median(axis=1,skipna=True)
        control['control_excess_rho']=control.rho-control.control_rho
        control['in_real_high_region']=[(a,b) in realkeys for a,b in zip(control.source_bin,control.tail_bin)]
        excess.append(control);h=control[control.in_real_high_region];both=h[['pseudo_left_rho','pseudo_right_rho']].notna().all(axis=1)
        neg=meta|dict(real_map_median=real.rho.median(),real_high_median=realhigh.rho.median(),real_peak=real.rho.max(),
           real_map_specificity=real.zero_lag_specificity.median(),real_high_specificity=realhigh.zero_lag_specificity.median(),
           real_peak_pm2_count=int(real.peak_within_pm2.sum()),control_map_median=control.control_rho.median(),
           control_high_median=h.control_rho.median(),excess_map_median=control.control_excess_rho.median(),
           excess_high_median=h.control_excess_rho.median(),high_cells_both_sides=int(both.sum()),real_high_cell_count=len(h),
           high_fraction_real_gt_both=float(((h.loc[both,'rho']>h.loc[both,'pseudo_left_rho'])&(h.loc[both,'rho']>h.loc[both,'pseudo_right_rho'])).mean()))
        for side,sm in control_maps.items():
            mask=np.array([(a,b) in realkeys for a,b in zip(sm.source_bin,sm.tail_bin)])
            neg[side+'_map_median']=sm.rho.median();neg[side+'_high_at_real_region']=sm.loc[mask,'rho'].median();neg[side+'_specificity']=sm.zero_lag_specificity.median();neg[side+'_similarity']=rho(real.rho,sm.rho)
        negative.append(neg)
        peak=real.sort_values(['rho','source_bin','tail_bin'],ascending=[False,True,True]).iloc[0]
        vr.update(negative_control_rho=neg['control_map_median'],real_high_rho=neg['real_high_median'],real_map_rho=neg['real_map_median'],
           real_peak_rho=neg['real_peak'],real_specificity=neg['real_map_specificity'],real_high_specificity=neg['real_high_specificity'],excess_high_rho=neg['excess_high_median'],
           high_source_u=realhigh.u_mid.median(),high_tail_depth_um=realhigh.tail_mid_um.median(),source_peak_coupling_bin=int(peak.source_bin),
           source_peak_detectability_bin=int(max(range(1,7),key=lambda i:metrics[f's{i}']['z_median'])),
           high_tail_z_median=np.median([metrics[f't{int(i):02}']['z_median'] for i in realhigh.tail_bin]),
           high_source_z_median=np.median([metrics[f's{int(i)}']['z_median'] for i in realhigh.source_bin]))
        if 'structural_source_sbr' in v:vr['structural_source_sbr']=v.structural_source_sbr.median()
        volumes.append(vr)
    metrics=pd.DataFrame(allmetrics);vol=pd.DataFrame(volumes);neg=pd.DataFrame(negative)
    csv('volume_region_detectability_summary.csv',metrics);csv('volume_qc_summary.csv',vol)
    csv('source_depth_detectability_summary.csv',metrics[metrics.region.str.fullmatch('s[1-6]')]);csv('tail_depth_detectability_summary.csv',metrics[metrics.region.str.fullmatch('t[0-9]{2}')])
    csv('proximal_tail_detectability_summary.csv',metrics[metrics.region.isin(['t01','t02','t03','t04','tail100'])])
    for n,rows in [('volume_geometry_qc.csv',geometry),('geometry_jitter_summary.csv',jitters),('background_roi_coverage.csv',coverage),
        ('volume_background_stability.csv',stability),('roi_effective_support.csv',support),('background_real_tail_correlation.csv',bgcorr),
        ('pseudo_source_tail_coupling_summary.csv',pseudo_summ),('negative_control_volume_summary.csv',negative),('denominator_special_cases.csv',special),
        ('structural_qc_summary.csv',struct),('independent_lag_spotchecks.csv.gz',spotchecks)]:csv(n,pd.DataFrame(rows))
    csv('background_slow_axis_curves.csv.gz',pd.concat(slow,ignore_index=True));csv('pseudo_source_tail_coupling_maps.csv',pd.concat(pseudo,ignore_index=True))
    csv('control_excess_coupling_summary.csv',pd.concat(excess,ignore_index=True))
    # Seven predeclared predictors, four outcomes, 5 independent volumes per diameter; no p-values.
    predictors=['source_z','proximal_z','background_rcv','z_top_jitter','X1_jitter','valid_fraction','negative_control_rho']
    outcomes=['real_high_rho','real_map_rho','real_peak_rho','real_specificity'];assoc=[];common=vol[vol.grid_scope.eq('common_grid')]
    for d in [128,235,285,500,'pooled_20']:
        a=common if d=='pooled_20' else common[common.diameter_um.eq(d)]
        for p in predictors:
            for o in outcomes:assoc.append(dict(diameter=d,predictor=p,outcome=o,n_volumes=len(a),spearman=rho(a[p],a[o]),
                        interpretation='potentially confounded by diameter; exploratory only' if d=='pooled_20' else 'descriptive n=5; no p-value'))
    csv('qc_vs_coupling_summary.csv',pd.DataFrame(assoc))
    comparisons=[]
    for d in [128,235,285]:
        a=common[common.diameter_um.eq(500)].set_index('flow_mm_s');b=common[common.diameter_um.eq(d)].set_index('flow_mm_s')
        for metric in predictors+['X4_jitter','background_drift','real_high_rho','excess_high_rho']:
            av,bv=a[metric],b[metric];delta=av-bv
            comparisons.append(dict(reference_diameter=d,metric=metric,d500_median=av.median(),reference_median=bv.median(),
                d500_min=av.min(),d500_max=av.max(),reference_min=bv.min(),reference_max=bv.max(),
                matched_difference_median=delta.median(),n_d500_lower=int((delta<0).sum()),n_d500_higher=int((delta>0).sum()),
                five_flow_d500=json.dumps(av.to_dict()),five_flow_reference=json.dumps(bv.to_dict())))
    csv('d500_quality_vs_other_diameters_summary.csv',pd.DataFrame(comparisons))
    summaries=[]
    for d,a in common.groupby('diameter_um'):
        for col in common.select_dtypes('number'):
            if col in ['diameter_um','flow_mm_s']:continue
            summaries.append(dict(diameter_um=d,metric=col,volume_values=json.dumps(dict(zip(a.scan_id,a[col]))))|summarize(a[col]))
    csv('diameter_qc_summary.csv',pd.DataFrame(summaries))
    # Structural relationship uses same retained cohort only.
    if ex['structural_complete']:
        aa=[]
        for d,g in common.groupby('diameter_um'):
            for p in ['source_z','proximal_z','real_high_rho']:
                aa.append(dict(diameter_um=d,metric=p,structural_metric='structural_source_sbr',spearman=rho(g[p],g.structural_source_sbr),n_volumes=5))
        csv('structural_qc_relationship.csv',pd.DataFrame(aa))
    before=json.loads((OUT/'frozen_files_before.json').read_text());unchanged=all(digest(ROOT/p)==h for p,h in before.items());assert unchanged
    effective=pd.DataFrame(support);assert (effective['min'].dropna()>0).all()
    val=ex|dict(status='passed',frozen_files_unchanged=unchanged,n_frozen_files_checked=len(before),
       independent_lag_spotcheck_count=len(spotchecks),independent_lag_max_absolute_error=pd.DataFrame(spotchecks).absolute_error.max(),
       independent_original_index_pair_counts_passed=True,background_masks_nonoverlapping=ex['max_overlap_weight']==0,
       effective_support_positive=True,no_volume_or_frame_excluded_by_qc=True,nominal_rows=len(data),valid_rows=int(data.valid_geometry.eq(True).sum()),
       background_coverage=pd.DataFrame(coverage).to_dict('records'),formal_coupling_recomputed=False,p_values_computed=False)
    js('validation.json',val)
    print('Analysis gates passed',flush=True)

if __name__=='__main__':main()
