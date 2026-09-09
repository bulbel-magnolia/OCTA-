#!/usr/bin/env python3
"""Within-volume depth maps, lag curves, controlled maps and volume summaries."""
from __future__ import annotations
import itertools
import json
import math
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from derive_depth_metrics import OUT,ROOT,START,digest,csv,js

ID=['scan_id','diameter_um','flow_mm_s','grid_scope']
MAPKEY=ID+['coordinate','source_bins','mode','metric']
CELL=['source_bin','tail_bin']
MODES=['raw','detrended_31','detrended_51','detrended_101']
COORDS=['absolute','normalized']
VALIDATION={}


def rankunit(x):
    a=rankdata(x,method='average',axis=0)
    a=a-a.mean(axis=0)
    norm=np.linalg.norm(a,axis=0)
    return np.divide(a,norm,out=np.full_like(a,np.nan),where=norm>0)


def block_rho(x,y):
    valid=np.isfinite(x).all(axis=1)&np.isfinite(y).all(axis=1)
    n=int(valid.sum())
    return np.clip(rankunit(x[valid]).T@rankunit(y[valid]),-1,1),n


def partial_block(x,y,controls):
    valid=np.isfinite(x).all(axis=1)&np.isfinite(y).all(axis=1)&np.isfinite(controls).all(axis=1)
    x,y,c=[rankdata(z[valid],axis=0,method='average') for z in [x,y,controls]]
    design=np.column_stack([np.ones(len(c)),c])
    def residual(a):
        r=a-design@np.linalg.lstsq(design,a,rcond=None)[0]
        r=r-r.mean(axis=0)
        norms=np.linalg.norm(r,axis=0)
        allowed=norms>1e-12*np.maximum(np.linalg.norm(a,axis=0),1.)
        return np.divide(r,norms,out=np.full_like(r,np.nan),where=allowed)
    return np.clip(residual(x).T@residual(y),-1,1),int(valid.sum()),int(np.linalg.matrix_rank(design))


def similarity(a,b):
    x,y=np.asarray(a,float).ravel(),np.asarray(b,float).ravel()
    m=np.isfinite(x)&np.isfinite(y)
    if m.sum()<3:return np.nan
    return float(block_rho(x[m,None],y[m,None])[0][0,0])


def src_cols(n,kind):
    return [f'source_s{j}_{kind}' if n==6 else f'source_b{n}_s{j}_{kind}' for j in range(1,n+1)]


def tail_cols(coord,kind):
    prefix,n=('t',20) if coord=='absolute' else ('n',5)
    return [f'tail_{prefix}{j:02d}_{kind}' for j in range(1,n+1)]


def cell_template(meta,coord,n):
    d=meta['diameter_um'];nt=20 if coord=='absolute' else 5
    j=np.repeat(np.arange(n),nt);k=np.tile(np.arange(nt),n)
    width=25. if coord=='absolute' else d/5
    data=pd.DataFrame(dict(source_bin=j+1,tail_bin=k+1,u_lo=j/n,u_hi=(j+1)/n,u_mid=(j+.5)/n,
        tail_lo_um=k*width,tail_hi_um=(k+1)*width,tail_mid_um=(k+.5)*width,
        eta_lo=k*width/d,eta_hi=(k+1)*width/d,eta_mid=(k+.5)*width/d))
    for key,value in meta.items():data[key]=value
    data['coordinate']=coord;data['source_bins']=n
    return data


def lag_cube(x,y):
    cubes=[];ns=[]
    for lag in range(-50,51):
        if lag>0:a,b=x[:-lag],y[lag:]
        elif lag<0:a,b=x[-lag:],y[:lag]
        else:a,b=x,y
        r,n=block_rho(a,b);cubes.append(r);ns.append(n)
    return np.asarray(cubes),np.asarray(ns)


def curve_summary(cube,ns,template,mode):
    values=cube.reshape(101,-1)
    rho0=values[50]
    peak=np.nanmax(values,axis=0)
    # First index in this ordering implements smallest |lag|, then smaller signed lag.
    order=sorted(range(101),key=lambda i:(abs(i-50),i-50))
    at=np.array([next(i for i in order if values[i,j]==peak[j]) for j in range(values.shape[1])])
    far=np.r_[0:21,80:101]
    remote=np.nanmedian(values[far],axis=0)
    s=template.copy();s['mode']=mode;s['metric']='mean'
    s['rho_lag0']=rho0;s['rho_peak']=peak;s['lag_at_peak']=at-50
    s['n_pairs_lag0']=ns[50];s['n_pairs_at_peak']=ns[at]
    s['zero_lag_rank']=1+(values>rho0[None,:]).sum(axis=0)
    s['peak_within_pm2']=np.abs(at-50)<=2;s['peak_within_pm5']=np.abs(at-50)<=5
    s['rho_remote_median']=remote;s['zero_lag_specificity']=rho0-remote
    s['peak_tie_count']=(values==peak[None,:]).sum(axis=0)
    return s


def high_region(g):
    positive=g[g.rho>0].sort_values(['rho','source_bin','tail_bin'],ascending=[False,True,True])
    target=int(math.ceil(.1*len(positive)))
    if target==0:return positive,0
    cutoff=positive.rho.iloc[target-1]
    return positive[positive.rho>=cutoff],target


def describe_map(g):
    peak=g.sort_values(['rho','source_bin','tail_bin'],ascending=[False,True,True]).iloc[0]
    h,target=high_region(g)
    def med(col):return float(h[col].median()) if len(h) and col in h else np.nan
    def minimum(col):return float(h[col].min()) if len(h) else np.nan
    def maximum(col):return float(h[col].max()) if len(h) else np.nan
    return dict(rho_map_median=g.rho.median(),rho_map_min=g.rho.min(),rho_map_max=g.rho.max(),
        n_cells=len(g),n_positive_cells=int((g.rho>0).sum()),n_negative_cells=int((g.rho<0).sum()),
        u_peak=peak.u_mid,source_bin_peak=int(peak.source_bin),tail_bin_peak=int(peak.tail_bin),
        tail_depth_peak_um=peak.tail_mid_um,eta_peak=peak.eta_mid,rho_peak_cell=peak.rho,
        specificity_at_peak=peak.get('zero_lag_specificity',np.nan),lag_at_peak_cell=peak.get('lag_at_peak',np.nan),
        peak_cell_lag_within_pm2=peak.get('peak_within_pm2',np.nan),peak_cell_lag_within_pm5=peak.get('peak_within_pm5',np.nan),
        high_positive_cells_available=len(g[g.rho>0]),high_target_count=target,high_actual_count=len(h),
        u_high_region_median=med('u_mid'),u_high_region_min=minimum('u_mid'),u_high_region_max=maximum('u_mid'),
        u_high_region_lo=minimum('u_lo'),u_high_region_hi=maximum('u_hi'),
        tail_high_region_median_um=med('tail_mid_um'),tail_high_region_min_um=minimum('tail_mid_um'),tail_high_region_max_um=maximum('tail_mid_um'),
        tail_high_region_lo_um=minimum('tail_lo_um'),tail_high_region_hi_um=maximum('tail_hi_um'),
        eta_high_region_median=med('eta_mid'),eta_high_region_lo=minimum('eta_lo'),eta_high_region_hi=maximum('eta_hi'),
        high_region_median_rho=med('rho'),high_region_median_specificity=med('zero_lag_specificity'),
        high_region_peak_pm2_count=int(h.peak_within_pm2.sum()) if 'peak_within_pm2' in h else np.nan,
        high_region_peak_pm5_count=int(h.peak_within_pm5.sum()) if 'peak_within_pm5' in h else np.nan,
        high_region_positive_specificity_count=int((h.zero_lag_specificity>0).sum()) if 'zero_lag_specificity' in h else np.nan)


def compute_maps(data):
    maps=[];summaries=[];partials=[];diags=[];full_lag_buffers={}
    lag_min=500;lag_max=0;zero_maxerr=0.;rowcount=0
    for key,g in data.groupby(ID,sort=True):
        meta=dict(zip(ID,key));print('Coupling maps '+meta['scan_id'],flush=True)
        columns=[c for c in g if c.endswith(('_mean_raw','_q_raw','_area_um2'))]+['source_area','X1','z_top']
        raw=g.set_index('frame_index')[columns].reindex(range(500))
        modes={'raw':raw,**{f'detrended_{w}':raw-raw.rolling(w,center=True,min_periods=1).median() for w in [31,51,101]}}
        for mode,d in modes.items():
            for coord in COORDS:
                nt=20 if coord=='absolute' else 5
                for n in [4,6,8]:
                    template=cell_template(meta,coord,n)
                    for metric in ['mean','q']:
                        x=d[src_cols(n,metric+'_raw')].to_numpy();y=d[tail_cols(coord,metric+'_raw')].to_numpy()
                        r,npairs=block_rho(x,y)
                        z=template.copy();z['mode']=mode;z['metric']=metric;z['rho']=r.ravel();z['n_pairs']=npairs
                        maps.append(z)
                        if n==6 and metric=='mean':
                            cube,ns=lag_cube(x,y)
                            zero_maxerr=max(zero_maxerr,float(np.max(np.abs(cube[50]-r))))
                            if not np.isfinite(cube).all():raise AssertionError('undefined primary mean lag cell')
                            if not np.all(ns<=500-np.abs(np.arange(-50,51))):raise AssertionError('lag counts')
                            lag_min=min(lag_min,int(ns.min()));lag_max=max(lag_max,int(ns.max()))
                            summaries.append(curve_summary(cube,ns,template,mode))
                            # Compact full lag payload; coordinates are defined in explicit lookup CSVs.
                            nl=101*n*nt
                            lf=pd.DataFrame(dict(scan_id=np.repeat(meta['scan_id'],nl),mode=mode,
                                source_bin=np.tile(np.repeat(np.arange(1,n+1),nt),101),
                                tail_bin=np.tile(np.arange(1,nt+1),101*n),
                                lag=np.repeat(np.arange(-50,51),n*nt),rho=cube.ravel(),n_pairs=np.repeat(ns,n*nt)))
                            filekey=(int(meta['diameter_um']),meta['grid_scope'],coord)
                            full_lag_buffers.setdefault(filekey,[]).append(lf);rowcount+=len(lf)
                            for control,cc in [('area_ztop',['source_area','z_top']),('x1_ztop',['X1','z_top'])]:
                                pr,pn,rank=partial_block(x,y,d[cc].to_numpy())
                                pp=template.copy();pp['mode']=mode;pp['metric']='mean';pp['control']=control
                                pp['rho']=pr.ravel();pp['n_pairs']=pn;pp['control_design_rank']=rank
                                partials.append(pp)
                    if n==6:
                        tq=d[tail_cols(coord,'q_raw')].to_numpy()
                        for label,cols in [('source_band_area',src_cols(6,'area_um2')),('X1',['X1']),('source_area',['source_area'])]:
                            rr,npairs=block_rho(d[cols].to_numpy(),tq)
                            dd=template.copy() if len(cols)==6 else template[template.source_bin.eq(1)].copy()
                            if len(cols)==1:
                                dd['source_bin']=0;dd[['u_lo','u_hi','u_mid']]=np.nan
                            dd['mode']=mode;dd['predictor']=label;dd['rho']=rr.ravel();dd['n_pairs']=npairs
                            diags.append(dd)
        # Flush completed diameter scopes later; buffers store only 7 compact columns.
    for (diameter,scope,coord),frames in full_lag_buffers.items():
        suffix='_extension' if scope=='d500_extension' else ''
        csv(f'lag_{coord}_D{diameter}{suffix}.csv.gz',pd.concat(frames,ignore_index=True))
    maps=pd.concat(maps,ignore_index=True);lag=pd.concat(summaries,ignore_index=True)
    part=pd.concat(partials,ignore_index=True);diag=pd.concat(diags,ignore_index=True)
    csv('coupling_maps_all_resolutions_modes.csv.gz',maps)
    csv('lag_peak_and_specificity_all_modes.csv.gz',lag)
    csv('lag_peak_and_specificity_summary.csv',lag[lag['mode'].eq('detrended_51')])
    csv('partial_spearman_all_coordinates_modes.csv.gz',part)
    csv('q_geometry_diagnostic_cells.csv.gz',diag)
    for control in ['area_ztop','x1_ztop']:
        csv(f'partial_spearman_{control}_6x20.csv',part[part.coordinate.eq('absolute')&part['mode'].eq('detrended_51')&part.control.eq(control)])
    VALIDATION.update(full_lag_rows=rowcount,full_lag_expected_rows=23*6*25*101*4,
                      lag_pair_count_min=lag_min,lag_pair_count_max=lag_max,
                      lag0_vs_direct_map_max_absolute_error=zero_maxerr,
                      primary_maps_finite=bool(np.isfinite(maps.rho).all()),partial_maps_finite=bool(np.isfinite(part.rho).all()))
    assert rowcount==23*6*25*101*4 and zero_maxerr<1e-12
    assert VALIDATION['primary_maps_finite'] and VALIDATION['partial_maps_finite']
    return maps,lag,part,diag


def attach_lag(maps,lag):
    fields=['rho_peak','lag_at_peak','n_pairs_lag0','zero_lag_rank','peak_within_pm2','peak_within_pm5',
            'rho_remote_median','zero_lag_specificity']
    return maps.merge(lag[MAPKEY+CELL+fields],on=MAPKEY+CELL,how='left',validate='one_to_one')


def consensus(maps,with_lag=False):
    dims=['diameter_um','coordinate','source_bins','mode','metric']+CELL
    coords=['u_lo','u_hi','u_mid','tail_lo_um','tail_hi_um','tail_mid_um','eta_lo','eta_hi','eta_mid']
    rows=[]
    for key,g in maps[maps.grid_scope.eq('common_grid')].groupby(dims):
        g=g.sort_values('flow_mm_s');r=dict(zip(dims,key))|{c:g[c].iloc[0] for c in coords}
        r.update(rho=g.rho.median(),rho_min=g.rho.min(),rho_max=g.rho.max(),number_positive=int((g.rho>0).sum()),
                 n_volumes=len(g),positive_5_of_5=bool((g.rho>0).sum()==5),positive_at_least_4_of_5=bool((g.rho>0).sum()>=4),
                 flow_rho_values_json=json.dumps([dict(flow_mm_s=float(x.flow_mm_s),rho=float(x.rho)) for x in g.itertuples()]))
        if with_lag:
            r.update(zero_lag_specificity=g.zero_lag_specificity.median(),specificity_min=g.zero_lag_specificity.min(),specificity_max=g.zero_lag_specificity.max(),
                number_positive_specificity=int((g.zero_lag_specificity>0).sum()),number_peak_within_pm2=int(g.peak_within_pm2.sum()),
                number_peak_within_pm5=int(g.peak_within_pm5.sum()))
        assert len(g)==5
        rows.append(r)
    return pd.DataFrame(rows)


def region_records(maps,keys):
    rows=[];members=[]
    for key,g in maps.groupby(keys):
        info=dict(zip(keys,key));rows.append(info|describe_map(g))
        h,target=high_region(g)
        members.append(h.assign(high_region_target_count=target))
    return pd.DataFrame(rows),pd.concat(members,ignore_index=True)


def matched_vector(g,n=None):
    return g.sort_values(CELL).rho.to_numpy().reshape(int(g.source_bins.iloc[0]),-1)


def resolution_similarity(g,ref):
    a,b=matched_vector(g),matched_vector(ref)
    common=math.lcm(a.shape[0],b.shape[0])
    # Only a map-shape diagnostic: exact constant-cell replication on common source intervals.
    return similarity(np.repeat(a,common//a.shape[0],axis=0),np.repeat(b,common//b.shape[0],axis=0))


def summarize(maps,lag,part,diag):
    enriched=attach_lag(maps,lag)
    main=enriched[enriched.source_bins.eq(6)&enriched['mode'].eq('detrended_51')&enriched.metric.eq('mean')].copy()
    maincons=consensus(main,True)
    allcons=consensus(maps)
    csv('diameter_consensus_all_resolutions_modes.csv.gz',allcons)
    regions,members=region_records(main,MAPKEY)
    csv('volume_high_coupling_region_summary.csv',regions)
    csv('volume_high_coupling_region_cells.csv',members)
    consensus_regions,consensus_members=region_records(maincons,['diameter_um','coordinate','source_bins','mode','metric'])
    csv('diameter_consensus_high_region_summary.csv',consensus_regions)
    csv('diameter_consensus_high_region_cells.csv',consensus_members)
    for coord in COORDS:
        mc=maincons[maincons.coordinate.eq(coord)]
        csv(f'diameter_{coord}_consensus_map.csv',mc)
        csv(f'diameter_{coord}_sign_consistency.csv',mc[['diameter_um']+CELL+['number_positive','positive_5_of_5','positive_at_least_4_of_5','number_positive_specificity','number_peak_within_pm2','number_peak_within_pm5']])
        csv(f'volume_{coord}_depth_coupling_summary.csv',main[main.coordinate.eq(coord)])
        csv(f'volume_{coord}_depth_peak_summary.csv',regions[regions.coordinate.eq(coord)])
    # Flow map shape: ten pairs per diameter; never pool spatial frames.
    flow=[]
    for (d,coord),g in main[main.grid_scope.eq('common_grid')].groupby(['diameter_um','coordinate']):
        for a,b in itertools.combinations(sorted(g.flow_mm_s.unique()),2):
            flow.append(dict(diameter_um=d,coordinate=coord,flow_a=a,flow_b=b,
                map_spearman=similarity(matched_vector(g[g.flow_mm_s.eq(a)]),matched_vector(g[g.flow_mm_s.eq(b)])),n_cells=len(g)//5))
    flow=pd.DataFrame(flow);csv('flow_map_similarity.csv',flow)
    cross=[]
    for coord,g in maincons.groupby('coordinate'):
        for a,b in itertools.combinations([128,235,285,500],2):
            cross.append(dict(coordinate=coord,diameter_a=a,diameter_b=b,map_spearman=similarity(
                matched_vector(g[g.diameter_um.eq(a)]),matched_vector(g[g.diameter_um.eq(b)])),n_cells=120 if coord=='absolute' else 30))
    cross=pd.DataFrame(cross)
    for coord in COORDS:csv(f'diameter_{coord}_map_similarity.csv',cross[cross.coordinate.eq(coord)])
    crosscompare=cross[cross.coordinate.eq('absolute')].merge(cross[cross.coordinate.eq('normalized')],on=['diameter_a','diameter_b'],suffixes=('_absolute','_normalized'))
    crosscompare['normalized_minus_absolute']=crosscompare.map_spearman_normalized-crosscompare.map_spearman_absolute
    csv('absolute_vs_normalized_similarity.csv',crosscompare)
    csv('cross_diameter_similarity_summary.csv',cross.groupby('coordinate').map_spearman.agg(['median','min','max','count']).reset_index())
    # Flow region spread, signs and specificity are objective quantities, no subjective grades.
    robustness=[]
    fields=['u_peak','tail_depth_peak_um','eta_peak','u_high_region_median','tail_high_region_median_um',
            'u_high_region_lo','u_high_region_hi','tail_high_region_lo_um','tail_high_region_hi_um',
            'rho_peak_cell','high_region_median_rho','high_region_median_specificity']
    for (d,coord),g in regions[regions.grid_scope.eq('common_grid')].groupby(['diameter_um','coordinate']):
        c=maincons[maincons.diameter_um.eq(d)&maincons.coordinate.eq(coord)]
        f=flow[flow.diameter_um.eq(d)&flow.coordinate.eq(coord)]
        r=dict(diameter_um=d,coordinate=coord,n_volumes=5,
            median_flow_map_similarity=f.map_spearman.median(),min_flow_map_similarity=f.map_spearman.min(),max_flow_map_similarity=f.map_spearman.max(),
            positive_5_of_5_cells=int(c.positive_5_of_5.sum()),positive_at_least_4_of_5_cells=int(c.positive_at_least_4_of_5.sum()),
            specificity_positive_5_of_5_cells=int(c.number_positive_specificity.eq(5).sum()),
            peak_cells_near_pm2_volumes=int(g.peak_cell_lag_within_pm2.sum()),peak_cells_near_pm5_volumes=int(g.peak_cell_lag_within_pm5.sum()),
            peak_cells_positive_specificity_volumes=int((g.specificity_at_peak>0).sum()))
        for field in fields:
            r[field+'_median']=g[field].median();r[field+'_min']=g[field].min();r[field+'_max']=g[field].max()
        robustness.append(r)
    csv('flow_robustness_summary.csv',pd.DataFrame(robustness))
    # Sensitivity summaries keep every volume's location; consensus summaries are separate.
    sensitivity=[]
    for key,g in maps[(maps.metric.eq('mean'))].groupby(MAPKEY):
        meta=dict(zip(MAPKEY,key));ref=main[(main.scan_id.eq(meta['scan_id']))&(main.coordinate.eq(meta['coordinate']))]
        if meta['mode']=='detrended_51':
            sensitivity.append(meta|describe_map(g)|dict(comparison='source_bins',map_similarity_to_6bin_51=resolution_similarity(g,ref)))
        elif meta['source_bins']==6:
            sensitivity.append(meta|describe_map(g)|dict(comparison='detrend_window',map_similarity_to_6bin_51=similarity(matched_vector(g),matched_vector(ref))))
    sensitivity=pd.DataFrame(sensitivity)
    csv('source_bin_sensitivity_summary.csv',sensitivity[sensitivity.comparison.eq('source_bins')])
    csv('detrend_sensitivity_summary.csv',sensitivity[sensitivity.comparison.eq('detrend_window')])
    cr,_=region_records(allcons[allcons['mode'].eq('detrended_51')&allcons.metric.eq('mean')],['diameter_um','coordinate','source_bins','mode','metric'])
    csv('source_bin_consensus_region_sensitivity.csv',cr)
    # Both geometry controls, both coordinate systems; raw and residual results are retained.
    gc=[];pc=[]
    for control,g in part.groupby('control'):
        c=consensus(g);c['control']=control;pc.append(c)
        for key,p in g[g['mode'].eq('detrended_51')].groupby(MAPKEY):
            meta=dict(zip(MAPKEY,key));ref=main[main.scan_id.eq(meta['scan_id'])&main.coordinate.eq(meta['coordinate'])]
            ph,_=high_region(p);rh,_=high_region(ref)
            a=set(map(tuple,ph[CELL].values));b=set(map(tuple,rh[CELL].values))
            gc.append(meta|dict(control=control,map_similarity_to_uncontrolled=similarity(matched_vector(p),matched_vector(ref)),
                high_region_jaccard=len(a&b)/len(a|b) if a|b else np.nan)|describe_map(p))
    pc=pd.concat(pc,ignore_index=True);gc=pd.DataFrame(gc)
    csv('geometry_control_summary.csv',gc)
    csv('diameter_geometry_control_consensus.csv',pc[pc['mode'].eq('detrended_51')])
    csv('diameter_geometry_control_all_modes.csv.gz',pc)
    pr,_=region_records(pc[pc['mode'].eq('detrended_51')],['diameter_um','coordinate','source_bins','mode','metric','control'])
    csv('geometry_control_consensus_region_summary.csv',pr)
    # Geometry/Q summaries are medians of map coefficients, not pooled-frame coefficients.
    qg=[]
    for key,g in diag.groupby(ID+['coordinate','mode','predictor']):
        qg.append(dict(zip(ID+['coordinate','mode','predictor'],key))|
                  dict(rho_median=g.rho.median(),rho_min=g.rho.min(),rho_max=g.rho.max(),positive_cells=int((g.rho>0).sum()),n_cells=len(g)))
    for key,g in maps[maps.source_bins.eq(6)&maps.metric.eq('q')].groupby(MAPKEY):
        qg.append(dict(zip(MAPKEY,key))|dict(predictor='source_band_Q_vs_tail_Q',rho_median=g.rho.median(),rho_min=g.rho.min(),rho_max=g.rho.max(),positive_cells=int((g.rho>0).sum()),n_cells=len(g)))
    csv('q_geometry_diagnostic_summary.csv',pd.DataFrame(qg))
    # A separate extension receipt prevents accidental entry into the 5-flow consensuses.
    csv('d500_extension_volume_summary.csv',regions[regions.grid_scope.eq('d500_extension')])
    VALIDATION.update(flow_similarity_pairs=len(flow),cross_diameter_similarity_pairs=len(cross),
        consensus_volumes_per_cell=int(maincons.n_volumes.min()),common_only_consensus=True,
        primary_absolute_volume_cells=int(len(main[main.coordinate.eq('absolute')])),
        primary_normalized_volume_cells=int(len(main[main.coordinate.eq('normalized')])))
    assert len(flow)==80 and len(cross)==12 and maincons.n_volumes.eq(5).all()


def main():
    v=json.loads((OUT/'data_validation.json').read_text(encoding='utf-8'))
    assert v['status']=='data_gate_passed_analysis_pending'
    files=[]
    for name,h in v['derived_file_sha256'].items():
        assert digest(OUT/name)==h
        files.append(pd.read_csv(OUT/name,float_precision='round_trip'))
    data=pd.concat(files,ignore_index=True)
    maps,lag,part,diag=compute_maps(data)
    summarize(maps,lag,part,diag)
    p=json.loads((OUT/'provenance.json').read_text(encoding='utf-8'))
    p.update(analysis_script_sha256=digest(__file__),analysis_parameters=dict(
        primary='6 source bins x 20 nonoverlapping 25um bands; mean-mean, 51-frame residual, lag0',
        normalized='6x5, direct 0.2D bins over 0-1D',source_bin_sensitivity=[4,8],optional_10_bins='not performed',
        modes=MODES,lag_range=[-50,50],full_lag_scope='6 bins, mean-mean, raw/31/51/101, both tail coordinates',
        detrend='centered median, original 0..499 grid, min_periods=1, truncated edges; invalid center remains NaN',
        spearman='average ranks per pairwise-valid set; no p-values',lag_sign='Source(i) versus Tail(i+lag); no wrap',
        controls='source_area+z_top primary; X1+z_top sensitivity; controls detrended with signals in residual modes',
        high_region='ceil(0.10 * number of positive cells), top signed rho, includes ties at cutoff; empty if no positive cells',
        peak_ties='cell: smallest source_bin then tail_bin; lag: smallest |lag| then smaller signed lag',
        specificity='rho(0) minus median rho at absolute lag 30..50 inclusive',
        flow_similarity='Spearman of aligned map cells, ten pairs per diameter; map cells not independent replicates',
        source_resolution_map_comparison='piecewise constant replication on LCM source rows only for shape diagnostic; no new quantitative SV and no tail-coordinate interpolation',
        Q_status='Q–Q analysis includes signal-intensity and ROI-geometry contributions and is secondary to mean–mean analysis.'))
    js('provenance.json',p)
    v.update(status='passed',coupling_checks=VALIDATION,missing_required_data=[],unfinished_computation=[],
             invalid_frames_zero_filled=False,invalid_frames_interpolated=False,p_values_computed=False,
             experimental_unit='scan volume')
    js('validation.json',v)
    print(json.dumps(VALIDATION,indent=2),flush=True)


if __name__=='__main__':main()
