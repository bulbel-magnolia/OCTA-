#!/usr/bin/env python3
"""Independent scalar/coordinate-join and raw-array checks for the depth analysis."""
import hashlib
import importlib.util
import io
import json
import subprocess
import sys
import zipfile
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.io import loadmat

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[1]
sys.path.insert(0,str(ROOT/'src'))
from svrecttail.geometry import VesselGeometry,rectangle_weights
read=lambda name:pd.read_csv(OUT/name,float_precision='round_trip')


def digest(path):
    with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def spearman(x,y):
    a,b=pd.Series(np.asarray(x,float)),pd.Series(np.asarray(y,float))
    m=a.notna()&b.notna()
    return float(a[m].rank(method='average').corr(b[m].rank(method='average')))


def main():
    validation=json.loads((OUT/'validation.json').read_text(encoding='utf-8'))
    assert validation['status']=='passed'
    dfs=[]
    for name,h in validation['derived_file_sha256'].items():
        assert digest(OUT/name)==h
        dfs.append(read(name))
    data=pd.concat(dfs,ignore_index=True)
    lag=read('lag_peak_and_specificity_all_modes.csv.gz')
    # Independently replay 8 formal images: one early and one late valid frame per diameter.
    # Full-image historical band/rectangle weights are independent of cropped optimized integration.
    identities=read('array_identity_audit.csv.gz').set_index(['scan_id','frame_index'])
    spec=importlib.util.spec_from_file_location('frozenbands',ROOT/'analysis/sv_diameter_stage2_scientific_v1/d500_source_axial_local_array_audit.py')
    old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
    qs=[];areas=[];replayed=0
    for diameter in [128,235,285,500]:
        g=data[data.diameter_um.eq(diameter)&data.flow_mm_s.eq(5)].sort_values('frame_index')
        for r in [next(g.head(1).itertuples()),next(g.tail(1).itertuples())]:
            identity=identities.loc[(r.scan_id,r.frame_index)]
            p=Path(identity.input_path)
            if diameter==128:
                with zipfile.ZipFile(p) as z:blob=z.read(identity.member)
                assert hashlib.sha256(blob).hexdigest()==identity.sha256
                with np.load(io.BytesIO(blob),allow_pickle=False) as f:sv=np.asarray(f['sv_raw'],float)
            else:
                assert digest(p)==identity.sha256
                sv=np.asarray(loadmat(p,variable_names=['sv_raw'],simplify_cells=True)['sv_raw'],float)
            geom=VesselGeometry(r.x_left_edge_px,r.x_right_edge_px,r.z_top,float(diameter),12.7,6.7)
            for n in [4,6,8]:
                for j in range(n):
                    w=old.band_weights(sv.shape,geom,j/n,(j+1)/n)
                    area,q,mean=old.weighted_stats(sv,w,12.7*6.7)
                    prefix=f'source_s{j+1}' if n==6 else f'source_b{n}_s{j+1}'
                    refq=getattr(r,prefix+'_q_raw');refa=getattr(r,prefix+'_area_um2')
                    qs.append(abs(q-refq)/abs(refq));areas.append(abs(area-refa)/abs(refa))
            for coord,n,width in [('t',20,25.),('n',5,diameter/5)]:
                for j in range(n):
                    w=rectangle_weights(sv.shape,x_left_edge_px=r.x_left_edge_px,x_right_edge_px=r.x_right_edge_px,
                        z_top_edge_px=geom.z_bottom_edge_px+j*width/6.7,z_bottom_edge_px=geom.z_bottom_edge_px+(j+1)*width/6.7)
                    area,q,mean=old.weighted_stats(sv,w,12.7*6.7)
                    refq=getattr(r,f'tail_{coord}{j+1:02d}_q_raw');refa=getattr(r,f'tail_{coord}{j+1:02d}_area_um2')
                    qs.append(abs(q-refq)/abs(refq));areas.append(abs(area-refa)/abs(refa))
            replayed+=1
    assert max(qs)<1e-12 and max(areas)<1e-12
    # Every volume and both coordinates: preselected cells and lag offsets; all four modes.
    errors=[];checked=0;cache={}
    for diameter in [128,235,285,500]:
        for scope in ['common_grid']+(['d500_extension'] if diameter==500 else []):
            suffix='_extension' if scope=='d500_extension' else ''
            for coord,nt,prefix in [('absolute',20,'t'),('normalized',5,'n')]:
                full=read(f'lag_{coord}_D{diameter}{suffix}.csv.gz')
                selected=full[full.source_bin.isin([1,3,6])&full.tail_bin.isin([1,nt])&full.lag.isin([-50,-5,0,5,50])]
                for r in selected.itertuples():
                    xcol=f'source_s{r.source_bin}_mean_raw';ycol=f'tail_{prefix}{r.tail_bin:02d}_mean_raw'
                    key=(r.scan_id,r.mode,xcol,ycol)
                    if key not in cache:
                        g=data[data.scan_id.eq(r.scan_id)][['frame_index',xcol,ycol]].copy()
                        if r.mode!='raw':
                            half=int(r.mode.split('_')[-1])//2
                            for col in [xcol,ycol]:
                                original=g[col].copy()
                                g[col]=original-np.array([np.median(original[g.frame_index.between(i-half,i+half)]) for i in g.frame_index])
                        cache[key]=g
                    g=cache[key];a=g[['frame_index',xcol]].copy();a['paired_frame']=a.frame_index+r.lag
                    b=g[['frame_index',ycol]].rename(columns={'frame_index':'paired_frame'})
                    pairs=a.merge(b,on='paired_frame',validate='one_to_one')
                    assert len(pairs)==r.n_pairs
                    errors.append(abs(spearman(pairs[xcol],pairs[ycol])-r.rho));checked+=1
                # Independently check every saved lag-summary row against complete curves.
                summary=lag[lag.diameter_um.eq(diameter)&lag.grid_scope.eq(scope)&lag.coordinate.eq(coord)]
                keyed=full.set_index(['scan_id','mode','source_bin','tail_bin']).sort_index()
                for r in summary.itertuples():
                    c=keyed.loc[(r.scan_id,r.mode,r.source_bin,r.tail_bin)].sort_values('lag')
                    zero=c[c.lag.eq(0)].iloc[0];remote=c[c.lag.abs().between(30,50)].rho.median()
                    expected=c[c.rho.eq(c.rho.max())].assign(abs_lag=lambda x:x.lag.abs()).sort_values(['abs_lag','lag']).iloc[0]
                    assert abs(zero.rho-r.rho_lag0)<1e-12
                    assert abs(zero.rho-remote-r.zero_lag_specificity)<1e-12
                    assert expected.lag==r.lag_at_peak
                    assert int(1+(c.rho>zero.rho).sum())==r.zero_lag_rank
                    assert zero.n_pairs==r.n_pairs_lag0
    assert max(errors)<1e-12
    # Primary partial cells: scalar residualization verifies the block implementation.
    part=read('partial_spearman_all_coordinates_modes.csv.gz')
    part=part[part.source_bin.eq(3)&part.tail_bin.eq(1)&part['mode'].eq('detrended_51')]
    perrors=[]
    for r in part.itertuples():
        pre='t' if r.coordinate=='absolute' else 'n'
        xc='source_s3_mean_raw';yc=f'tail_{pre}01_mean_raw'
        controls=['source_area','z_top'] if r.control=='area_ztop' else ['X1','z_top']
        g=data[data.scan_id.eq(r.scan_id)].set_index('frame_index')[[xc,yc]+controls].reindex(range(500))
        rr=(g-g.rolling(51,center=True,min_periods=1).median()).dropna().rank(method='average')
        c=np.column_stack([np.ones(len(rr)),rr[controls]])
        a=rr[xc].to_numpy();b=rr[yc].to_numpy()
        a=a-c@np.linalg.lstsq(c,a,rcond=None)[0];b=b-c@np.linalg.lstsq(c,b,rcond=None)[0]
        expected=np.corrcoef(a,b)[0,1];perrors.append(abs(expected-r.rho))
    assert max(perrors)<1e-12
    primary=read('volume_absolute_depth_coupling_summary.csv')
    reg=read('volume_high_coupling_region_summary.csv')
    for r in reg[reg.coordinate.eq('absolute')].itertuples():
        g=primary[primary.scan_id.eq(r.scan_id)]
        positive=g[g.rho>0].rho.sort_values(ascending=False)
        target=int(np.ceil(.1*len(positive)));h=g[g.rho>=positive.iloc[target-1]]
        assert len(h)==r.high_actual_count and target==r.high_target_count
        assert abs(h.u_mid.median()-r.u_high_region_median)<1e-12
        assert abs(h.tail_mid_um.median()-r.tail_high_region_median_um)<1e-12
    for coord in ['absolute','normalized']:
        con=read(f'diameter_{coord}_consensus_map.csv')
        volume=read(f'volume_{coord}_depth_coupling_summary.csv')
        volume=volume[volume.grid_scope.eq('common_grid')]
        for r in con.itertuples():
            g=volume[volume.diameter_um.eq(r.diameter_um)&volume.source_bin.eq(r.source_bin)&volume.tail_bin.eq(r.tail_bin)]
            assert len(g)==5 and abs(g.rho.median()-r.rho)<1e-12
            assert (g.rho>0).sum()==r.number_positive
        cross=read(f'diameter_{coord}_map_similarity.csv')
        for r in cross.itertuples():
            a=con[con.diameter_um.eq(r.diameter_a)].sort_values(['source_bin','tail_bin']).rho
            b=con[con.diameter_um.eq(r.diameter_b)].sort_values(['source_bin','tail_bin']).rho
            assert abs(spearman(a,b)-r.map_spearman)<1e-12
    prefix='analysis/sv_source_tail_depth_coupling_v1/'
    changed=subprocess.check_output(['git','diff','--name-only',START,'--'],cwd=ROOT,text=True).splitlines()
    assert all(p.startswith(prefix) for p in changed)
    report=dict(status='passed',raw_array_replays=replayed,raw_Q_comparisons=len(qs),
        raw_Q_max_relative_error=max(qs),raw_area_max_relative_error=max(areas),
        independently_recomputed_lag_rows=checked,lag_max_absolute_error=max(errors),
        all_lag_summary_rows_verified=len(lag),partial_scalar_checks=len(perrors),partial_max_absolute_error=max(perrors),
        high_region_selection_checked=True,common_flow_consensus_and_cross_diameter_similarity_checked=True,
        original_coordinate_pairing_and_edge_medians_checked=True,existing_frozen_files_modified=0)
    (OUT/'independent_verification.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(report,indent=2),flush=True)


START='98af95d4e16f3e117362ae2ae8a3c2600cc36ddc'
if __name__=='__main__':main()
