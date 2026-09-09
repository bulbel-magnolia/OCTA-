#!/usr/bin/env python3
"""Independent coordinate-join audit of published spatial results and file scope."""
import hashlib
import json
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[1]


def digest(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f,'sha256').hexdigest()


def main():
    read=lambda name: pd.read_csv(OUT/name,float_precision='round_trip')
    fw=read('framewise_source_tail_metrics.csv')
    zero=read('volume_detrended_coupling.csv')
    lag=pd.concat([read('lag_correlation_raw.csv'),read('lag_correlation_detrended_51.csv'),
                   read('lag_correlation_detrended_sensitivity.csv')],ignore_index=True)
    validation=json.loads((OUT/'validation.json').read_text())
    assert validation['status']=='passed'
    manifest=read('input_manifest.csv')
    assert all(digest(ROOT/r.input_path)==r.sha256 for r in manifest.itertuples())
    assert digest(OUT/'framewise_source_tail_metrics.csv')==validation['framewise_sha256']
    a=pd.read_csv(ROOT/'analysis/sv_diameter_stage2_scientific_v1/d128_d235_source_axial_audit/d128_d235_source_axial_framewise.csv')
    b=pd.read_csv(ROOT/'analysis/sv_diameter_stage2_scientific_v1/d500_source_axial_audit/source_axial_band_framewise.csv')
    a=a.rename(columns={'input_sha256':'band_array_sha256'})
    b=b.rename(columns={'input_mat_sha256':'band_array_sha256'})
    hashes=pd.concat([a[['scan_id','frame_index_0based','band_array_sha256']],
                      b[['scan_id','frame_index_0based','band_array_sha256']]]).rename(columns={'frame_index_0based':'frame_index'})
    identity=fw.merge(hashes,on=['scan_id','frame_index'],validate='one_to_one')
    assert len(identity)==10931 and identity.input_sha256.eq(identity.band_array_sha256).all()
    # Independent calculation on actual frame coordinates, including edge windows.
    # A stable, predetermined subset checks every volume, mode, tail depth and signal family.
    selected=lag[lag.lag.isin([-50,-31,-5,0,5,31,50]) &
                 lag.predictor.isin(['source_mean_raw','lower_q_raw'])]
    cache={}
    errors=[]
    for r in selected.itertuples():
        key=(r.scan_id,r.mode,r.predictor,r.tail_variable)
        if key not in cache:
            g=fw[fw.scan_id.eq(r.scan_id)].copy()
            cols=[r.predictor,r.tail_variable]
            if r.mode!='raw':
                w=int(r.mode.split('_')[-1]); half=w//2
                for col in cols:
                    original=g[col].copy()
                    # Uses coordinate-distance predicates, not rolling or compressed row position.
                    baseline=np.array([np.median(original[g.frame_index.between(i-half,i+half)]) for i in g.frame_index])
                    g[col]=original-baseline
            cache[key]=g[['frame_index']+cols]
        g=cache[key]
        left=g[['frame_index',r.predictor]].copy()
        left['tail_frame']=left.frame_index+r.lag
        right=g[['frame_index',r.tail_variable]].rename(columns={'frame_index':'tail_frame'})
        matched=left.merge(right,on='tail_frame',validate='one_to_one')
        expected=matched[r.predictor].rank(method='average').corr(matched[r.tail_variable].rank(method='average'))
        assert len(matched)==r.n_pairs
        errors.append(abs(expected-r.rho))
    assert max(errors)<1e-12
    # Every first difference must be an actual adjacent spatial pair.
    diff=read('volume_first_difference_coupling.csv')
    diff_errors=[]
    for r in diff.itertuples():
        g=fw[fw.scan_id.eq(r.scan_id)].set_index('frame_index')
        ids=g.index[g.index.to_series().sub(1).isin(g.index)]
        dx=g.loc[ids,r.predictor].to_numpy()-g.loc[ids-1,r.predictor].to_numpy()
        dy=g.loc[ids,r.tail_variable].to_numpy()-g.loc[ids-1,r.tail_variable].to_numpy()
        expected=pd.Series(dx).rank().corr(pd.Series(dy).rank())
        assert len(ids)==r.n_pairs
        diff_errors.append(abs(expected-r.rho))
    assert max(diff_errors)<1e-12
    peaks=read('lag_peak_summary.csv')
    assert len(peaks)==1472
    assert len(zero)==1380
    assert len(lag)==148672
    assert set(fw[fw.grid_scope.eq('common_grid')].flow_mm_s)=={1,3,5,7,10}
    assert fw[fw.grid_scope.eq('common_grid')].groupby('diameter_um').scan_id.nunique().eq(5).all()
    assert set(fw[fw.grid_scope.eq('d500_extension')].diameter_um)=={500}
    # All input edits must be confined to this newly added analysis directory.
    base='0d1c428d901447f17d0fa600b37dfdfa41e776dd'
    diffpaths=subprocess.check_output(['git','diff','--name-only',base,'--'],cwd=ROOT,text=True).splitlines()
    prefix='analysis/sv_source_tail_framewise_coupling_v1/'
    assert all(p.startswith(prefix) for p in diffpaths),diffpaths
    tracked_before=subprocess.check_output(['git','ls-tree','-r','--name-only',base,'--',prefix],cwd=ROOT,text=True)
    assert not tracked_before.strip()
    report=dict(status='passed',input_manifest_files_rehashed=len(manifest),
        band_array_identities_matched_to_formal_frames=len(identity),
        independently_recomputed_lag_rows=len(selected),lag_max_absolute_error=max(errors),
        independently_recomputed_first_difference_rows=len(diff),first_difference_max_absolute_error=max(diff_errors),
        coordinate_join_not_compressed_frames=True,edge_median_coordinate_neighborhood_verified=True,
        only_new_analysis_directory_changed=True,preexisting_files_modified=0,
        notes='No p-values computed. Additional publication integrity is verified by Git after commit/push.')
    (OUT/'independent_verification.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
