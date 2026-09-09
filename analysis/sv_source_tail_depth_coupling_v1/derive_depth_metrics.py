#!/usr/bin/env python3
"""Direct, SHA-verified subpixel depth integration from retained formal raw SV."""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import io
import json
import platform
import sys
import zipfile
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
from scipy.io import loadmat

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[1]
PREV=ROOT/'analysis/sv_source_tail_framewise_coupling_v1'
START='98af95d4e16f3e117362ae2ae8a3c2600cc36ddc'
sys.path.insert(0,str(ROOT/'src'))
from svrecttail.geometry import VesselGeometry, interval_overlap_weights, ellipse_weights
SOURCE_BINS=[4,6,8]
DX,DZ=12.7,6.7
TOL=1e-10
MANIFEST=[]
CHECKS={}


def digest(path):
    with Path(path).open('rb') as f: return hashlib.file_digest(f,'sha256').hexdigest()


def js(name,obj):
    (OUT/name).write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8',newline='\n')


def csv(name,data):
    opts={'method':'gzip','mtime':0} if name.endswith('.gz') else None
    data.to_csv(OUT/name,index=False,float_format='%.17g',lineterminator='\n',compression=opts)


def record(path,role,h=None):
    path=Path(path).resolve()
    try: name=path.relative_to(ROOT).as_posix()
    except ValueError: name=path.as_posix()
    MANIFEST.append(dict(input_path=name,sha256=h or digest(path),role=role,bytes=path.stat().st_size))


def check(name,condition,**details):
    CHECKS[name]=dict(status='passed' if condition else 'failed',**details)
    if not condition:
        js('validation.json',dict(status='failed',checks=CHECKS))
        raise AssertionError(name)


def error(name,actual,expected):
    a,b=np.asarray(actual,float),np.asarray(expected,float)
    absolute=float(np.max(np.abs(a-b)))
    relative=float(np.max(np.abs(a-b)/np.maximum(np.abs(b),np.finfo(float).tiny)))
    good=np.isfinite(a).all() and np.isfinite(b).all() and relative<TOL
    old=CHECKS.get(name,dict(max_absolute_error=0.,max_relative_error=0.,n_values=0))
    check(name,bool(good),max_absolute_error=max(absolute,old['max_absolute_error']),
          max_relative_error=max(relative,old['max_relative_error']),n_values=old['n_values']+a.size,
          tolerance_relative=TOL)


def source_local_weights(shape,g,nbins):
    """Partition exactly the frozen 16x16 sample lattice; no geometry estimation."""
    xi=np.flatnonzero(interval_overlap_weights(shape[1],g.x_left_edge_px,g.x_right_edge_px)>0)
    zi=np.flatnonzero(interval_overlap_weights(shape[0],g.z_top_edge_px,g.z_bottom_edge_px)>0)
    offsets=(np.arange(16,dtype=float)+.5)/16-.5
    xs=xi[:,None]+offsets
    zs=zi[:,None]+offsets
    xt=((xs-g.x_center_px)*DX/(g.lateral_width_um/2))**2
    zt=((zs-g.z_center_px)*DZ/(g.diameter_um/2))**2
    inside=zt[:,None,:,None]+xt[None,:,None,:]<=1
    counts=inside.sum(axis=3)
    u=(zs-g.z_top_edge_px)*DZ/g.diameter_um
    weights=[]
    for j in range(nbins):
        m=(u>=j/nbins)&((u<(j+1)/nbins) if j<nbins-1 else (u<=1))
        weights.append((counts*m[:,None,:]).sum(axis=2)/256.)
    weights=np.asarray(weights)
    full=inside.mean(axis=(2,3),dtype=np.float64)
    if not np.array_equal(weights.sum(axis=0),full): raise AssertionError('source sample partition')
    return zi,xi,weights,full


def source_stats(sv,g,nbins):
    zi,xi,w,full=source_local_weights(sv.shape,g,nbins)
    pixels=sv[np.ix_(zi,xi)]
    area=w.sum(axis=(1,2))*DX*DZ
    q=(w*pixels[None,:,:]).sum(axis=(1,2))*DX*DZ
    return area,q,q/area


def tail_stats(sv,g,bounds):
    xw=interval_overlap_weights(sv.shape[1],g.x_left_edge_px,g.x_right_edge_px)
    xi=np.flatnonzero(xw>0)
    areas,qs=[],[]
    for lo,hi in zip(bounds[:-1],bounds[1:]):
        zw=interval_overlap_weights(sv.shape[0],g.z_bottom_edge_px+lo/DZ,g.z_bottom_edge_px+hi/DZ)
        zi=np.flatnonzero(zw>0)
        w=np.outer(zw[zi],xw[xi])
        area=w.sum()*DX*DZ
        q=(sv[np.ix_(zi,xi)]*w).sum()*DX*DZ
        areas.append(area);qs.append(q)
    areas,qs=np.asarray(areas),np.asarray(qs)
    return areas,qs,qs/areas


def derive(sv,r):
    g=VesselGeometry(r.x_left_edge_px,r.x_right_edge_px,r.z_top,float(r.diameter_um),DX,DZ)
    rec={c:getattr(r,c) for c in ['scan_id','diameter_um','flow_mm_s','grid_scope','frame_index',
            'valid_geometry','X1','X1_px','X4','z_top','z_bottom_edge_px','x_left_edge_px',
            'x_right_edge_px','source_area','source_mean_raw','source_q_raw','input_sha256']}
    for n in SOURCE_BINS:
        a,q,m=source_stats(sv,g,n)
        error(f'source_{n}_q_reconstruction',q.sum(),r.source_q_raw)
        error(f'source_{n}_area_reconstruction',a.sum(),r.source_area)
        error(f'source_{n}_mean_reconstruction',q.sum()/a.sum(),r.source_mean_raw)
        if n==6:
            for j,band in enumerate(['upper','middle','lower']):
                error(f'source_6_reconstruct_{band}_q',q[2*j:2*j+2].sum(),getattr(r,f'{band}_q_raw'))
                error(f'source_6_reconstruct_{band}_area',a[2*j:2*j+2].sum(),getattr(r,f'{band}_area_um2'))
        for j in range(n):
            prefix=f'source_s{j+1}' if n==6 else f'source_b{n}_s{j+1}'
            rec.update({prefix+'_area_um2':a[j],prefix+'_q_raw':q[j],prefix+'_mean_raw':m[j]})
    a,q,m=tail_stats(sv,g,np.arange(21)*25.)
    for d,count in [(100,4),(500,20)]:
        error(f'tail{d}_q_reconstruction',q[:count].sum(),getattr(r,f'tail_q_raw_{d}um'))
        error(f'tail{d}_area_reconstruction',a[:count].sum(),getattr(r,f'tail_area_um2_{d}um'))
        error(f'tail{d}_mean_reconstruction',q[:count].sum()/a[:count].sum(),getattr(r,f'tail_mean_raw_{d}um'))
    for j in range(20):
        rec.update({f'tail_t{j+1:02d}_area_um2':a[j],f'tail_t{j+1:02d}_q_raw':q[j],f'tail_t{j+1:02d}_mean_raw':m[j]})
    error('absolute_tail_band_physical_area',a,np.full(20,r.X1*25))
    a,q,m=tail_stats(sv,g,np.arange(6)*float(r.diameter_um)/5)
    fa,fq,fm=tail_stats(sv,g,[0,float(r.diameter_um)])
    error('normalized_tail_0_to_D_q_reconstruction',q.sum(),fq[0])
    error('normalized_tail_0_to_D_area_reconstruction',a.sum(),fa[0])
    error('normalized_tail_band_physical_area',a,np.full(5,r.X1*r.diameter_um/5))
    for j in range(5):
        rec.update({f'tail_n{j+1:02d}_area_um2':a[j],f'tail_n{j+1:02d}_q_raw':q[j],f'tail_n{j+1:02d}_mean_raw':m[j]})
    return rec


def implementation_checks():
    spec=importlib.util.spec_from_file_location('frozen_axial',ROOT/'analysis/sv_diameter_stage2_scientific_v1/d500_source_axial_local_array_audit.py')
    old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
    maxerr=0.
    for diameter in [128.,235.,285.,500.]:
        g=VesselGeometry(31.125,57.75,30.5,diameter,DX,DZ)
        shape=(160,80)
        for n in SOURCE_BINS:
            zi,xi,w,full=source_local_weights(shape,g,n)
            ref=ellipse_weights(shape,g,supersample=16)
            assert np.array_equal(full,ref[np.ix_(zi,xi)])
            for j in range(n):
                bw=old.band_weights(shape,g,j/n,(j+1)/n)
                maxerr=max(maxerr,float(np.max(np.abs(w[j]-bw[np.ix_(zi,xi)]))))
    check('source_partition_exact_match_frozen_sample_implementation',maxerr==0.,max_weight_absolute_error=maxerr,
          tested_diameters=[128,235,285,500],tested_bins=SOURCE_BINS)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--mat-root',type=Path,required=True)
    ap.add_argument('--d235-f01-mat-root',type=Path,required=True)
    ap.add_argument('--d128-package-root',type=Path,required=True)
    args=ap.parse_args()
    js('validation.json',dict(status='running_input_gate'))
    implementation_checks()
    fw=pd.read_csv(PREV/'framewise_source_tail_metrics.csv',float_precision='round_trip')
    previous=json.loads((PREV/'validation.json').read_text(encoding='utf-8'))
    check('previous_framewise_identity',digest(PREV/'framewise_source_tail_metrics.csv')==previous['framewise_sha256'])
    check('scope',len(fw)==10931 and fw.scan_id.nunique()==23 and
          len(fw[fw.grid_scope.eq('common_grid')])==9456 and fw[fw.grid_scope.eq('common_grid')].scan_id.nunique()==20)
    for p in [PREV/'framewise_source_tail_metrics.csv',PREV/'validation.json',PREV/'provenance.json',
              ROOT/'src/svrecttail/geometry.py',ROOT/'analysis/sv_diameter_stage2_scientific_v1/d500_source_axial_local_array_audit.py',
              ROOT/'analysis/sv_diameter_stage2_scientific_v1/all_volume_metrics_with_d500_extension.csv']:
        record(p,'frozen metrics, geometry implementation, or validation reference')
    # Resolve every required frame before integration; a missing array is a hard gate.
    mat_paths={}
    for r in fw[fw.diameter_um.ne(128)].itertuples():
        p=(args.d235_f01_mat_root if r.scan_id=='D235_F01_V01' else args.mat_root/r.scan_id)/f'frame_{r.frame_index:03d}.mat'
        if not p.is_file(): raise FileNotFoundError(p)
        mat_paths[(r.scan_id,r.frame_index)]=p
    check('retained_mat_coverage',len(mat_paths)==8509)
    formal=ROOT/'results/formal_sv_d128_v21_full2500_run001'
    pk=pd.read_csv(formal/'download_packages.csv')
    for p in [formal/'download_packages.csv',formal/'arrays_sha256.csv']:record(p,'frozen release identity list')
    for r in pk.itertuples():
        if not (args.d128_package_root/r.file).is_file():raise FileNotFoundError(r.file)
    rows=[];identities=[]
    for pkg in pk.itertuples():
        p=args.d128_package_root/pkg.file
        h=digest(p)
        check('package_'+pkg.file,h==pkg.sha256 and p.stat().st_size==pkg.bytes)
        record(p,'D128 retained formal release ZIP',h)
        parts=p.stem.split('_');flow=int(parts[-3].replace('flow',''));lo,hi=map(int,parts[-2:])
        g=fw[fw.diameter_um.eq(128)&fw.flow_mm_s.eq(flow)&fw.frame_index.between(lo,hi)]
        print(f'Derive D128 flow {flow}: {lo}-{hi} ({len(g)} frames)',flush=True)
        with zipfile.ZipFile(p) as z:
            for r in g.itertuples():
                member=f'arrays/flow{flow:02d}/frame_{r.frame_index:03d}.npz'
                blob=z.read(member);h=hashlib.sha256(blob).hexdigest()
                if h!=r.input_sha256:raise AssertionError(member)
                with np.load(io.BytesIO(blob),allow_pickle=False) as data: sv=np.asarray(data['sv_raw'],dtype=float)
                if sv.shape!=(351,500) or not np.isfinite(sv).all():raise AssertionError(member)
                rows.append(derive(sv,r))
                identities.append(dict(scan_id=r.scan_id,frame_index=r.frame_index,input_path=p.as_posix(),member=member,
                                       sha256=h,identity_matched=True))
    for scan,g in fw[fw.diameter_um.ne(128)].groupby('scan_id'):
        print(f'Derive {scan}: {len(g)} retained formal MAT frames',flush=True)
        for r in g.itertuples():
            p=mat_paths[(scan,r.frame_index)];h=digest(p)
            if h!=r.input_sha256:raise AssertionError(f'MAT SHA mismatch {p}')
            record(p,'retained formal MAT sv_raw',h)
            sv=np.asarray(loadmat(p,variable_names=['sv_raw'],simplify_cells=True)['sv_raw'],float)
            if sv.shape!=(351,500) or not np.isfinite(sv).all():raise AssertionError(p)
            rows.append(derive(sv,r))
            identities.append(dict(scan_id=scan,frame_index=r.frame_index,input_path=p.as_posix(),member='',sha256=h,identity_matched=True))
    data=pd.DataFrame(rows).sort_values(['diameter_um','flow_mm_s','frame_index']).reset_index(drop=True)
    check('all_array_identities_verified',len(data)==10931 and len(identities)==10931)
    check('all_depth_values_finite',np.isfinite(data.select_dtypes(include='number')).all().all())
    check('frozen_frame_keys_preserved',set(zip(data.scan_id,data.frame_index))==set(zip(fw.scan_id,fw.frame_index)))
    merged=data.merge(fw,on=['scan_id','frame_index'],validate='one_to_one',suffixes=('','_ref'))
    for c in ['source_mean_raw','source_q_raw','source_area','X1','X4','z_top']:
        error('frozen_'+c,merged[c],merged[c+'_ref'])
    # Rebuild source/tail medians from bins and compare all prior Stage 2 volume values.
    ref=pd.read_csv(ROOT/'analysis/sv_diameter_stage2_scientific_v1/all_volume_metrics_with_d500_extension.csv')
    vr=[]
    for r in ref.itertuples():
        g=data[data.scan_id.eq(r.scan_id)]
        source=g[[f'source_s{i}_q_raw' for i in range(1,7)]].sum(axis=1)/g[[f'source_s{i}_area_um2' for i in range(1,7)]].sum(axis=1)
        vals={'source_mean_raw_median':source.median()}
        for depth,n in [(100,4),(500,20)]:
            tail=g[[f'tail_t{i:02d}_q_raw' for i in range(1,n+1)]].sum(axis=1)/g[[f'tail_t{i:02d}_area_um2' for i in range(1,n+1)]].sum(axis=1)
            vals[f'tail_mean_raw_{depth}um_median']=tail.median()
            if depth==500:vals['ri_tail_median']=(tail/source).median()
        for k,v in vals.items():vr.append(dict(scan_id=r.scan_id,metric=k,recomputed=v,reference=getattr(r,k)))
    vr=pd.DataFrame(vr);error('stage2_volume_medians',vr.recomputed,vr.reference)
    csv('stage2_volume_validation.csv',vr)
    outputs={}
    for (d,scope),g in data.groupby(['diameter_um','grid_scope']):
        name=f'framewise_depth_metrics_D{int(d)}'+('_extension' if scope=='d500_extension' else '')+'.csv.gz'
        csv(name,g);outputs[name]=digest(OUT/name)
    csv('array_identity_audit.csv.gz',pd.DataFrame(identities))
    csv('input_manifest.csv',pd.DataFrame(MANIFEST).drop_duplicates('input_path'))
    source=[];tail=[]
    for n in SOURCE_BINS:
        for j in range(n):source.append(dict(source_bins=n,source_bin=j+1,u_lo=j/n,u_hi=(j+1)/n,u_mid=(j+.5)/n))
    for d in [128,235,285,500]:
        for coord,bounds in [('absolute',np.arange(21)*25),('normalized',np.arange(6)*d/5)]:
            for k,(lo,hi) in enumerate(zip(bounds[:-1],bounds[1:])):
                tail.append(dict(diameter_um=d,coordinate=coord,tail_bin=k+1,tail_lo_um=lo,tail_hi_um=hi,tail_mid_um=(lo+hi)/2,
                                 eta_lo=lo/d,eta_hi=hi/d,eta_mid=(lo+hi)/2/d))
    csv('source_bin_coordinates.csv',pd.DataFrame(source));csv('tail_band_coordinates.csv',pd.DataFrame(tail))
    coverage=data.groupby('grid_scope').agg(volumes=('scan_id','nunique'),frames=('scan_id','size')).to_dict('index')
    val=dict(status='data_gate_passed_analysis_pending',checks=CHECKS,coverage=coverage,derived_file_sha256=outputs,
             source_bins=SOURCE_BINS,optional_10_bins='not performed; 4/6/8 are fixed primary and sensitivity',
             identities_verified=dict(d128_npz=2422,formal_mat=8509,release_zip=25),
             no_background_subtraction=True,no_oct_reconstruction=True)
    js('data_validation.json',val);js('validation.json',val)
    js('provenance.json',dict(starting_head=START,analysis_code_commit='Git commit containing this new directory',
       repository='bulbel-magnolia/OCTA-',branch='analysis/sv-diameter-stage2',
       extraction_script_sha256=digest(__file__),input_manifest_sha256=digest(OUT/'input_manifest.csv'),
       software=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,matlab='not used; read retained formal MAT'),
       signal='var(abs(E),1,3), denominator N; linear raw SV',
       geometry=dict(frozen='continuity-first v2.1',source='unchanged X1 x physical diameter ellipse, centered X4',dx_um=DX,dz_um=DZ,supersample=16,
                     tail='frozen X1 rectangle below true physical bottom, guard 0'),
       source_bins=SOURCE_BINS,primary_source_bins=6,absolute_tail_edges_um=list(range(0,501,25)),normalized_tail_edges_D=[0,.2,.4,.6,.8,1],
       normalized_tail_calculation='direct integration at physical boundaries; no map interpolation',
       experimental_unit='scan volume; no p-values, no pooled frames'))
    print(json.dumps(dict(status='data_gate_passed',coverage=coverage,source6=CHECKS['source_6_q_reconstruction'],tail500=CHECKS['tail500_q_reconstruction']),indent=2),flush=True)


if __name__=='__main__':main()
