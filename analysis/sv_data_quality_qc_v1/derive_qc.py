"""Read-only retained-array audit and deterministic spatially matched QC."""
from __future__ import annotations
import hashlib, io, json, platform, sys, zipfile, subprocess
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
from scipy.io import loadmat

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[1]
DEP=ROOT/'analysis/sv_source_tail_depth_coupling_v1'
PREV=ROOT/'analysis/sv_source_tail_framewise_coupling_v1'
sys.dont_write_bytecode=True
sys.path.insert(0,str(DEP))
from derive_depth_metrics import source_local_weights, VesselGeometry, interval_overlap_weights, DX, DZ
from svrecttail.geometry import interval_is_complete, ellipse_weights, rectangle_weights
ID=['scan_id','diameter_um','flow_mm_s','grid_scope']
REGIONS=['source']+[f's{i}' for i in range(1,7)]+[f't{i:02}' for i in range(1,21)]+['tail100','tail500']
STATS=['mean','median','mad','sd','area','q','neff']

def digest(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def clean(x):
    if isinstance(x,dict):return {str(k):clean(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [clean(v) for v in x]
    if isinstance(x,np.generic):return clean(x.item())
    if isinstance(x,float) and not np.isfinite(x):return None
    return x
def js(n,x):(OUT/n).write_text(json.dumps(clean(x),ensure_ascii=False,indent=2)+'\n',encoding='utf-8',newline='\n')
def csv(n,x):x.to_csv(OUT/n,index=False,float_format='%.17g',lineterminator='\n',compression={'method':'gzip','mtime':0} if n.endswith('.gz') else None)
def read(p):return pd.read_csv(p,float_precision='round_trip',low_memory=False)
def wmedian(v,w):
    order=np.argsort(v,kind='stable');v=v[order];w=w[order];c=np.cumsum(w)
    k=np.searchsorted(c,c[-1]/2,side='left')
    return float((v[k]+v[k+1])/2 if c[k]==c[-1]/2 and k+1<len(v) else v[k])
def stats(v,w):
    m=(w>0);v=v[m];w=w[m]
    if not len(w):return dict.fromkeys(STATS,np.nan)
    assert np.isfinite(v).all()
    sw=w.sum();mu=np.sum(v*w)/sw;med=wmedian(v,w)
    return dict(mean=mu,median=med,mad=wmedian(np.abs(v-med),w),sd=np.sqrt(np.sum(w*(v-mu)**2)/sw),
                area=sw*DX*DZ,q=np.sum(v*w)*DX*DZ,neff=sw**2/np.sum(w*w))
def contrast(real,bg):
    m,b=real['mean'],bg['mean'];den=1.4826*bg['mad']
    sbr=m/b if b>0 else np.nan
    return dict(z=(real['median']-bg['median'])/den if den>0 else np.nan,
       cnr=(m-b)/bg['sd'] if bg['sd']>0 else np.nan,sbr=sbr,
       sbr_db=10*np.log10(sbr) if sbr>0 else np.nan,contrast=(m-b)/(m+b) if m+b!=0 else np.nan,
       robust_den_zero=bool(den==0),sd_den_zero=bool(bg['sd']==0))
def rois(shape,g):
    zi,xi,w,full=source_local_weights(shape,g,6)
    out={'source':(zi,xi,full)}|{f's{i+1}':(zi,xi,w[i]) for i in range(6)}
    xw=interval_overlap_weights(shape[1],g.x_left_edge_px,g.x_right_edge_px);xx=np.flatnonzero(xw)
    for name,lo,hi in [(f't{k+1:02}',k*25,(k+1)*25) for k in range(20)]+[('tail100',0,100),('tail500',0,500)]:
        zw=interval_overlap_weights(shape[0],g.z_bottom_edge_px+lo/DZ,g.z_bottom_edge_px+hi/DZ);zz=np.flatnonzero(zw)
        out[name]=(zz,xx,np.outer(zw[zz],xw[xx]))
    return out
def pixels(a,m):
    z,x,w=m;v=a[np.ix_(z,x)];ok=w>0
    return v[ok],w[ok]
def overlap(a,b):
    az,ax,aw=a;bz,bx,bw=b
    z,ai,bi=np.intersect1d(az,bz,return_indices=True);x,aj,bj=np.intersect1d(ax,bx,return_indices=True)
    return float(np.minimum(aw[np.ix_(ai,aj)],bw[np.ix_(bi,bj)]).sum())
def audit_raw(sv,meta):
    f=np.isfinite(sv);v=sv[f];q=np.quantile(v,[0,.5,.95,.99,.999,1])
    return meta|dict(shape=str(list(sv.shape)),n_pixels=sv.size,finite_fraction=f.mean(),nan_count=np.isnan(sv).sum(),
        inf_count=np.isinf(sv).sum(),negative_count=(sv<0).sum(),zero_fraction=(sv==0).mean(),
        **dict(zip(['min','median','p95','p99','p999','max'],q)))
def main():
    start=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    frozen={str(p.relative_to(ROOT)):digest(p) for folder in [PREV,DEP] for p in folder.rglob('*') if p.is_file() and '__pycache__' not in str(p)}
    js('frozen_files_before.json',frozen)
    js('validation.json',dict(status='running',starting_head=start))
    fw=read(PREV/'framewise_source_tail_metrics.csv')
    assert digest(PREV/'framewise_source_tail_metrics.csv')==json.loads((PREV/'validation.json').read_text())['framewise_sha256']
    identities=read(DEP/'array_identity_audit.csv.gz').fillna('')
    prior=read(DEP/'input_manifest.csv');manifest=[]
    # Verify all frozen table and implementation identities used by this layer.
    for path,h in frozen.items():manifest.append(dict(input_path=path,sha256=h,role='frozen coupling reference',bytes=(ROOT/path).stat().st_size))
    for r in prior.itertuples():
        p=Path(r.input_path);p=p if p.is_absolute() else ROOT/p
        if p.suffix=='.zip':
            assert digest(p)==r.sha256
            manifest.append(dict(input_path=str(p.resolve()),sha256=r.sha256,role='release ZIP',bytes=p.stat().st_size))
    expected=fw.set_index(['scan_id','frame_index'])
    rows=[];raws=[];struct=[];identity=[];errs={};maxov=0.;maxtrans=0.;spot=0.;negative=[]
    for scan,ii in identities.groupby('scan_id',sort=True):
        print('Extract QC '+scan,flush=True);cache={};volume_rows=[]
        for ir in ii.itertuples():
            r=expected.loc[(scan,ir.frame_index)];meta={k:r[k] if k!='scan_id' else scan for k in ID}|dict(frame_index=int(ir.frame_index),valid_geometry=True)
            p=Path(ir.input_path);p=p if p.is_absolute() else ROOT/p
            if ir.member:
                if str(p) not in cache:cache[str(p)]=zipfile.ZipFile(p)
                blob=cache[str(p)].read(ir.member)
                with np.load(io.BytesIO(blob),allow_pickle=False) as d:arrays={k:d[k] for k in ['sv_raw','stru_amp'] if k in d}
                h=hashlib.sha256(blob).hexdigest()
            else:
                blob=p.read_bytes();h=hashlib.sha256(blob).hexdigest();arrays=loadmat(io.BytesIO(blob),variable_names=['sv_raw','stru_amp'],simplify_cells=True)
                manifest.append(dict(input_path=str(p.resolve()),sha256=h,role='retained formal MAT',bytes=len(blob)))
            assert h==ir.sha256==r.input_sha256
            identity.append(meta|dict(input_path=str(p.resolve()),member=ir.member,sha256=h,identity_matched=True))
            sv=np.asarray(arrays['sv_raw'],float);raws.append(audit_raw(sv,meta))
            if (sv<0).any():
                for z,x in np.argwhere(sv<0):negative.append(meta|dict(z=int(z),x=int(x),value=sv[z,x]))
            assert sv.shape==(351,500) and np.isfinite(sv).all(), 'raw integrity failure'
            st=arrays.get('stru_amp');available=st is not None and st.shape==sv.shape and np.isfinite(st).all()
            struct.append(meta|dict(has_stru_amp=st is not None,shape_matches=st is not None and st.shape==sv.shape,finite_all=available))
            g=VesselGeometry(r.x_left_edge_px,r.x_right_edge_px,r.z_top,float(r.diameter_um),DX,DZ)
            assert abs(g.x_center_px-r.X4)<1e-10
            real=rois(sv.shape,g);masks={'real':real};rec=meta|{k:r[k] for k in ['X1','X1_px','X4','z_top']}
            for side,sign in [('left',-1),('right',1)]:
                shift=sign*1.5*r.X1_px
                bg=VesselGeometry(g.x_left_edge_px+shift,g.x_right_edge_px+shift,g.z_top_edge_px,g.diameter_um,DX,DZ)
                trans=max(abs(bg.lateral_width_um-g.lateral_width_um),abs(bg.z_bottom_edge_px-g.z_bottom_edge_px),abs(bg.x_center_px-(r.X4+shift)))
                maxtrans=max(maxtrans,trans);assert trans<1e-10
                valid=interval_is_complete(sv.shape[1],bg.x_left_edge_px,bg.x_right_edge_px) and interval_is_complete(sv.shape[0],bg.z_top_edge_px,bg.z_bottom_edge_px+500/DZ)
                rec[side+'_valid']=valid;rec[side+'_center_px']=bg.x_center_px
                if valid:
                    masks[side]=rois(sv.shape,bg)
                    for name in ['source','tail500']:
                        ov=overlap(real[name],masks[side][name]);maxov=max(maxov,ov);assert ov<1e-12
                    if ir.frame_index==ii.frame_index.min():
                        z,x,w=masks[side]['source'];ref=ellipse_weights(sv.shape,bg,supersample=16)
                        spot=max(spot,float(np.max(np.abs(w-ref[np.ix_(z,x)]))));assert spot==0
            rec['any_valid']=rec['left_valid'] or rec['right_valid'];rec['both_valid']=rec['left_valid'] and rec['right_valid']
            for name in REGIONS:
                packs={s:pixels(sv,mm[name]) for s,mm in masks.items()}
                sides=[packs[s] for s in ['left','right'] if s in packs]
                packs['bg']=(np.concatenate([a[0] for a in sides]),np.concatenate([a[1] for a in sides])) if sides else (np.array([]),np.array([]))
                ss={s:stats(*v) for s,v in packs.items()}
                for side in ['real','left','right','bg']:
                    for k,v in ss.get(side,dict.fromkeys(STATS,np.nan)).items():rec[f'{name}_{side}_{k}']=v
                for k,v in contrast(ss['real'],ss['bg']).items():rec[f'{name}_{k}']=v
                for t in [0,1,2,3]:rec[f'{name}_z_gt{t}']=float(rec[f'{name}_z']>t) if np.isfinite(rec[f'{name}_z']) else np.nan
                if len(sides)==2:
                    l,rr=ss['left']['mean'],ss['right']['mean'];eps=np.finfo(float).tiny
                    ratio=np.log(l+eps)-np.log(rr+eps)
                    rec[f'{name}_bg_asym_abs']=abs(l-rr);rec[f'{name}_bg_asym_relative']=2*abs(l-rr)/(l+rr) if l+rr>0 else np.nan
                    rec[f'{name}_bg_logratio']=ratio;rec[f'{name}_bg_asym_log']=abs(ratio)
                if name in ['source','tail100','tail500']:
                    for k in ['mean','q']:
                        ref=r[f'source_{k}_raw'] if name=='source' else r[f'tail_{k}_raw_{name[4:]}um']
                        val=ss['real'][k];e=abs(val-ref);rel=e/abs(ref) if ref else e
                        old=errs.setdefault(name+'_'+k,dict(max_absolute_error=0.,max_relative_error=0.))
                        old['max_absolute_error']=max(old['max_absolute_error'],e);old['max_relative_error']=max(old['max_relative_error'],rel)
                        assert rel<1e-10,'formal reintegration mismatch'
                if available and name in ['source','tail100']:
                    sp={s:pixels(st,mm[name]) for s,mm in masks.items()};sb=[sp[s] for s in ['left','right'] if s in sp]
                    bgstats=stats(np.concatenate([a[0] for a in sb]),np.concatenate([a[1] for a in sb])) if sb else dict.fromkeys(STATS,np.nan)
                    sr=stats(*sp['real'])
                    for side,ss0 in [('real',sr),('bg',bgstats)]:
                        for k,v in ss0.items():rec[f'structural_{name}_{side}_{k}']=v
                    for k,v in contrast(sr,bgstats).items():
                        if k!='sbr_db':rec[f'structural_{name}_{k}']=v
            rows.append(rec);volume_rows.append(rec)
        for z in cache.values():z.close()
        csv('checkpoint_'+scan+'.csv.gz',pd.DataFrame(volume_rows))
        js('extraction_progress.json',dict(last_completed_volume=scan,completed_frames=len(rows),formal_errors=errs,max_overlap_weight=maxov))
    data=pd.DataFrame(rows)
    assert len(data)==10931 and data.scan_id.nunique()==23
    for (d,scope),gg in data.groupby(['diameter_um','grid_scope']):
        nominal=[]
        for key,g0 in gg.groupby(ID):
            g0=g0.set_index('frame_index').reindex(range(500)).rename_axis('frame_index').reset_index()
            for k,v in zip(ID,key):g0[k]=v
            g0['valid_geometry']=g0.valid_geometry.fillna(False).astype(bool);nominal.append(g0)
        csv(f'framewise_qc_D{int(d)}'+('_extension' if scope=='d500_extension' else '')+'.csv.gz',pd.concat(nominal,ignore_index=True))
    csv('raw_sv_array_audit.csv.gz',pd.DataFrame(raws));csv('array_identity_audit.csv.gz',pd.DataFrame(identity))
    csv('negative_values.csv.gz',pd.DataFrame(negative,columns=ID+['frame_index','z','x','value']))
    csv('structural_availability.csv',pd.DataFrame(struct));csv('input_manifest.csv',pd.DataFrame(manifest).drop_duplicates('input_path'))
    js('extraction_validation.json',dict(status='passed',starting_head=start,formal_reintegration=errs,max_overlap_weight=maxov,
       translation_parameter_max_error=maxtrans,full_image_ellipse_spotcheck_max_weight_error=spot,
       common_grid=dict(volumes=20,valid_frames=9456,nominal_frames=10000),extension=dict(volumes=3,valid_frames=1475,nominal_frames=1500),
       all_array_sha_passed=True,raw_finite_all=all(a['finite_fraction']==1 for a in raws),negative_count=sum(a['negative_count'] for a in raws),
       structural_complete=all(a['finite_all'] for a in struct),original_invalid_frames_preserved_as_missing=True))
    js('provenance.json',dict(starting_head=start,repository='bulbel-magnolia/OCTA-',branch='analysis/sv-diameter-stage2',
       software=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__),
       formal_signal='var(abs(IMG),1,3), linear raw SV; denominator N; unchanged',
       geometry='frozen continuity-first v2.1; X4 center, X1 span, physical D; dx12.7 dz6.7, supersample16',
       background='centers X4 +/- 1.5 X1_px; source ellipse translated in physical geometry and rasterized on same lattice; tail same boundaries',
       weighted_statistics='mean, median, MAD and population SD use fractional weights; weighted median inverse cumulative mass with midpoint for exact half ties; BG pools side pixels by weight',
       fov='side invalid if full source plus 0-500um tail outside image edges; no recenter or clipping',
       asymmetry='signed log ratio uses float64 tiny epsilon; absolute log, absolute difference, 2abs(L-R)/(L+R)',
       descriptive_thresholds=[0,1,2,3],denominator_zero='NaN, explicit flags; no epsilon for Z/CNR',
       physical_inference='none',calibrated_noise_reference='none documented in retained cohort provenance; no calibrated instrument SNR computed'))
    print('Extraction passed '+str(errs),flush=True)

if __name__=='__main__':main()
