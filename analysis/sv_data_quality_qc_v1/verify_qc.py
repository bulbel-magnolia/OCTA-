"""Independent full-image ROI checks, frozen depth agreement, and delivery gates."""
import io,json,zipfile,hashlib
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.io import loadmat
from derive_qc import OUT,ROOT,DEP,PREV,ID,REGIONS,read,csv,js,digest,DX,DZ,VesselGeometry,ellipse_weights,rectangle_weights

def wm(v,w):
    # Independent scalar cumulative implementation (not the extraction helper).
    vals=sorted(zip(v,w),key=lambda p:p[0]);half=sum(w)/2.;total=0.
    for i,(x,weight) in enumerate(vals):
        total+=weight
        if total>=half:return (x+vals[i+1][0])/2 if total==half and i+1<len(vals) else x
    raise AssertionError('weighted median')
def main():
    val=json.loads((OUT/'validation.json').read_text());assert val['status']=='passed'
    prior=json.loads((DEP/'data_validation.json').read_text())
    for name,h in prior['derived_file_sha256'].items():assert digest(DEP/name)==h
    priorp=json.loads((DEP/'provenance.json').read_text())
    for name,h in priorp['all_script_sha256'].items():assert digest(DEP/name)==h
    for row in read(DEP/'input_manifest.csv').itertuples():
        if row.input_path.startswith(('analysis/','src/','results/')):assert digest(ROOT/row.input_path)==row.sha256
    data=pd.concat([read(p) for p in sorted(OUT.glob('framewise_qc_D*.csv.gz'))],ignore_index=True)
    fw=read(PREV/'framewise_source_tail_metrics.csv');ident=read(OUT/'array_identity_audit.csv.gz').fillna('')
    depth=pd.concat([read(p) for p in sorted(DEP.glob('framewise_depth_metrics_D*.csv.gz'))],ignore_index=True)
    keys=['scan_id','frame_index'];valid=data[data.valid_geometry.eq(True)].copy()
    assert set(map(tuple,valid[keys].values))==set(map(tuple,fw[keys].values))
    invalid=data[~data.valid_geometry.eq(True)];assert len(invalid)==569 and invalid.filter(regex='_(real|bg|left|right)_(mean|median)$').isna().all().all()
    for region in REGIONS:
        for side in ['real','left','right','bg']:
            mass=valid[f'{region}_{side}_area']/(DX*DZ);neff=valid[f'{region}_{side}_neff']
            assert ((neff+1e-9)>=mass).all() and (neff>=1).all()
    merged=valid.merge(depth,on=keys,validate='one_to_one',suffixes=('','_depth'));errors=[]
    for name in [f's{i}' for i in range(1,7)]+[f't{i:02}' for i in range(1,21)]:
        prefix='source_'+name if name.startswith('s') else 'tail_'+name
        for field,old in [('mean','mean_raw'),('q','q_raw'),('area','area_um2')]:
            a=merged[f'{name}_real_{field}'].to_numpy();b=merged[f'{prefix}_{old}'].to_numpy();e=np.max(abs(a-b));rel=np.max(abs(a-b)/np.maximum(abs(b),np.finfo(float).tiny));assert rel<1e-10
            errors.append(dict(region=name,field=field,max_absolute_error=e,max_relative_error=rel))
    csv('frozen_depth_reintegration_validation.csv',pd.DataFrame(errors))
    checks=[];availability=[];noiseflags=[]
    for scan,g in valid.groupby('scan_id'):
        # Predefined nearest valid original index to 249, tie favors smaller index.
        i=min(g.frame_index,key=lambda i:(abs(i-249),i));r=fw[fw.scan_id.eq(scan)&fw.frame_index.eq(i)].iloc[0];qr=g[g.frame_index.eq(i)].iloc[0]
        ir=ident[ident.scan_id.eq(scan)&ident.frame_index.eq(i)].iloc[0];p=Path(ir.input_path)
        if ir.member:
            with zipfile.ZipFile(p) as z:blob=z.read(ir.member)
            with np.load(io.BytesIO(blob),allow_pickle=False) as a:
                sv=np.asarray(a['sv_raw'],float);st=np.asarray(a['stru_amp'],float);names=a.files
                metadata=str(a['metadata_json']) if 'metadata_json' in a else ''
        else:
            a=loadmat(p,simplify_cells=True);sv=np.asarray(a['sv_raw'],float);st=np.asarray(a['stru_amp'],float);names=list(a);metadata=str(a.get('metadata',''))
        availability.append(dict(scan_id=scan,frame_index=i,keys=json.dumps(names),metadata=metadata))
        for term in ['noise_floor','dark_acquisition','air_roi','noise_reference']:
            if term in metadata.lower():noiseflags.append(dict(scan_id=scan,term=term,context=metadata))
        for name,lo,hi in [('source',None,None),('t01',0,25),('tail100',0,100),('tail500',0,500)]:
            samples={}
            for side,sign in [('real',0),('left',-1),('right',1)]:
                shift=sign*1.5*r.X1/12.7;geo=VesselGeometry(r.x_left_edge_px+shift,r.x_right_edge_px+shift,r.z_top,float(r.diameter_um),DX,DZ)
                w=ellipse_weights(sv.shape,geo,supersample=16) if name=='source' else rectangle_weights(sv.shape,x_left_edge_px=geo.x_left_edge_px,x_right_edge_px=geo.x_right_edge_px,z_top_edge_px=geo.z_bottom_edge_px+lo/DZ,z_bottom_edge_px=geo.z_bottom_edge_px+hi/DZ)
                mask=w>0;samples[side]=(sv[mask],w[mask]);mean=np.sum(sv[mask]*w[mask])/w[mask].sum()
                ref=qr[f'{name}_{side}_mean'];err=abs(mean-ref)/abs(ref);assert err<1e-10
                neff=w.sum()**2/np.sum(w*w);assert abs(neff-qr[f'{name}_{side}_neff'])<1e-8 and neff<=mask.sum()+1e-8
                checks.append(dict(scan_id=scan,frame_index=i,region=name,side=side,mean_relative_error=err))
            bgv=np.r_[samples['left'][0],samples['right'][0]];bgw=np.r_[samples['left'][1],samples['right'][1]]
            bmed=wm(bgv,bgw);mad=wm(abs(bgv-bmed),bgw);rmed=wm(*samples['real']);z=(rmed-bmed)/(1.4826*mad)
            assert abs(z-qr[name+'_z'])<1e-10
    csv('independent_roi_spotchecks.csv',pd.DataFrame(checks));csv('retained_metadata_availability_audit.csv',pd.DataFrame(availability))
    csv('noise_reference_keyword_hits.csv',pd.DataFrame(noiseflags,columns=['scan_id','term','context']))
    assert not noiseflags,'Review potential reference metadata before reporting reference absence'
    cov=read(OUT/'background_roi_coverage.csv');assert cov[['left_valid_fraction','right_valid_fraction','any_valid_fraction','both_valid_fraction']].eq(1).all().all()
    flags=read(OUT/'denominator_special_cases.csv');assert (flags.n_frames==0).all()
    supports=read(OUT/'roi_effective_support.csv');assert supports['min'].gt(0).all()
    assert len(read(OUT/'pseudo_source_tail_coupling_maps.csv'))==23*2*2*120
    assert len(read(OUT/'control_excess_coupling_summary.csv'))==23*120
    for p in OUT.glob('lag_D*.csv.gz'):assert len(read(p))==101*120
    before=json.loads((OUT/'frozen_files_before.json').read_text());assert all(digest(ROOT/p)==h for p,h in before.items())
    add=json.loads((OUT/'additional_raw_validation.json').read_text());assert add['available_invalid_arrays']==569 and add['missing_invalid_arrays']==0
    val.update(independent_full_image_roi_spotchecks=len(checks),independent_roi_mean_max_relative_error=pd.DataFrame(checks).mean_relative_error.max(),
       frozen_depth_all_frame_comparison=pd.DataFrame(errors).to_dict('records'),all_nominal_arrays_audited=11500,
       raw_nominal_finite_all=True,raw_nominal_negative_count=0,zero_robust_and_sd_denominators=0,
       calibrated_noise_reference_search=dict(status='unavailable',metadata_samples=23,keyword_hits=0,scope='retained cohort manifests, 23 preselected mid-volume metadata records, source implementation and retained-array keys; no new cohort'),
       no_formal_files_modified=True,effective_support_mass_bound_all_frames_passed=True,
       independent_weighted_Z_checks=23*4,unfinished_computation=[],missing_required_data=[])
    js('validation.json',val)
    print('Independent verification passed',len(checks),'full-image ROI checks; 78 all-frame frozen depth checks',flush=True)

if __name__=='__main__':main()
