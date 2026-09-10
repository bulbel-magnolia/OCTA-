"""New fixed-D weighted distributions, matched pseudo ROIs and rigid offsets.

Official real means/Q/area are read from v2 tables; real pixel integrals are
computed only on the frozen 45 frames for independent numerical checks.
"""
from common import *
import io,zipfile
from scipy.io import loadmat
sys.path.insert(0,str(DEP))
from derive_depth_metrics import source_local_weights,VesselGeometry,interval_overlap_weights
from svrecttail.geometry import interval_is_complete
FIELDS=['median','mad','mean','q','area','neff','n_pixels','sum_w2','z_span_pixels']
def geometry(r,dx=0,dz=0):
    return VesselGeometry(r.X4-r.diameter_um/(2*DX)+dx,r.X4+r.diameter_um/(2*DX)+dx,r.z_top+dz,r.diameter_um,DX,DZ)
def complete(g,shape=(351,500)):
    return interval_is_complete(shape[1],g.x_left_edge_px,g.x_right_edge_px) and interval_is_complete(shape[0],g.z_top_edge_px,g.z_bottom_edge_px+500/DZ)
def masks(g,shape=(351,500),limited=False):
    z,x,w,full=source_local_weights(shape,g,6)
    out={'source':(z,x,full),'s36':(z,x,w[2:].sum(axis=0))}|{f's{i+1}':(z,x,w[i]) for i in range(6)}
    xw=interval_overlap_weights(shape[1],g.x_left_edge_px,g.x_right_edge_px);xx=np.flatnonzero(xw>0)
    for name,lo,hi in [(f't{i:02}',(i-1)*25,i*25) for i in range(1,5 if limited else 21)]+[('tail100',0,100),('tail500',0,500)]:
        zw=interval_overlap_weights(shape[0],g.z_bottom_edge_px+lo/DZ,g.z_bottom_edge_px+hi/DZ);zz=np.flatnonzero(zw>0)
        out[name]=(zz,xx,np.outer(zw[zz],xw[xx]))
    return out
def pixels(a,m):
    z,x,w=m;v=a[np.ix_(z,x)];ok=w>0;return v[ok],w[ok]
def wmedian(v,w):
    m=(w>0)&np.isfinite(w)&np.isfinite(v);v,w=np.asarray(v)[m],np.asarray(w)[m]
    if not len(w):return np.nan
    order=np.argsort(v,kind='stable');v,w=v[order],w[order];c=np.cumsum(w);half=c[-1]/2;k=np.searchsorted(c,half,side='left')
    # Near a half-weight boundary, floating summation order can select either
    # side or the midpoint. Resolve that numerical uncertainty by exact integer
    # arithmetic on the stored binary floating weights, without an epsilon in
    # the statistic. The error bound only selects the exact-arithmetic path.
    bound=np.finfo(float).eps*len(w)*c[-1]*2
    if min(abs(c[k]-half),abs(c[k-1]-half) if k else np.inf)<=bound:
        ratios={float(x):float(x).as_integer_ratio() for x in np.unique(w)};den=max(b for a,b in ratios.values())
        integers={x:a*(den//b) for x,(a,b) in ratios.items()};iw=[integers[float(x)] for x in w];total=sum(iw);cum=0
        for j,weight in enumerate(iw):
            cum+=weight
            if 2*cum==total:return float((v[j]+v[j+1])/2) if j+1<len(v) else float(v[j])
            if 2*cum>total:return float(v[j])
    return float((v[k]+v[k+1])/2 if c[k]==half and k+1<len(v) else v[k])
def stats(v,w,integrate=True):
    m=w>0;v,w=v[m],w[m]
    if not len(v):return {k:np.nan for k in FIELDS}|dict(status='empty_ROI')
    med=wmedian(v,w);mad=wmedian(np.abs(v-med),w);sw=w.sum();out=dict(median=med,mad=mad,area=sw*DX*DZ,neff=sw*sw/np.dot(w,w),n_pixels=len(w),sum_w2=np.dot(w,w),status='ok')
    if integrate:out.update(mean=np.dot(v,w)/sw,q=np.dot(v,w)*DX*DZ)
    return out
def measure(a,m):
    v,w=pixels(a,m);sw=w.sum();return dict(q=np.dot(v,w)*DX*DZ,area=sw*DX*DZ,mean=np.dot(v,w)/sw)
def shared(a,b):
    az,ax,aw=a;bz,bx,bw=b
    _,ai,bi=np.intersect1d(az,bz,return_indices=True);_,aj,bj=np.intersect1d(ax,bx,return_indices=True)
    x,y=aw[np.ix_(ai,aj)],bw[np.ix_(bi,bj)];den=np.sqrt(np.sum(aw**2)*np.sum(bw**2))
    return dict(shared_pixels=int(((x>0)&(y>0)).sum()),overlap_index=float(np.sum(x*y)/den) if den else np.nan)
def load_array(ir,cache):
    p=Path(ir.input_path)
    if ir.member:
        if str(p) not in cache:cache[str(p)]=zipfile.ZipFile(p)
        blob=cache[str(p)].read(ir.member);h=hashlib.sha256(blob).hexdigest();assert h==ir.sha256
        with np.load(io.BytesIO(blob),allow_pickle=False) as z:arrays={k:np.asarray(z[k],float) for k in ['sv_raw','stru_amp'] if k in z}
    else:
        blob=p.read_bytes();h=hashlib.sha256(blob).hexdigest();assert h==ir.sha256
        arrays={k:np.asarray(v,float) for k,v in loadmat(io.BytesIO(blob),variable_names=['sv_raw','stru_amp']).items() if k in ['sv_raw','stru_amp']}
    assert arrays['sv_raw'].shape==(351,500) and np.isfinite(arrays['sv_raw']).all()
    return arrays,h
def independent_masks(g,shape=(351,500),limited=False):
    """Different construction: flattened sample coordinates, bin assignment,
    then bincount to pixel+bin, versus the production tensor masks."""
    xs=np.arange(max(0,int(np.ceil(g.x_left_edge_px-.5))),min(shape[1],int(np.floor(g.x_right_edge_px+.5))+1))
    zs=np.arange(max(0,int(np.ceil(g.z_top_edge_px-.5))),min(shape[0],int(np.floor(g.z_bottom_edge_px+.5))+1))
    xx,zz=np.meshgrid(xs,zs);flat_x=xx.ravel();flat_z=zz.ravel();counts=np.zeros((6,len(flat_x)))
    for iz in range(16):
        z=flat_z+(iz+.5)/16-.5;u=(z-g.z_top_edge_px)*DZ/g.diameter_um
        binid=np.minimum(np.floor(u*6).astype(int),5)
        for ix in range(16):
            x=flat_x+(ix+.5)/16-.5;inside=((x-g.x_center_px)*DX/(g.diameter_um/2))**2+((z-g.z_center_px)*DZ/(g.diameter_um/2))**2<=1
            pos=np.flatnonzero(inside&(u>=0)&(u<=1));counts[binid[pos],pos]+=1
    weights=counts.reshape(6,len(zs),len(xs))/256
    out={f's{i+1}':(zs,xs,weights[i]) for i in range(6)}|{'source':(zs,xs,weights.sum(axis=0)),'s36':(zs,xs,weights[2:].sum(axis=0))}
    for name,lo,hi in [(f't{i:02}',25*(i-1),25*i) for i in range(1,5 if limited else 21)]+[('tail100',0,100),('tail500',0,500)]:
        bot=g.z_bottom_edge_px;zs2=np.arange(max(0,int(np.ceil(bot+lo/DZ-.5))),min(shape[0],int(np.floor(bot+hi/DZ+.5))+1))
        wx=np.array([max(0,min(x+.5,g.x_right_edge_px)-max(x-.5,g.x_left_edge_px)) for x in xs])
        wz=np.array([max(0,min(z+.5,bot+hi/DZ)-max(z-.5,bot+lo/DZ)) for z in zs2]);out[name]=(zs2,xs,np.outer(wz,wx))
    return out
def independent_wmedian(v,w):
    from fractions import Fraction
    pairs=sorted((float(x),float(y)) for x,y in zip(v,w) if y>0 and np.isfinite(x) and np.isfinite(y))
    if not pairs:return np.nan
    pairs=[(x,Fraction.from_float(y)) for x,y in pairs];half=sum(y for x,y in pairs)/2;c=0
    for i,(x,y) in enumerate(pairs):
        c+=y
        if c==half and i+1<len(pairs):return (x+pairs[i+1][0])/2
        if c>half:return x
    return pairs[-1][0]
def weighted_tests():
    tests=[([1,3,9],[1,2,1],3),([1,3,9],[1,1,2],6),([1,1,1],[1,3,4],1),([],[],np.nan),([1,2,8,9],[0,1,0,1],5.5),([1,2,3,4],[.1,.2,.1,.2],2.5)]
    rows=[]
    for i,(v,w,expected) in enumerate(tests):
        v,w=np.array(v,float),np.array(w,float);got=wmedian(v,w);ref=independent_wmedian(v,w)
        assert (np.isnan(got) and np.isnan(expected) and np.isnan(ref)) or got==expected==ref
        mad=wmedian(abs(v-got),w);rmad=independent_wmedian(abs(v-ref),w)
        assert (np.isnan(mad) and np.isnan(rmad)) or mad==rmad
        rows.append(dict(case=i,median=got,mad=mad,status='passed'))
    csv('weighted_statistics_tests.csv',rows)
def verify_integrals(a,gg,mm,meta,rows,reference=None,limited=False):
    refm=independent_masks(gg,limited=limited)
    for region,m in mm.items():
        vals=measure(a,m);refs=measure(a,refm[region]);v,w=pixels(a,m)
        for field in ['q','area','mean']:
            ref=refs[field];err=abs(vals[field]-ref);rel=err/abs(ref) if ref!=0 else np.nan
            assert (rel<1e-10) if ref!=0 else err==0
            rows.append(meta|dict(region=region,field=field,reference_value=ref,actual_value=vals[field],relative_error=rel,zero_reference=ref==0,absolute_error=err,status='passed'))
            if reference is not None:
                vr=reference[region][field];err=abs(vals[field]-vr);rel=err/abs(vr) if vr!=0 else np.nan
                assert (rel<1e-10) if vr!=0 else err==0
                rows.append(meta|dict(region=region,field='formal_table_'+field,reference_value=vr,actual_value=vals[field],relative_error=rel,zero_reference=vr==0,absolute_error=err,status='passed'))
        # Every region's weighted median/MAD is checked against a separate sorting implementation.
        med=wmedian(v,w);refmed=independent_wmedian(v,w);mad=wmedian(abs(v-med),w);refmad=independent_wmedian(abs(v-refmed),w)
        assert med==refmed and mad==refmad
    for target,parts in [('source',[f's{s}' for s in range(1,7)]),('tail100',[f't{t:02}' for t in range(1,5)])]+([] if limited else [('tail500',[f't{t:02}' for t in range(1,21)])]):
        for field in ['q','area']:
            val=sum(measure(a,mm[p])[field] for p in parts);ref=measure(a,mm[target])[field]
            assert abs((val-ref)/ref)<1e-10 if ref else val==0
def draw_qc(scan,shots,limits,folder='localization'):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse,Rectangle
    fig,axes=plt.subplots(len(shots),2,figsize=(13,3.33*len(shots)),layout='constrained',squeeze=False)
    for row,(r,arrays) in enumerate(shots):
        for col0,name in enumerate(['stru_amp','sv_raw']):
            ax=axes[row,col0];arr=arrays.get(name);g=geometry(r);d=r.diameter_um
            if arr is None:ax.text(.5,.5,'Structural array unavailable',ha='center');continue
            ax.imshow(arr,origin='upper',extent=[-.5*DX,499.5*DX,350.5*DZ,-.5*DZ],cmap='gray' if name=='stru_amp' else 'magma',vmin=0,vmax=limits[name],aspect='equal')
            for side,shift,color in [('real',0,'#00ffff'),('L',-1.5*d/DX,'#00ff00'),('R',1.5*d/DX,'#ff66ff')]:
                q=geometry(r,shift,0);ax.add_patch(Ellipse((q.x_center_px*DX,q.z_center_px*DZ),d,d,fill=False,edgecolor=color,lw=1.2))
                ax.add_patch(Rectangle((q.x_left_edge_px*DX,q.z_bottom_edge_px*DZ),d,500,fill=False,edgecolor=color,lw=1))
            ax.axvline(r.X4*DX,color='#00ffff',ls=':',lw=.6);ax.plot([g.x_left_edge_px*DX,g.x_right_edge_px*DX],[g.z_top_edge_px*DZ]*2,color='white',lw=1)
            ax.plot([g.x_left_edge_px*DX,g.x_right_edge_px*DX],[g.z_bottom_edge_px*DZ]*2,color='yellow',ls='--',lw=1)
            ax.set(xlim=((r.X4*DX)-2.8*d,(r.X4*DX)+2.8*d),ylim=(min(350.5*DZ,g.z_bottom_edge_px*DZ+550),max(-.5*DZ,r.z_top*DZ-120)),xlabel='Lateral position (um)',ylabel='Axial position (um)',title=f'frame {int(r.frame_index)} | {name} | linear [0,{limits[name]:.3g}]')
    fig.suptitle(scan+' | physical aspect; cyan real / green L / magenta R; yellow = prior lower boundary\nDisplay range only; numerical SV unchanged')
    (OUT/folder).mkdir(exist_ok=True);fig.savefig(OUT/folder/f'{scan}.png',dpi=150);plt.close(fig)
def main():
    weighted_tests();fw=frame_table();ids=read(INPUT/'array_identity_audit.csv.gz').fillna('').set_index(KEY);fixed=read(OUT/'localization_fixed_frames.csv')
    allident=[];allchecks=[];coverage=[];localization=[];support_files=[];overlap_files=[]
    for key,g in fw.groupby(ID,sort=True):
        meta=dict(zip(ID,key));scan=meta['scan_id'];print('New fixed-D pixels',scan,flush=True)
        outfile=OUT/f'framewise_local_background_qc_{scan}.csv.gz';assert not outfile.exists(),'No checkpoint overwrite'
        cache={};records=[];perturb=[];support=[];overlaps=[];shots=[];limits={'sv_raw':0.,'stru_amp':0.}
        spotids=set(fixed.loc[fixed.scan_id.eq(scan),'frame_index']);sig=signals(g);qs=signals(g,'q');ars=signals(g,'area')
        for r in g.itertuples():
            i=int(r.frame_index);rec=meta|dict(frame_index=i,valid_geometry=bool(r.valid_geometry),left_valid=False,right_valid=False,both_valid=False,
                left_status='invalid_geometry',right_status='invalid_geometry',combined_status='invalid_geometry')
            pr=meta|dict(frame_index=i,valid_geometry=bool(r.valid_geometry))
            if not r.valid_geometry:
                for name in OFFSETS:pr[name+'_status']='invalid_geometry';pr[name+'_valid']=False
                records.append(rec);perturb.append(pr);continue
            ir=ids.loc[(scan,i)];arrays,h=load_array(ir,cache);assert h==r.input_sha256;sv=arrays['sv_raw']
            allident.append(meta|dict(frame_index=i,input_path=ir.input_path,member=ir.member,sha256=h,status='matched_on_read'))
            for name in limits:
                if name in arrays and np.isfinite(arrays[name]).all():limits[name]=max(limits[name],float(arrays[name].max()))
            if i in spotids:shots.append((r,arrays))
            gm={'real':geometry(r)}|{side:geometry(r,sign*1.5*r.diameter_um/DX) for side,sign in [('left',-1),('right',1)]}
            mm={'real':masks(gm['real'])}
            for side in ['left','right']:
                ok=complete(gm[side]);rec[side+'_valid']=ok;rec[side+'_status']='ok' if ok else 'out_of_FOV';rec[side+'_center_px']=gm[side].x_center_px
                assert abs(gm[side].lateral_width_um-r.diameter_um)<1e-10 and abs(gm[side].z_top_edge_px-r.z_top)<1e-12
                if ok:mm[side]=masks(gm[side])
            both=rec['left_valid'] and rec['right_valid'];rec['both_valid']=both;rec['combined_status']='ok' if both else 'requires_both_complete_sides'
            for region in REGIONS:
                packs={s:pixels(sv,mm0[region]) for s,mm0 in mm.items()}
                if both:packs['bg']=(np.concatenate([packs[s][0] for s in ['left','right']]),np.concatenate([packs[s][1] for s in ['left','right']]))
                outstats={s:stats(v,w,integrate=s!='real') for s,(v,w) in packs.items()}
                outstats['real'].update(mean=sig[region][i],q=qs[region][i],area=ars[region][i])
                for side in ['real','left','right','bg']:
                    for field in FIELDS:rec[f'{region}_{side}_{field}']=outstats.get(side,{}).get(field,np.nan)
                    if side in mm:
                        z,x,w=mm[side][region];span=int(np.count_nonzero((w>0).any(axis=1)));rec[f'{region}_{side}_z_span_pixels']=span
                        # Source table Neff is reused and independently compared with geometric weight Neff.
                        ss=outstats[side];given=getattr(r,'source_effective_pixels' if region=='source' else f'source_{region}_effective_pixels',np.nan)
                        if side=='real' and np.isfinite(given):assert abs(given-ss['neff'])/given<1e-10;rec[f'{region}_{side}_neff']=given
                        support.append(meta|dict(frame_index=i,region=region,side=side,neff=rec[f'{region}_{side}_neff'],nonzero_pixels=ss['n_pixels'],sum_w2=ss['sum_w2'],actual_area_um2=ss['area'],z_span_pixels=span))
                if both:
                    real,bg=outstats['real'],outstats['bg'];den=1.4826*bg['mad'];rec[region+'_z_local']=(real['median']-bg['median'])/den if den!=0 else np.nan
                    rec[region+'_z_status']='ok' if den!=0 else 'zero_background_MAD';rec[region+'_real_minus_bg_mean']=real['mean']-bg['mean']
                    l,rr=outstats['left']['mean'],outstats['right']['mean'];rec[region+'_signed_asymmetry']=2*(l-rr)/(l+rr) if l+rr!=0 else np.nan
                    rec[region+'_asymmetry_status']='ok' if l+rr!=0 else 'zero_sum_means'
                    # Validate combined mean is the Q/area pool, never median-of-side-means.
                    pooled=(outstats['left']['q']+outstats['right']['q'])/(outstats['left']['area']+outstats['right']['area'])
                    assert abs(pooled-bg['mean'])<=abs(bg['mean'])*1e-10 if bg['mean'] else pooled==0
                else:
                    rec[region+'_z_local']=np.nan;rec[region+'_z_status']='requires_both_complete_sides';rec[region+'_real_minus_bg_mean']=np.nan;rec[region+'_signed_asymmetry']=np.nan;rec[region+'_asymmetry_status']='requires_both_complete_sides'
            for s in range(1,7):
                for t in range(1,21):overlaps.append(meta|dict(frame_index=i,source_bin=s,tail_bin=t)|shared(mm['real'][f's{s}'],mm['real'][f't{t:02}']))
            if i in spotids:
                reference={region:{field:d[region][i] for field,d in [('mean',sig),('q',qs),('area',ars)]} for region in REGIONS}
                for side,m in mm.items():verify_integrals(sv,gm[side],m,meta|dict(frame_index=i,side=side),allchecks,reference if side=='real' else None)
                localization.append(meta|dict(frame_index=i,figure=f'localization/{scan}.png',structural_available='stru_amp' in arrays,review_status='pending_visual_review',source_boundary='prior circle and lower boundary; no truth claim'))
            for name,(dx,dz) in OFFSETS.items():
                gg=geometry(r,dx,dz);ok=complete(gg);pr[name+'_valid']=ok;pr[name+'_status']='ok' if ok else 'out_of_FOV'
                if not ok:continue
                pm=masks(gg,limited=True)
                for region,m in pm.items():
                    for field,value in measure(sv,m).items():pr[f'{name}_{region}_{field}']=value
                if i in spotids:verify_integrals(sv,gg,pm,meta|dict(frame_index=i,side=name),allchecks,limited=True)
            records.append(rec);perturb.append(pr)
        for z in cache.values():z.close()
        out=pd.DataFrame(records);csv(outfile.name,out);csv(f'position_integrals_{scan}.csv.gz',perturb)
        sp=f'geometry_support_{scan}.csv.gz';ov=f'shared_pixel_support_{scan}.csv.gz';csv(sp,support);csv(ov,overlaps);support_files.append(sp);overlap_files.append(ov)
        coverage.append(meta|dict(n_geometry_valid=int(g.valid_geometry.sum()),left_valid=int(out.left_valid.sum()),right_valid=int(out.right_valid.sum()),both_valid=int(out.both_valid.sum()),n_nominal=500))
        draw_qc(scan,shots,limits)
        status('pixel_extraction',last_completed_scan=scan,arrays_verified=len(allident))
    csv('array_identity_on_read.csv.gz',allident);csv('independent_integration_checks.csv.gz',allchecks);csv('background_coverage.csv',coverage);csv('localization_review.csv',localization)
    js('pixel_validation.json',dict(status='passed',arrays_verified=len(allident),geometry_valid_frames=7018,checks=len(allchecks),fixed_check_frames=45,
        weighted_median_mad_tests='passed',all_combined_background_requires_both=True,official_real_integrals_reused=True,pseudo_v2_recomputed=True,offsets_all_four=True,
        support_files=support_files,overlap_files=overlap_files,localization_visual_review='pending'))
    status('pixel_extraction_complete')
if __name__=='__main__':main()
