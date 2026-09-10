"""Frozen definitions for fixed physical diameter selective analysis v2."""
from pathlib import Path
import sys, json, hashlib, platform, datetime
sys.dont_write_bytecode=True
import numpy as np
import pandas as pd
from scipy.stats import rankdata
OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[1]
INPUT=ROOT/'analysis/sv_physical_diameter_geometry_v2'
DEP=ROOT/'analysis/sv_source_tail_depth_coupling_v1'
PREV=ROOT/'analysis/sv_source_tail_framewise_coupling_v1'
QC=ROOT/'analysis/sv_data_quality_qc_v1'
DX,DZ=12.7,6.7
ID=['scan_id','diameter_um','flow_mm_s']
KEY=['scan_id','frame_index']
REGIONS=['source']+[f's{i}' for i in range(1,7)]+['s36']+[f't{i:02}' for i in range(1,21)]+['tail100','tail500']
REP={'W100':('source','tail100'),'W500':('source','tail500'),'C100':('s36','tail100')}
OFFSETS={'xp':(1,0),'xm':(-1,0),'zp':(0,1),'zm':(0,-1)}
LAGS=np.arange(-50,51)
FAR=np.abs(LAGS)>=30
C=[(s,t) for s in range(3,7) for t in range(1,5)]
def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def clean(x):
    if isinstance(x,dict):return {str(k):clean(v) for k,v in x.items()}
    if isinstance(x,(tuple,list,np.ndarray)):return [clean(v) for v in x]
    if isinstance(x,np.generic):return clean(x.item())
    if isinstance(x,float) and not np.isfinite(x):return None
    return x
def js(n,x):(OUT/n).write_text(json.dumps(clean(x),ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def csv(n,x):
    if not isinstance(x,pd.DataFrame):x=pd.DataFrame(x)
    x.to_csv(OUT/n,index=False,float_format='%.17g',lineterminator='\n',compression={'method':'gzip','mtime':0} if n.endswith('.gz') else None)
def read(p):return pd.read_csv(p,float_precision='round_trip',low_memory=False)
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def status(stage,**kw):
    p=OUT/'run_status.json';x=json.loads(p.read_text()) if p.exists() else {}
    x.update(stage=stage,updated_utc=now(),**kw);js(p.name,x)
def frame_table():
    geo=read(INPUT/'geometry_manifest.csv.gz').sort_values(KEY).reset_index(drop=True)
    fw=pd.concat([read(INPUT/f'framewise_fixed_diameter_D{d}.csv.gz') for d in [128,235,285]],ignore_index=True)
    assert not geo.duplicated(KEY).any() and not fw.duplicated(KEY).any()
    assert set(map(tuple,geo[KEY].values))==set(map(tuple,fw[KEY].values))
    data=geo.merge(fw.drop(columns=[c for c in geo if c not in KEY]),on=KEY,how='left',validate='one_to_one')
    assert len(data)==7500 and data.valid_geometry.sum()==7018 and data.scan_id.nunique()==15
    for scan,g in data.groupby('scan_id'):assert np.array_equal(g.frame_index,np.arange(500))
    return data
def col(region,measure='mean'):
    suffix={'mean':'mean_raw','q':'q_raw','area':'area_um2'}[measure]
    if region=='source':return 'source_'+suffix
    if region.startswith('s') and region!='s36':return 'source_'+region+'_'+suffix
    if region.startswith('tail'):return 'tail_'+suffix+'_'+region[4:]+'um'
    return 'tail_'+region+'_'+suffix
def signals(g,measure='mean'):
    out={r:g[col(r,measure)].to_numpy(float) for r in REGIONS if r!='s36'}
    q=sum(g[col(f's{s}','q')].to_numpy(float) for s in range(3,7))
    a=sum(g[col(f's{s}','area')].to_numpy(float) for s in range(3,7))
    out['s36']=q/a if measure=='mean' else q if measure=='q' else a
    return out
def ratio(a,b):
    a,b=np.broadcast_arrays(np.asarray(a,float),np.asarray(b,float));out=np.full(a.shape,np.nan)
    np.divide(a,b,out=out,where=np.isfinite(a)&np.isfinite(b)&(b!=0));return out
def summary(v):
    a=np.asarray(v,float);b=a[np.isfinite(a)]
    return dict(n_defined=len(b),n_missing=len(a)-len(b),median=float(np.median(b)) if len(b) else np.nan,
                min=float(b.min()) if len(b) else np.nan,max=float(b.max()) if len(b) else np.nan,
                n_positive=int((b>0).sum()),n_negative=int((b<0).sum()),n_zero=int((b==0).sum()))
def strict(v,n):
    z=summary(v);z['value']=z['median'] if z['n_defined']==n and len(v)==n else np.nan
    z['status']='ok' if z['n_defined']==n and len(v)==n else 'incomplete_components'
    return z
def rho(a,b):
    a,b=np.asarray(a,float),np.asarray(b,float);m=np.isfinite(a)&np.isfinite(b);a,b=a[m],b[m];n=len(a)
    if n<3:return np.nan,n,'fewer_than_3_pairs'
    ca=np.all(a==a[0]);cb=np.all(b==b[0])
    if ca or cb:return np.nan,n,'constant_both' if ca and cb else 'constant_source' if ca else 'constant_tail'
    x,y=rankdata(a,method='average'),rankdata(b,method='average');x-=x.mean();y-=y.mean()
    return float(np.dot(x,y)/np.sqrt(np.dot(x,x)*np.dot(y,y))),n,'ok'
def detrend(a,w):
    a=np.asarray(a,float);s=pd.Series(a);roll=s.rolling(w,center=True,min_periods=1)
    trend=roll.median().to_numpy();support=roll.count().to_numpy().astype(int)
    return a-trend,trend,support
def lag_curve(a,b):
    rows=[]
    for lag in LAGS:
        i=np.arange(max(0,-lag),min(500,500-lag));rr,n,reason=rho(a[i],b[i+lag])
        rows.append(dict(lag=int(lag),rho=rr,n_pairs=n,status=reason))
    return rows
def curve_summary(rows):
    a=np.array([r['rho'] for r in rows]);ok=np.isfinite(a);r0=a[50];far=a[FAR]
    spec=r0-np.median(far) if np.isfinite(r0) and np.isfinite(far).all() else np.nan
    peak=sorted([int(x) for x in LAGS[a==a.max()]],key=lambda x:(abs(x),x))[0] if ok.all() else np.nan
    return dict(rho0=r0,n_pairs0=rows[50]['n_pairs'],rho0_status=rows[50]['status'],
                specificity=spec,specificity_status='ok' if np.isfinite(spec) else 'undefined_zero_or_far_lag',n_far_defined=int(np.isfinite(far).sum()),
                peak_lag=peak,peak_rho=float(a.max()) if ok.all() else np.nan,
                zero_lag_rank=int(1+(a>r0).sum()) if ok.all() else np.nan,
                peak_rank_status='ok' if ok.all() else 'incomplete_101_lags',n_lags_defined=int(ok.sum()),boundary_peak=bool(abs(peak)==50) if np.isfinite(peak) else None)
def matched_specificity(a,b):
    rows=[]
    for lag in LAGS[FAR]:
        i=np.arange(max(0,-lag),min(500,500-lag));m=np.isfinite(a[i])&np.isfinite(b[i])&np.isfinite(b[i+lag]);i=i[m]
        r0,n,s0=rho(a[i],b[i]);rl,_,sl=rho(a[i],b[i+lag]);delta=r0-rl
        rows.append(dict(lag=int(lag),n_pairs=n,rho0_matched=r0,rho_lag_matched=rl,delta=delta,status='ok' if s0==sl=='ok' else s0+';'+sl))
    return strict([r['delta'] for r in rows],42),rows
def cell_meta(s,t):return dict(source_bin=s,u_lo=(s-1)/6,u_hi=s/6,u_mid=(s-.5)/6,tail_bin=t,tail_lo_um=(t-1)*25,tail_hi_um=t*25,tail_mid_um=(t-.5)*25,in_C16=(s,t) in C)
def map_rows(x,meta,cells=None):
    rows=[]
    for s,t in cells or [(s,t) for s in range(1,7) for t in range(1,21)]:
        rr,n,reason=rho(x[f's{s}'],x[f't{t:02}']);rows.append(meta|cell_meta(s,t)|dict(rho=rr,n_pairs=n,status=reason))
    return rows
def region_rows(rows,meta):
    df=pd.DataFrame(rows);ans=[meta|dict(region='C16')|strict(df[df.in_C16].rho.to_numpy(),16)]
    for s in range(1,7):ans.append(meta|dict(region=f'S{s}_T1T4')|strict(df[(df.source_bin==s)&(df.tail_bin<=4)].rho.to_numpy(),4))
    return ans
def detrended_signals(sig,w,meta,support_rows):
    res={}
    for region,a in sig.items():
        res[region],trend,n=detrend(a,w)
        support_rows.extend(meta|dict(region=region,window=w,frame_index=i,center_available=bool(np.isfinite(a[i])),n_window=int(n[i]),trend=trend[i],residual=res[region][i]) for i in range(500))
    return res
