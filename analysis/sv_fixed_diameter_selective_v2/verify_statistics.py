"""Independent implementations: no scipy correlation test or p-value interface."""
from common import *
def reference_rho(a,b):
    pairs=[(float(x),float(y)) for x,y in zip(a,b) if np.isfinite(x) and np.isfinite(y)]
    if len(pairs)<3:return np.nan,len(pairs)
    x,y=map(list,zip(*pairs))
    if min(x)==max(x) or min(y)==max(y):return np.nan,len(pairs)
    def ranks(v):
        order=sorted(range(len(v)),key=lambda i:v[i]);r=[0.]*len(v);i=0
        while i<len(v):
            j=i+1
            while j<len(v) and v[order[j]]==v[order[i]]:j+=1
            for k in range(i,j):r[order[k]]=(i+1+j)/2
            i=j
        return np.array(r)
    return float(np.corrcoef(ranks(x),ranks(y))[0,1]),len(x)
def reference_trend(a,w):
    trend=[];n=[]
    for i in range(len(a)):
        vals=[x for x in a[max(0,i-w//2):min(len(a),i+w//2+1)] if np.isfinite(x)]
        trend.append(np.median(vals) if vals else np.nan);n.append(len(vals))
    return np.asarray(trend),np.asarray(n)
def assert_close(a,b):
    if np.isnan(a) or np.isnan(b):assert np.isnan(a) and np.isnan(b)
    else:assert abs(a-b)<=1e-12,(a,b)
def verify_pair(a,b,name,rows):
    curve=lag_curve(a,b)
    for r in curve:
        lag=r['lag'];pairs=[(a[i],b[i+lag]) for i in range(500) if 0<=i+lag<500]
        rr,n=reference_rho(*map(np.array,zip(*pairs)));assert_close(r['rho'],rr);assert n==r['n_pairs']
        rows.append(dict(case=name,module='lag',lag=lag,absolute_error=abs(rr-r['rho']) if np.isfinite(rr) else np.nan,n_pairs=n,status='passed'))
    ms,details=matched_specificity(a,b);ds=[]
    for lag in LAGS[FAR]:
        ids=[i for i in range(500) if 0<=i+lag<500 and np.isfinite(a[i]) and np.isfinite(b[i]) and np.isfinite(b[i+lag])]
        rr,_=reference_rho(a[ids],b[ids]);rl,_=reference_rho(a[ids],b[np.array(ids,dtype=int)+lag]);ds.append(rr-rl)
    expected=np.median(ds) if np.isfinite(ds).all() else np.nan;assert_close(ms['value'],expected)
    for got,exp in zip(details,ds):assert_close(got['delta'],exp)
    da,db=np.diff(a),np.diff(b);rr,n,s=rho(da,db);ref,n0=reference_rho(da,db);assert_close(rr,ref);assert n==n0
    rows.append(dict(case=name,module='matched_specificity_and_adjacent_difference',absolute_error=0,status='passed'))
def main():
    rows=[];rng=np.random.default_rng(98434);a=rng.normal(size=500);b=rng.normal(size=500)
    cases={'random':(a,b),'ties':(np.round(a),np.round(b)),'constant':(np.ones(500),b),'empty':(np.full(500,np.nan),b),'negative':(a,-a)}
    x=a.copy();y=b.copy();x[[0,1,6,80,250,251,499]]=np.nan;y[[1,3,80,252,498]]=np.nan;cases['missing_internal_endpoints']=(x,y)
    sh=np.full(500,np.nan);sh[7:]=a[:-7];cases['known_shift_7']=(a,sh)
    for name,(x,y) in cases.items():
        verify_pair(x,y,name,rows)
        for w in [31,51,101]:
            res,t,n=detrend(x,w);rt,rn=reference_trend(x,w)
            np.testing.assert_allclose(t,rt,rtol=0,atol=0,equal_nan=True);assert np.array_equal(n,rn)
            assert np.array_equal(np.isnan(res),np.isnan(x));rows.append(dict(case=name,module='detrend',window=w,status='passed',absolute_error=0))
    cs=curve_summary(lag_curve(a,sh));assert cs['peak_lag']==7 and cs['peak_rho']>=1-1e-12
    tie=[dict(lag=int(l),rho=0.5 if l in [-1,1] else 0.,n_pairs=500-abs(int(l)),status='ok') for l in LAGS]
    assert curve_summary(tie)['peak_lag']==-1
    tie[0]['rho']=np.nan;cs=curve_summary(tie);assert np.isnan(cs['peak_lag']) and np.isnan(cs['specificity'])
    assert np.isnan(strict([.1]*15+[np.nan],16)['value']) and strict([-1.]*16,16)['value']==-1
    fw=frame_table()
    for d,g in fw[fw.flow_mm_s.eq(1)].groupby('diameter_um'):
        sig=signals(g)
        for rep,(s,t) in REP.items():
            for w in [0,51]:
                x,y=sig[s],sig[t]
                if w:x=detrend(x,w)[0];y=detrend(y,w)[0]
                verify_pair(x,y,f'D{d}_{rep}_{w}',rows)
                rr,_,_=rho(x,y);scaled,_,_=rho(x*123.7,y*997.3);assert_close(rr,scaled)
        # Common support must be imposed before trend; test differing missing supports.
        left=sig['source'].copy();right=sig['tail100'].copy();left[70:85]=np.nan;right[400:420]=np.nan
        m=np.isfinite(sig['source'])&np.isfinite(left)&np.isfinite(right)
        counts=[]
        for z in [sig['source'],left,right]:
            z=np.where(m,z,np.nan);r,t,n=detrend(z,51);rt,rn=reference_trend(z,51);np.testing.assert_allclose(r,z-rt,rtol=0,atol=0,equal_nan=True);assert np.array_equal(n,rn);counts.append(n)
        assert np.array_equal(counts[0],counts[1]) and np.array_equal(counts[1],counts[2])
    csv('statistical_independent_checks.csv.gz',rows)
    js('statistical_validation.json',dict(status='passed',cases=len(cases),production_volumes=3,checks=len(rows),rho_tolerance=1e-12,
        covered=['lag original coordinates and pairs','missing endpoints and internal','ties','constant','empty','signed correlation','known shift','positive scaling','strict summaries','peak tie and undefined far lag','common mask before trend','raw adjacent difference','matched specificity','window counts and center missingness']))
    print('Independent statistics passed:',len(rows),flush=True)
if __name__=='__main__':main()
