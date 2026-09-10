"""Table-driven v2 real statistics; old values used only for paired geometry comparison."""
from common import *
def metrics(sig):return {r:sig[r] for r in ['source','tail100','tail500']}|{'RI100':ratio(sig['tail100'],sig['source']),'RI500':ratio(sig['tail500'],sig['source'])}
def main():
    assert json.loads((OUT/'statistical_validation.json').read_text())['status']=='passed'
    fw=frame_table();vol=[];profiles=[];blocks=[];maps=[];regions=[];lags=[];whole=[];matched=[];diff=[];supports=[];sensitivity=[];senscells=[];qchecks=[];recon=[];areas=[];rirows=[];geom=[];jitters=[]
    oldwhole=pd.concat([read(PREV/'volume_zero_lag_coupling.csv'),read(PREV/'volume_detrended_coupling.csv')])
    oldmaps=read(DEP/'coupling_maps_all_resolutions_modes.csv.gz')
    oldmaps=oldmaps[oldmaps.scan_id.isin(fw.scan_id.unique())&oldmaps.coordinate.eq('absolute')&oldmaps.source_bins.eq(6)&oldmaps.metric.eq('mean')&oldmaps['mode'].isin(['raw','detrended_51'])]
    oldframe=read(PREV/'framewise_source_tail_metrics.csv');oldframe=oldframe[oldframe.scan_id.isin(fw.scan_id.unique())]
    assert set(map(tuple,oldframe[KEY].values))==set(map(tuple,fw.loc[fw.valid_geometry,KEY].values))
    for key,g in fw.groupby(ID,sort=True):
        meta=dict(zip(ID,key));scan=meta['scan_id'];g=g.sort_values('frame_index');sig=signals(g);qs=signals(g,'q');ar=signals(g,'area');valid=g.valid_geometry.to_numpy()
        print('Table statistics',scan,flush=True)
        for target,parts in [('source',[f's{i}' for i in range(1,7)]),('tail100',[f't{i:02}' for i in range(1,5)]),('tail500',[f't{i:02}' for i in range(1,21)])]:
            for name,vals in [('q',qs),('area',ar)]:
                a=sum(vals[p] for p in parts)[valid];b=vals[target][valid];err=np.max(np.abs((a-b)/b));assert err<1e-10
                recon.append(meta|dict(region=target,measure=name,max_relative_error=err,status='passed'))
        mm=metrics(sig)
        for name,a in mm.items():vol.append(meta|dict(metric=name,unit='dimensionless' if name.startswith('RI') else 'raw_SV_instrument_units',n_geometry_valid=int(valid.sum()))|summary(a))
        for i in range(500):rirows.append(meta|dict(frame_index=i,valid_geometry=bool(valid[i]),RI100=mm['RI100'][i],RI500=mm['RI500'][i],ri_status='invalid_geometry' if not valid[i] else 'zero_source' if sig['source'][i]==0 else 'ok'))
        for t in range(1,21):
            name=f't{t:02}'
            for metric,a in [('tail_mean',sig[name]),('RI',ratio(sig[name],sig['source']))]:profiles.append(meta|dict(tail_bin=t,tail_lo_um=25*(t-1),tail_hi_um=25*t,tail_mid_um=25*(t-.5),metric=metric)|summary(a))
        for lo in range(0,500,100):
            for name,a in mm.items():blocks.append(meta|dict(frame_lo=lo,frame_hi=lo+99,metric=name,n_geometry_valid=int(valid[lo:lo+100].sum()))|summary(a[lo:lo+100]))
        for region,a in ar.items():
            d=meta['diameter_um']
            if region=='source':theory=np.pi*d*d/4
            elif region.startswith('tail'):theory=d*int(region[4:])
            elif region.startswith('t'):theory=d*25
            else:
                lo,hi=(2/6,1) if region=='s36' else ((int(region[1:])-1)/6,int(region[1:])/6)
                def primitive(u):
                    z=2*u-1;return (z*np.sqrt(max(0,1-z*z))+np.arcsin(z))/2
                theory=d*d/2*(primitive(hi)-primitive(lo))
            for i in range(500):areas.append(meta|dict(frame_index=i,valid_geometry=bool(valid[i]),region=region,analytic_area_um2=theory,actual_area_um2=a[i],relative_area_error=ratio(a[i]-theory,theory).item()))
        for parameter in ['X4','z_top']:
            z=g[parameter].to_numpy();a=np.diff(z);m=np.isfinite(a)
            for i in np.flatnonzero(m):jitters.append(meta|dict(parameter=parameter,frame_i=int(i),frame_next=int(i+1),signed_delta_px=a[i],absolute_delta_px=abs(a[i]),signed_delta_um=a[i]*(DX if parameter=='X4' else DZ)))
        for mode,w in [('raw',0),('detrended_51',51)]:
            context=meta|dict(mode=mode,support='real_all')
            x=sig if not w else detrended_signals(sig,w,context,supports)
            mr=map_rows(x,context);maps.extend(mr);regions.extend(region_rows(mr,context))
            for rep,(s,t) in REP.items():
                rs=lag_curve(x[s],x[t]);base=context|dict(pair=rep);lags.extend(base|r for r in rs);cs=curve_summary(rs)
                if w:
                    ms,det=matched_specificity(x[s],x[t]);cs.update(matched_specificity=ms['value'],matched_status=ms['status'],n_matched_far_defined=ms['n_defined']);matched.extend(base|r for r in det)
                whole.append(base|cs)
                qx,qy=qs[s],qs[t]
                if w:qx=detrend(qx,w)[0];qy=detrend(qy,w)[0]
                qr,n,qstatus=rho(qx,qy);rr,_,_=rho(x[s],x[t]);theory_x=float(np.median(ar[s][valid]));theory_y=meta['diameter_um']*(100 if t=='tail100' else 500)
                # Use analytic source/pooled source constant from the already recorded analytic support, not actual fluctuating area.
                theory_x=next(r['analytic_area_um2'] for r in areas if r['scan_id']==scan and r['region']==s)
                tx=sig[s]*theory_x;ty=sig[t]*theory_y
                if w:tx=detrend(tx,w)[0];ty=detrend(ty,w)[0]
                tr,_,_=rho(tx,ty);assert abs(tr-rr)<1e-12 or (np.isnan(tr) and np.isnan(rr))
                qchecks.append(base|dict(mean_mean_rho=rr,actual_q_q_rho=qr,actual_q_minus_mean=qr-rr,theoretical_q_q_rho=tr,theoretical_scale_error=tr-rr,n_pairs=n,status=qstatus))
            old=oldmaps[(oldmaps.scan_id==scan)&(oldmaps['mode']==mode)];assert len(old)==120
            oldc=old[(old.source_bin>=3)&(old.tail_bin<=4)].rho.to_numpy();newc=np.array([r['rho'] for r in mr if r['in_C16']]);os=strict(oldc,16);ns=strict(newc,16)
            geom.append(context|dict(metric='C16',old_value=os['value'],new_value=ns['value'],new_minus_old=ns['value']-os['value'],old_status=os['status'],new_status=ns['status']))
            for rep,h in [('W100',100),('W500',500)]:
                row=oldwhole[(oldwhole.scan_id==scan)&(oldwhole['mode']==mode)&oldwhole.predictor.eq('source_mean_raw')&oldwhole.tail_depth_um.eq(h)&oldwhole.metric_family.eq('mean')];assert len(row)==1
                nv=next(r['rho0'] for r in whole if r['scan_id']==scan and r['mode']==mode and r['pair']==rep)
                geom.append(context|dict(metric=rep,old_value=row.rho.iloc[0],new_value=nv,new_minus_old=nv-row.rho.iloc[0],old_status='saved_v1',new_status='ok' if np.isfinite(nv) else 'undefined'))
        for rep,(s,t) in REP.items():
            rr,n,reason=rho(np.diff(sig[s]),np.diff(sig[t]));diff.append(meta|dict(pair=rep,rho=rr,n_pairs=n,status=reason,definition='raw differences, original i/i+1 only'))
        for w in [31,101]:
            need=['source','s36','tail100','tail500']+[f's{s}' for s in range(3,7)]+[f't{t:02}' for t in range(1,5)]
            context=meta|dict(mode=f'detrended_{w}',support='real_all');x=detrended_signals({r:sig[r] for r in need},w,context,supports)
            rows=map_rows(x,context,C);senscells.extend(rows);sensitivity.append(context|dict(pair='C16')|strict([r['rho'] for r in rows],16))
            for rep,(s,t) in REP.items():
                rs=lag_curve(x[s],x[t]);lags.extend(context|dict(pair=rep)|r for r in rs);sensitivity.append(context|dict(pair=rep)|curve_summary(rs))
        for h in [100,500]:
            oldri=ratio(g[f'old_tail_mean_raw_{h}um'],g.old_source_mean_raw);newri=mm[f'RI{h}'];vr=ratio(newri,oldri);right=ratio(ratio(sig[f'tail{h}'],g[f'old_tail_mean_raw_{h}um']),ratio(sig['source'],g.old_source_mean_raw))
            m=np.isfinite(vr)&np.isfinite(right);assert np.max(abs(vr[m]-right[m]))<1e-10
            geom.append(meta|dict(mode='per_frame_then_median',metric=f'RI{h}',old_value=summary(oldri)['median'],new_value=summary(newri)['median'],new_minus_old=summary(newri)['median']-summary(oldri)['median'],paired_percent_change_median=summary((vr-1)*100)['median'],identity_max_absolute_error=np.max(abs(vr[m]-right[m])),old_status='derived_from_saved_same_identity_means',new_status='ok'))
    vol=pd.DataFrame(vol);contr=[]
    for (f,metric),g in vol.groupby(['flow_mm_s','metric']):
        a=g.set_index('diameter_um')
        for d1,d2 in [(128,235),(128,285),(235,285)]:
            v1,v2=a.loc[d1,'median'],a.loc[d2,'median'];contr.append(dict(flow_mm_s=f,metric=metric,reference_D=d1,comparison_D=d2,reference_value=v1,comparison_value=v2,difference=v2-v1,percent_change=ratio((v2-v1)*100,v1).item(),status='ok' if np.isfinite(v1) and v1!=0 and np.isfinite(v2) else 'undefined_denominator_or_metric'))
    mapsdf=pd.DataFrame(maps);cons=[]
    for (d,mode,s,t),g in mapsdf.groupby(['diameter_um','mode','source_bin','tail_bin']):
        assert len(g)==5
        cons.append(dict(diameter_um=d,mode=mode)|cell_meta(s,t)|summary(g.rho)|{f'flow{int(r.flow_mm_s):02}_rho':r.rho for r in g.itertuples()})
    for name,rows in [('volume_metrics.csv',vol),('diameter_contrasts_by_flow.csv',contr),('raw_depth_profiles.csv',profiles),('spatial_block_summary.csv',blocks),('depth_maps_raw51.csv',maps),('fixed_region_summary.csv',regions),('representative_lag.csv.gz',lags),('whole_coupling.csv',whole),('representative_matched_lag.csv.gz',matched),('adjacent_difference.csv',diff),('detrend_support_real.csv.gz',supports),('detrend_sensitivity.csv',sensitivity),('detrend_sensitivity_cells.csv',senscells),('q_mean_numerical_check.csv',qchecks),('reconstruction_checks.csv',recon),('geometry_area.csv.gz',areas),('framewise_RI.csv.gz',rirows),('diameter_consensus_maps.csv',cons),('geometry_version_comparison.csv',geom),('geometry_jitter_adjacent.csv.gz',jitters)]:csv(name,rows)
    csv('geometry_version_existing_migration_summary.csv',read(INPUT/'volume_geometry_change_summary.csv'))
    js('table_validation.json',dict(status='passed',main_map_rows=len(maps),expected_map_rows=3600,nominal_rows=len(fw),valid_rows=int(fw.valid_geometry.sum()),reconstruction_checks=len(recon),no_old_arrays_read=True))
    status('table_statistics_complete')
if __name__=='__main__':main()
