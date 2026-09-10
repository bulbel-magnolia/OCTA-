"""Matched support controls and prespecified small rigid localization sensitivity."""
from common import *
from analyze_tables import metrics
from verify_statistics import reference_rho,assert_close,verify_pair
def main():
    assert json.loads((OUT/'pixel_validation.json').read_text())['status']=='passed'
    fw=frame_table();pmaps=[];cells=[];excess=[];curves=[];coupling=[];msrows=[];cross=[];bgcorr=[];drift=[];driftcurves=[];profiles=[];support=[];posmetrics=[];pospairs=[];poscells=[];posregions=[];poscurves=[];possupports=[];masks=[];checks=[];extensions=[];extrows=[];extpairs=[];supsummary=[];ovsummary=[]
    for key,g in fw.groupby(ID,sort=True):
        meta=dict(zip(ID,key));scan=meta['scan_id'];print('Controls and sensitivity',scan,flush=True)
        bg=read(OUT/f'framewise_local_background_qc_{scan}.csv.gz');p=read(OUT/f'position_integrals_{scan}.csv.gz')
        assert np.array_equal(bg.frame_index,np.arange(500)) and np.array_equal(p.frame_index,np.arange(500))
        assert np.array_equal(bg.valid_geometry,g.valid_geometry) and np.array_equal(p.valid_geometry,g.valid_geometry)
        sig=signals(g);M=g.valid_geometry.to_numpy()&bg.both_valid.to_numpy();allraw={'real_matched':sig}|{side:{r:bg[f'{r}_{side}_mean'].to_numpy() for r in REGIONS} for side in ['left','right']}
        allraw={side:{r:np.where(M,a,np.nan) for r,a in ss.items()} for side,ss in allraw.items()}
        det={side:detrended_signals(ss,51,meta|dict(side=side,support='real_L_R_common'),support) for side,ss in allraw.items()}
        mm={}
        for side,ss in det.items():
            context=meta|dict(side=side,mode='detrended_51',support='real_L_R_common',n_common=int(M.sum()))
            mr=map_rows(ss,context);pmaps.extend(mr);mm[side]=pd.DataFrame(mr)
            for rep,(s,t) in REP.items():
                rs=lag_curve(ss[s],ss[t]);cs=curve_summary(rs);msp,msdet=matched_specificity(ss[s],ss[t]);cs.update(matched_specificity=msp['value'],matched_status=msp['status'],n_matched_far_defined=msp['n_defined'])
                coupling.append(context|dict(pair=rep)|cs);curves.extend(context|dict(pair=rep)|r for r in rs);msrows.extend(context|dict(pair=rep)|r for r in msdet)
                if meta['flow_mm_s']==1:verify_pair(ss[s],ss[t],scan+'_'+side+'_'+rep,checks)
            if meta['flow_mm_s']==1:
                for r in mr:
                    rr,n=reference_rho(ss[f's{r["source_bin"]}'],ss[f't{r["tail_bin"]:02}']);assert_close(rr,r['rho']);assert n==r['n_pairs']
        rmap,lmap,rmap0=mm['real_matched'],mm['left'],mm['right']
        c=rmap[ID+['source_bin','tail_bin','in_C16']].copy();c['real_matched_rho']=rmap.rho;c['left_rho']=lmap.rho;c['right_rho']=rmap0.rho
        c['control_rho']=(c.left_rho+c.right_rho)/2;c['excess_cell']=c.real_matched_rho-c.control_rho;c['real_minus_L']=c.real_matched_rho-c.left_rho;c['real_minus_R']=c.real_matched_rho-c.right_rho;c['n_common']=int(M.sum())
        c['status']=np.where(c[['real_matched_rho','left_rho','right_rho']].notna().all(axis=1),'ok','undefined_real_or_side_rho');cells.extend(c.to_dict('records'))
        reg=c[c.in_C16];er=meta|dict(n_common=int(M.sum()))
        for field in ['real_matched_rho','left_rho','right_rho','control_rho','excess_cell','real_minus_L','real_minus_R']:
            z=strict(reg[field].to_numpy(),16)
            for k,v in z.items():er[field+'_'+k]=v
        er['side_dependent_sign']=bool(er['real_minus_L_value']*er['real_minus_R_value']<0);excess.append(er)
        if er['side_dependent_sign']:
            extensions.append(meta|dict(extension='C',trigger='C16 median cellwise real-minus-L and real-minus-R have opposite signs',evidence_table='fixed_region_control_excess.csv',evidence=f"real-L={er['real_minus_L_value']:.17g}; real-R={er['real_minus_R_value']:.17g}",scope='once only +/-2.5D same shape and depth; fixed region / three representative pairs / all-depth background',status='triggered_pending_FOV_check'))
        for rep,(s,t) in REP.items():
            for side in ['left','right']:
                rr,n,reason=rho(det['real_matched'][s],det[side][t]);cross.append(meta|dict(pair=rep,side=side,rho=rr,n_pairs=n,status=reason,support='real_L_R_common'))
        for t in range(1,21):
            region=f't{t:02}';b=bg[f'{region}_bg_mean'].to_numpy();a=sig[region];a,b=np.where(M,a,np.nan),np.where(M,b,np.nan)
            for mode,w in [('raw',0),('detrended_51',51)]:
                aa,bb=(detrend(a,w)[0],detrend(b,w)[0]) if w else (a,b);rr,n,reason=rho(aa,bb)
                bgcorr.append(meta|dict(tail_bin=t,tail_lo_um=(t-1)*25,tail_hi_um=t*25,tail_mid_um=(t-.5)*25,mode=mode,rho=rr,n_pairs=n,status=reason,support='real_L_R_common'))
        for region in REGIONS:
            for metric in ['z_local','real_minus_bg_mean','signed_asymmetry']:
                a=bg[f'{region}_{metric}'];profiles.append(meta|dict(region=region,metric=metric,n_geometry_valid=int(g.valid_geometry.sum()))|summary(a))
            for side in ['real','left','right','bg']:
                for field in ['mean','median','mad']:
                    profiles.append(meta|dict(region=region,metric=f'{side}_{field}',n_geometry_valid=int(g.valid_geometry.sum()))|summary(bg[f'{region}_{side}_{field}']))
        for region in ['source','tail100','tail500']:
            for side in ['left','right','bg']:
                b=bg[f'{region}_{side}_median'].to_numpy();sm=summary(b);med=sm['median'];mad=np.nanmedian(abs(b-med));first=b[:51];last=b[449:];a0=summary(first)['median'];a1=summary(last)['median'];diff=a1-a0
                row=meta|dict(region=region,side=side,n_defined=sm['n_defined'],curve_median=med,curve_MAD=mad,RCV=ratio(1.4826*mad,med).item(),first_median=a0,last_median=a1,first_coverage=int(np.isfinite(first).sum()),last_coverage=int(np.isfinite(last).sum()),signed_difference=diff,absolute_difference=abs(diff),relative_difference=ratio(diff,med).item(),status='ok' if np.isfinite(med) and med!=0 and np.isfinite(a0) and np.isfinite(a1) else 'zero_denominator_or_empty_endpoint')
                drift.append(row);res,tr,n=detrend(b,51)
                driftcurves.extend(meta|dict(region=region,side=side,frame_index=i,raw_pixel_median=b[i],trend_51=tr[i],residual_51=res[i],n_window=int(n[i]),center_available=bool(np.isfinite(b[i]))) for i in range(500))
        sp=read(OUT/f'geometry_support_{scan}.csv.gz');ov=read(OUT/f'shared_pixel_support_{scan}.csv.gz')
        for (region,side),gg in sp.groupby(['region','side']):
            for metric in ['neff','nonzero_pixels','sum_w2','z_span_pixels']:supsummary.append(meta|dict(region=region,side=side,metric=metric)|summary(gg[metric]))
        for (s,t),gg in ov.groupby(['source_bin','tail_bin']):ovsummary.append(meta|cell_meta(s,t)|dict(n_shared_positive=int((gg.shared_pixels>0).sum()),max_shared_pixels=int(gg.shared_pixels.max()),max_overlap_index=gg.overlap_index.max(),median_overlap_index=gg.overlap_index.median(),n_frames=len(gg)))
        hit=ov[(ov.tail_bin==1)&(ov.overlap_index>0)]
        if len(hit):
            extensions.append(meta|dict(extension='A',trigger='source-T1 shared pixels detected',evidence_table=f'shared_pixel_support_{scan}.csv.gz',evidence=f'{len(hit)} source-bin/frame records; first frame {int(hit.frame_index.min())}',scope='fixed S3-S6 x T2-T4 12-cell median and pooled S3-S6 vs tail25-100 51-frame',status='triggered_before_execution'))
        # All four positions intersect with baseline; primary valid flag remains untouched.
        PM=g.valid_geometry.to_numpy().copy()
        for name in OFFSETS:PM &= p[name+'_valid'].to_numpy()
        masks.extend(meta|dict(frame_index=i,valid_geometry=bool(g.valid_geometry.iloc[i]),background_common=bool(M[i]),position_common=bool(PM[i])) for i in range(500))
        need=['source','s36']+[f's{s}' for s in range(1,7)]+[f't{t:02}' for t in range(1,5)]+['tail100','tail500']
        poss={'baseline':{r:sig[r] for r in need}}|{name:{r:p.get(f'{name}_{r}_mean',pd.Series(np.nan,index=range(500))).to_numpy() for r in need} for name in OFFSETS}
        pcm={};pcmrows={};pairrows={}
        for name,ss in poss.items():
            for selection in ['all_available','five_position_common']:
                xx=ss if selection=='all_available' else {r:np.where(PM,a,np.nan) for r,a in ss.items()}
                for metric,a in metrics(xx).items():posmetrics.append(meta|dict(offset=name,dx_pixel=0 if name=='baseline' else OFFSETS[name][0],dz_pixel=0 if name=='baseline' else OFFSETS[name][1],support=selection,metric=metric)|summary(a))
            xx={r:np.where(PM,a,np.nan) for r,a in ss.items()};ctx=meta|dict(offset=name,support='five_position_common',mode='detrended_51',n_common=int(PM.sum()))
            d=detrended_signals(xx,51,ctx,possupports);pcm[name]=d;mr=map_rows(d,ctx,C);pcmrows[name]=mr
            for rep,(s,t) in REP.items():
                rs=lag_curve(d[s],d[t]);row=ctx|dict(pair=rep)|curve_summary(rs);pairrows[(name,rep)]=row;poscurves.extend(ctx|dict(pair=rep)|r for r in rs)
                if meta['flow_mm_s']==1:verify_pair(d[s],d[t],scan+'_offset_'+name+'_'+rep,checks)
            if meta['flow_mm_s']==1:
                for row in mr:
                    rr,n=reference_rho(d[f's{row["source_bin"]}'],d[f't{row["tail_bin"]:02}']);assert_close(rr,row['rho']);assert n==row['n_pairs']
        for name,mr in pcmrows.items():
            for row,base in zip(mr,pcmrows['baseline']):row['baseline_common_rho']=base['rho'];row['rho_change']=row['rho']-base['rho'];poscells.append(row)
            st=strict([r['rho'] for r in mr],16);bst=strict([r['rho'] for r in pcmrows['baseline']],16)
            posregions.append(meta|dict(offset=name,n_common=int(PM.sum()),baseline_common_C16=bst['value'],C16_change=st['value']-bst['value'])|st)
        for (name,rep),row in pairrows.items():
            base=pairrows[('baseline',rep)];row.update(baseline_common_rho0=base['rho0'],rho0_change=row['rho0']-base['rho0'],baseline_common_specificity=base['specificity'],specificity_change=row['specificity']-base['specificity']);pospairs.append(row)
    # Persist triggers before any extension results are computed.
    for letter in ['B','D','E']:extensions.append(dict(extension=letter,scan_id='ALL',trigger='awaiting visual review or reporting dependency',scope='exact frozen conditional scope in analysis_specification.md',status='pending_trigger_assessment'))
    csv('conditional_extension_log.csv',extensions)
    for key,g in fw.groupby(ID,sort=True):
        meta=dict(zip(ID,key));scan=meta['scan_id']
        if not any(r.get('extension')=='A' and r.get('scan_id')==scan for r in extensions):continue
        sig=signals(g);d={r:detrend(a,51)[0] for r,a in sig.items()};ctx=meta|dict(extension='A',mode='detrended_51',support='real_all');rr=map_rows(d,ctx,[(s,t) for s in range(3,7) for t in range(2,5)]);extrows.extend(rr)
        st=strict([r['rho'] for r in rr],12);q=sum(g[col(f't{t:02}','q')].to_numpy() for t in range(2,5));a=sum(g[col(f't{t:02}','area')].to_numpy() for t in range(2,5));tail=detrend(q/a,51)[0];rho0,n,reason=rho(d['s36'],tail)
        extpairs.append(ctx|dict(C12=st['value'],C12_status=st['status'],C12_defined=st['n_defined'],pooled_S36_tail25_100_rho=rho0,n_pairs=n,status=reason))
        for r in extensions:
            if r.get('extension')=='A' and r.get('scan_id')==scan:r['status']='complete'
    for name,rows in [('pseudo_maps.csv',pmaps),('fixed_region_control_cells.csv',cells),('fixed_region_control_excess.csv',excess),('pseudo_representative_lag.csv.gz',curves),('pseudo_whole_coupling.csv',coupling),('pseudo_matched_lag.csv.gz',msrows),('crossed_controls.csv',cross),('background_real_tail_covariation.csv',bgcorr),('background_drift.csv',drift),('background_drift_curves.csv.gz',driftcurves),('local_contrast_profiles.csv',profiles),('detrend_support_controls.csv.gz',support),('position_sensitivity_metrics.csv',posmetrics),('position_sensitivity_pairs.csv',pospairs),('position_sensitivity_cells.csv',poscells),('position_sensitivity_fixed_region.csv',posregions),('position_sensitivity_lag.csv.gz',poscurves),('detrend_support_positions.csv.gz',possupports),('common_support_masks.csv.gz',masks),('control_position_independent_checks.csv.gz',checks),('geometry_support_summary.csv',supsummary),('shared_pixel_support_summary.csv',ovsummary),('conditional_A_cells.csv',extrows),('conditional_A_summary.csv',extpairs),('conditional_extension_log.csv',extensions)]:csv(name,rows)
    js('controls_validation.json',dict(status='passed',pseudo_and_real_matched_map_rows=len(pmaps),expected_map_rows=5400,excess_cell_rows=len(cells),n_C16_volumes=len(excess),production_independent_checks=len(checks),comparison_support='mask before detrend verified',conditional_C_triggered=sum(r.get('extension')=='C' for r in extensions),conditional_A_triggered=len(extpairs)))
    status('controls_and_position_statistics_complete',conditional_extensions='visual assessment and triggered C pending')
if __name__=='__main__':main()
