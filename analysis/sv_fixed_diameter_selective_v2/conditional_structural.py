from extract_pixels import *
from verify_statistics import reference_rho,assert_close
def main():
    logs=read(OUT/'conditional_extension_log.csv');todo=logs[logs.extension.eq('D')&logs.status.eq('triggered_before_execution')].scan_id.tolist()
    fw=frame_table();ids=read(INPUT/'array_identity_audit.csv.gz').fillna('').set_index(KEY);fixed=read(OUT/'localization_fixed_frames.csv');rows=[];receipts=[];checks=[];support=[];summaries=[]
    for key,g in fw[fw.scan_id.isin(todo)].groupby(ID,sort=True):
        meta=dict(zip(ID,key));scan=meta['scan_id'];print('Triggered structural diagnostic',scan,flush=True);cache={};rec=[]
        for r in g.itertuples():
            row=meta|dict(frame_index=r.frame_index,valid_geometry=r.valid_geometry,structural_status='invalid_geometry')
            if not r.valid_geometry:rec.append(row);continue
            arrays,h=load_array(ids.loc[(scan,r.frame_index)],cache);receipts.append(meta|dict(frame_index=r.frame_index,sha256=h,status='matched_on_read_for_trigger_D'))
            if 'stru_amp' not in arrays or not np.isfinite(arrays['stru_amp']).all():row['structural_status']='unavailable_or_nonfinite';rec.append(row);continue
            a=arrays['stru_amp'];row['structural_status']='ok'
            for side,dx in [('real',0),('left',-1.5*r.diameter_um/DX),('right',1.5*r.diameter_um/DX)]:
                gg=geometry(r,dx);ok=complete(gg);row[side+'_valid']=ok
                if not ok:continue
                mm=masks(gg,limited=True)
                if r.frame_index in set(fixed[fixed.scan_id==scan].frame_index):ref=independent_masks(gg,limited=True)
                for region in ['source','tail100','tail500']:
                    mu=measure(a,mm[region])['mean'];row[f'{side}_{region}_structural_mean']=mu
                    if r.frame_index in set(fixed[fixed.scan_id==scan].frame_index):
                        z,x,w=ref[region];refmu=np.sum(a[np.ix_(z,x)]*w)/np.sum(w);err=abs(mu-refmu)/abs(refmu) if refmu else abs(mu);assert err<1e-10
                        checks.append(meta|dict(frame_index=r.frame_index,side=side,region=region,relative_error=err,status='passed'))
            rec.append(row)
        for z in cache.values():z.close()
        data=pd.DataFrame(rec);rows.extend(rec);bg=read(OUT/f'framewise_local_background_qc_{scan}.csv.gz');sig=signals(g)
        M=g.valid_geometry.to_numpy()&bg.both_valid.to_numpy()&data.structural_status.eq('ok').to_numpy()
        for side in ['real','left','right']:
            rawst={r:np.where(M,data[f'{side}_{r}_structural_mean'].to_numpy(),np.nan) for r in ['source','tail100','tail500']}
            rawsv={r:np.where(M,sig[r] if side=='real' else bg[f'{r}_{side}_mean'].to_numpy(),np.nan) for r in rawst}
            for mode,w in [('raw',0),('detrended_51',51)]:
                if w:
                    st=detrended_signals(rawst,w,meta|dict(side=side,channel='structural',support='real_L_R_structural_common'),support)
                    sv=detrended_signals(rawsv,w,meta|dict(side=side,channel='SV',support='real_L_R_structural_common'),support)
                else:st,sv=rawst,rawsv
                for region in rawst:
                    rr,n,reason=rho(st[region],sv[region]);summaries.append(meta|dict(side=side,mode=mode,kind='same_ROI_structural_vs_SV',region=region,rho=rr,n_pairs=n,status=reason))
                    if meta['flow_mm_s']==1:ref,n0=reference_rho(st[region],sv[region]);assert_close(rr,ref);assert n==n0
                for region in ['tail100','tail500']:
                    rr,n,reason=rho(st['source'],st[region]);summaries.append(meta|dict(side=side,mode=mode,kind='structural_source_vs_tail',region=region,rho=rr,n_pairs=n,status=reason))
        logs.loc[logs.extension.eq('D')&logs.scan_id.eq(scan),'status']='complete'
    csv('conditional_D_structural_framewise.csv.gz',rows);csv('conditional_D_structural_correlations.csv',summaries);csv('conditional_D_identity_receipts.csv.gz',receipts);csv('conditional_D_integration_checks.csv',checks);csv('conditional_D_detrend_support.csv.gz',support)
    # B uses fixed maxima pairs; use the same scan display range as original panels.
    pairs=read(OUT/'conditional_B_fixed_pairs.csv');settings=read(OUT/'localization_display_settings.csv');breview=[];cache={}
    for r in pairs.itertuples():
        shots=[]
        for rr in fw[(fw.scan_id==r.scan_id)&fw.frame_index.isin([r.frame_i,r.frame_next])].itertuples():
            arrays,h=load_array(ids.loc[(rr.scan_id,rr.frame_index)],cache);shots.append((rr,arrays))
        limits=dict(zip(settings[settings.scan_id==r.scan_id].channel,settings[settings.scan_id==r.scan_id].linear_display_upper));name=r.scan_id+'_'+r.parameter
        draw_qc(name,shots,limits,folder='localization_extension_B');breview.append(dict(scan_id=r.scan_id,parameter=r.parameter,frame_i=r.frame_i,frame_next=r.frame_next,figure='localization_extension_B/'+name+'.png',review_status='pending_visual_review'))
    for z in cache.values():z.close()
    csv('conditional_B_review.csv',breview);csv('conditional_extension_log.csv',logs)
    js('conditional_validation.json',dict(status='statistics_passed_visual_B_pending',structural_volumes=len(todo),structural_rows=len(rows),structural_valid_arrays=len(receipts),structural_integral_checks=len(checks),structural_correlations=len(summaries),B_pairs=len(breview)))
    status('conditional_statistics_complete',conditional_B_visual_review='pending')
if __name__=='__main__':main()
