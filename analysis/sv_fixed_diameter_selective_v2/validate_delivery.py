"""Final definition/coverage gates; same-input checks are not independent experiments."""
from common import *
from verify_statistics import reference_rho,assert_close
def main():
    checks=[]
    def gate(name,ok,**detail):
        checks.append(dict(check=name,status='passed' if ok else 'failed',**detail))
        if not ok:js('validation.json',dict(status='failed',checks=checks));raise AssertionError(name)
    fw=frame_table();ids=read(OUT/'array_identity_on_read.csv.gz');original=read(INPUT/'array_identity_audit.csv.gz');joined=ids.merge(original[KEY+['sha256']],on=KEY,suffixes=('_read','_input'),validate='one_to_one')
    gate('array identity on all actual reads',len(joined)==7018 and joined.sha256_read.eq(joined.sha256_input).all(),arrays=len(joined))
    gate('cohort and nominal/valid identity',len(fw)==7500 and fw.valid_geometry.sum()==7018 and set(fw.diameter_um)=={128,235,285})
    coverage=[('depth_maps_raw51.csv',3600,ID+['mode','source_bin','tail_bin']),('pseudo_maps.csv',5400,ID+['side','source_bin','tail_bin']),('fixed_region_control_cells.csv',1800,ID+['source_bin','tail_bin']),('position_sensitivity_cells.csv',1200,ID+['offset','source_bin','tail_bin']),('volume_metrics.csv',75,ID+['metric']),('raw_depth_profiles.csv',600,ID+['tail_bin','metric']),('spatial_block_summary.csv',375,ID+['frame_lo','metric']),('whole_coupling.csv',90,ID+['mode','pair']),('representative_lag.csv.gz',18180,ID+['mode','pair','lag']),('pseudo_representative_lag.csv.gz',13635,ID+['side','pair','lag']),('position_sensitivity_lag.csv.gz',22725,ID+['offset','pair','lag']),('adjacent_difference.csv',45,ID+['pair']),('crossed_controls.csv',90,ID+['pair','side']),('background_real_tail_covariation.csv',600,ID+['mode','tail_bin']),('conditional_A_cells.csv',180,ID+['source_bin','tail_bin']),('conditional_D_structural_correlations.csv',450,ID+['side','mode','kind','region'])]
    for file,n,keys in coverage:
        a=read(OUT/file);gate(file+' complete keys',len(a)==n and not a.duplicated(keys).any() and set(a.scan_id)==set(fw.scan_id),expected_rows=n,actual_rows=len(a))
        if 'rho' in a:gate(file+' finite coefficients have reasons',a.loc[a.rho.notna(),'rho'].between(-1-1e-12,1+1e-12).all() and a.status.notna().all())
    mainmaps=read(OUT/'depth_maps_raw51.csv');controls=read(OUT/'fixed_region_control_cells.csv');exc=read(OUT/'fixed_region_control_excess.csv');reg=read(OUT/'fixed_region_summary.csv');vm=read(OUT/'volume_metrics.csv')
    ncheck=0
    for key,g in fw.groupby(ID):
        meta=dict(zip(ID,key));scan=meta['scan_id'];sig=signals(g);m=g.valid_geometry.to_numpy();bg=read(OUT/f'framewise_local_background_qc_{scan}.csv.gz');p=read(OUT/f'position_integrals_{scan}.csv.gz')
        gate(scan+' pixel and offset nominal skeleton',len(bg)==len(p)==500 and np.array_equal(bg.frame_index,g.frame_index) and np.array_equal(bg.valid_geometry,m) and np.array_equal(p.valid_geometry,m))
        gate(scan+' missing centers not zero-filled',bg.loc[~m,'source_real_mean'].isna().all() and bg.loc[~m,'tail500_bg_median'].isna().all())
        for region in REGIONS:
            both=bg.left_valid&bg.right_valid;z=ratio(bg[f'{region}_real_median']-bg[f'{region}_bg_median'],1.4826*bg[f'{region}_bg_mad']);np.testing.assert_allclose(z,bg[f'{region}_z_local'],rtol=1e-14,atol=1e-14,equal_nan=True)
            pool=ratio(bg[f'{region}_left_q']+bg[f'{region}_right_q'],bg[f'{region}_left_area']+bg[f'{region}_right_area']);np.testing.assert_allclose(pool[both],bg.loc[both,f'{region}_bg_mean'],rtol=1e-10,atol=0,equal_nan=True)
        for h in [100,500]:
            v=np.median([t/s for s,t in zip(sig['source'],sig[f'tail{h}']) if np.isfinite(s) and np.isfinite(t) and s!=0]);v0=vm[(vm.scan_id==scan)&vm.metric.eq(f'RI{h}')]['median'].iloc[0];gate(scan+f' median-of-frame RI{h}',abs(v-v0)<1e-14)
        c=controls[controls.scan_id==scan];np.testing.assert_allclose(c.control_rho,(c.left_rho+c.right_rho)/2,rtol=0,atol=0);np.testing.assert_allclose(c.excess_cell,c.real_matched_rho-c.control_rho,rtol=0,atol=0)
        e=exc[exc.scan_id==scan].iloc[0]
        for field in ['real_matched_rho','left_rho','right_rho','excess_cell','real_minus_L','real_minus_R']:assert_close(e[field+'_value'],strict(c.loc[c.in_C16,field].to_numpy(),16)['value'])
        if meta['flow_mm_s']==1:
            for mode,w in [('raw',0),('detrended_51',51)]:
                x={r:detrend(a,w)[0] for r,a in sig.items()} if w else sig
                a=mainmaps[(mainmaps.scan_id==scan)&mainmaps['mode'].eq(mode)]
                for row in a.itertuples():rr,n=reference_rho(x[f's{row.source_bin}'],x[f't{row.tail_bin:02}']);assert_close(row.rho,rr);assert n==row.n_pairs;ncheck+=1
            a=read(OUT/'conditional_A_summary.csv');row=a[a.scan_id==scan].iloc[0];q=sum(g[col(f't{t:02}','q')].to_numpy() for t in range(2,5));ar=sum(g[col(f't{t:02}','area')].to_numpy() for t in range(2,5));rr,n=reference_rho(detrend(sig['s36'],51)[0],detrend(q/ar,51)[0]);assert_close(rr,row.pooled_S36_tail25_100_rho);assert n==row.n_pairs
    gate('independent final production map and conditional A checks',True,cell_checks=ncheck,representative_volumes=3)
    for file,group in [('detrend_support_controls.csv.gz',ID+['frame_index','window']),('detrend_support_positions.csv.gz',ID+['frame_index','mode'])]:
        a=read(OUT/file);v=a.groupby(group).n_window.nunique();gate(file+' all compared signals share support',v.eq(1).all())
    cc=read(OUT/'independent_integration_checks.csv.gz');nonzero=~cc.zero_reference;gate('independent 45-frame integrations',cc.loc[nonzero,'relative_error'].lt(1e-10).all() and cc.loc[~nonzero,'absolute_error'].eq(0).all(),checks=len(cc),max_relative_error=cc.relative_error.max(),zero_reference_checks=int((~nonzero).sum()),fixed_frames=cc[KEY].drop_duplicates().shape[0])
    for name in ['statistical_validation.json','pixel_validation.json','table_validation.json','controls_validation.json','conditional_validation.json']:
        gate(name,json.loads((OUT/name).read_text())['status']=='passed')
    logs=read(OUT/'conditional_extension_log.csv');gate('all conditional triggers resolved',logs.status.isin(['complete','not_triggered']).all())
    loc=read(OUT/'localization_review.csv');b=read(OUT/'conditional_B_review.csv');gate('all fixed and triggered images reviewed',len(loc)==45 and len(b)==4 and not loc.review_status.str.contains('pending').any() and not b.review_status.str.contains('pending').any())
    unchanged=[]
    for r in read(OUT/'input_manifest.csv').itertuples():
        h=sha(ROOT/r.input_path);unchanged.append(dict(input_path=r.input_path,before_sha256=r.sha256,after_sha256=h,unchanged=h==r.sha256))
    csv('input_hashes_after.csv',unchanged);gate('all inputs and historical results unchanged',all(x['unchanged'] for x in unchanged),files=len(unchanged))
    gate('frozen contract unchanged',sha(OUT/'analysis_plan.json')==json.loads((OUT/'provenance.json').read_text())['contract_sha256'])
    csv('cohort_primary.csv',fw[ID].drop_duplicates().merge(fw.groupby(ID).valid_geometry.sum().reset_index(name='valid_frames'),on=ID));csv('d500_exclusion_record.csv',read(INPUT/'d500_exclusion_record.csv'))
    coverage_rows=[]
    for file,n,keys in coverage:coverage_rows.append(dict(file=file,expected_rows=n,actual_rows=len(read(OUT/file)),key=' + '.join(keys),status='passed'))
    csv('coverage_ledger.csv',coverage_rows)
    js('validation.json',dict(status='passed',checks=checks,all_required_modules_complete=True,nominal_frames=7500,valid_frames=7018,original_invalid_frames=482,
        scientific_independent_validation=False,numerical_independent_implementations=True,anatomical_localization_truth='unresolved; explicit result, not computational blocker',p_values=False,D500_arrays_read=0,
        sign_and_undefined_results_preserved=True,local_output_review_pending=True))
    status('validation_complete',overall='running',remaining='final narrative and figure visual QA')
    print('Final numerical/coverage/identity gates passed',len(checks),flush=True)
if __name__=='__main__':main()
