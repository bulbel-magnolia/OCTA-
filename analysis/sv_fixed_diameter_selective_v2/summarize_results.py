"""Explicit five-flow summaries, signed completeness and traceable comparisons."""
from common import *
def five_flow(df,groupcols,value,source,output):
    rows=[]
    for key,g in df.groupby(groupcols,dropna=False,sort=True):
        key=(key,) if not isinstance(key,tuple) else key;meta=dict(zip(groupcols,key));assert len(g)==5 and set(g.flow_mm_s)=={1,3,5,7,10}
        rows.append(meta|dict(input_table=source,value_column=value)|summary(g[value])|{f'flow{int(r.flow_mm_s):02}':getattr(r,value) for r in g.itertuples()})
    csv(output,rows);return pd.DataFrame(rows)
def main():
    fw=frame_table();absolute=[];position_change=[];reasons=[];diagnostic=[]
    for key,g in fw.groupby(ID):
        meta=dict(zip(ID,key));scan=meta['scan_id'];q=signals(g,'q');area=signals(g,'area')
        for region in ['source','tail100','tail500']:
            for metric,vals,unit in [('Q',q[region],'raw_SV_instrument_units_um2'),('area',area[region],'um2')]:absolute.append(meta|dict(region=region,metric=metric,unit=unit)|summary(vals))
        bg=read(OUT/f'framewise_local_background_qc_{scan}.csv.gz');p=read(OUT/f'position_integrals_{scan}.csv.gz')
        for region in REGIONS:
            for metric,statuscol in [('z_local',region+'_z_status'),('signed_asymmetry',region+'_asymmetry_status')]:
                s=bg[statuscol].fillna('invalid_geometry');assert s[~g.valid_geometry.to_numpy()].eq('invalid_geometry').all()
                for reason,n in s.value_counts().items():reasons.append(meta|dict(module='local_background',region=region,metric=metric,reason=reason,n_positions=int(n)))
        for name in OFFSETS:
            for reason,n in p[name+'_status'].value_counts().items():reasons.append(meta|dict(module='position_integral',region=name,metric='validity',reason=reason,n_positions=int(n)))
        sp=read(OUT/f'geometry_support_{scan}.csv.gz');ar=read(OUT/'geometry_area.csv.gz');ar=ar[ar.scan_id==scan]
        theory=ar[['region','analytic_area_um2']].drop_duplicates();assert len(theory)==30
        sp=sp.merge(theory,on='region',validate='many_to_one');sp['relative_analytic_area_error']=(sp.actual_area_um2-sp.analytic_area_um2)/sp.analytic_area_um2
        for (region,side),gg in sp.groupby(['region','side']):diagnostic.append(meta|dict(region=region,side=side)|summary(gg.relative_analytic_area_error))
    csv('absolute_integral_metrics.csv',absolute);csv('diagnostic_status_counts.csv',reasons);csv('geometry_support_area_error.csv',diagnostic)
    for source,groups,value,out in [
        ('volume_metrics.csv',['diameter_um','metric'],'median','volume_metrics_across_flow.csv'),
        ('raw_depth_profiles.csv',['diameter_um','tail_bin','tail_lo_um','tail_hi_um','tail_mid_um','metric'],'median','depth_profiles_across_flow.csv'),
        ('fixed_region_summary.csv',['diameter_um','mode','support','region'],'value','fixed_region_across_flow.csv'),
        ('local_contrast_profiles.csv',['diameter_um','region','metric'],'median','local_contrast_across_flow.csv'),
        ('adjacent_difference.csv',['diameter_um','pair'],'rho','adjacent_difference_across_flow.csv'),
        ('background_real_tail_covariation.csv',['diameter_um','mode','tail_bin','tail_lo_um','tail_hi_um','tail_mid_um'],'rho','background_covariation_across_flow.csv'),
        ('crossed_controls.csv',['diameter_um','pair','side'],'rho','crossed_controls_across_flow.csv')]:five_flow(read(OUT/source),groups,value,source,out)
    joint=[]
    for source,groupcols,cols in [
        ('whole_coupling.csv',['diameter_um','mode','pair'],['rho0','specificity','matched_specificity','peak_lag','zero_lag_rank']),
        ('fixed_region_control_excess.csv',['diameter_um'],['real_matched_rho_value','left_rho_value','right_rho_value','excess_cell_value','real_minus_L_value','real_minus_R_value']),
        ('background_drift.csv',['diameter_um','region','side'],['RCV','signed_difference','absolute_difference','relative_difference']),
        ('position_sensitivity_pairs.csv',['diameter_um','offset','pair'],['rho0','specificity','rho0_change','specificity_change']),
        ('position_sensitivity_fixed_region.csv',['diameter_um','offset'],['value','C16_change'])]:
        df=read(OUT/source)
        for key,g in df.groupby(groupcols,sort=True):
            key=(key,) if not isinstance(key,tuple) else key;meta=dict(zip(groupcols,key));assert len(g)==5
            for c in cols:
                if c=='matched_specificity' and meta.get('mode')=='raw':continue
                joint.append(meta|dict(source_table=source,metric=c)|summary(g[c])|{f'flow{int(r.flow_mm_s):02}':getattr(r,c) for r in g.itertuples()})
    csv('five_flow_diagnostic_summaries.csv',joint)
    maps=read(OUT/'depth_maps_raw51.csv');signs=[]
    for key,g in maps.groupby(ID+['mode']):
        meta=dict(zip(ID+['mode'],key))
        for scope,gg in [('all_120',g),('C16',g[g.in_C16]),('outside_C16',g[~g.in_C16]),('tail_100_to_500',g[g.tail_bin>4])]:signs.append(meta|dict(scope=scope)|summary(gg.rho))
    csv('map_signed_coverage.csv',signs)
    # No new cutoffs for weak association; every coefficient is delivered.
    csv('negative_map_cells.csv',maps[maps.rho<0]);csv('undefined_map_cells.csv',maps[maps.rho.isna()])
    c=read(OUT/'fixed_region_control_cells.csv');csv('nonpositive_control_excess_cells.csv',c[c.excess_cell<=0])
    pos=read(OUT/'position_sensitivity_metrics.csv');pos=pos[pos.support.eq('five_position_common')]
    base=pos[pos.offset.eq('baseline')].set_index(['scan_id','metric'])
    for r in pos[~pos.offset.eq('baseline')].itertuples():
        b=base.loc[(r.scan_id,r.metric),'median'];position_change.append({k:getattr(r,k) for k in ID}|dict(offset=r.offset,metric=r.metric,baseline_common_median=b,offset_common_median=r.median,difference=r.median-b,percent_change=ratio((r.median-b)*100,b).item(),n_common=r.n_defined))
    csv('position_sensitivity_metric_changes.csv',position_change)
    sensitivity=read(OUT/'detrend_sensitivity.csv');baseline=read(OUT/'whole_coupling.csv');baseline=baseline[baseline['mode']=='detrended_51'];basereg=read(OUT/'fixed_region_summary.csv');basereg=basereg[(basereg['mode']=='detrended_51')&basereg.region.eq('C16')]
    sens=[]
    for r in sensitivity.itertuples():
        if r.pair=='C16':b=basereg[basereg.scan_id==r.scan_id].value.iloc[0];sens.append({k:getattr(r,k) for k in ID}|dict(mode=r.mode,pair='C16',metric='rho_median',baseline_51=b,value=r.value,change=r.value-b))
        else:
            b=baseline[(baseline.scan_id==r.scan_id)&baseline.pair.eq(r.pair)].iloc[0]
            for metric in ['rho0','specificity']:sens.append({k:getattr(r,k) for k in ID}|dict(mode=r.mode,pair=r.pair,metric=metric,baseline_51=b[metric],value=getattr(r,metric),change=getattr(r,metric)-b[metric]))
    csv('detrend_sensitivity_changes.csv',sens)
    geom=read(OUT/'geometry_version_comparison.csv');geom['direction_status']=np.where(geom.old_value.isna()|geom.new_value.isna(),'undefined',np.where(geom.old_value*geom.new_value<0,'reversed',np.where(abs(geom.new_value)<abs(geom.old_value),'same_sign_reduced',np.where(abs(geom.new_value)>abs(geom.old_value),'same_sign_increased','unchanged'))))
    csv('geometry_version_direction_ledger.csv',geom)
    # Exact literal units and unavailable physical interpretations are part of coverage.
    js('undefined_physical_interpretations.json',dict(instrument_setting_consistency='not verified across scans; raw recorded units retained',slow_axis_um_per_frame='not recorded; no conversion',independent_lower_wall='not resolved by fixed retained structural panels; geometry prior only',PSF='not recorded',noise_floor='not calibrated',independent_speckle_count='not inferred from Neff',detection_depth='not defined by local robust Z',causal_effect='not identified by coefficient subtraction'))
    status('summaries_complete')
if __name__=='__main__':main()
