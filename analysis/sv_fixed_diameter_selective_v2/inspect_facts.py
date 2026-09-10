from common import *
def main():
    for file,cols,filters in [
        ('fixed_region_across_flow.csv',['diameter_um','mode','region','median','min','max'],{'region':'C16'}),
        ('geometry_version_direction_ledger.csv',['scan_id','mode','metric','old_value','new_value','new_minus_old','direction_status'],{'metric':'C16'}),
        ('local_contrast_across_flow.csv',['diameter_um','region','metric','median','min','max'],{'metric':'z_local'}),
        ('adjacent_difference_across_flow.csv',None,{}),
        ('conditional_A_summary.csv',None,{})]:
        df=read(OUT/file)
        for k,v in filters.items():df=df[df[k]==v]
        if file=='local_contrast_across_flow.csv':df=df[df.region.isin(['source','s36','tail100','tail500','t01','t04','t10','t20'])]
        print('\n'+file+'\n'+(df[cols] if cols else df).to_string(index=False),flush=True)
    ex=read(OUT/'fixed_region_control_cells.csv');print('nonpositive excess',len(ex[ex.excess_cell<=0]),'of',len(ex),'C16',len(ex[ex.in_C16&(ex.excess_cell<=0)]),'negative',ex.excess_cell.min())
    maps=read(OUT/'depth_maps_raw51.csv');print('minimum map',maps.loc[maps.rho.idxmin()].to_dict())
    for filename,col0 in [('detrend_sensitivity_changes.csv','change'),('position_sensitivity_metric_changes.csv','percent_change')]:
        df=read(OUT/filename);print(filename,df.groupby('metric')[col0].agg(['min','max']).to_string())
    st=read(OUT/'conditional_D_structural_correlations.csv');print('structural',st.groupby(['diameter_um','side','mode','kind']).rho.agg(['min','max']).to_string())
    dr=read(OUT/'background_drift.csv');print('drift combined',dr[dr.side.eq('bg')].groupby('region')[['RCV','relative_difference']].agg(['min','max']).to_string())
if __name__=='__main__':main()
