from extract_pixels import *
def main():
    fw=frame_table();fixed=read(OUT/'localization_fixed_frames.csv');ids=read(INPUT/'array_identity_audit.csv.gz').fillna('').set_index(KEY);cache={};rows=[]
    for r in fw[(fw.scan_id=='D285_F03_V01')&fw.frame_index.isin(fixed[fixed.scan_id=='D285_F03_V01'].frame_index)].itertuples():
        arrays,h=load_array(ids.loc[(r.scan_id,r.frame_index)],cache)
        for side,dx in [('real',0),('left',-1.5*r.diameter_um/DX),('right',1.5*r.diameter_um/DX)]:
            for region,m in masks(geometry(r,dx)).items():
                v,w=pixels(arrays['sv_raw'],m);med=wmedian(v,w);ref=independent_wmedian(v,w);mad=wmedian(abs(v-med),w);rmd=independent_wmedian(abs(v-ref),w)
                if med!=ref or mad!=rmd:
                    row=dict(scan_id=r.scan_id,frame_index=r.frame_index,side=side,region=region,median=med,reference_median=ref,mad=mad,reference_mad=rmd)
                    order=np.argsort(v,kind='stable');c=np.cumsum(w[order]);k=np.searchsorted(c,c[-1]/2);row.update(half_weight=c[-1]/2,cumulative=c[k],at_k=v[order[k]],previous=v[order[k-1]],next=v[order[min(k+1,len(v)-1)]])
                    rows.append(row);print(row,flush=True)
    csv('weighted_discrepancy_diagnosis.csv',rows);status('weighted_statistic_diagnosis',overall='running',pixel_validation='failed_pending_fix')
if __name__=='__main__':main()
