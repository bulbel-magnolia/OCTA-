from extract_pixels import *
def main():
    fw=frame_table();fixed=read(OUT/'localization_fixed_frames.csv');ids=read(INPUT/'array_identity_audit.csv.gz').fillna('').set_index(KEY);cache={};records=[]
    for scan,gg in fixed.groupby('scan_id'):
        shots=[];samples={k:[] for k in ['sv_raw','stru_amp']}
        for r in fw[(fw.scan_id==scan)&fw.frame_index.isin(gg.frame_index)].itertuples():
            arrays,h=load_array(ids.loc[(scan,r.frame_index)],cache);shots.append((r,arrays));g=geometry(r)
            x0=max(0,int(np.floor(r.X4-2.8*r.diameter_um/DX)));x1=min(500,int(np.ceil(r.X4+2.8*r.diameter_um/DX)))
            z0=max(0,int(np.floor(r.z_top-120/DZ)));z1=min(351,int(np.ceil(g.z_bottom_edge_px+550/DZ)))
            for name in samples:
                if name in arrays:samples[name].append(arrays[name][z0:z1,x0:x1].ravel())
        limits={k:float(np.quantile(np.concatenate(v),.995)) if v else 1 for k,v in samples.items()}
        draw_qc(scan,shots,limits,folder='localization_review_display')
        records.extend(dict(scan_id=scan,channel=k,linear_display_lower=0,linear_display_upper=v,rule='99.5th percentile of pooled displayed crops in the same fixed three frames; common range within volume; display only',figure=f'localization_review_display/{scan}.png') for k,v in limits.items())
    for z in cache.values():z.close()
    csv('localization_display_settings.csv',records)
    (OUT/'localization_display_note.md').write_text('Initial full-volume maximum scaling rendered ROI pixels visually dark because maxima occurred elsewhere. All initial figures are retained. The same 45 fixed IDs are additionally shown with a single linear range per scan/channel, 0 to the 99.5th percentile pooled across the three displayed fixed-frame crops. This uniform display-only rule was adopted before detailed contour assessment, selects no frames and changes no analytical values.\n',encoding='utf-8')
    print('45 fixed frames rendered with common per-volume linear display range',flush=True)
if __name__=='__main__':main()
