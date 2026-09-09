"""Audit available raw arrays at invalid geometry positions, without constructing ROIs."""
import io,json,zipfile,hashlib
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.io import loadmat
from derive_qc import OUT,ROOT,DEP,ID,read,csv,js,audit_raw,digest

def main():
    ids=read(DEP/'array_identity_audit.csv.gz').fillna('');fw=read(ROOT/'analysis/sv_source_tail_framewise_coupling_v1/framewise_source_tail_metrics.csv')
    pkg=read(ROOT/'results/formal_sv_d128_v21_full2500_run001/download_packages.csv')
    hashes=read(ROOT/'results/formal_sv_d128_v21_full2500_run001/arrays_sha256.csv').set_index(['scan_id','frame_index_0based'])
    extra=[];ident=[];missing=[];metadata=[]
    for scan,g in ids.groupby('scan_id'):
        m=fw[fw.scan_id.eq(scan)].iloc[0];meta={k:m[k] for k in ID};cache={};root=Path(g.iloc[0].input_path).parent
        root=root if root.is_absolute() else ROOT/root
        for i in sorted(set(range(500))-set(g.frame_index)):
            base=meta|dict(frame_index=i,valid_geometry=False)
            if m.diameter_um==128:
                flow=int(m.flow_mm_s);row=pkg[pkg.file.str.contains(f'flow{flow:02}_{i//100*100:03}_{i//100*100+99:03}')].iloc[0]
                p=root/row.file;member=f'arrays/flow{flow:02}/frame_{i:03}.npz'
                if str(p) not in cache:cache[str(p)]=zipfile.ZipFile(p)
                blob=cache[str(p)].read(member);h=hashlib.sha256(blob).hexdigest();assert h==hashes.loc[(f'flow{flow:02}',i),'sha256']
                with np.load(io.BytesIO(blob),allow_pickle=False) as ar:a=np.asarray(ar['sv_raw'],float)
                status='matched frozen release NPZ SHA256'
            else:
                p=root/f'frame_{i:03}.mat';member=''
                if not p.exists():missing.append(base|dict(expected_path=str(p.resolve())));continue
                blob=p.read_bytes();h=hashlib.sha256(blob).hexdigest();a=np.asarray(loadmat(io.BytesIO(blob),variable_names=['sv_raw'])['sv_raw'],float)
                status='new read-only SHA receipt; invalid frame not present in prior coupling array identity table'
            extra.append(audit_raw(a,base));ident.append(base|dict(input_path=str(p.resolve()),member=member,sha256=h,identity_status=status))
        for z in cache.values():z.close()
    csv('additional_invalid_geometry_raw_audit.csv.gz',pd.DataFrame(extra))
    csv('additional_invalid_array_identities.csv',pd.DataFrame(ident));csv('missing_invalid_geometry_arrays.csv',pd.DataFrame(missing,columns=ID+['frame_index','expected_path']))
    js('additional_raw_validation.json',dict(available_invalid_arrays=len(extra),missing_invalid_arrays=len(missing),
        all_finite=all(x['finite_fraction']==1 for x in extra),negative_count=sum(int(x['negative_count']) for x in extra),
        no_geometry_or_roi_constructed=True))
    print('Additional invalid geometry arrays',len(extra),'missing',len(missing),flush=True)

if __name__=='__main__':main()
