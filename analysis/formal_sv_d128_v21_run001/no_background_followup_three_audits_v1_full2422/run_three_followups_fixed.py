#!/usr/bin/env python3
"""Corrected runner for analyze_three_followups.py.

The primary CSV contains both frame_index and frame_index_0based. The original
exact-window helper renames frame_index_0based to frame_index, so this runner
passes a canonical primary table without the redundant legacy frame_index column.
Scientific definitions and all analysis parameters are unchanged.
"""
from __future__ import annotations
import argparse, importlib.util, json
from pathlib import Path
import numpy as np
import pandas as pd


def load(path: Path):
    spec=importlib.util.spec_from_file_location('three_followups_core',path)
    mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod); return mod


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--root',default='.'); ap.add_argument('--output-dir',default='analysis/formal_sv_d128_v21_run001/no_background_followup_three_audits_v1_full2422'); ap.add_argument('--source-sha',required=True); ap.add_argument('--workflow-source-sha',required=True); a=ap.parse_args()
    root=Path(a.root).resolve(); out=(root/a.output_dir).resolve(); out.mkdir(parents=True,exist_ok=True)
    core_path=out/'analyze_three_followups.py'; C=load(core_path)
    primary_path=root/'analysis/formal_sv_d128_v21_run001/no_background_relative_tail_v1_full2422/framewise_primary.csv'
    primary=pd.read_csv(primary_path)
    if 'frame_index' in primary.columns: primary=primary.drop(columns=['frame_index'])
    if len(primary)!=C.EXPECTED_VALID or primary.groupby('scan_id').size().to_dict()!=C.EXPECTED_COUNTS: raise AssertionError('primary coverage mismatch')
    if primary.duplicated(['scan_id','frame_index_0based']).any(): raise AssertionError('duplicate primary keys')
    for col in ['source_mean_raw','tail_mean_raw','ri_tail']:
        if not np.isfinite(primary[col].to_numpy(float)).all(): raise AssertionError(f'nonfinite {col}')
    coupling,deciles,ctx=C.coupling_audit(primary,out)
    wf,ws,wctx,wval=C.exact_windows(root,primary,out)
    seg,st,sctx=C.segment_audit(primary,out)
    validation={'valid_frames':C.EXPECTED_VALID,'scan_valid_counts':C.EXPECTED_COUNTS,'window_frame_rows':len(wf),'window_lengths_um':list(C.WINDOWS_UM),
                'segments':[[x,y] for x,y in C.SEGMENTS],**wval,'p_values_computed':False,'background_subtracted':False,
                'invalid_frames_interpolated':False,'invalid_frames_zero_filled':False,'bscans_treated_as_independent_replicates':False,
                'implementation_note':'canonical primary input drops redundant legacy frame_index before exact-window merge'}
    (out/'validation.json').write_text(json.dumps(validation,indent=2)+'\n',encoding='utf-8')
    provenance={'analysis':'three descriptive follow-ups for no-background RI_tail','source_sha':a.source_sha,'workflow_source_sha':a.workflow_source_sha,
                'primary_definition':'mean raw SV tail / mean raw SV vessel','audits':['source-tail-ratio coupling','exact 100/200/300/500 um windows','five 100-frame spatial segments'],
                'statistical_scope':'within-volume descriptive only; no flow p-values; B-scans not independent experimental units'}
    (out/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n',encoding='utf-8')
    inputs=[primary_path,root/'results/formal_sv_d128_v21_full2500_run001/download_packages.csv',root/'results/formal_sv_d128_v21_full2500_run001/arrays_sha256.csv',root/'results/formal_sv_d128_v21_full2500_run001/run_config.json',
            root/'analysis/formal_sv_d128_v21_run001/fixed128_width_sensitivity_full2422/audit_fixed128_width.py',root/'src/svrecttail/geometry.py',core_path]
    pd.DataFrame([{'input_path':p.relative_to(root).as_posix(),'sha256':C.sha256(p),'size_bytes':p.stat().st_size} for p in inputs]).to_csv(out/'input_sha256.csv',index=False)
    C.write_readme(out,coupling,ctx,ws,wctx,seg,st,sctx,validation)
    print(json.dumps({'valid_frames':C.EXPECTED_VALID,'coupling':coupling[['scan_id','spearman_source_vs_ri','spearman_tail_vs_ri']].to_dict('records'),
                      'window_flow07':wctx.to_dict('records'),'segment_flow07_ri':sctx[sctx.metric.eq('ri_median')].to_dict('records'),'validation':wval},indent=2))

if __name__=='__main__': main()
