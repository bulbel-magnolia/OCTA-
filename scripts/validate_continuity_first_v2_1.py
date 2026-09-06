"""Reproduce continuity-first v2.1 localization audits from external Flow DICOMs.

Only tables, hashes and rendered figures are published. All indices are zero-based.
SV intensity arrays and tail integrals are outside this localization validation.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import platform
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from svrecttail.mentor import tracking_core as core
from svrecttail.mentor_tracking import (
    write_mentor_tracking_bundle, load_flow_dicom, build_localization_from_tracking, build_mentor_tracking_tables,
)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False)+'\n', encoding='utf-8')


def geometry_valid(row):
    if not np.isfinite(row.z_upper_px):
        return False
    try:
        return bool(build_localization_from_tracking(row, diameter_um=128, dx_um=12.7, dz_um=6.7).source_qc_valid)
    except ValueError:
        return False


def max_adjacent(values):
    delta = np.abs(np.diff(np.asarray(values, float)))
    return float(np.max(delta[np.isfinite(delta)])) if np.isfinite(delta).any() else None


def frame_panel(ax, image, row, old, title):
    x = float(row.x_center_px)
    positive = image[image > 0]
    ax.imshow(image, cmap='gray', vmin=np.percentile(positive, 1), vmax=np.percentile(positive, 99), aspect='auto')
    ax.set_xlim(x-40, x+40); ax.set_ylim(280, 165)
    ax.axhline(old.z_upper_px, color='#e9ae32', ls='--', lw=1.5, label='v2')
    if np.isfinite(row.z_upper_px):
        ax.axhline(row.z_upper_px, color='#00e5ff', lw=1.5, label='v2.1')
        ax.axhline(row.z_upper_px+19, color='#00e5ff', ls=':', lw=1)
    if np.isfinite(row.seed_z_upper_px):
        ax.scatter(x, row.seed_z_upper_px, marker='o' if row.z_candidate_accepted else 'x', color='#ff6ad5', s=32, zorder=5)
    z = f'{row.z_upper_px:.0f}' if np.isfinite(row.z_upper_px) else 'NA'
    ax.set_title(f'{title} | v2 {old.z_upper_px:.0f} / v2.1 {z}\n{row.z_continuity_status}', fontsize=9)
    ax.set_xlabel('x (px)'); ax.set_ylabel('z (px)')


def backward_compatibility(raw_root, out):
    volume = load_flow_dicom(raw_root / '3.oct-Flow_ed.dcm')
    checks = []
    for name in ['pilot_mentor_tracking_15', 'pilot_mentor_tracking_upper_edge_v2_15']:
        path = ROOT / 'results' / name / 'tracking/flow03/flow03_mentor_tracking.csv'
        metadata = json.loads(path.with_suffix('.metadata.json').read_text(encoding='utf-8-sig'))
        current, _, _ = build_mentor_tracking_tables(volume, scan_id='flow03', diameter_um=128,
                                                     tracking_config=metadata['effective_tracking_config'])
        old = pd.read_csv(path)
        columns = ['x_center_px', 'seed_z_upper_px', 'z_upper_px',
                   'x4_centroid_isolated_jump_corrected_px', 'assessability_score']
        differences = {col:float(np.nanmax(abs(current[col]-old[col]))) for col in columns}
        exact = {col:bool((current[col]==old[col]).all()) for col in
                 ['qc_valid', 'new_tracking_class', 'vessel_presence_prediction']}
        assert max(differences.values()) < 1e-12 and all(exact.values())
        checks.append(dict(baseline=name, scan_id='flow03', frames=500,
                           max_abs_difference=differences, categorical_exact=exact, status='passed'))
    write_json(out/'backward_compatibility.json', checks)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    out = args.output.resolve()
    if out.exists() and any(out.iterdir()):
        raise FileExistsError('Use a new output directory; historical results are immutable.')
    out.mkdir(parents=True, exist_ok=True)
    cfg_path = ROOT/'config/tracking_config.continuity_first_v2_1.a020_n4.json'
    cfg = json.loads(cfg_path.read_text(encoding='utf-8'))
    baseline = ROOT/'results/pilot_mentor_tracking_upper_edge_v2_15'
    summaries, gaps, relocks, selected, comparison, keys = [], [], [], [], [], []
    fig, axes = plt.subplots(5, 1, figsize=(14, 14), layout='constrained')
    sample_fig, sample_axes = plt.subplots(5, 3, figsize=(15, 17), layout='constrained')
    exceptions = []
    source_hashes = []
    for speed in [3, 1, 5, 7, 10]:
        scan = f'flow{speed:02d}'
        source = args.raw_root/f'{speed}.oct-Flow_ed.dcm'
        bundle = write_mentor_tracking_bundle(source, scan_id=scan, diameter_um=128,
            output_dir=out/'tracking'/scan, tracking_config=cfg)
        metadata = json.loads(bundle.metadata_path.read_text(encoding='utf-8'))
        metadata['flow_dicom_path'] = source.name
        metadata['mentor_tracking_core_path'] = 'src/svrecttail/mentor/tracking_core.py'
        metadata['z_core_feedback_method'] = 'global_viterbi_bootstrap_then_upper_slab_causal_penalized_restart_v2_1'
        write_json(bundle.metadata_path, metadata)
        source_hashes.append({'file':source.name, 'bytes':source.stat().st_size, 'sha256':metadata['flow_dicom_sha256']})
        d = pd.read_csv(bundle.primary_table_path)
        old_path = baseline/'tracking'/scan/f'{scan}_mentor_tracking.csv'
        old = pd.read_csv(old_path)
        old_meta = json.loads((old_path.with_suffix('.metadata.json')).read_text(encoding='utf-8-sig'))
        assert old_meta['flow_dicom_sha256'] == metadata['flow_dicom_sha256'], 'Baseline input hash differs'
        valid = d.z_upper_px.notna().to_numpy()
        geom = d.apply(geometry_valid, axis=1).to_numpy(bool)
        n = len(d)
        missing_runs = core._boolean_runs(~valid)
        for start, stop in missing_runs:
            gaps.append(dict(scan_id=scan, start_frame=start, end_frame=stop-1, length_frames=stop-start,
                             terminal=start==0 or stop==n))
        for frame in np.flatnonzero(d.z_relock_start):
            relocks.append(dict(scan_id=scan,start_frame=int(frame),confirmation_frame=int(d.loc[frame,'z_relock_confirmation_frame']),
                                seed_z_upper_px=float(d.loc[frame,'seed_z_upper_px']),z_upper_px=float(d.loc[frame,'z_upper_px'])))
        filled_lengths = [b-a for a,b in core._boolean_runs(d.z_short_gap_filled.to_numpy(bool))]
        core_steps = d.z_core_path_px.diff().abs()
        uninterrupted = ~d.z_core_restart.astype(bool)
        assert core_steps[uninterrupted].max() <= 2
        assert max_adjacent(d.z_upper_px) <= 2
        assert not filled_lengths or max(filled_lengths) <= 5
        assert (d.loc[~valid,'new_tracking_class'] == 'failed').all()
        assert (d.loc[~valid,'vessel_presence_prediction'] == 'not_assessable').all()
        accepted = d.z_candidate_accepted.to_numpy(bool)
        adjacent = accepted[1:] & accepted[:-1]
        assert (np.abs(np.diff(d.seed_z_upper_px))[adjacent] <= 2).all()
        checked = accepted & d.z_local_prediction_px.notna().to_numpy()
        assert ((d.seed_z_upper_px-d.z_local_prediction_px).abs()[checked] <= 3).all()
        summary = dict(scan_id=scan,total_frames=n,z_valid_frames=int(valid.sum()),z_excluded_frames=int((~valid).sum()),
            z_coverage=float(valid.mean()),tracking_qc_valid_frames=int(d.qc_valid.sum()),
            geometry_qc_valid_frames=int(geom.sum()),geometry_qc_excluded_frames=int((~geom).sum()),
            accepted_candidates=int(accepted.sum()),short_gap_filled_frames=int(d.z_short_gap_filled.sum()),
            missing_segments=len(missing_runs),max_missing_length_frames=max([b-a for a,b in missing_runs],default=0),
            relock_count=int(d.z_relock_start.sum()),max_adjacent_z_jump_px=max_adjacent(d.z_upper_px),
            v2_max_adjacent_z_jump_px=max_adjacent(old.z_upper_px),max_continuous_core_jump_px=float(core_steps[uninterrupted].max()))
        summaries.append(summary)
        comp = d[['scan_id','frame_index','seed_z_upper_px','z_upper_px','z_candidate_accepted','z_continuity_status',
                  'z_relock_start','z_relock_confirmation_frame','z_valid_segment_id','z_core_path_px','z_core_restart',
                  'qc_valid','new_tracking_class','vessel_presence_prediction']].copy()
        comp['v2_z_upper_px'] = old.z_upper_px
        comp['delta_z_upper_px'] = d.z_upper_px-old.z_upper_px
        comp['geometry_qc_valid'] = geom
        comparison.append(comp)
        row_index = [1,3,5,7,10].index(speed)
        ax=axes[row_index]
        ax.plot(old.frame_index,old.z_upper_px,color='#bb8a22',label='v2',lw=1)
        ax.scatter(d.frame_index,d.seed_z_upper_px,color='#c5c9cc',s=5,label='raw single candidate')
        ax.plot(d.frame_index,d.z_upper_px,color='#007b9e',lw=1.6,label='v2.1')
        for start,stop in missing_runs:
            ax.axvspan(start-.5,stop-.5,color='#d05b5b',alpha=.16)
        ax.set_title(f'{scan}: z coverage {valid.sum()}/{n} ({valid.mean():.1%}); missing {sum(~valid)}')
        ax.set_ylabel('Upper edge z (px)');ax.invert_yaxis();ax.grid(alpha=.15)
        if row_index==0: ax.legend(ncol=3,fontsize=8)
        volume=load_flow_dicom(source)
        for column,frame in enumerate([0,249,499]):
            record=comp.loc[frame].to_dict()
            record['within_2px'] = bool(np.isfinite(record['delta_z_upper_px']) and abs(record['delta_z_upper_px'])<=2)
            selected.append(record)
            frame_panel(sample_axes[row_index,column],volume[frame],d.loc[frame],old.loc[frame],f'{scan} / {frame}')
            if not record['within_2px']:
                exceptions.append((scan,frame,volume[frame].copy(),d.loc[frame].copy(),old.loc[frame].copy()))
        if speed==3:
            frame_list=[420,423,436,447,469,471]
            key_fig,key_axes=plt.subplots(2,3,figsize=(15,10),layout='constrained')
            for ax,frame in zip(key_axes.flat,frame_list):
                frame_panel(ax,volume[frame],d.loc[frame],old.loc[frame],f'flow03 / {frame}')
                record=d.loc[frame].to_dict();record['v2_z_upper_px']=float(old.loc[frame,'z_upper_px']);keys.append(record)
            key_fig.savefig(out/'flow03_key_frames.png',dpi=150);plt.close(key_fig)
            # Show the same unmodified v2 detector evidence at the recovered shallow core.
            effective=core.merge_tracking_config(cfg)
            xt=core.locate_lateral_track(volume,128,effective)
            profiles=core.extract_profiles(volume,xt['x_center'],128,effective)
            support=core.extract_persistent_row_support(volume,xt['x_center'],128,effective)
            evidence_fig,evidence_axes=plt.subplots(2,3,figsize=(15,9),layout='constrained')
            evidence_rows=[]
            for ax,frame in zip(evidence_axes.flat,frame_list):
                r=d.loc[frame];z=np.arange(volume.shape[1])
                ax.plot(z,profiles.excess[frame]/r.background_sigma,color='#333333',lw=1,label='excess / noise sigma')
                ax.axvline(r.seed_z_upper_px,color='#007b9e',label='raw upper candidate')
                ax.axvline(r.z_core_path_px,color='#d25999',ls=':',label='core hint')
                ax2=ax.twinx();ax2.plot(z,support[frame],color='#b8860b',alpha=.6,lw=1);ax2.set_ylim(0,1);ax2.set_ylabel('row support')
                ax.set_xlim(175,245);ax.set_xlabel('z (px)');ax.set_ylabel('excess / sigma')
                ax.set_title(f'{frame}: top contrast {r.upper_top_contrast_snr:.2f}; persistent {r.upper_persistent_fraction:.2f}',fontsize=10)
                if frame==420:ax.legend(fontsize=7)
                for depth in range(165,281):
                    evidence_rows.append(dict(frame_index=frame,z_px=depth,excess=float(profiles.excess[frame,depth]),
                                              noise_sigma=float(r.background_sigma),row_support=float(support[frame,depth])))
            evidence_fig.savefig(out/'flow03_candidate_evidence.png',dpi=140);plt.close(evidence_fig)
            pd.DataFrame(evidence_rows).to_csv(out/'flow03_candidate_profiles.csv',index=False)
            assert d.loc[420,'z_upper_px'] in [190,191]
            assert not d.loc[425:465,'z_upper_px'].between(225,236).any()
            assert abs(d.loc[469,'z_upper_px']-199)<=2
        del volume
        print(json.dumps(summary),flush=True)
    axes[-1].set_xlabel('B-scan frame index (zero-based)')
    fig.savefig(out/'v2_vs_v2_1_trajectories_5x500.png',dpi=140);plt.close(fig)
    sample_fig.suptitle('Yellow dashed: v2 | cyan: v2.1 upper / operational lower | pink: raw candidate (x = rejected)',fontsize=12)
    sample_fig.savefig(out/'v2_vs_v2_1_sampled_15.png',dpi=140);plt.close(sample_fig)
    if exceptions:
        ef,ea=plt.subplots(1,len(exceptions),figsize=(6*len(exceptions),5),layout='constrained',squeeze=False)
        for ax,(scan,frame,img,row,old) in zip(ea.flat,exceptions):
            frame_panel(ax,img,row,old,f'{scan} / {frame}: exception')
        ef.savefig(out/'sampled_exceptions.png',dpi=150);plt.close(ef)
    pd.DataFrame(summaries).sort_values('scan_id').to_csv(out/'full_volume_tracking_audit.csv',index=False)
    pd.DataFrame(gaps,columns=['scan_id','start_frame','end_frame','length_frames','terminal']).to_csv(out/'missing_segments.csv',index=False)
    pd.DataFrame(relocks,columns=['scan_id','start_frame','confirmation_frame','seed_z_upper_px','z_upper_px']).to_csv(out/'relock_events.csv',index=False)
    pd.concat(comparison).sort_values(['scan_id','frame_index']).to_csv(out/'comparison_v2_vs_v2_1_2500.csv',index=False)
    pd.DataFrame(selected).sort_values(['scan_id','frame_index']).to_csv(out/'sampled_localization_15.csv',index=False)
    pd.DataFrame(keys).to_csv(out/'flow03_key_frames.csv',index=False)
    pd.DataFrame(source_hashes).to_csv(out/'input_sha256.csv',index=False)
    method_files=['src/svrecttail/mentor/tracking_core.py','src/svrecttail/mentor_tracking.py','src/svrecttail/mentor/xroi.py',
                  'config/tracking_config.continuity_first_v2_1.a020_n4.json','scripts/validate_continuity_first_v2_1.py']
    write_json(out/'validation.json',dict(status='complete_with_documented_sample_exceptions',
        scope='five full Flow volumes and 15 fixed B-scan localization/QC checks; SV quantification not rerun',
        method=core.UPPER_EDGE_CONTINUITY,manual_clicks=0,code_sha256={p:digest(ROOT/p) for p in method_files},
        environment=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__),
        all_2500_frames_audited=True,sampled_total=15,sampled_within_2px=sum(r['within_2px'] for r in selected),sampled_exceptions=len(exceptions),
        total_z_valid_frames=sum(r['z_valid_frames'] for r in summaries),total_z_excluded_frames=sum(r['z_excluded_frames'] for r in summaries),
        flow03_wrong_branch_removed=True,final_adjacent_jump_limit_passed=True,short_gap_limit_passed=True,
        raw_candidate_continuity_limits_passed=True,missing_z_invalid_and_not_assessable_passed=True))
    backward_compatibility(args.raw_root, out)
    # Public artifact guard: reject local absolute paths and raw/array files.
    for path in out.rglob('*'):
        if not path.is_file():continue
        assert path.suffix.lower() not in {'.oct','.dcm','.mat','.npy','.npz','.h5'}
        if path.suffix in {'.csv','.json','.md'}:
            text=path.read_text(encoding='utf-8-sig')
            assert 'C:\\Users' not in text and 'C:/Users' not in text and 'C:\\\\Users' not in text
    return 0


if __name__=='__main__':
    raise SystemExit(main())
