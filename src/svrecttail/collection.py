"""Frozen v2.1 SV collection, QC and portable numerical handoff.

This layer schedules and records existing algorithms. It never calls the tracker
or moves a recorded vessel to an image maximum. All published paths are relative.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
from dataclasses import asdict
import hashlib
import json
import os
from pathlib import Path
import shutil
import struct
import subprocess
import time
import zipfile

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from .config import load_config
from .io import FrameMaps, load_frame_maps
from .mentor_tracking import build_localization_from_tracking, load_flow_dicom
from .pipeline import _quantify, _sensitivity_records, _validate_map_metadata, _json_compatible
from .qc import save_qc_figure, _display_log, _overlay_geometry

BASELINE = 'ee3012c04be0d02289520a34409d0855a3e5ca58'
ROOT = Path(__file__).resolve().parents[2]
TRACKING = ROOT / 'results/pilot_mentor_tracking_continuity_first_v2_1_15'
QUANT_CONFIG = 'config/run_config.mentor_tracking.continuity_first_v2_1.pilot.json'
TRACK_CONFIG = 'config/tracking_config.continuity_first_v2_1.a020_n4.json'
SPEEDS = (1, 3, 5, 7, 10)
PROTECTED = [QUANT_CONFIG, TRACK_CONFIG, 'src/svrecttail/mentor/tracking_core.py',
             'src/svrecttail/mentor_tracking.py', 'src/svrecttail/quantification.py',
             'src/svrecttail/geometry.py', 'src/svrecttail/background.py',
             'src/svrecttail/pipeline.py', 'src/svrecttail/localization.py'] + [
                 str(p.relative_to(ROOT)).replace('\\', '/') for p in sorted((ROOT/'matlab').glob('*'))
                 if p.suffix in {'.m', '.ini'} and p.name != 'export_sv_collection.m']
TRACK_FIELDS = ['seed_z_upper_px', 'z_upper_px', 'z_candidate_accepted', 'z_short_gap_filled',
                'z_continuity_status', 'z_valid_segment_id', 'new_tracking_class',
                'vessel_presence_prediction', 'x4_centroid_isolated_jump_corrected_px',
                'x1_local_geometry_px', 'local_body_run_width_px', 'x_center_px',
                'z_relock_start', 'z_relock_confirmation_frame']


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(4*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def json_write(path: Path, data) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(json.dumps(_json_compatible(data), ensure_ascii=False, indent=2,
                                    allow_nan=False)+'\n', encoding='utf-8', newline='\n')
    temporary.replace(path)


def portable(value):
    if isinstance(value, dict):
        return {key: portable(item) for key, item in value.items()}
    if isinstance(value, list):
        return [portable(item) for item in value]
    if isinstance(value, str) and (':\\' in value or ':/' in value):
        return value.replace('\\', '/').rsplit('/', 1)[-1]
    return value


def freeze_evidence() -> dict:
    evidence = {}
    for relative in PROTECTED:
        reference = subprocess.check_output(['git', 'show', f'{BASELINE}:{relative}'], cwd=ROOT)
        actual = (ROOT/relative).read_bytes()
        if reference.replace(b'\r\n', b'\n') != actual.replace(b'\r\n', b'\n'):
            raise ValueError(f'Frozen algorithm/config differs from baseline: {relative}')
        evidence[relative] = dict(executed_sha256=hashlib.sha256(actual).hexdigest(),
                                  git_lf_sha256=hashlib.sha256(reference).hexdigest())
    return evidence


def oct_header(path: Path) -> dict:
    with path.open('rb') as f:
        fmt = '<IdIIIIddddIII'
        values = struct.unpack(fmt, f.read(struct.calcsize(fmt)))
    keys = ['bob', 'SPL', 'nX', 'nY_raw', 'Boffset', 'Blength_minus_1', 'Xcenter',
            'Xspan', 'Ycenter', 'Yspan', 'frame_per_pos', 'reserved', 'sizeBck']
    result = dict(zip(keys, values))
    result['n_bscan_positions'] = result['nY_raw']//result['frame_per_pos']
    result['minimum_file_bytes'] = int(result['bob'] + result['nY_raw']*(result['SPL']*result['nX']+2)*2)
    if path.stat().st_size < result['minimum_file_bytes']:
        raise ValueError(f'Truncated acquisition: {path.name}')
    return result


def prepare(raw_root: Path, output: Path, archive: Path, stage: str) -> None:
    frozen = freeze_evidence()
    if stage == 'full':
        gate = output.parent/'formal_sv_d128_v21_bridge15_run001/bridge_validation.json'
        if not gate.exists() or not json.loads(gate.read_text())['stage_b_authorized_by_checks']:
            raise ValueError('Stage A input, numerical and visual checks must pass first')
    expected_run = f'formal_sv_d128_v21_{"bridge15" if stage == "bridge" else "full2500"}_run001'
    if output.exists() and any(output.iterdir()):
        raise FileExistsError('Prepare requires a new output directory; use process for resume.')
    output.mkdir(parents=True, exist_ok=True)
    archive.mkdir(parents=True, exist_ok=True)
    inventory = pd.read_csv(ROOT/'data/raw_inventory.csv')
    sources, rows, jobs = [], [], []
    frame_indices = [0, 249, 499] if stage == 'bridge' else list(range(500))
    for speed in SPEEDS:
        scan = f'flow{speed:02d}'
        header = oct_header(raw_root/f'{speed}.oct')
        if (header['n_bscan_positions'], header['nX'], header['frame_per_pos']) != (500, 500, 3):
            raise ValueError(f'Unexpected OCT acquisition dimensions for {scan}')
        for suffix, kind in [('.oct', 'raw_oct'), ('.oct-Flow_ed.dcm', 'localization_flow_dicom')]:
            source = raw_root/f'{speed}{suffix}'
            expected = inventory.loc[(inventory.scan_id == scan) & (inventory.asset_type == kind)].iloc[0]
            digest = sha256(source)
            if digest != expected.sha256:
                raise ValueError(f'Raw source differs from frozen inventory: {source.name}')
            sources.append(dict(scan_id=scan, asset_type=kind, file=source.name,
                                bytes=source.stat().st_size, sha256=digest))
        target = output/'tracking'/scan
        target.mkdir(parents=True)
        for suffix in ['.csv', '.metadata.json']:
            source = TRACKING/'tracking'/scan/f'{scan}_mentor_tracking{suffix}'
            expected = subprocess.check_output(['git', 'show', f'{BASELINE}:{source.relative_to(ROOT).as_posix()}'], cwd=ROOT)
            if source.read_bytes().replace(b'\r\n', b'\n') != expected:
                raise ValueError(f'Frozen localization changed: {scan}')
            shutil.copy2(source, target/source.name)
        for frame in frame_indices:
            rows.append(dict(scan_id=scan, frame_index_0based=frame, bscan_index_matlab_1based=frame+1,
                source_file=f'arrays/{scan}/frame_{frame:03d}.npz', raw_oct_file=f'{speed}.oct',
                flow_dicom_file=f'{speed}.oct-Flow_ed.dcm', tracking_file=f'tracking/{scan}/{scan}_mentor_tracking.csv',
                diameter_um=128, flow_speed_mm_s=speed, dx_um=12.7, dz_um=6.7,
                vessel_id=None, phantom_id=None, session_id=None, independent_acquisition_id=None,
                temporal_repeat_id=None, temporal_repeat_count=header['frame_per_pos'],
                scan_time_interval_s=None, acquisition_order=None, slow_axis_position_um=None,
                slow_axis_spacing_um=None, reconstruction_version=BASELINE,
                position_label={0:'front',249:'middle',499:'rear'}.get(frame, 'whole_volume'),
                manual_adjustment_source='', manual_adjustment_count=0))
        jobs.append(dict(scan_id=scan, source_file=str((raw_root/f'{speed}.oct').resolve()),
                         interim_dir=str((archive.parent/'interim'/scan).resolve()),
                         frame_indices_0based=frame_indices))
        print(f'PREPARED {scan}: OCT/DICOM hashes match; 500 positions x 3 repeats', flush=True)
    shutil.copy2(ROOT/QUANT_CONFIG, output/'run_config.json')
    shutil.copy2(ROOT/TRACK_CONFIG, output/'tracking_config.json')
    pd.DataFrame(rows).to_csv(output/'manifest.csv', index=False)
    pd.DataFrame(sources).to_csv(output/'input_sha256.csv', index=False)
    state = dict(run_id=expected_run, stage=stage, baseline_commit=BASELINE, algorithm_frozen=frozen,
                 planned_frames=len(rows), status='prepared', blank_profiles='not_available',
                 identity_note='Spatial B-scans are not independent vessels; within-position repeats are not independent acquisitions.')
    json_write(output/'collection_run.json', state)
    json_write(archive/'local_jobs.json', jobs)
    json_write(archive/'local_paths.json', dict(raw_root=str(raw_root.resolve()), output=str(output.resolve()),
                                             archive=str(archive.resolve())))


def validate_maps(maps: FrameMaps, row: pd.Series, flow_shape: tuple) -> None:
    _validate_map_metadata(maps, pd.Series({'bscan_index':row.frame_index_0based}))
    expected_name = row.raw_oct_file
    actual_name = str(maps.metadata.get('source_file', '')).replace('\\', '/').rsplit('/', 1)[-1]
    if actual_name != expected_name:
        raise ValueError('interim_source_file_mismatch')
    if maps.sv_raw.shape != flow_shape:
        raise ValueError('frame_shape_mismatch')
    if maps.stru_amp is None:
        raise ValueError('structural_amplitude_missing')
    metadata = maps.metadata.get('export_metadata', {})
    if metadata.get('variance_denominator') != 'N':
        raise ValueError('variance_denominator_mismatch')
    reconstruction = metadata.get('reconstruction', {})
    if reconstruction.get('IMcropRg') != list(range(50, 401)):
        raise ValueError('reconstruction_crop_mismatch')
    if reconstruction.get('IMG_size') != [flow_shape[0], flow_shape[1], 3]:
        raise ValueError('reconstruction_repeat_or_axis_mismatch')


def alignment_record(maps: FrameMaps, flow: np.ndarray, scan: str, frame: int) -> dict:
    # Diagnostic comparisons only; never transform the measurement or geometry.
    variants = dict(identity=maps.omag_raw, flip_z=maps.omag_raw[::-1],
                    flip_x=maps.omag_raw[:, ::-1], flip_zx=maps.omag_raw[::-1, ::-1])
    scores = {key:float(spearmanr(value.ravel(), flow.ravel()).statistic) for key,value in variants.items()}
    return dict(scan_id=scan, frame_index_0based=frame, shape_z=flow.shape[0], shape_x=flow.shape[1],
                crop_first_matlab=50, crop_last_matlab=400, expected_transform='identity',
                best_flip_diagnostic=max(scores, key=scores.get), **{f'spearman_{k}':v for k,v in scores.items()})


def qc_selected() -> dict[str, set[int]]:
    selected = {f'flow{s:02d}':set([0,49,99,149,199,249,299,349,399,449,499]) for s in SPEEDS}
    selected['flow03'].update([420,423,432,434,436,447,469,470,471])
    gaps = pd.read_csv(TRACKING/'missing_segments.csv')
    relocks = pd.read_csv(TRACKING/'relock_events.csv')
    for row in gaps.itertuples():
        selected[row.scan_id].update([max(0,row.start_frame-1), row.start_frame, row.end_frame, min(499,row.end_frame+1)])
    for row in relocks.itertuples():
        selected[row.scan_id].update([max(0,row.start_frame-1), row.start_frame, row.confirmation_frame,
                                     min(499,row.confirmation_frame+1)])
    return selected


def mapping_figure(path: Path, maps: FrameMaps, flow: np.ndarray, tracking: pd.Series,
                   result=None) -> None:
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1,4,figsize=(16,6),layout='constrained')
    for ax,img,title in zip(axes,[flow,maps.omag_raw,maps.sv_raw,maps.stru_amp],
                            ['Localization Flow DICOM','Rebuilt OMAG raw','Rebuilt SV raw','Rebuilt OCT amplitude']):
        shown,lo,hi=_display_log(img)
        ax.imshow(shown,cmap='gray',vmin=lo,vmax=hi,aspect=6.7/12.7)
        ax.set_title(title)
        x=float(tracking.x_center_px)
        ax.set_xlim(x-35,x+35);ax.set_ylim(310,160)
        if result is not None:
            _overlay_geometry(ax,result)
        else:
            ax.text(.5,.04,'NA geometry: excluded',transform=ax.transAxes,color='red',ha='center')
        ax.set_xlabel('x (zero-based)');ax.set_ylabel('z (zero-based)')
    fig.suptitle(f'{tracking.scan_id}, frame {tracking.frame_index}: shared coordinates; display transforms only')
    path.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(path,dpi=125);plt.close(fig)


def quantitative_checks(result, sv_raw: np.ndarray) -> dict:
    b = result.background
    mask = (b.left_valid_count > 0) & (b.right_valid_count > 0)
    combined = (b.left*b.left_valid_count + b.right*b.right_valid_count)/(b.left_valid_count+b.right_valid_count)
    assert np.allclose(b.combined[mask],combined[mask],rtol=2e-12,atol=1e-12)
    corrected = sv_raw-b.combined[:,None]
    assert np.array_equal(result.corrected_sv,corrected,equal_nan=True)
    width=result.vessel_x_weights.sum()*result.geometry.dx_um
    mask=np.isfinite(result.tail_linear_density)&np.isfinite(result.tail_contrast_profile)
    assert np.allclose(result.tail_linear_density[mask],width*result.tail_contrast_profile[mask],rtol=2e-12,atol=1e-12)
    support=np.multiply.outer(result.tail_z_weights,result.vessel_x_weights)
    direct=float(np.sum(corrected[support>0]*support[support>0])*result.geometry.dx_um*result.geometry.dz_um)
    if np.isfinite(result.q_tail):
        assert np.isclose(direct,result.q_tail,rtol=2e-12,atol=1e-12)
    if result.valid:
        assert np.isclose(result.ratio_tail_to_vessel,result.q_tail/result.q_vessel,rtol=2e-12,atol=1e-12)
    return dict(background_combination_check=True,negative_residual_retention_check=True,
                negative_corrected_pixels=int((corrected<0).sum()),p_width_t_check=True,
                tail_2d_profile_check=True,ratio_identity_check=bool(result.valid),
                raw_sv_min=float(np.min(sv_raw)),raw_sv_max=float(np.max(sv_raw)))


def archive_frame(path: Path, maps: FrameMaps, flow: np.ndarray, record: dict,
                  result, frozen_tracking: dict) -> str:
    content = dict(sv_raw=maps.sv_raw,stru_amp=maps.stru_amp,omag_raw=maps.omag_raw,
                   localization_flow_dicom=flow.astype(np.uint16),
                   metadata_json=np.array(json.dumps(_json_compatible(dict(
                       frame=record,tracking=frozen_tracking,reconstruction=portable(maps.metadata),
                       coordinate_base=0,array_axis_order='z,x',crop_origin_zx=[0,0],
                       full_bscan=True,negative_background_residuals='retain')),ensure_ascii=False,allow_nan=False)))
    if result is not None:
        content.update(source_weights=result.vessel_ellipse_weights,
                       tail_x_weights=result.vessel_x_weights,tail_z_weights=result.tail_z_weights,
                       background_left_columns=result.background.left_columns,
                       background_right_columns=result.background.right_columns,
                       B_left=result.background.left,B_right=result.background.right,B=result.background.combined)
    path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_suffix('.tmp.npz')
    np.savez_compressed(temporary,**content)
    with np.load(temporary,allow_pickle=False) as loaded:
        for name in ['sv_raw','stru_amp','omag_raw']:
            if not np.array_equal(loaded[name],content[name],equal_nan=True):
                raise AssertionError('Numerical archive round-trip failed')
    temporary.replace(path)
    return sha256(path)


@contextmanager
def scan_lock(archive: Path, scan: str):
    """OS-owned lock releases after crashes; two runners cannot write one frame."""
    path = archive/'locks'/f'{scan}.lock'
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('a+b') as handle:
        if path.stat().st_size == 0:
            handle.write(b'0'); handle.flush()
        handle.seek(0)
        if os.name == 'nt':
            import msvcrt
            while True:
                try:
                    msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
                    break
                except OSError:
                    time.sleep(1)
            try:
                yield
            finally:
                handle.seek(0)
                msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
        else:
            import fcntl
            fcntl.flock(handle, fcntl.LOCK_EX)
            try:
                yield
            finally:
                fcntl.flock(handle, fcntl.LOCK_UN)


def process_scan(output: Path, archive: Path, scan: str, remove_interim: bool = False,
                 write_summary: bool = True) -> None:
    with scan_lock(archive, scan):
        _process_scan(output, archive, scan, remove_interim, write_summary)


def _process_scan(output: Path, archive: Path, scan: str, remove_interim: bool,
                  write_summary: bool) -> None:
    freeze_evidence()
    local=json.loads((archive/'local_paths.json').read_text(encoding='utf-8'))
    run=json.loads((output/'collection_run.json').read_text(encoding='utf-8'))
    config=load_config(output/'run_config.json')
    manifest=pd.read_csv(output/'manifest.csv')
    planned=manifest.loc[manifest.scan_id==scan]
    tracking=pd.read_csv(output/f'tracking/{scan}/{scan}_mentor_tracking.csv').set_index('frame_index',drop=False)
    source=Path(local['raw_root'])/str(planned.iloc[0].flow_dicom_file)
    flow=load_flow_dicom(source)
    checkpoints=archive/'checkpoints'/scan
    checkpoints.mkdir(parents=True,exist_ok=True)
    selected=qc_selected()[scan]
    code_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    for row in planned.itertuples(index=False):
        row=pd.Series(row._asdict());frame=int(row.frame_index_0based);key=f'frame_{frame:03d}'
        checkpoint=checkpoints/f'{key}.json';array=archive/'arrays'/scan/f'{key}.npz'
        profile=checkpoints/f'{key}.profile.csv'
        if checkpoint.exists():
            prior=json.loads(checkpoint.read_text(encoding='utf-8'))
            if (prior.get('status')=='complete' and array.exists() and sha256(array)==prior.get('array_sha256')
                    and prior.get('config_sha256')==sha256(output/'run_config.json')
                    and prior.get('tracking_sha256')==sha256(output/f'tracking/{scan}/{scan}_mentor_tracking.csv')):
                continue
            if prior.get('status') != 'input_missing':
                raise ValueError(f'Resume checkpoint conflicts with inputs: {scan}/{frame}')
        tr=tracking.loc[frame];maps=None;result=None;loc=None;sensitivity=[];alignment={}
        record={**row.to_dict(),**{k:tr[k] for k in TRACK_FIELDS},'bscan_index':frame,
                'run_id':run['run_id'],'code_commit':code_commit,'raw_input_available':False,
                'input_mapping_valid':False,'geometry_qc_valid':False,'source_qc_valid':False,
                'background_qc_valid':False,'window_qc_valid':False,'valid':False,
                'q_vessel':np.nan,'q_tail':np.nan,'ratio_tail_to_vessel':np.nan,
                'source_area_um2':np.nan,'source_mean':np.nan,
                'detection_status':'not_evaluated_no_matched_blank','detected':None,
                'detectable_length_um':np.nan,'detection_evaluable':False,
                'localization_source':('missing' if not np.isfinite(tr.z_upper_px) else
                                       'direct_candidate_supported' if tr.z_candidate_accepted else 'short_gap_assisted')}
        source_mat=archive.parent/'interim'/scan/f'{key}.mat'
        try:
            maps=load_frame_maps(source_mat);record['raw_input_available']=True
            validate_maps(maps,row,flow.shape[1:]);record['input_mapping_valid']=True
            if not np.isfinite(tr.z_upper_px):
                record['invalid_reason']='no_reliable_z_geometry'
            else:
                loc=build_localization_from_tracking(tr,diameter_um=128,dx_um=12.7,dz_um=6.7)
                record['geometry_qc_valid']=bool(loc.source_qc_valid)
                record['geometry']=asdict(loc.geometry)
                record['localization_invalid_reason']=loc.mentor_tracking.invalid_reason
                result=_quantify(maps.sv_raw,loc,config)
                record.update(result.summary())
                record.update(quantitative_checks(result,maps.sv_raw))
                record.update(background_qc_valid=result.background_complete,
                              window_qc_valid=result.tail_window_complete and result.source_window_complete)
                record.update({f'diagnostic_{k}':record[k] for k in ['q_vessel','q_tail','ratio_tail_to_vessel']})
                if not result.valid:
                    for k in ['q_vessel','q_tail','ratio_tail_to_vessel']:record[k]=np.nan
                sensitivity=_sensitivity_records(scan,f'{scan}_Bscan{frame:03d}',maps.sv_raw,loc,config,None)
                for r in sensitivity:
                    r.update(frame_index_0based=frame,bscan_index_matlab_1based=frame+1,
                             primary_valid=result.valid,localization_source=record['localization_source'])
                profiles=pd.DataFrame(result.profile_records(scan))
                profiles.insert(1,'frame_index_0based',frame)
                profiles.insert(2,'bscan_index_matlab_1based',frame+1)
                profiles['frame_valid']=result.valid
                profiles['localization_source']=record['localization_source']
                profiles.to_csv(profile,index=False)
            if frame in selected or run['stage']=='bridge':
                alignment=alignment_record(maps,flow[frame],scan,frame)
                mapping_figure(output/f'qc/{scan}_Bscan{frame:03d}_mapping.png',maps,flow[frame],tr,result)
                if result is not None:
                    save_qc_figure(output/f'qc/{scan}_Bscan{frame:03d}_quantification.png',scan_id=scan,
                                   maps=maps,localization=loc,result=result)
        except FileNotFoundError:
            record['invalid_reason']='raw_interim_file_missing'
        except ValueError as error:
            record['invalid_reason']=('input_mapping_error:' if not record['input_mapping_valid'] else 'geometry_unavailable:')+str(error)
        if maps is None:
            # Record the missing planned position, but do not claim a complete archive.
            json_write(checkpoint,dict(status='input_missing',record=record,sensitivity=[],alignment={}))
            continue
        if not profile.exists():
            pd.DataFrame(dict(scan_id=scan,frame_index_0based=frame,bscan_index_matlab_1based=frame+1,
                z_index_0based=np.arange(maps.sv_raw.shape[0]),z_um=np.arange(maps.sv_raw.shape[0])*6.7,
                r_um=np.nan,V=np.nan,B_left=np.nan,B_right=np.nan,B=np.nan,T=np.nan,P=np.nan,
                tail_z_fraction=np.nan,validity='unavailable_geometry',frame_valid=False,
                localization_source=record['localization_source'])).to_csv(profile,index=False)
        array_hash=archive_frame(array,maps,flow[frame],record,result,tr.to_dict())
        json_write(checkpoint,dict(status='complete',record=record,sensitivity=sensitivity,alignment=alignment,
                                  array_sha256=array_hash,config_sha256=sha256(output/'run_config.json'),
                                  tracking_sha256=sha256(output/f'tracking/{scan}/{scan}_mentor_tracking.csv')))
        if remove_interim and source_mat.resolve().is_relative_to((archive.parent/'interim').resolve()):
            source_mat.unlink()
        if frame%25==0 or frame==int(planned.iloc[-1].frame_index_0based):
            print(f'PROCESSED {scan} {frame}: valid={record["valid"]} {record.get("invalid_reason")}',flush=True)
    if write_summary:
        summarize(output,archive)


def summarize(output: Path, archive: Path) -> None:
    manifest=pd.read_csv(output/'manifest.csv');records=[];sensitivity=[];alignment=[];arrays=[]
    for row in manifest.itertuples(index=False):
        checkpoint=archive/f'checkpoints/{row.scan_id}/frame_{row.frame_index_0based:03d}.json'
        if not checkpoint.exists():continue
        c=json.loads(checkpoint.read_text(encoding='utf-8'));records.append(c['record']);sensitivity.extend(c['sensitivity'])
        if c['alignment']:alignment.append(c['alignment'])
        if c.get('array_sha256'):
            path=archive/row.source_file
            arrays.append(dict(scan_id=row.scan_id,frame_index_0based=row.frame_index_0based,
                               file=row.source_file,bytes=path.stat().st_size,sha256=c['array_sha256']))
    if not records:return
    table=pd.DataFrame(records).sort_values(['scan_id','frame_index_0based'])
    geometries=table.pop('geometry') if 'geometry' in table else pd.Series([None]*len(table),index=table.index)
    geometry_table=pd.DataFrame([r if isinstance(r, dict) else {} for r in geometries])
    pd.concat([table[['scan_id','frame_index_0based','bscan_index_matlab_1based',*TRACK_FIELDS]].reset_index(drop=True),
               geometry_table],axis=1).to_csv(output/'localization.csv',index=False)
    table['assessability'] = table.vessel_presence_prediction
    table.to_csv(output/'frame_results.csv',index=False)
    table.loc[~table.valid.astype(bool)].to_csv(output/'exclusions.csv',index=False)
    pd.DataFrame(sensitivity).to_csv(output/'sensitivity_results.csv',index=False)
    pd.DataFrame(alignment).to_csv(output/'coordinate_mapping_checks.csv',index=False)
    pd.DataFrame(arrays).to_csv(output/'arrays_sha256.csv',index=False)
    summaries=[]
    for scan,group in table.groupby('scan_id',sort=True):
        valid=group.valid.astype(bool)
        summary=dict(scan_id=scan,planned_frames=int((manifest.scan_id==scan).sum()),processed_frames=len(group),
            input_available=int(group.raw_input_available.sum()),z_present=int(group.z_upper_px.notna().sum()),
            geometry_qc_valid=int(group.geometry_qc_valid.sum()),quantification_valid=int(valid.sum()),
            quantification_excluded=int((~valid).sum()),
            direct_supported_total=int(group.localization_source.eq('direct_candidate_supported').sum()),
            assisted_total=int(group.localization_source.eq('short_gap_assisted').sum()),
            direct_supported_included=int((valid&group.localization_source.eq('direct_candidate_supported')).sum()),
            assisted_included=int((valid&group.localization_source.eq('short_gap_assisted')).sum()))
        for col in ['q_vessel','q_tail','ratio_tail_to_vessel','source_area_um2','source_mean']:
            values=pd.to_numeric(group.loc[valid,col],errors='coerce')
            summary[col+'_median']=values.median();summary[col+'_q25']=values.quantile(.25);summary[col+'_q75']=values.quantile(.75)
        summaries.append(summary)
    pd.DataFrame(summaries).to_csv(output/'scan_summary.csv',index=False)
    reason_counts=table.loc[~table.valid.astype(bool)].groupby(['scan_id','invalid_reason'],dropna=False).size().reset_index(name='frames')
    reason_counts.to_csv(output/'exclusion_counts.csv',index=False)
    profile_dir=output/'profiles';profile_dir.mkdir(exist_ok=True)
    for scan,group in manifest.groupby('scan_id',sort=True):
        for chunk_start in range(0,500,100):
            frames=group.loc[group.frame_index_0based.between(chunk_start,chunk_start+99)]
            parts=[archive/f'checkpoints/{scan}/frame_{f:03d}.profile.csv' for f in frames.frame_index_0based]
            parts=[p for p in parts if p.exists()]
            if parts:
                pd.concat([pd.read_csv(p) for p in parts],ignore_index=True).to_csv(
                    profile_dir/f'{scan}_frames_{chunk_start:03d}_{chunk_start+99:03d}.csv.gz',index=False,
                    compression={'method':'gzip','compresslevel':6,'mtime':0})
    json_write(output/'execution_status.json',dict(planned=len(manifest),processed=len(table),
        inputs_archived=len(arrays),quantification_valid=int(table.valid.sum()),
        status='numerical_extraction_complete' if len(table)==len(manifest)==len(arrays) else 'incomplete',
        coordinate_mapping_review='pending_visual_review',detection='not_evaluated_no_matched_blank'))


def package(output: Path, archive: Path, destination: Path) -> None:
    status=json.loads((output/'execution_status.json').read_text(encoding='utf-8'))
    if status['status']!='numerical_extraction_complete':raise ValueError('Incomplete extraction cannot be packaged as complete')
    manifest=pd.read_csv(output/'manifest.csv');index=pd.read_csv(output/'arrays_sha256.csv').set_index('file')
    destination.mkdir(parents=True,exist_ok=True);products=[]
    run=json.loads((output/'collection_run.json').read_text(encoding='utf-8'))
    groups=[('bridge15',manifest)] if run['stage']=='bridge' else [
        (f'{scan}_{start:03d}_{start+99:03d}',group.loc[group.frame_index_0based.between(start,start+99)])
        for scan,group in manifest.groupby('scan_id') for start in range(0,500,100)]
    for name,group in groups:
        path=destination/f'{run["run_id"]}_{name}.zip'
        if path.exists():raise FileExistsError(path)
        with zipfile.ZipFile(path,'w',compression=zipfile.ZIP_STORED,allowZip64=True) as z:
            for file in sorted(output.rglob('*')):
                relative=file.relative_to(output).as_posix()
                if not file.is_file() or file.suffix=='.zip':
                    continue
                if run['stage'] != 'bridge':
                    if relative.startswith('profiles/') and name not in file.name.replace('frames_', ''):
                        continue
                    if relative.startswith('qc/'):
                        if not any(file.name.startswith(f'{r.scan_id}_Bscan{r.frame_index_0based:03d}_')
                                   for r in group.itertuples(index=False)):
                            continue
                z.write(file,relative)
            for row in group.itertuples(index=False):
                array=archive/row.source_file
                if sha256(array)!=index.loc[row.source_file,'sha256']:raise ValueError('Array changed before packaging')
                z.write(array,row.source_file)
        with zipfile.ZipFile(path) as z:
            assert z.testzip() is None
        products.append(dict(file=path.name,bytes=path.stat().st_size,sha256=sha256(path),frames=len(group)))
        print(f'PACKAGED {path.name}',flush=True)
    pd.DataFrame(products).to_csv(output/'download_packages.csv',index=False)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=['prepare','process','summarize','package'])
    parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--archive',required=True,type=Path)
    parser.add_argument('--raw-root',type=Path)
    parser.add_argument('--stage',choices=['bridge','full'],default='bridge')
    parser.add_argument('--scan')
    parser.add_argument('--remove-interim',action='store_true')
    parser.add_argument('--no-summary',action='store_true',
                        help='Write only this scan checkpoints/QC; summarize after independent scans finish')
    parser.add_argument('--packages',type=Path)
    args=parser.parse_args(argv)
    if args.command=='prepare':prepare(args.raw_root,args.output,args.archive,args.stage)
    elif args.command=='process':process_scan(args.output,args.archive,args.scan,args.remove_interim,
                                            write_summary=not args.no_summary)
    elif args.command=='summarize':summarize(args.output,args.archive)
    else:package(args.output,args.archive,args.packages)
    return 0
