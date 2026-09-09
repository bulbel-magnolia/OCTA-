#!/usr/bin/env python3
"""Frozen raw-SV within-volume spatial coupling. No inferential tests."""
from __future__ import annotations
import argparse
import hashlib
import io
import json
import platform
import subprocess
import sys
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy.stats import rankdata

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'src'))
from svrecttail.geometry import VesselGeometry
from svrecttail.raw_sv_quantification import quantify_raw_sv

START_HEAD = '0d1c428d901447f17d0fa600b37dfdfa41e776dd'
STAGE = 'analysis/sv_diameter_stage2_scientific_v1/'
D128 = 'analysis/formal_sv_d128_v21_run001/'
FORMAL128 = 'results/formal_sv_d128_v21_full2500_run001/'
COMMON = [1, 3, 5, 7, 10]
BANDS = ['upper', 'middle', 'lower']
PREDICTORS = [f'{s}_{m}_raw' for s in ['source'] + BANDS for m in ['mean', 'q']]
TAILS = [f'tail_{m}_raw_{d}um' for d in [100, 500] for m in ['mean', 'q']]
VARIABLES = PREDICTORS + TAILS + ['source_area', 'z_top']
KEYS = ['scan_id', 'frame_index_0based']
ID = ['scan_id', 'diameter_um', 'flow_mm_s', 'grid_scope']
MANIFEST = []
CHECKS = {}
TOL = 1e-10


def sha(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def record(path, role, digest=None):
    path = Path(path).resolve()
    try:
        name = path.relative_to(ROOT).as_posix()
    except ValueError:
        name = path.as_posix()
    MANIFEST.append(dict(input_path=name, sha256=digest or sha(path), role=role,
                         size_bytes=path.stat().st_size))


def read(name, role='frozen derived metrics'):
    path = ROOT / name
    record(path, role)
    return pd.read_csv(path, float_precision='round_trip')


def check(name, condition, **details):
    CHECKS[name] = dict(status='passed' if bool(condition) else 'failed', **details)
    if not condition:
        write_json('validation.json', dict(status='failed', checks=CHECKS))
        raise AssertionError(name)


def compare(name, actual, expected):
    a, b = np.asarray(actual, float), np.asarray(expected, float)
    error = np.abs(a-b)
    relative = error / np.maximum(np.abs(b), np.finfo(float).tiny)
    check(name, a.shape == b.shape and np.isfinite(a).all() and
          np.isfinite(b).all() and float(relative.max()) < TOL,
          max_absolute_error=float(error.max()), max_relative_error=float(relative.max()),
          tolerance_relative=TOL, n_values=int(a.size))


def csv(name, frame):
    frame.to_csv(OUT / name, index=False, float_format='%.17g', lineterminator='\n')


def write_json(name, obj):
    (OUT / name).write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False)+'\n', encoding='utf-8')


def pearson(a, b):
    a, b = a-a.mean(), b-b.mean()
    den = np.linalg.norm(a)*np.linalg.norm(b)
    return float(np.clip(np.dot(a, b)/den, -1, 1)) if den > 0 else np.nan


def rho(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    valid = np.isfinite(x) & np.isfinite(y)
    n = int(valid.sum())
    if n < 3:
        return np.nan, n
    return pearson(rankdata(x[valid], method='average'), rankdata(y[valid], method='average')), n


def detrend(frame, window):
    """Integer spatial grid, edge-truncated centered median, ignore missing values."""
    return frame - frame.rolling(window, center=True, min_periods=1).median()


def lag_pair(x, y, lag):
    # Positive lag = tail at higher original frame index than source.
    if lag > 0:
        return x[:-lag], y[lag:]
    if lag < 0:
        return x[-lag:], y[:lag]
    return x, y


def partial(x, y, controls):
    matrix = np.column_stack([x, y, controls])
    matrix = matrix[np.isfinite(matrix).all(axis=1)]
    n = len(matrix)
    if n < 3:
        return np.nan, n, 0, 'insufficient_pairs'
    ranks = np.column_stack([rankdata(v) for v in matrix.T])
    design = np.column_stack([np.ones(n), ranks[:, 2:]])
    rank = int(np.linalg.matrix_rank(design))
    residual = ranks[:, :2] - design @ np.linalg.lstsq(design, ranks[:, :2], rcond=None)[0]
    # Exact/nearly exact rank collinearity has undefined partial correlation.
    if any(np.linalg.norm(residual[:, i]) <= 1e-12 * max(np.linalg.norm(ranks[:, i]), 1) for i in [0, 1]):
        return np.nan, n, rank, 'constant_or_fully_explained_rank_signal'
    return pearson(residual[:, 0], residual[:, 1]), n, rank, 'ok'


def self_checks():
    x = np.array([1., 2, np.nan, 3, 4])
    y = np.array([8., 1, 2, np.nan, 3])
    r, n = rho(*lag_pair(x, y, 1))
    check('lag_uses_original_indices_and_positive_tail_shift', abs(r-1) < 1e-14 and n == 3)
    f = pd.DataFrame({'x': [1., 10, np.nan, 100, 7]})
    med = f.rolling(3, center=True, min_periods=1).median().x.to_numpy()
    check('moving_median_missing_and_edges', np.array_equal(med, [5.5, 5.5, 55, 53.5, 53.5])
          and np.isnan(detrend(f, 3).x.iloc[2]))
    check('difference_does_not_bridge_gap', int(np.isfinite(pd.Series(x).diff()).sum()) == 2)
    r, n = rho([1, 1, 3, np.nan, 5], [2, 2, 6, 9, 10])
    check('average_ties_pairwise_ranking', abs(r-1) < 1e-14 and n == 4)
    rng = np.random.default_rng(20260909)
    a, b, c = rng.normal(size=(3, 70))
    p, n, rank, status = partial(a, b, c[:, None])
    ra, rb, rc = [rankdata(z) for z in [a, b, c]]
    ab, ac, bc = pearson(ra, rb), pearson(ra, rc), pearson(rb, rc)
    expected = (ab-ac*bc)/np.sqrt((1-ac*ac)*(1-bc*bc))
    check('partial_matches_single_control_identity', abs(p-expected) < 1e-12)
    p, _, _, status = partial(c, b, c[:, None])
    check('partial_full_collinearity_is_undefined', np.isnan(p) and status != 'ok')
    check('constant_spearman_is_undefined', np.isnan(rho([1, 1, 1], [1, 2, 3])[0]))


def load_framewise(package_root):
    new = read('analysis/formal_sv_diameter_v1/framewise_primary.csv')
    base128 = read(D128+'no_background_relative_tail_v1_full2422/framewise_primary.csv')
    obs128 = read(D128+'observed_tail_intensity_full2422/observed_tail_intensity_framewise.csv',
                  'D128 observed raw Q and area only; corrected columns never used')
    loc128 = read(FORMAL128+'localization.csv', 'frozen continuity-first v2.1 geometry')
    window128 = read(D128+'no_background_followup_three_audits_v1_full2422/window_framewise.csv.gz')
    packages = read(FORMAL128+'download_packages.csv', 'frozen release package identities')
    arrays = read(FORMAL128+'arrays_sha256.csv', 'frozen per-frame array identities')
    bands_a = read(STAGE+'d128_d235_source_axial_audit/d128_d235_source_axial_framewise.csv')
    bands_b = read(STAGE+'d500_source_axial_audit/source_axial_band_framewise.csv')
    for name in [STAGE+'provenance.json', STAGE+'validation.json',
                 STAGE+'d128_d235_source_axial_audit/provenance.json',
                 STAGE+'d128_d235_source_axial_audit/validation.json',
                 STAGE+'d500_source_axial_audit/provenance.json',
                 STAGE+'d500_source_axial_audit/validation.json',
                 'src/svrecttail/geometry.py', 'src/svrecttail/raw_sv_quantification.py',
                 STAGE+'d500_source_axial_local_array_audit.py']:
        record(ROOT/name, 'inherited validation, definitions or frozen implementation')
    # Verify the predecessor's explicitly recorded identities, not just file presence.
    for provname, hashkey in [(STAGE+'provenance.json', 'frozen_input_sha256'),
                             (STAGE+'d128_d235_source_axial_audit/provenance.json', 'source_file_sha256')]:
        expected = json.loads((ROOT/provname).read_text(encoding='utf-8'))[hashkey]
        check('frozen_identity_'+hashkey, all(sha(ROOT/p) == h for p, h in expected.items()),
              files_verified=len(expected))
    check('formal_frame_counts', len(new) == 8509 and len(base128) == 2422)
    b = base128.merge(loc128, on=KEYS, validate='one_to_one')
    check('d128_geometry_mapping', len(b) == 2422 and b.z_top_edge_px.notna().all())
    b = b.set_index(KEYS)
    array_hash = arrays.set_index(KEYS).sha256.to_dict()
    replay = []
    array_audit = []
    for pkg in packages.itertuples(index=False):
        p = package_root / pkg.file
        digest = sha(p)
        check('package_'+pkg.file, digest == pkg.sha256 and p.stat().st_size == pkg.bytes)
        record(p, 'retained formal D128 release package; direct raw-SV integration', digest)
        bits = pkg.file.removesuffix('.zip').split('_')
        scan, lo, hi = bits[-3], int(bits[-2]), int(bits[-1])
        wanted = b.loc[scan].loc[lambda g: g.index.to_series().between(lo, hi)]
        print(f'Replay {scan} {lo:03d}-{hi:03d}: {len(wanted)} frozen valid frames', flush=True)
        with zipfile.ZipFile(p) as zf:
            for fi, r in wanted.iterrows():
                member = f'arrays/{scan}/frame_{fi:03d}.npz'
                blob = zf.read(member)
                h = hashlib.sha256(blob).hexdigest()
                if h != array_hash[(scan, fi)]:
                    raise AssertionError(f'NPZ identity mismatch: {member}')
                with np.load(io.BytesIO(blob), allow_pickle=False) as d:
                    sv = np.asarray(d['sv_raw'], dtype=np.float64)
                if sv.shape != (351, 500) or not np.isfinite(sv).all():
                    raise AssertionError(member)
                geom = VesselGeometry(float(r.x_left_edge_px), float(r.x_right_edge_px),
                                      float(r.z_top_edge_px), 128., 12.7, 6.7)
                metrics, _ = quantify_raw_sv(sv, geom, geometry_qc_valid=True)
                replay.append(dict(scan_id=scan, frame_index_0based=fi, diameter_um=128,
                                   flow_mm_s=r.flow_mm_s, geometry_qc_valid=True,
                                   x4_px=r.x4_centroid_isolated_jump_corrected_px,
                                   x_left_edge_px=r.x_left_edge_px, x_right_edge_px=r.x_right_edge_px,
                                   apparent_width_um=geom.lateral_width_um,
                                   apparent_width_px=geom.lateral_width_um/12.7,
                                   z_top_edge_px=r.z_top_edge_px, z_bottom_edge_px=geom.z_bottom_edge_px,
                                   input_sha256=h, **metrics))
                array_audit.append(dict(scan_id=scan, frame_index=fi, package=p.name, member=member,
                                        sha256=h, matches_frozen_identity=True))
    d128 = pd.DataFrame(replay)
    check('d128_array_sha_coverage', len(array_audit) == 2422)
    csv('d128_array_identity_audit.csv', pd.DataFrame(array_audit))
    cmp = d128.merge(base128, on=KEYS, validate='one_to_one', suffixes=('', '_ref'))
    compare('d128_source_mean_vs_formal', cmp.source_mean_raw, cmp.source_mean_raw_ref)
    compare('d128_tail500_vs_formal', cmp.tail_mean_raw_500um, cmp.tail_mean_raw)
    compare('d128_ri_vs_formal', cmp.ri_tail, cmp.ri_tail_ref)
    cmp = d128.merge(obs128, on=KEYS, validate='one_to_one', suffixes=('', '_ref'))
    compare('d128_source_q_vs_formal_observed', cmp.source_q_raw, cmp.q_vessel_observed)
    compare('d128_source_area_vs_formal', cmp.source_area_um2, cmp.source_area_um2_ref)
    compare('d128_tail500_q_vs_formal_observed', cmp.tail_q_raw_500um, cmp.q_tail_observed)
    for depth in [100, 500]:
        cmp = d128.merge(window128[window128.window_um.eq(depth)], on=KEYS, validate='one_to_one', suffixes=('', '_ref'))
        compare(f'd128_tail{depth}_vs_exact_window_mean', cmp[f'tail_mean_raw_{depth}um'], cmp.tail_mean_raw_window)
        compare(f'd128_tail{depth}_vs_exact_window_area', cmp[f'tail_area_um2_{depth}um'], cmp.tail_area_um2)
    d128.scan_id = d128.flow_mm_s.map(lambda f: f'D128_F{int(f):02d}_V01')
    new['input_sha256'] = new.input_mat_sha256
    columns = ID[:3]+KEYS[1:]+['geometry_qc_valid', 'x4_px', 'x_left_edge_px', 'x_right_edge_px',
             'apparent_width_um', 'apparent_width_px', 'z_top_edge_px', 'z_bottom_edge_px',
             'source_area_um2', 'source_q_raw', 'source_mean_raw', 'input_sha256']+TAILS+[
             'tail_area_um2_100um', 'tail_area_um2_500um']
    fw = pd.concat([d128[columns], new[columns]], ignore_index=True)
    bandcols = [f'{s}_{m}' for s in BANDS for m in ['mean_raw', 'q_raw', 'area_um2', 'q_fraction']]
    bands = pd.concat([bands_a[KEYS+bandcols+['source_mean_raw_frozen']],
                       bands_b[KEYS+bandcols+['source_mean_raw_frozen']]], ignore_index=True)
    check('band_key_sets_equal_formal', set(map(tuple, fw[KEYS].values)) == set(map(tuple, bands[KEYS].values)))
    fw = fw.merge(bands, on=KEYS, validate='one_to_one')
    compare('all_source_mean_vs_band_formal', fw.source_mean_raw, fw.source_mean_raw_frozen)
    fw = fw.drop(columns='source_mean_raw_frozen').rename(columns={
        'frame_index_0based': 'frame_index', 'geometry_qc_valid': 'valid_geometry',
        'apparent_width_um': 'X1', 'apparent_width_px': 'X1_px', 'x4_px': 'X4',
        'z_top_edge_px': 'z_top', 'source_area_um2': 'source_area'})
    fw['grid_scope'] = np.where(fw.flow_mm_s.isin(COMMON), 'common_grid', 'd500_extension')
    fw = fw.sort_values(['diameter_um', 'flow_mm_s', 'frame_index']).reset_index(drop=True)
    check('unique_frame_keys', not fw.duplicated(['scan_id', 'frame_index']).any())
    check('valid_frame_complete_metrics', fw.valid_geometry.all() and np.isfinite(fw[VARIABLES]).all().all())
    compare('band_q_reconstruction', fw[[f'{s}_q_raw' for s in BANDS]].sum(axis=1), fw.source_q_raw)
    compare('band_area_reconstruction', fw[[f'{s}_area_um2' for s in BANDS]].sum(axis=1), fw.source_area)
    compare('band_fraction_sum', fw[[f'{s}_q_fraction' for s in BANDS]].sum(axis=1), np.ones(len(fw)))
    compare('source_q_equals_per_frame_mean_times_area', fw.source_q_raw, fw.source_mean_raw*fw.source_area)
    compare('frozen_x4_is_source_center', fw.X4, (fw.x_left_edge_px+fw.x_right_edge_px)/2)
    compare('physical_bottom_frozen', fw.z_bottom_edge_px, fw.z_top+fw.diameter_um/6.7)
    for d in [100, 500]:
        compare(f'tail{d}_q_identity', fw[f'tail_q_raw_{d}um'], fw[f'tail_mean_raw_{d}um']*fw[f'tail_area_um2_{d}um'])
        compare(f'tail{d}_area_identity', fw[f'tail_area_um2_{d}um'], fw.X1*d)
    ref = read(STAGE+'all_volume_metrics_with_d500_extension.csv', 'Stage 2 formal volume medians and coverage')
    counts = fw.groupby('scan_id').size()
    compare('per_volume_counts_vs_stage2', counts.reindex(ref.scan_id), ref.geometry_valid_frames)
    vrows = []
    for r in ref.itertuples(index=False):
        g = fw[fw.scan_id.eq(r.scan_id)]
        vals = {'source_mean_raw_median': g.source_mean_raw.median(),
                'tail_mean_raw_100um_median': g.tail_mean_raw_100um.median(),
                'tail_mean_raw_500um_median': g.tail_mean_raw_500um.median(),
                'ri_tail_median': (g.tail_mean_raw_500um/g.source_mean_raw).median()}
        for metric, value in vals.items():
            vrows.append(dict(scan_id=r.scan_id, metric=metric, recomputed=value,
                              stage2_reference=getattr(r, metric)))
    vt = pd.DataFrame(vrows)
    compare('volume_medians_vs_stage2', vt.recomputed, vt.stage2_reference)
    csv('volume_median_validation.csv', vt)
    check('analysis_scope', fw.scan_id.nunique() == 23 and fw[fw.grid_scope.eq('common_grid')].scan_id.nunique() == 20
          and len(fw) == 10931 and set(fw[fw.grid_scope.eq('d500_extension')].flow_mm_s) == {2, 9, 12})
    csv('framewise_source_tail_metrics.csv', fw)
    coverage = fw.groupby(ID).size().rename('valid_frames').reset_index()
    coverage['nominal_frames'] = 500
    coverage['invalid_geometry_frames'] = 500-coverage.valid_frames
    csv('volume_frame_coverage.csv', coverage)
    return fw


def pair_definitions(include_area=True):
    pairs = []
    for p in PREDICTORS:
        family = 'mean' if '_mean_' in p else 'q'
        for depth in [100, 500]:
            pairs.append((p, f'tail_{family}_raw_{depth}um', depth, family))
    if include_area:
        for family in ['mean', 'q']:
            for depth in [100, 500]:
                pairs.append(('source_area', f'tail_{family}_raw_{depth}um', depth, 'area_vs_'+family))
    return pairs


def analyze(fw):
    zero, diff, lagrows, partialrows = [], [], [], []
    diagnostic_rows = []
    for meta, group in fw.groupby(ID, sort=True):
        info = dict(zip(ID, meta))
        print(f'Coupling {info["scan_id"]}', flush=True)
        raw = group.set_index('frame_index')[VARIABLES].reindex(range(500))
        modes = {'raw': raw, **{f'detrended_{w}': detrend(raw, w) for w in [31, 51, 101]}}
        dd = raw.diff()
        diagnostic_rows.append(info | dict(n_valid=len(group), n_adjacent_valid_pairs=int(dd.source_mean_raw.notna().sum())))
        for mode, data in modes.items():
            for predictor, tail, depth, family in pair_definitions():
                r, n = rho(data[predictor], data[tail])
                zero.append(info | dict(mode=mode, predictor=predictor, tail_variable=tail,
                                       tail_depth_um=depth, metric_family=family, rho=r, n_pairs=n))
            for predictor, tail, depth, family in pair_definitions(False):
                x, y = data[predictor].to_numpy(), data[tail].to_numpy()
                for lag in range(-50, 51):
                    r, n = rho(*lag_pair(x, y, lag))
                    lagrows.append(info | dict(mode=mode, predictor=predictor, tail_variable=tail,
                           tail_depth_um=depth, metric_family=family, lag=lag, rho=r, n_pairs=n))
            # Also report detrended controls, using the same residual transform for each control.
            for predictor in [f'{s}_mean_raw' for s in ['source']+BANDS]:
                for depth in [100, 500]:
                    tail = f'tail_mean_raw_{depth}um'
                    for controls in [['source_area'], ['z_top'], ['source_area', 'z_top']]:
                        r, n, design_rank, status = partial(data[predictor], data[tail], data[controls].to_numpy())
                        partialrows.append(info | dict(mode=mode, predictor=predictor,
                            tail_variable=tail, tail_depth_um=depth, controls='+'.join(controls),
                            partial_rho=r, n_pairs=n, control_design_rank=design_rank, status=status))
        for predictor, tail, depth, family in pair_definitions():
            r, n = rho(dd[predictor], dd[tail])
            diff.append(info | dict(mode='first_difference', predictor=predictor, tail_variable=tail,
                                   tail_depth_um=depth, metric_family=family, rho=r, n_pairs=n))
    zero, diff, lag, part = map(pd.DataFrame, [zero, diff, lagrows, partialrows])
    csv('volume_zero_lag_coupling.csv', zero[zero['mode'].eq('raw') & ~zero.predictor.str.startswith(tuple(BANDS))])
    csv('band_zero_lag_coupling.csv', zero[zero['mode'].eq('raw') & zero.predictor.str.startswith(tuple(BANDS))])
    csv('volume_detrended_coupling.csv', zero[~zero['mode'].eq('raw')])
    csv('volume_first_difference_coupling.csv', diff)
    csv('partial_spearman_area_ztop.csv', part)
    csv('volume_pair_diagnostics.csv', pd.DataFrame(diagnostic_rows))
    csv('lag_correlation_raw.csv', lag[lag['mode'].eq('raw')])
    csv('lag_correlation_detrended_51.csv', lag[lag['mode'].eq('detrended_51')])
    csv('lag_correlation_detrended_sensitivity.csv', lag[lag['mode'].isin(['detrended_31', 'detrended_101'])])
    peaks = []
    grouping = ID+['mode', 'predictor', 'tail_variable', 'tail_depth_um', 'metric_family']
    for key, g in lag.groupby(grouping):
        g = g.sort_values('lag')
        z = g[g.lag.eq(0)].iloc[0]
        maxrho = g.rho.max()
        ties = g[g.rho.eq(maxrho)].copy()
        ties['abs_lag'] = ties.lag.abs()
        peak = ties.sort_values(['abs_lag', 'lag']).iloc[0]
        remote = g[g.lag.abs().between(30, 50)].rho.median()
        peaks.append(dict(zip(grouping, key)) | dict(rho_lag0=z.rho, rho_peak=peak.rho,
            lag_at_peak=int(peak.lag), n_pairs_at_lag0=int(z.n_pairs), n_pairs_at_peak=int(peak.n_pairs),
            rho_remote_median=remote, zero_lag_specificity=z.rho-remote,
            zero_lag_rank=int(1+(g.rho > z.rho).sum()), peak_within_5=bool(abs(peak.lag) <= 5),
            peak_tie_count=len(ties), n_finite_lags=int(g.rho.notna().sum())))
    peaks = pd.DataFrame(peaks)
    csv('lag_peak_summary.csv', peaks)
    check('all_lag_curve_counts', len(lag) == 23*16*4*101 and len(peaks) == 23*16*4)
    check('all_requested_coefficients_defined', np.isfinite(lag.rho).all() and np.isfinite(zero.rho).all()
          and np.isfinite(diff.rho).all() and np.isfinite(part.partial_rho).all())
    z = lag[lag.lag.eq(0)].merge(zero, on=grouping, suffixes=('_lag', '_zero'), validate='one_to_one')
    compare('lag0_matches_zero_lag_tables', z.rho_lag, z.rho_zero)
    check('lag_pair_counts_bounded', (lag.n_pairs <= 500-lag.lag.abs()).all() and (lag.n_pairs >= 3).all())
    # Matched, per-volume 100 versus 500 comparisons, for all predictors and modes.
    keys = ID+['mode', 'predictor', 'metric_family']
    comp = peaks[peaks.tail_depth_um.eq(100)].merge(peaks[peaks.tail_depth_um.eq(500)],
                                                 on=keys, suffixes=('_100', '_500'), validate='one_to_one')
    for metric in ['rho_lag0', 'rho_peak', 'zero_lag_specificity']:
        comp[metric+'_100_minus_500'] = comp[metric+'_100']-comp[metric+'_500']
    comp['abs_peak_lag_100_minus_500'] = comp.lag_at_peak_100.abs()-comp.lag_at_peak_500.abs()
    csv('tail100_vs_tail500_volume_comparison.csv', comp)
    bandcomp = []
    for key, g in zero[zero.predictor.str.startswith(tuple(BANDS))].groupby(ID+['mode','tail_depth_um','metric_family']):
        vals = {r.predictor.split('_')[0]: r.rho for r in g.itertuples()}
        bandcomp.append(dict(zip(ID+['mode','tail_depth_um','metric_family'], key)) |
                        {f'{s}_rho': vals[s] for s in BANDS} |
                        dict(middle_minus_upper=vals['middle']-vals['upper'],
                             lower_minus_upper=vals['lower']-vals['upper'],
                             upper_minus_lower=vals['upper']-vals['lower'],
                             middle_minus_lower=vals['middle']-vals['lower']))
    bandcomp = pd.DataFrame(bandcomp)
    csv('band_volume_comparison.csv', bandcomp)
    summarize(zero, diff, part, peaks, comp, bandcomp)
    return zero, diff, part, lag, peaks, comp, bandcomp


def summarize(zero, diff, part, peaks, comp, bandcomp):
    records = []
    def add(table, metric, analysis, dims):
        for key, g in table.groupby(['grid_scope', 'diameter_um']+dims, dropna=False):
            g = g.sort_values('flow_mm_s')
            vals = g[metric].to_numpy(float)
            finite = np.isfinite(vals)
            records.append(dict(zip(['grid_scope','diameter_um']+dims, key)) |
                dict(analysis=analysis, metric=metric, n_volumes=len(g), n_finite=int(finite.sum()),
                     volume_values_json=json.dumps([dict(scan_id=r.scan_id, flow_mm_s=float(r.flow_mm_s),
                         value=float(getattr(r, metric)) if np.isfinite(getattr(r, metric)) else None) for r in g.itertuples()]),
                     median=float(np.nanmedian(vals)), min=float(np.nanmin(vals)), max=float(np.nanmax(vals)),
                     positive_count=int((vals > 0).sum()), negative_count=int((vals < 0).sum()),
                     zero_count=int((vals == 0).sum()),
                     peak_within_5_count=int(g.peak_within_5.sum()) if 'peak_within_5' in g else np.nan))
    dims = ['mode', 'predictor', 'tail_depth_um', 'metric_family']
    add(zero, 'rho', 'zero_lag', dims)
    add(diff, 'rho', 'first_difference', dims)
    add(part, 'partial_rho', 'partial_spearman', ['mode','predictor','tail_depth_um','controls'])
    for m in ['rho_lag0','rho_peak','lag_at_peak','zero_lag_specificity','zero_lag_rank']:
        add(peaks, m, 'lag_curve', dims)
    for m in ['rho_lag0_100_minus_500','rho_peak_100_minus_500','zero_lag_specificity_100_minus_500','abs_peak_lag_100_minus_500']:
        add(comp, m, 'tail_depth_paired_difference', ['mode','predictor','metric_family'])
    for m in ['middle_minus_upper','lower_minus_upper','upper_minus_lower','middle_minus_lower']:
        add(bandcomp, m, 'band_paired_difference', ['mode','tail_depth_um','metric_family'])
    s = pd.DataFrame(records)
    check('summary_units_and_separation', s[s.grid_scope.eq('common_grid')].n_volumes.eq(5).all()
          and s[s.grid_scope.eq('d500_extension')].n_volumes.eq(3).all())
    csv('diameter_common_flow_summary.csv', s[s.grid_scope.eq('common_grid')])
    csv('d500_extension_summary.csv', s[s.grid_scope.eq('d500_extension')])


def figures(fw, lag, peaks, zero, comp):
    plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False})
    colors = {100: '#2374ab', 500: '#d36d32'}
    def save(fig, name):
        fig.savefig(OUT/(name+'.png'), dpi=160, bbox_inches='tight')
        fig.savefig(OUT/(name+'.pdf'), bbox_inches='tight')
        plt.close(fig)
    fig, axes = plt.subplots(4, 1, figsize=(11, 10), sharex=True)
    for ax, d in zip(axes, [128,235,285,500]):
        g = fw[fw.diameter_um.eq(d) & fw.flow_mm_s.eq(5)].set_index('frame_index').reindex(range(500))
        for col, label, color in [('source_mean_raw','source','#333333'),('tail_mean_raw_100um','tail100',colors[100]),('tail_mean_raw_500um','tail500',colors[500])]:
            ax.plot(g.index, g[col]/g[col].median(), lw=.85, label=label, color=color)
        ax.set_title(f'D{d}, flow 5 mm/s', loc='left'); ax.set_ylabel('Value / median'); ax.legend(ncol=3)
    axes[-1].set_xlabel('Original B-scan index (0-based)')
    fig.suptitle('Figure 1. Display normalization only; quantitative input is linear raw SV')
    fig.tight_layout(); save(fig, 'figure1_framewise_flow5')
    for mode, num in [('raw', 2),('detrended_51',3)]:
        fig, axes = plt.subplots(2,2, figsize=(10,7), sharex=True, sharey=True)
        for ax,d in zip(axes.flat,[128,235,285,500]):
            g = lag[lag.diameter_um.eq(d) & lag.flow_mm_s.eq(5) & lag['mode'].eq(mode) & lag.predictor.eq('source_mean_raw')]
            for dep in [100,500]:
                h=g[g.tail_depth_um.eq(dep)]; ax.plot(h.lag,h.rho,label=f'tail{dep}', color=colors[dep])
            ax.axvline(0,color='gray',lw=.8); ax.axhline(0,color='gray',lw=.5)
            ax.set_title(f'D{d}, flow 5'); ax.set_xlabel('Lag (B-scans)'); ax.set_ylabel('Descriptive Spearman'); ax.legend()
        fig.suptitle(f'Figure {num}. Source mean / tail mean: {mode}')
        fig.tight_layout(); save(fig,f'figure{num}_lag_{mode}_flow5')
    primary = peaks[peaks.grid_scope.eq('common_grid') & peaks['mode'].eq('detrended_51') & peaks.predictor.eq('source_mean_raw')]
    fig,axes=plt.subplots(2,2,figsize=(12,7),sharex=True)
    for row,metric in enumerate(['rho_lag0','zero_lag_specificity']):
        for col,dep in enumerate([100,500]):
            ax=axes[row,col]
            for j,d in enumerate([128,235,285,500]):
                g=primary[primary.diameter_um.eq(d)&primary.tail_depth_um.eq(dep)].sort_values('flow_mm_s')
                xs=j+np.linspace(-.18,.18,5); ax.scatter(xs,g[metric],color=colors[dep])
                for xx,r in zip(xs,g.itertuples()): ax.annotate(str(int(r.flow_mm_s)),(xx,getattr(r,metric)),xytext=(3,3),textcoords='offset points',fontsize=7)
            ax.axhline(0,color='gray',lw=.5);ax.set_xticks(range(4),['D128','D235','D285','D500']);ax.set_ylabel(metric);ax.set_title(f'Tail{dep}')
    fig.suptitle('Figure 4. 51-frame residuals; every common-flow volume (labels = mm/s)')
    fig.tight_layout();save(fig,'figure4_all_common_volumes')
    fig,axes=plt.subplots(2,4,figsize=(15,7),sharey=True)
    for col,d in enumerate([128,235,285,500]):
        for row,family in enumerate(['mean','q']):
            ax=axes[row,col]
            for dep in [100,500]:
                for flow in COMMON:
                    vals=[]
                    for band in BANDS:
                        h=zero[zero.diameter_um.eq(d)&zero.flow_mm_s.eq(flow)&zero['mode'].eq('detrended_51')&zero.tail_depth_um.eq(dep)&zero.predictor.eq(f'{band}_{family}_raw')]
                        vals.append(float(h.rho.iloc[0]))
                    ax.plot(range(3),vals,color=colors[dep],alpha=.5,lw=.8,marker='o',ms=3,label=f'tail{dep}' if flow==1 else None)
            ax.set_xticks(range(3),['Upper','Middle','Lower']);ax.axhline(0,color='gray',lw=.5);ax.set_title(f'D{d}, {family}');ax.set_ylabel('Descriptive Spearman');ax.legend()
    fig.suptitle('Figure 5. Axial-band 51-frame detrended zero-lag coupling; all common flows')
    fig.tight_layout();save(fig,'figure5_axial_bands')
    fig,axes=plt.subplots(1,2,figsize=(11,5))
    for ax,family in zip(axes,['mean','q']):
        for d in [128,235,285,500]:
            g=comp[comp.grid_scope.eq('common_grid')&comp['mode'].eq('detrended_51')&comp.predictor.eq(f'source_{family}_raw')&comp.diameter_um.eq(d)]
            ax.scatter(g.zero_lag_specificity_500,g.zero_lag_specificity_100,label=f'D{d}')
        lim=[min(ax.get_xlim()[0],ax.get_ylim()[0]),max(ax.get_xlim()[1],ax.get_ylim()[1])]
        ax.plot(lim,lim,color='gray',ls='--');ax.set_xlim(lim);ax.set_ylim(lim);ax.set_xlabel('Tail500 zero-lag specificity');ax.set_ylabel('Tail100 zero-lag specificity');ax.set_title(f'Source {family}, 51-frame residuals');ax.legend()
    fig.suptitle('Figure 6. Matched tail-depth comparison; all 20 common-flow volumes')
    fig.tight_layout();save(fig,'figure6_tail_depth_specificity')


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--d128-package-root',type=Path,required=True)
    ap.add_argument('--reuse-verified-framewise',action='store_true',help='Reuse this run\'s metrics only after manifest identity verification')
    args=ap.parse_args()
    self_checks()
    if args.reuse_verified_framewise:
        prior=json.loads((OUT/'validation.json').read_text(encoding='utf-8'))
        check('reused_validation_previously_passed', prior['status']=='passed')
        previous=pd.read_csv(OUT/'input_manifest.csv')
        check('reused_input_identities_unchanged',all(sha(ROOT/r.input_path)==r.sha256 for r in previous.itertuples()))
        check('reused_framewise_identity',sha(OUT/'framewise_source_tail_metrics.csv')==prior['framewise_sha256'])
        MANIFEST.extend(previous.to_dict('records'))
        CHECKS.update(prior['checks'])
        fw=pd.read_csv(OUT/'framewise_source_tail_metrics.csv',float_precision='round_trip')
    else:
        fw=load_framewise(args.d128_package_root.resolve())
        # Persist a validated data checkpoint before correlations, without claiming whole-run completion.
        csv('input_manifest.csv',pd.DataFrame(MANIFEST).drop_duplicates('input_path'))
        write_json('validation.json',dict(status='data_gate_passed_analysis_pending',checks=CHECKS,
                                         framewise_sha256=sha(OUT/'framewise_source_tail_metrics.csv')))
    zero,diff,part,lag,peaks,comp,bandcomp=analyze(fw)
    figures(fw,lag,peaks,zero,comp)
    csv('input_manifest.csv',pd.DataFrame(MANIFEST).drop_duplicates('input_path'))
    provenance=dict(repository='bulbel-magnolia/OCTA-',branch='analysis/sv-diameter-stage2',
        starting_head=START_HEAD,result_commit='The Git commit containing this directory; reported after commit/push.',
        script='analysis/sv_source_tail_framewise_coupling_v1/run_analysis.py',script_sha256=sha(__file__),
        completed_utc=datetime.now(timezone.utc).isoformat(),
        software=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__),
        signal='SV_raw = Var_t(abs(E)), denominator N; linear raw SV',
        geometry=dict(version='continuity-first v2.1',source='X1 apparent lateral span x physical diameter ellipse',
                      lateral_center='frozen X4',dx_um=12.7,dz_um=6.7,supersample=16,
                      bands='physical axial thirds with frozen subpixel weights',tail='X1 rectangle at true physical bottom',guard_um=0,tail_windows_um=[100,500]),
        analysis=dict(experimental_unit='scan volume',frame_role='within-volume spatial position',
            common_flows_mm_s=COMMON,d500_extension_flows_mm_s=[2,9,12],frame_index_base=0,nominal_frame_grid=[0,499],
            detrend_windows=[31,51,101],primary_detrend_window=51,moving_median='centered on original integer grid, min_periods=1; edge truncated; NaN ignored for baseline, invalid center stays NaN',
            lag_definition='Source(i) versus Tail(i+lag)',lags=[-50,50],spearman='Pearson of pairwise-valid average ranks; min 3 pairs',
            specificity='rho(0) - median rho at absolute lag 30 through 50, inclusive',
            peak='maximum signed rho; exact ties resolved by smallest absolute lag then smaller signed lag',
            zero_rank='competition rank: 1 + number of strictly greater rho values',
            partial='average ranks then OLS with intercept and indicated controls; Pearson of residuals; rank is re-computed in each pairwise-complete subset',
            partial_modes='raw plus residual signals and residual controls at 31/51/101',
            no_p_values=True,no_frame_pooled_correlation=True,no_background_subtraction=True,no_raw_signal_normalization=True),
        d128_q='Direct quantify_raw_sv replay from SHA-verified retained formal NPZ; never median-times-area',
        reuse_verified_framewise=bool(args.reuse_verified_framewise))
    write_json('provenance.json',provenance)
    coverage=fw.groupby('grid_scope').agg(volumes=('scan_id','nunique'),frames=('scan_id','size')).to_dict('index')
    write_json('validation.json',dict(status='passed',checks=CHECKS,coverage=coverage,
        counts_by_diameter={str(int(k)):int(v) for k,v in fw.groupby('diameter_um').size().items()},
        framewise_sha256=sha(OUT/'framewise_source_tail_metrics.csv'),
        missing_required_data=[],unfinished_analysis_items=[],p_values_computed=False,
        raw_sv_reconstruction_from_oct=False,background_subtraction=False,
        inferential_statistics=False,invalid_frames_zero_filled=False,
        rows=dict(framewise=len(fw),zero_lag=len(zero),first_difference=len(diff),partial=len(part),lag=len(lag),lag_peak=len(peaks))))
    print(json.dumps(dict(status='passed',coverage=coverage,band_q=CHECKS['band_q_reconstruction']),indent=2),flush=True)


if __name__=='__main__':
    main()
