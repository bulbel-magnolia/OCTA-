import numpy as np
from svrecttail.mentor import tracking_core as core


def trajectory(values):
    values = np.array(values, float)
    cfg = core.merge_tracking_config({'upper_edge_method': core.UPPER_EDGE_CONTINUITY,
        'max_model_assisted_gap_frames': 5})
    return core.continuity_first_trajectory(values, np.isfinite(values), cfg)


def test_large_jump_rejected_before_smoothing():
    model, accepted, _, _, audit = trajectory([190]*10 + [235] + [190]*10)
    assert not accepted[10]
    assert audit['z_continuity_status'][10] == 'rejected_continuity'
    np.testing.assert_allclose(model, 190)


def test_reliable_gap_up_to_five_frames_is_filled():
    for gap in range(1, 6):
        model, accepted, _, _, audit = trajectory([190]*10 + [np.nan]*gap + [192]*10)
        assert np.isfinite(model).all()
        assert audit['z_short_gap_filled'][10:10+gap].all()
        assert not accepted[10:10+gap].any()
        assert np.max(np.abs(np.diff(np.rint(model)))) <= 2


def test_long_gap_stays_na_and_needs_three_consecutive_candidates():
    model, accepted, _, _, audit = trajectory([190]*10 + [np.nan]*6 + [199, 197, 197])
    assert np.isnan(model[10:16]).all()
    assert accepted[16:19].all()
    assert audit['z_relock_start'][16]
    np.testing.assert_array_equal(audit['z_relock_confirmation_frame'][16:19], [18]*3)
    short, accepted_short, *_ = trajectory([190]*10 + [np.nan]*6 + [199, 197])
    assert np.isnan(short[10:]).all()
    assert not accepted_short[16:].any()


def test_relock_counter_resets_at_missing_or_discontinuous_candidate():
    model, accepted, *_ = trajectory([190]*10 + [np.nan]*6 + [199, 197, np.nan, 197, 235, 235])
    assert np.isnan(model[10:]).all()
    assert not accepted[16:].any()


def test_smoothing_does_not_cross_invalid_segments():
    first, *_ = trajectory([190]*10 + [np.nan]*6 + [200]*10)
    second, *_ = trajectory([190]*10 + [np.nan]*6 + [240]*10)
    np.testing.assert_allclose(first[:10], second[:10])
    assert np.isnan(first[10:16]).all()
    np.testing.assert_allclose(first[16:], 200)
    np.testing.assert_allclose(second[16:], 240)


def test_no_endpoint_extrapolation_or_all_missing_crash():
    model, *_ = trajectory([np.nan]*3 + [190]*10 + [np.nan]*3)
    assert np.isnan(model[:3]).all() and np.isnan(model[-3:]).all()
    empty, accepted, *_ = trajectory([np.nan]*20)
    assert np.isnan(empty).all() and not accepted.any()


def test_prediction_gate_is_separate_from_per_frame_slope_gate():
    model, accepted, _, _, audit = trajectory([190]*10 + [np.nan]*2 + [194] + [190]*10)
    assert abs(194 - 190) <= 2 * 3
    assert not accepted[12]
    assert audit['z_local_prediction_px'][12] == 190
    np.testing.assert_allclose(model, 190)


def test_pending_restart_clears_prediction_from_abandoned_run():
    _, accepted, _, _, audit = trajectory([190]*10 + [np.nan]*6 + [210, 240, 240, 240])
    assert accepted[17:20].all()
    assert np.isnan(audit['z_local_prediction_px'][17])
    assert audit['z_relock_confirmation_frame'][17] == 19


def test_initial_candidate_is_preserved_without_endpoint_extrapolation():
    model, accepted, *_ = trajectory([190, 190, 190, np.nan])
    assert accepted[0]
    assert model[0] == 190
    assert np.isnan(model[-1])


def _synthetic_volume():
    rng = np.random.default_rng(19)
    volume = rng.normal(20, 0.3, (35, 100, 60)).astype(np.float32)
    volume[:, 30:50, 27:34] += 15
    cfg = core.merge_tracking_config({'upper_edge_method':core.UPPER_EDGE_CONTINUITY,
        'upper_edge_alphas':[0.2], 'primary_alpha':0.2, 'axial_um_per_px':1,
        'lateral_um_per_px':5, 'z_search_px':[10,90], 'z_viterbi_max_jump_px':2,
        'z_viterbi_jump_penalty':0.6,'max_model_assisted_gap_frames':5})
    return volume, cfg


def test_core_feedback_resists_later_bright_deep_tail(monkeypatch):
    volume, cfg = _synthetic_volume()
    volume[12:, 65:80, 27:34] += 90
    # Isolate the feedback stage from the documented global bootstrap limitation.
    monkeypatch.setattr(core, 'locate_peak_path', lambda profiles,cfg,diameter: np.full(len(volume),40))
    outputs, _, _ = core.track_volume(volume, scan_id='synthetic', diameter_um=20, config=cfg)
    d = outputs[0.2]
    assert d.z_upper_px.notna().sum() >= 30
    assert d.z_upper_px.max() < 40
    steps = d.z_core_path_px.diff().abs()
    assert steps[~d.z_core_restart].max() <= 2
    assert d.z_core_path_px[12:].max() < 55


def test_completely_unavailable_volume_keeps_invalid_schema(monkeypatch):
    from dataclasses import fields
    from svrecttail.mentor_tracking import build_mentor_tracking_tables, REQUIRED_TRACKING_COLUMNS
    volume, cfg = _synthetic_volume()
    invalid = core.UpperEdgeCandidate(**{f.name: np.nan for f in fields(core.UpperEdgeCandidate)} |
                                      {'valid':False,'high_component_width_px':0})
    monkeypatch.setattr(core, 'persistent_upper_edge_candidate', lambda *args, **kwargs: invalid)
    d, _, _ = build_mentor_tracking_tables(volume, scan_id='missing', diameter_um=20, tracking_config=cfg)
    assert REQUIRED_TRACKING_COLUMNS <= set(d.columns)
    assert d.z_upper_px.isna().all()
    assert (d.new_tracking_class == 'failed').all()
    assert (d.vessel_presence_prediction == 'not_assessable').all()
    assert not d.qc_valid.any()
