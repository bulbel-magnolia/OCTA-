import numpy as np

from svrecttail.detection import detect_tail_extent


def _blanks() -> np.ndarray:
    offsets = np.array([-0.2, -0.1, 0.0, 0.1, 0.2])
    return np.repeat(offsets[:, None], 20, axis=1)


def test_detects_furthest_sustained_pair_of_bins() -> None:
    signal = np.r_[np.ones(15), np.zeros(5)]
    result = detect_tail_extent(signal, _blanks(), dz_um=2.0)
    assert result.detected
    assert result.detectable_length_um == 30.0
    assert not result.right_censored
    np.testing.assert_array_equal(result.sustained_detection, [True, True, True, False])


def test_single_above_threshold_bin_is_rejected() -> None:
    signal = np.r_[np.ones(5), np.zeros(15)]
    result = detect_tail_extent(signal, _blanks(), dz_um=2.0)
    assert not result.detected
    assert np.isnan(result.detectable_length_um)
    assert result.invalid_reason == "not_detected"


def test_detection_at_final_bin_is_right_censored() -> None:
    signal = np.ones(20)
    result = detect_tail_extent(signal, _blanks(), dz_um=2.0)
    assert result.detected
    assert result.detectable_length_um == 40.0
    assert result.right_censored


def test_zero_blank_mad_is_not_evaluable() -> None:
    signal = np.ones(20)
    blanks = np.zeros((5, 20))
    result = detect_tail_extent(signal, blanks, dz_um=2.0)
    assert not result.detected
    assert np.isnan(result.detectable_length_um)
    assert result.invalid_reason == "degenerate_blank_scale"
