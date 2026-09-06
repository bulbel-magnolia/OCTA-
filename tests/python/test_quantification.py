import numpy as np

from svrecttail.background import estimate_background
from svrecttail.geometry import VesselGeometry, ellipse_weights
from svrecttail.quantification import quantify_frame


def _geometry() -> VesselGeometry:
    return VesselGeometry.from_inclusive_centres(
        x_left_center_px=30,
        x_right_center_px=39,
        z_top_center_px=20,
        diameter_um=20,
        dx_um=1,
        dz_um=1,
    )


def test_background_is_pixel_count_weighted_bilateral_mean() -> None:
    geometry = _geometry()
    image = np.zeros((5, 80), dtype=float)
    image[:, 22:27] = 2.0
    image[:, 43:48] = 10.0
    image[0, 22] = np.nan
    background = estimate_background(image, geometry)
    assert background.left_valid_count[0] == 4
    assert background.right_valid_count[0] == 5
    assert np.isclose(background.combined[0], 58.0 / 9.0)
    assert np.isclose(background.combined[1], 6.0)


def test_constant_image_background_correction_is_exactly_zero() -> None:
    geometry = _geometry()
    result = quantify_frame(
        np.full((120, 80), 7.25), geometry, tail_length_um=10.0
    )
    np.testing.assert_array_equal(result.corrected_sv, 0.0)
    assert result.q_vessel == 0.0
    assert result.q_tail == 0.0
    assert np.isnan(result.ratio_tail_to_vessel)
    assert "q_vessel_nonpositive" in result.invalid_reason


def test_formal_integrals_and_profile_identity() -> None:
    geometry = _geometry()
    rows = np.arange(160, dtype=float)
    baseline = 5.0 + 0.01 * rows
    image = np.repeat(baseline[:, None], 80, axis=1)
    source_weight = ellipse_weights(image.shape, geometry, supersample=32)
    source_amplitude = 4.0
    image[source_weight > 0] += source_amplitude
    tail_amplitude = 2.5
    image[40:80, 30:40] += tail_amplitude

    result = quantify_frame(
        image,
        geometry,
        tail_length_um=40.0,
        ellipse_supersample=32,
    )
    assert result.valid
    assert np.isclose(result.q_vessel, source_amplitude * source_weight.sum())
    assert np.isclose(result.source_mean, source_amplitude)
    assert np.isclose(result.source_area_um2, source_weight.sum())
    assert np.isclose(result.q_tail, tail_amplitude * 10 * 40)
    assert result.tail_area_um2 == result.requested_tail_area_um2 == 400.0
    assert result.summary()["lateral_width_um"] == 10.0
    assert np.isclose(result.q_tail, result.q_tail_direct_check)
    np.testing.assert_allclose(result.tail_linear_density[40:80], tail_amplitude * 10)
    np.testing.assert_allclose(
        result.tail_linear_density[40:80],
        result.tail_contrast_profile[40:80] * geometry.lateral_width_um,
    )
    assert np.isclose(result.ratio_tail_to_vessel, result.q_tail / result.q_vessel)


def test_negative_background_residual_is_retained() -> None:
    geometry = _geometry()
    image = np.full((120, 80), 5.0)
    source_weight = ellipse_weights(image.shape, geometry, supersample=16)
    image[source_weight > 0] += 3.0
    image[40:50, 30:40] += 2.0
    image[45, 30] = 4.0
    result = quantify_frame(image, geometry, tail_length_um=10.0)
    assert result.corrected_sv[45, 30] == -1.0
    expected = (10 * 10 * 2.0) - 3.0
    assert np.isclose(result.q_tail, expected)


def test_incomplete_tail_window_is_na_not_truncated_result() -> None:
    geometry = _geometry()
    image = np.full((80, 80), 5.0)
    source_weight = ellipse_weights(image.shape, geometry)
    image[source_weight > 0] += 3.0
    result = quantify_frame(image, geometry, tail_length_um=500.0)
    assert not result.valid
    assert not result.tail_window_complete
    assert np.isnan(result.q_tail)
    assert np.isnan(result.ratio_tail_to_vessel)
    assert "tail_window_out_of_frame" in result.invalid_reason


def test_failed_source_qc_blocks_ratio() -> None:
    geometry = _geometry()
    image = np.full((120, 80), 5.0)
    source_weight = ellipse_weights(image.shape, geometry)
    image[source_weight > 0] += 3.0
    result = quantify_frame(
        image,
        geometry,
        tail_length_um=10.0,
        source_qc_valid=False,
    )
    assert np.isfinite(result.q_vessel)
    assert np.isnan(result.ratio_tail_to_vessel)
    assert "source_qc_failed" in result.invalid_reason


def test_unapproved_single_side_background_is_qc_invalid() -> None:
    geometry = _geometry()
    image = np.full((120, 80), 5.0)
    source_weight = ellipse_weights(image.shape, geometry)
    image[source_weight > 0] += 3.0
    image[:, 22:27] = np.nan
    result = quantify_frame(image, geometry, tail_length_um=10.0)
    assert not result.background_complete
    assert np.isnan(result.q_tail)
    assert "background_incomplete" in result.invalid_reason
    assert set(result.background.row_mode) == {"right_only_unapproved"}


def test_explicit_single_side_background_is_audited_and_valid() -> None:
    geometry = _geometry()
    image = np.full((120, 80), 5.0)
    source_weight = ellipse_weights(image.shape, geometry)
    image[source_weight > 0] += 3.0
    image[:, 22:27] = np.nan
    result = quantify_frame(
        image,
        geometry,
        tail_length_um=10.0,
        background_excluded_side="left",
    )
    assert result.background_complete
    assert result.valid
    assert result.background.excluded_side == "left"
    assert set(result.background.row_mode) == {"right_only_explicit"}
