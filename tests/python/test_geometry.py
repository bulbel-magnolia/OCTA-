import numpy as np

from svrecttail.background import background_columns
from svrecttail.geometry import (
    VesselGeometry,
    ellipse_weights,
    interval_overlap_weights,
    rectangle_weights,
)


def test_interval_overlap_uses_pixel_edges() -> None:
    weights = interval_overlap_weights(5, 0.25, 2.75)
    np.testing.assert_allclose(weights, [0.25, 1.0, 1.0, 0.25, 0.0])
    assert weights.sum() == 2.5


def test_rectangle_area_matches_continuous_area() -> None:
    weights = rectangle_weights(
        (20, 30),
        x_left_edge_px=3.2,
        x_right_edge_px=12.8,
        z_top_edge_px=2.1,
        z_bottom_edge_px=9.9,
    )
    assert np.isclose(weights.sum(), (12.8 - 3.2) * (9.9 - 2.1))


def test_physical_ellipse_area_is_accurate() -> None:
    geometry = VesselGeometry(
        x_left_edge_px=20.0,
        x_right_edge_px=30.0,
        z_top_edge_px=15.0,
        diameter_um=12.8,
        dx_um=1.0,
        dz_um=1.0,
    )
    weights = ellipse_weights((60, 60), geometry, supersample=64)
    expected_area_px2 = np.pi * 5.0 * 6.4
    assert np.isclose(weights.sum(), expected_area_px2, rtol=0.003)


def test_background_columns_follow_three_skip_five_take_rule() -> None:
    geometry = VesselGeometry.from_inclusive_centres(
        x_left_center_px=20,
        x_right_center_px=29,
        z_top_center_px=10,
        diameter_um=20,
        dx_um=1,
        dz_um=1,
    )
    left, right = background_columns(80, geometry)
    np.testing.assert_array_equal(left, np.arange(12, 17))
    np.testing.assert_array_equal(right, np.arange(33, 38))
