import numpy as np

from svrecttail.localization import localize_geometry, localize_lateral_body


def _synthetic_omag() -> np.ndarray:
    frame = np.ones((90, 70), dtype=float)
    frame[20:40, 25:35] = 11.0
    return frame


def test_x1_local_geometry_returns_body_boundaries_and_midpoint() -> None:
    result = localize_lateral_body(
        _synthetic_omag(),
        x_anchor_center_px=30.0,
        z_top_anchor_center_px=20.0,
        diameter_um=20.0,
        dx_um=2.0,
        dz_um=1.0,
    )
    assert result.valid
    assert not result.fallback
    assert result.x_left_center_px == 25
    assert result.x_right_center_px == 34
    assert result.x1_local_geometry_center_px == 29.5
    assert result.qc_valid


def test_geometry_reestablishes_top_on_x1_central_line() -> None:
    result = localize_geometry(
        _synthetic_omag(),
        x_anchor_center_px=30.0,
        z_anchor_center_px=21.0,
        diameter_um=20.0,
        dx_um=2.0,
        dz_um=1.0,
    )
    assert result.top_edge.z_top_center_px == 20.0
    assert result.geometry.z_top_edge_px == 19.5
    assert result.geometry.z_bottom_edge_px == 39.5
    assert result.source_qc_valid


def test_x1_strongest_run_fallback_is_qc_flagged() -> None:
    result = localize_lateral_body(
        _synthetic_omag(),
        x_anchor_center_px=15.0,
        z_top_anchor_center_px=20.0,
        diameter_um=20.0,
        dx_um=2.0,
        dz_um=1.0,
    )
    assert result.valid
    assert result.fallback
    assert not result.qc_valid
