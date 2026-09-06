import numpy as np

from svrecttail.localization import (
    find_surface_guided_anchor,
    localize_geometry,
    localize_geometry_from_surface,
    localize_lateral_body,
)


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


def test_surface_guided_search_finds_global_vessel_without_manifest_anchor() -> None:
    frame = np.ones((100, 100), dtype=float)
    frame[30:50, 45:55] = 11.0
    coarse = find_surface_guided_anchor(
        frame,
        surface_z_center_px=10.0,
        surface_to_vessel_top_um=20.0,
        effective_refractive_index=1.0,
        diameter_um=20.0,
        dx_um=2.0,
        dz_um=1.0,
        edge_exclusion_px=5,
    )
    assert coarse.valid
    assert coarse.predicted_z_top_center_px == 30.0
    assert 49.0 <= coarse.x_anchor_center_px <= 50.0

    result = localize_geometry_from_surface(
        frame,
        surface_z_center_px=10.0,
        surface_to_vessel_top_um=20.0,
        effective_refractive_index=1.0,
        diameter_um=20.0,
        dx_um=2.0,
        dz_um=1.0,
        edge_exclusion_px=5,
    )
    assert result.coarse_anchor is not None
    assert result.lateral.x_left_center_px == 45
    assert result.lateral.x_right_center_px == 54
    assert result.top_edge.z_top_center_px == 30.0
    assert result.source_qc_valid


def test_surface_guided_search_rejects_a_blank_frame() -> None:
    result = find_surface_guided_anchor(
        np.ones((100, 100), dtype=float),
        surface_z_center_px=10.0,
        surface_to_vessel_top_um=20.0,
        effective_refractive_index=1.0,
        diameter_um=20.0,
        dx_um=2.0,
        dz_um=1.0,
        edge_exclusion_px=5,
    )
    assert not result.valid
    assert result.invalid_reason == "coarse_peak_cnr_below_threshold"
