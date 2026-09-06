import numpy as np
import pytest

from svrecttail.mentor import tracking_core


def test_persistent_edge_ignores_unlinked_bright_rows() -> None:
    nz = 100
    profile = np.zeros(nz, dtype=float)
    profile[20:22] = 8.0
    profile[40:75] = 10.0
    background = 0.2 * np.sin(np.arange(nz, dtype=float))
    support = np.zeros(nz, dtype=float)
    support[20:22] = 0.8
    support[40:75] = 0.8
    cfg = tracking_core.merge_tracking_config(
        {
            "axial_um_per_px": 1.0,
            "z_search_px": [10, 90],
            "upper_edge_alphas": [0.20],
            "primary_alpha": 0.20,
            "upper_edge_noise_multiplier": 4.0,
        }
    )

    result = tracking_core.persistent_upper_edge_candidate(
        profile,
        background,
        support,
        peak_hint=50.0,
        diameter_um=20.0,
        alpha=0.20,
        cfg=cfg,
    )

    assert result.valid
    assert result.z_upper_px == 40.0
    assert result.core_start_px == 40.0
    assert result.high_component_width_px >= 20
    assert result.top_contrast_snr > 0


def test_row_support_uses_physical_lateral_scale() -> None:
    rng = np.random.default_rng(4)
    volume = rng.normal(20.0, 0.5, size=(2, 80, 60)).astype(np.float32)
    volume[:, 30:55, 27:34] += 12.0
    cfg = tracking_core.merge_tracking_config(
        {
            "lateral_um_per_px": 5.0,
            "upper_edge_support_width_fraction_of_lateral_diameter": 1.0,
        }
    )

    support = tracking_core.extract_persistent_row_support(
        volume,
        np.array([30.0, 30.0]),
        diameter_um=30.0,
        cfg=cfg,
    )

    assert support.shape == (2, 80)
    assert float(np.median(support[:, 36:49])) > 0.70
    assert float(np.median(support[:, 5:20])) < 0.20


def test_upper_edge_method_is_versioned_and_validated() -> None:
    assert (
        tracking_core.merge_tracking_config(None)["upper_edge_method"]
        == tracking_core.UPPER_EDGE_LEGACY
    )
    persistent = tracking_core.merge_tracking_config(
        {"upper_edge_method": tracking_core.UPPER_EDGE_PERSISTENT}
    )
    assert (
        persistent["upper_edge_method"]
        == tracking_core.UPPER_EDGE_PERSISTENT
    )
    with pytest.raises(tracking_core.TrackingError, match="upper_edge_method"):
        tracking_core.merge_tracking_config({"upper_edge_method": "unknown"})
