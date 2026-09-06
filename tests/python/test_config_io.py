import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import savemat

from svrecttail.config import load_config
from svrecttail.io import load_frame_maps, load_manifest


ROOT = Path(__file__).resolve().parents[2]


def test_pilot_config_freezes_protocol_critical_values() -> None:
    config = load_config(ROOT / "config" / "run_config.pilot.json")
    assert config.calibration.diameter_um == 128.0
    assert config.calibration.dx_um == 12.7
    assert config.calibration.dz_um == 6.7
    assert config.geometry.tail_length_um == 500.0
    assert config.background.skip_columns == 3
    assert config.background.strip_width_columns == 5


def test_surface_guided_config_freezes_confirmed_surface_geometry() -> None:
    config = load_config(
        ROOT / "config" / "run_config.surface_z176.pilot.json"
    )
    assert config.localization.mode == "fixed_surface_global_x"
    assert config.localization.fixed_surface_z_center_px == 176.0
    assert config.localization.surface_to_vessel_top_um == 200.0
    assert config.localization.effective_refractive_index == 1.12


def test_manifest_template_has_complete_15_frame_pilot_grid() -> None:
    table = load_manifest(ROOT / "data" / "manifest_template.csv", require_complete=False)
    assert len(table) == 15
    assert sorted(table["flow_speed_mm_s"].unique().tolist()) == [1, 3, 5, 7, 10]
    assert sorted(table["position_label"].unique().tolist()) == ["front", "middle", "rear"]
    assert table.groupby("scan_id").size().eq(3).all()


def test_manifest_allows_unknown_identity_and_acquisition_metadata(tmp_path: Path) -> None:
    (tmp_path / "frame.mat").touch()
    table = pd.DataFrame(
        [
            {
                "scan_id": "scan",
                "source_file": "frame.mat",
                "vessel_id": pd.NA,
                "phantom_id": pd.NA,
                "session_id": pd.NA,
                "diameter_um": 128,
                "flow_speed_mm_s": 1,
                "position_label": "front",
                "dx_um": 12.7,
                "dz_um": 6.7,
                "bscan_index": 0,
                "slow_axis_position_um": pd.NA,
                "temporal_repeat_id": pd.NA,
                "temporal_repeat_count": pd.NA,
                "scan_time_interval_s": pd.NA,
                "acquisition_order": pd.NA,
                "reconstruction_version": "commit",
                "x_anchor_center_px": 10,
                "z_anchor_center_px": 20,
                "geometry_source": "mentor_tracking",
                "background_excluded_side": pd.NA,
                "background_exclusion_reason": pd.NA,
            }
        ]
    )
    path = tmp_path / "manifest.csv"
    table.to_csv(path, index=False)
    loaded = load_manifest(path, require_complete=True)
    assert pd.isna(loaded.loc[0, "vessel_id"])


def test_surface_guided_manifest_allows_empty_legacy_anchor_columns(
    tmp_path: Path,
) -> None:
    (tmp_path / "frame.mat").touch()
    table = pd.DataFrame(
        [
            {
                "scan_id": "scan",
                "source_file": "frame.mat",
                "vessel_id": pd.NA,
                "phantom_id": pd.NA,
                "session_id": pd.NA,
                "diameter_um": 128,
                "flow_speed_mm_s": 1,
                "position_label": "front",
                "dx_um": 12.7,
                "dz_um": 6.7,
                "bscan_index": 0,
                "slow_axis_position_um": pd.NA,
                "temporal_repeat_id": pd.NA,
                "temporal_repeat_count": pd.NA,
                "scan_time_interval_s": pd.NA,
                "acquisition_order": pd.NA,
                "reconstruction_version": "commit",
                "x_anchor_center_px": pd.NA,
                "z_anchor_center_px": pd.NA,
                "geometry_source": "fixed_surface_global_x",
                "background_excluded_side": pd.NA,
                "background_exclusion_reason": pd.NA,
            }
        ]
    )
    path = tmp_path / "manifest.csv"
    table.to_csv(path, index=False)
    loaded = load_manifest(
        path,
        require_complete=True,
        require_localization_anchors=False,
    )
    assert pd.isna(loaded.loc[0, "x_anchor_center_px"])
    assert pd.isna(loaded.loc[0, "z_anchor_center_px"])


def test_npz_loader_preserves_linear_maps(tmp_path: Path) -> None:
    sv = np.arange(20, dtype=float).reshape(4, 5)
    omag = sv + 1.0
    structure = sv + 2.0
    path = tmp_path / "maps.npz"
    cv2 = sv / (structure**2 + np.finfo(float).eps)
    np.savez(path, sv_raw=sv, omag_raw=omag, stru_amp=structure, sv_cv2=cv2)
    maps = load_frame_maps(path)
    np.testing.assert_array_equal(maps.sv_raw, sv)
    np.testing.assert_array_equal(maps.omag_raw, omag)
    np.testing.assert_array_equal(maps.stru_amp, structure)
    np.testing.assert_array_equal(maps.sv_cv2, cv2)


def test_classic_mat_result_struct_loader(tmp_path: Path) -> None:
    sv = np.arange(20, dtype=float).reshape(4, 5)
    path = tmp_path / "maps.mat"
    savemat(
        path,
        {
            "result": {
                "sv_raw": sv,
                "omag_raw": sv + 1.0,
                "stru_amp": sv + 2.0,
                "sv_norm": sv / ((sv + 2.0) ** 2 + np.finfo(float).eps),
                "x_idx": 4,
                "z_top": 7,
            }
        },
    )
    maps = load_frame_maps(path)
    np.testing.assert_array_equal(maps.sv_raw, sv)
    assert maps.sv_cv2 is not None
    assert maps.metadata["stored_anchor_index_base"] == 1


def test_matlab_export_metadata_is_preserved(tmp_path: Path) -> None:
    sv = np.arange(20, dtype=float).reshape(4, 5)
    path = tmp_path / "export.mat"
    savemat(
        path,
        {
            "sv_raw": sv,
            "omag_raw": sv + 1.0,
            "metadata": {
                "source_file": "scan.oct",
                "bscan_index_matlab_1based": 3,
                "formal_signal_definition": "var(abs(E), 1, 3)",
                "dimension_order": "depth x A-line",
                "reconstruction": {"registration_Nsub": 10},
            },
        },
    )
    maps = load_frame_maps(path)
    assert maps.metadata["bscan_index_matlab_1based"] == 3
    assert maps.metadata["formal_signal_definition"] == "var(abs(E), 1, 3)"
    assert maps.metadata["export_metadata"]["reconstruction"]["registration_Nsub"] == 10


def test_config_rejects_normalized_formal_input(tmp_path: Path) -> None:
    source = json.loads((ROOT / "config" / "run_config.pilot.json").read_text(encoding="utf-8"))
    source["signal"]["formal_input"] = "sv_norm"
    path = tmp_path / "bad.json"
    path.write_text(json.dumps(source), encoding="utf-8")
    try:
        load_config(path)
    except ValueError as error:
        assert "sv_raw" in str(error)
    else:
        raise AssertionError("normalized input was accepted")


def test_config_rejects_fractional_integer_parameter(tmp_path: Path) -> None:
    source = json.loads((ROOT / "config" / "run_config.pilot.json").read_text(encoding="utf-8"))
    source["detection"]["bin_rows"] = 5.5
    path = tmp_path / "fractional.json"
    path.write_text(json.dumps(source), encoding="utf-8")
    try:
        load_config(path)
    except ValueError as error:
        assert "bin_rows" in str(error)
    else:
        raise AssertionError("fractional bin_rows was accepted")
