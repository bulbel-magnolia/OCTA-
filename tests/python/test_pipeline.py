import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.io import loadmat, savemat

from svrecttail.geometry import VesselGeometry, ellipse_weights
from svrecttail.pipeline import run_batch


ROOT = Path(__file__).resolve().parents[2]


def test_manifest_batch_writes_complete_audit_bundle(tmp_path: Path) -> None:
    shape = (200, 100)
    omag = np.ones(shape, dtype=float)
    omag[30:49, 45:55] = 11.0
    sv = np.repeat((5.0 + 0.001 * np.arange(shape[0]))[:, None], shape[1], axis=1)
    geometry = VesselGeometry.from_inclusive_centres(
        x_left_center_px=45,
        x_right_center_px=54,
        z_top_center_px=30,
        diameter_um=128.0,
        dx_um=12.7,
        dz_um=6.7,
    )
    source_weights = ellipse_weights(shape, geometry)
    sv[source_weights > 0] += 3.0
    z_stop = geometry.z_bottom_edge_px + 500.0 / geometry.dz_um
    rows = (np.arange(shape[0]) > geometry.z_bottom_edge_px) & (
        np.arange(shape[0]) <= z_stop
    )
    sv[np.ix_(rows, np.arange(45, 55))] += 1.0
    map_path = tmp_path / "frame.mat"
    savemat(
        map_path,
        {
            "sv_raw": sv,
            "omag_raw": omag,
            "metadata": {
                "bscan_index_matlab_1based": 1,
                "formal_signal_definition": "var(abs(E), 1, 3)",
                "dimension_order": "depth x A-line",
            },
        },
    )

    manifest = pd.DataFrame(
        [
            {
                "scan_id": "synthetic_001",
                "source_file": str(map_path),
                "vessel_id": "v1",
                "phantom_id": "p1",
                "session_id": "s1",
                "diameter_um": 128.0,
                "flow_speed_mm_s": 3.0,
                "position_label": "middle",
                "dx_um": 12.7,
                "dz_um": 6.7,
                "bscan_index": 0,
                "slow_axis_position_um": 0.0,
                "temporal_repeat_id": "r1",
                "temporal_repeat_count": 4,
                "scan_time_interval_s": 0.02,
                "acquisition_order": 1,
                "reconstruction_version": "synthetic",
                "x_anchor_center_px": 50.0,
                "z_anchor_center_px": 30.0,
                "geometry_source": "mentor_tracking",
                "background_excluded_side": pd.NA,
                "background_exclusion_reason": pd.NA,
                "notes": "test",
            }
        ]
    )
    bad_map_path = tmp_path / "bad_frame.npz"
    np.savez(
        bad_map_path,
        sv_raw=np.ones(shape, dtype=float),
        omag_raw=np.full(shape, np.nan, dtype=float),
    )
    bad_row = manifest.iloc[0].copy()
    bad_row["scan_id"] = "synthetic_bad"
    bad_row["source_file"] = str(bad_map_path)
    manifest = pd.concat([manifest, bad_row.to_frame().T], ignore_index=True)
    manifest_path = tmp_path / "manifest.csv"
    manifest.to_csv(manifest_path, index=False)
    blank_path = tmp_path / "blank_profiles.npz"
    blank_offsets = np.array([-0.2, -0.1, 0.0, 0.1, 0.2])[:, None]
    np.savez(blank_path, all=np.repeat(blank_offsets, 75, axis=1))
    output = tmp_path / "run"
    result = run_batch(
        config_path=ROOT / "config" / "run_config.pilot.json",
        manifest_path=manifest_path,
        output_dir=output,
        blank_profiles_path=blank_path,
    )
    assert result.frame_count == 2
    assert result.valid_frame_count == 1
    expected = {
        "run_config.json",
        "manifest.csv",
        "frame_results.csv",
        "scan_summary.csv",
        "profiles.csv",
        "localization.csv",
        "sensitivity_results.csv",
        "detection_results.csv",
        "detection_bins.csv",
        "run_complete.json",
        "arrays/synthetic_001_Bscan000_Tr1_SV_Rectangle_v1_pilot.mat",
        "arrays/synthetic_001_Bscan000_Tr1_SV_Rectangle_v1_pilot_detection.mat",
        "profiles/synthetic_001_Bscan000_Tr1_SV_Rectangle_v1_pilot.csv",
        "qc/synthetic_001_Bscan000_Tr1_SV_Rectangle_v1_pilot_QC01-03.png",
        "qc/synthetic_001_Bscan000_Tr1_SV_Rectangle_v1_pilot_QC03_detection.png",
        "logs/run_complete.json",
        "logs/manual_adjustments.csv",
    }
    assert all((output / name).exists() for name in expected)
    frame_table = pd.read_csv(output / "frame_results.csv")
    frame = frame_table.loc[frame_table["scan_id"] == "synthetic_001"].iloc[0]
    assert bool(frame["valid"])
    assert frame["input_metadata_qc_status"] == "verified_export_metadata"
    assert np.isfinite(frame["ratio_tail_to_vessel"])
    assert frame["source_area_um2"] > 0
    assert np.isclose(frame["tail_area_um2"], frame["requested_tail_area_um2"])
    assert bool(frame["detection_qc_valid"])
    assert bool(frame["tail_detected"])
    assert np.isfinite(frame["detectable_length_um"])
    assert frame["detectable_length_um"] <= 500.0
    invalid_frame = frame_table.loc[frame_table["scan_id"] == "synthetic_bad"].iloc[0]
    assert not bool(invalid_frame["valid"])
    assert str(invalid_frame["invalid_reason"]).startswith("frame_processing_error:")
    archive = loadmat(
        output
        / "arrays"
        / "synthetic_001_Bscan000_Tr1_SV_Rectangle_v1_pilot.mat"
    )
    required_arrays = {
        "sv_raw",
        "omag_raw",
        "corrected_sv",
        "source_ellipse_fraction",
        "source_mask",
        "tail_rectangle_fraction",
        "tail_rectangle_mask",
        "background_left_pixels",
        "background_right_pixels",
        "background_left_mask",
        "background_right_mask",
    }
    assert required_arrays <= set(archive)
    q_tail_recomputed = (
        archive["corrected_sv"] * archive["tail_rectangle_fraction"]
    ).sum() * 12.7 * 6.7
    assert np.isclose(q_tail_recomputed, frame["q_tail"])
    detection_archive = loadmat(
        output
        / "arrays"
        / "synthetic_001_Bscan000_Tr1_SV_Rectangle_v1_pilot_detection.mat"
    )
    assert detection_archive["matched_blank_T_profiles"].shape == (5, 75)
    assert detection_archive["signal_bin_mean"].size == 14
    detection = pd.read_csv(output / "detection_results.csv").iloc[0]
    assert bool(detection["detected"])
    assert bool(detection["right_censored"])
    assert len(pd.read_csv(output / "detection_bins.csv")) == 14
    summary_table = pd.read_csv(output / "scan_summary.csv")
    summary = summary_table.loc[summary_table["scan_id"] == "synthetic_001"].iloc[0]
    assert summary["frame_count"] == 1
    assert summary["valid_frame_count"] == 1
    run_metadata = json.loads((output / "run_complete.json").read_text(encoding="utf-8"))
    assert run_metadata["frame_processing_error_count"] == 1
    assert run_metadata["manual_adjustment_count"] == 0
    assert pd.read_csv(output / "logs" / "manual_adjustments.csv").empty
    with pytest.raises(FileExistsError):
        run_batch(
            config_path=ROOT / "config" / "run_config.pilot.json",
            manifest_path=manifest_path,
            output_dir=output,
            blank_profiles_path=blank_path,
        )


def test_mentor_tracking_batch_uses_tracking_table_without_surface(
    tmp_path: Path,
) -> None:
    shape = (160, 100)
    omag = np.ones(shape, dtype=float)
    sv = np.ones(shape, dtype=float)
    tracked_geometry = VesselGeometry(
        x_left_edge_px=45.25,
        x_right_edge_px=55.25,
        z_top_edge_px=29.5,
        diameter_um=128.0,
        dx_um=12.7,
        dz_um=6.7,
    )
    sv[ellipse_weights(shape, tracked_geometry) > 0] += 3.0
    map_path = tmp_path / "frame.mat"
    savemat(
        map_path,
        {
            "sv_raw": sv,
            "omag_raw": omag,
            "metadata": {
                "bscan_index_matlab_1based": 1,
                "formal_signal_definition": "var(abs(E), 1, 3)",
                "dimension_order": "depth x A-line",
            },
        },
    )
    tracking_path = tmp_path / "tracking.csv"
    pd.DataFrame(
        [
            {
                "scan_id": "tracked_scan",
                "frame_index": 0,
                "alpha": 0.15,
                "x_center_px": 51.0,
                "z_upper_px": 30.0,
                "new_tracking_class": "high_confidence",
                "x_path_confidence_class": "strong_evidence",
                "z_edge_confidence_class": "strong_evidence",
                "valid_local_body": True,
                "x1_local_geometry_px": 50.5,
                "x2_robust_centroid_px": 50.25,
                "x4_centroid_isolated_jump_corrected_px": 50.25,
                "x4_jump_corrected": False,
                "x1_fallback": False,
                "local_body_run_width_px": 10,
                "expected_lateral_width_px": 10,
                "local_body_background": 1.0,
                "local_body_sigma": 0.1,
                "local_body_peak_cnr": 8.0,
                "local_body_axial_completeness": 0.9,
                "assessability_score": 0.8,
                "vessel_presence_prediction": "assessable",
                "peak_snr": 6.0,
            }
        ]
    ).to_csv(tracking_path, index=False)
    manifest = pd.DataFrame(
        [
            {
                "scan_id": "tracked_scan",
                "source_file": str(map_path),
                "tracking_file": str(tracking_path),
                "vessel_id": "v1",
                "phantom_id": "p1",
                "session_id": "s1",
                "diameter_um": 128.0,
                "flow_speed_mm_s": 3.0,
                "position_label": "front",
                "dx_um": 12.7,
                "dz_um": 6.7,
                "bscan_index": 0,
                "slow_axis_position_um": 0.0,
                "temporal_repeat_id": "r1",
                "temporal_repeat_count": 4,
                "scan_time_interval_s": 0.02,
                "acquisition_order": 1,
                "reconstruction_version": "synthetic",
                "x_anchor_center_px": pd.NA,
                "z_anchor_center_px": pd.NA,
                "geometry_source": "mentor_tracking",
                "background_excluded_side": pd.NA,
                "background_exclusion_reason": pd.NA,
                "notes": "mentor pipeline integration",
            }
        ]
    )
    manifest_path = tmp_path / "manifest.csv"
    manifest.to_csv(manifest_path, index=False)
    output = tmp_path / "run"
    result = run_batch(
        config_path=ROOT / "config" / "run_config.mentor_tracking.pilot.json",
        manifest_path=manifest_path,
        output_dir=output,
    )
    assert result.valid_frame_count == 1
    localization = pd.read_csv(output / "localization.csv").iloc[0]
    assert localization["coarse_method"] == "mentor_full_volume_slow_axis_viterbi"
    assert localization["mentor_frame_index"] == 0
    assert localization["mentor_vessel_presence_prediction"] == "assessable"
    assert localization["x_left_edge_px"] == 45.25
    assert localization["x_right_edge_px"] == 55.25
    assert localization["z_top_edge_px"] == 29.5
    assert "coarse_surface_z_center_px" not in localization.index
