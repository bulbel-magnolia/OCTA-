"""Manifest-driven formal batch analysis."""

from __future__ import annotations

import hashlib
import json
import platform
import re
import shutil
import subprocess
import sys
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import h5py
import matplotlib
import numpy as np
import pandas as pd
import scipy
from scipy.io import savemat

from .config import RunConfig, load_config
from .detection import DetectionResult, detect_tail_extent
from .io import FrameMaps, load_frame_maps, load_manifest, resolve_source_path
from .localization import LocalizationResult, localize_geometry, shifted_geometry
from .qc import save_detection_qc, save_qc_figure
from .quantification import QuantificationResult, quantify_frame


@dataclass(frozen=True)
class BatchRunResult:
    output_dir: Path
    frame_count: int
    valid_frame_count: int
    frame_results_path: Path
    profiles_path: Path


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _safe_name(value: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("._")
    return safe or "frame"


def _frame_id(row: pd.Series) -> str:
    base = f"{row['scan_id']}_Bscan{int(row['bscan_index']):03d}"
    repeat = row.get("temporal_repeat_id")
    if repeat is not None and not pd.isna(repeat) and str(repeat).strip():
        base += f"_T{_safe_name(str(repeat))}"
    return base


def _frame_stem(row: pd.Series, config: RunConfig) -> str:
    return f"{_safe_name(_frame_id(row))}_{_safe_name(config.schema_version)}"


def _git_state() -> dict[str, Any]:
    repo_root = Path(__file__).resolve().parents[2]
    try:
        commit = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=repo_root,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
        porcelain = subprocess.run(
            ["git", "status", "--porcelain"],
            cwd=repo_root,
            check=True,
            capture_output=True,
            text=True,
        ).stdout
        return {"commit": commit, "worktree_dirty": bool(porcelain.strip())}
    except (OSError, subprocess.CalledProcessError):
        return {"commit": None, "worktree_dirty": None}


def _scalar_manifest_values(row: pd.Series) -> dict[str, Any]:
    output: dict[str, Any] = {}
    for key, value in row.items():
        if pd.isna(value):
            output[key] = None
        elif isinstance(value, np.generic):
            output[key] = value.item()
        else:
            output[key] = value
    return output


def _json_compatible(value: Any) -> Any:
    """Convert NumPy and non-finite metadata to strict JSON values."""

    if isinstance(value, dict):
        return {str(key): _json_compatible(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_compatible(item) for item in value]
    if isinstance(value, np.ndarray):
        return _json_compatible(value.tolist())
    if isinstance(value, np.generic):
        return _json_compatible(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


def _localization_record(localization: LocalizationResult) -> dict[str, Any]:
    lateral = asdict(localization.lateral)
    top = asdict(localization.top_edge)
    return {
        "localization_image": "omag_raw",
        "lateral_method": "mentor_X1_local_q75_continuous_run",
        "lateral_threshold_sigma": 1.5,
        "top_method": "X1_centreline_three_column_sustained_edge",
        "coordinate_convention": "zero_based_pixel_centres",
        **{f"lateral_{key}": value for key, value in lateral.items()},
        **{f"top_{key}": value for key, value in top.items()},
        "source_qc_valid": localization.source_qc_valid,
    }


def _manifest_completeness_record(row: pd.Series) -> dict[str, bool]:
    def present(name: str) -> bool:
        value = row.get(name)
        return value is not None and not pd.isna(value) and bool(str(value).strip())

    return {
        "manifest_identity_complete": all(
            present(name) for name in ("vessel_id", "session_id")
        ),
        "manifest_phantom_id_known": present("phantom_id"),
        "manifest_acquisition_metadata_complete": all(
            present(name)
            for name in (
                "acquisition_order",
                "temporal_repeat_count",
                "scan_time_interval_s",
            )
        ),
    }


def _validate_map_metadata(maps: FrameMaps, row: pd.Series) -> dict[str, Any]:
    metadata = maps.metadata
    exported_index = metadata.get("bscan_index_matlab_1based")
    if exported_index is not None:
        try:
            exported_index_value = float(exported_index)
        except (TypeError, ValueError) as error:
            raise ValueError("MATLAB export bscan index is not numeric") from error
        if (
            not np.isfinite(exported_index_value)
            or exported_index_value != round(exported_index_value)
            or int(exported_index_value) != int(row["bscan_index"]) + 1
        ):
            raise ValueError(
                "manifest bscan_index disagrees with MATLAB export metadata"
            )
    formal_definition = metadata.get("formal_signal_definition")
    if formal_definition is not None and formal_definition != "var(abs(E), 1, 3)":
        raise ValueError("interim map metadata declares a non-formal SV definition")
    dimension_order = metadata.get("dimension_order")
    if dimension_order is not None and dimension_order != "depth x A-line":
        raise ValueError("interim map metadata dimension_order must be depth x A-line")
    verified_fields = (exported_index, formal_definition, dimension_order)
    return {
        "input_map_format": metadata.get("format"),
        "input_metadata_qc_status": (
            "verified_export_metadata"
            if all(value is not None for value in verified_fields)
            else "legacy_or_unspecified"
        ),
        "input_bscan_index_matlab_1based": exported_index,
        "input_declared_sv_definition": (
            formal_definition
            if formal_definition is not None
            else metadata.get("SV_definition")
        ),
        "input_dimension_order": dimension_order,
    }


def _frame_processing_error_record(
    row: pd.Series,
    *,
    frame_id: str,
    source_path: Path,
    config: RunConfig,
    error: Exception,
) -> dict[str, Any]:
    reason = f"frame_processing_error:{type(error).__name__}:{error}"
    return {
        **_scalar_manifest_values(row),
        "frame_id": frame_id,
        "resolved_source_file": str(source_path),
        **_manifest_completeness_record(row),
        "localization_image": "omag_raw",
        "lateral_method": "mentor_X1_local_q75_continuous_run",
        "lateral_threshold_sigma": 1.5,
        "top_method": "X1_centreline_three_column_sustained_edge",
        "coordinate_convention": "zero_based_pixel_centres",
        "geometry_qc_valid": False,
        "source_qc_valid": False,
        "background_qc_valid": False,
        "window_qc_valid": False,
        "background_modes": "unavailable",
        "input_metadata_qc_status": "unavailable",
        "q_vessel": np.nan,
        "q_tail": np.nan,
        "ratio_tail_to_vessel": np.nan,
        "q_tail_direct_check": np.nan,
        "source_area_um2": np.nan,
        "tail_area_um2": np.nan,
        "requested_tail_area_um2": np.nan,
        "source_mean": np.nan,
        "tail_window_complete": False,
        "source_window_complete": False,
        "background_complete": False,
        "valid": False,
        "invalid_reason": reason,
        "tail_gap_um": config.geometry.tail_gap_um,
        "requested_tail_length_um": config.geometry.tail_length_um,
    }


def _save_frame_archive(
    path: Path,
    *,
    maps: FrameMaps,
    result: QuantificationResult,
    localization: LocalizationResult,
    row: pd.Series,
    config: RunConfig,
) -> None:
    background = result.background
    tail_rectangle = np.multiply.outer(
        result.tail_z_weights, result.vessel_x_weights
    )
    left_background_mask = np.zeros(maps.sv_raw.shape, dtype=bool)
    right_background_mask = np.zeros(maps.sv_raw.shape, dtype=bool)
    left_background_mask[:, background.left_columns] = True
    right_background_mask[:, background.right_columns] = True
    metadata = {
        "schema_version": config.schema_version,
        "formal_signal": "var(abs(E), 1, 3)",
        "coordinate_convention": (
            "zero-based pixel centres; image edges -0.5 and n-0.5"
        ),
        "manifest_json": json.dumps(
            _json_compatible(_scalar_manifest_values(row)),
            ensure_ascii=False,
            allow_nan=False,
        ),
        "geometry_json": json.dumps(
            _json_compatible(result.summary()), ensure_ascii=False, allow_nan=False
        ),
        "localization_json": json.dumps(
            _json_compatible(_localization_record(localization)),
            ensure_ascii=False,
            allow_nan=False,
        ),
        "input_map_metadata_json": json.dumps(
            _json_compatible(maps.metadata), ensure_ascii=False, allow_nan=False
        ),
    }
    content = {
        "sv_raw": maps.sv_raw,
        "sv_cv2": np.empty((0, 0)) if maps.sv_cv2 is None else maps.sv_cv2,
        "omag_raw": maps.omag_raw,
        "stru_amp": np.empty((0, 0)) if maps.stru_amp is None else maps.stru_amp,
        "corrected_sv": result.corrected_sv,
        "vessel_profile_V": result.vessel_profile,
        "background_B": background.combined,
        "background_B_left": background.left,
        "background_B_right": background.right,
        "background_standard_deviation": background.standard_deviation,
        "background_median_diagnostic": background.median,
        "background_scaled_mad_diagnostic": background.scaled_mad,
        "background_left_columns_0based": background.left_columns,
        "background_right_columns_0based": background.right_columns,
        "background_left_pixels": maps.sv_raw[:, background.left_columns],
        "background_right_pixels": maps.sv_raw[:, background.right_columns],
        "background_left_mask": left_background_mask,
        "background_right_mask": right_background_mask,
        "background_left_valid_count": background.left_valid_count,
        "background_right_valid_count": background.right_valid_count,
        "background_row_mode": background.row_mode,
        "tail_contrast_profile_T": result.tail_contrast_profile,
        "tail_linear_density_P": result.tail_linear_density,
        "vessel_x_weights": result.vessel_x_weights,
        "source_ellipse_fraction": result.vessel_ellipse_weights,
        "source_mask": result.vessel_ellipse_weights > 0,
        "tail_z_weights": result.tail_z_weights,
        "tail_rectangle_fraction": tail_rectangle,
        "tail_rectangle_mask": tail_rectangle > 0,
        "x_pixel_centres_0based": np.arange(maps.sv_raw.shape[1]),
        "z_pixel_centres_0based": np.arange(maps.sv_raw.shape[0]),
        "metadata": metadata,
    }
    savemat(path, content, do_compression=True, long_field_names=True)


def _save_detection_archive(
    path: Path,
    *,
    signal_profile: np.ndarray,
    blank_profiles: np.ndarray,
    detection: DetectionResult,
    blank_profile_key: str,
    blank_source_path: Path,
    blank_source_sha256: str,
    config: RunConfig,
) -> None:
    savemat(
        path,
        {
            "signal_T_profile": np.asarray(signal_profile, dtype=np.float64),
            "matched_blank_T_profiles": np.asarray(blank_profiles, dtype=np.float64),
            "signal_bin_mean": detection.signal_bin_mean,
            "blank_median": detection.blank_median,
            "blank_scaled_mad": detection.blank_scaled_mad,
            "threshold": detection.threshold,
            "blank_count": detection.blank_count,
            "exceeds_threshold": detection.exceeds_threshold,
            "sustained_detection": detection.sustained_detection,
            "bin_rows": detection.bin_rows,
            "dz_um": detection.dz_um,
            "threshold_mad_multiplier": config.detection.threshold_mad_multiplier,
            "minimum_consecutive_bins": config.detection.minimum_consecutive_bins,
            "minimum_blank_samples_per_bin": (
                config.detection.minimum_blank_samples_per_bin
            ),
            "blank_profile_key": blank_profile_key,
            "blank_source_path": str(blank_source_path),
            "blank_source_sha256": blank_source_sha256,
            "status": detection.invalid_reason,
            "detectable_length_um": detection.detectable_length_um,
            "right_censored": detection.right_censored,
        },
        do_compression=True,
        long_field_names=True,
    )


def _quantify(
    maps_sv: np.ndarray,
    localization: LocalizationResult,
    config: RunConfig,
    *,
    geometry=None,
    background_skip_columns: int | None = None,
    background_excluded_side: str | None = None,
) -> QuantificationResult:
    return quantify_frame(
        maps_sv,
        localization.geometry if geometry is None else geometry,
        tail_gap_um=config.geometry.tail_gap_um,
        tail_length_um=config.geometry.tail_length_um,
        background_skip_columns=(
            config.background.skip_columns
            if background_skip_columns is None
            else background_skip_columns
        ),
        background_strip_width_columns=config.background.strip_width_columns,
        background_excluded_side=background_excluded_side,
        ellipse_supersample=config.geometry.ellipse_supersample,
        source_qc_valid=localization.source_qc_valid,
    )


def _sensitivity_records(
    scan_id: str,
    frame_id: str,
    sv_raw: np.ndarray,
    localization: LocalizationResult,
    config: RunConfig,
    background_excluded_side: str | None,
) -> list[dict[str, Any]]:
    variants = [
        ("x_minus_1px", shifted_geometry(localization.geometry, x_shift_px=-1), None),
        ("x_plus_1px", shifted_geometry(localization.geometry, x_shift_px=1), None),
        ("z_minus_1px", shifted_geometry(localization.geometry, z_shift_px=-1), None),
        ("z_plus_1px", shifted_geometry(localization.geometry, z_shift_px=1), None),
        (
            "background_skip_plus_2px",
            localization.geometry,
            config.background.skip_columns + 2,
        ),
    ]
    records: list[dict[str, Any]] = []
    for variant, geometry, background_skip in variants:
        result = _quantify(
            sv_raw,
            localization,
            config,
            geometry=geometry,
            background_skip_columns=background_skip,
            background_excluded_side=background_excluded_side,
        )
        records.append(
            {
                "scan_id": scan_id,
                "frame_id": frame_id,
                "variant": variant,
                **result.summary(),
            }
        )
    return records


def _tail_centre_profile(result: QuantificationResult) -> np.ndarray:
    """Select T(r) rows without letting complete bins extend beyond L."""

    centres = np.arange(result.tail_linear_density.size, dtype=np.float64)
    start = result.geometry.z_bottom_edge_px + result.tail_gap_um / result.geometry.dz_um
    stop = start + result.requested_tail_length_um / result.geometry.dz_um
    selected = result.tail_contrast_profile[(centres > start) & (centres <= stop)]
    maximum_complete_rows = int(
        np.floor(result.requested_tail_length_um / result.geometry.dz_um)
    )
    return selected[:maximum_complete_rows]


def _check_calibration(row: pd.Series, config: RunConfig) -> None:
    expected = {
        "diameter_um": config.calibration.diameter_um,
        "dx_um": config.calibration.dx_um,
        "dz_um": config.calibration.dz_um,
    }
    mismatches = [
        name
        for name, value in expected.items()
        if not np.isclose(float(row[name]), value, rtol=0.0, atol=1e-9)
    ]
    if mismatches:
        raise ValueError(
            f"scan {row['scan_id']} differs from frozen config: {', '.join(mismatches)}"
        )


def _write_json(path: Path, content: dict[str, Any]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(content, handle, ensure_ascii=False, indent=2, allow_nan=False)
        handle.write("\n")


def _median_iqr(values: pd.Series) -> tuple[float, float]:
    finite = pd.to_numeric(values, errors="coerce").dropna()
    if finite.empty:
        return float("nan"), float("nan")
    quartiles = finite.quantile([0.25, 0.75])
    return float(finite.median()), float(quartiles.loc[0.75] - quartiles.loc[0.25])


def _build_scan_summary(
    frame_table: pd.DataFrame, detection_table: pd.DataFrame
) -> pd.DataFrame:
    records: list[dict[str, Any]] = []
    for scan_id, group in frame_table.groupby("scan_id", sort=False):
        valid_group = group[group["valid"].astype(bool)]
        record: dict[str, Any] = {
            "scan_id": scan_id,
            "frame_count": len(group),
            "valid_frame_count": int(group["valid"].sum()),
            "invalid_frame_count": int((~group["valid"].astype(bool)).sum()),
            "geometry_qc_failure_count": int(
                (~group["geometry_qc_valid"].astype(bool)).sum()
            ),
            "source_qc_failure_count": int(
                (~group["source_qc_valid"].astype(bool)).sum()
            ),
            "background_qc_failure_count": int(
                (~group["background_qc_valid"].astype(bool)).sum()
            ),
            "window_qc_failure_count": int(
                (~group["window_qc_valid"].astype(bool)).sum()
            ),
        }
        speeds = sorted(pd.to_numeric(group["flow_speed_mm_s"]).unique())
        record["flow_speed_mm_s"] = speeds[0] if len(speeds) == 1 else ";".join(
            str(value) for value in speeds
        )
        for column in ("q_vessel", "q_tail", "ratio_tail_to_vessel", "source_mean"):
            median, iqr = _median_iqr(valid_group[column])
            record[f"{column}_median"] = median
            record[f"{column}_iqr"] = iqr
        scan_detection = detection_table[detection_table["scan_id"] == scan_id]
        record["detection_detected_count"] = int(scan_detection["detected"].sum())
        record["detection_right_censored_count"] = int(
            scan_detection["right_censored"].sum()
        )
        record["detection_not_detected_count"] = int(
            (scan_detection["status"] == "not_detected").sum()
        )
        record["detection_not_evaluable_count"] = int(
            (~scan_detection["status"].isin(["ok", "right_censored", "not_detected"])).sum()
        )
        records.append(record)
    return pd.DataFrame.from_records(records)


def run_batch(
    *,
    config_path: str | Path,
    manifest_path: str | Path,
    output_dir: str | Path,
    blank_profiles_path: str | Path | None = None,
) -> BatchRunResult:
    """Execute a frozen one-row-per-frame batch and write auditable products."""

    config_source = Path(config_path).resolve()
    manifest_source = Path(manifest_path).resolve()
    output = Path(output_dir).resolve()
    config = load_config(config_source)
    manifest = load_manifest(manifest_source, require_complete=True)
    for _, calibration_row in manifest.iterrows():
        _check_calibration(calibration_row, config)
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"output directory is not empty: {output}")
    blank_data: dict[str, np.ndarray] | None = None
    blank_source: Path | None = None
    blank_source_sha256: str | None = None
    if blank_profiles_path is not None:
        blank_source = Path(blank_profiles_path).resolve()
        if blank_source.suffix.lower() != ".npz":
            raise ValueError("blank profiles input must be an NPZ file")
        blank_source_sha256 = _sha256(blank_source)
        with np.load(blank_source, allow_pickle=False) as blank_content:
            blank_data = {
                key: np.asarray(blank_content[key], dtype=np.float64)
                for key in blank_content.files
            }
    output.mkdir(parents=True, exist_ok=True)
    arrays_dir = output / "arrays"
    qc_dir = output / "qc"
    profiles_dir = output / "profiles"
    logs_dir = output / "logs"
    arrays_dir.mkdir(exist_ok=True)
    qc_dir.mkdir(exist_ok=True)
    profiles_dir.mkdir(exist_ok=True)
    logs_dir.mkdir(exist_ok=True)
    shutil.copy2(config_source, output / "run_config.json")
    shutil.copy2(manifest_source, output / "manifest.csv")

    frame_records: list[dict[str, Any]] = []
    profile_records: list[dict[str, Any]] = []
    localization_records: list[dict[str, Any]] = []
    sensitivity_records: list[dict[str, Any]] = []
    detection_records: list[dict[str, Any]] = []
    detection_bin_records: list[dict[str, Any]] = []
    manual_adjustment_records: list[dict[str, Any]] = []

    try:
        for _, row in manifest.iterrows():
            scan_id = str(row["scan_id"])
            frame_id = _frame_id(row)
            frame_stem = _frame_stem(row, config)
            excluded_value = row.get("background_excluded_side")
            background_excluded_side = (
                None
                if excluded_value is None or pd.isna(excluded_value)
                else str(excluded_value)
            )
            if background_excluded_side is not None:
                manual_adjustment_records.append(
                    {
                        "scan_id": scan_id,
                        "frame_id": frame_id,
                        "adjustment_type": "background_side_exclusion",
                        "value": background_excluded_side,
                        "reason": str(row["background_exclusion_reason"]),
                        "source": "manifest_predeclared",
                    }
                )
            source_path = resolve_source_path(manifest_source, str(row["source_file"]))
            try:
                maps = load_frame_maps(source_path)
                map_metadata_record = _validate_map_metadata(maps, row)
                localization = localize_geometry(
                    maps.omag_raw,
                    x_anchor_center_px=float(row["x_anchor_center_px"]),
                    z_anchor_center_px=float(row["z_anchor_center_px"]),
                    diameter_um=float(row["diameter_um"]),
                    dx_um=float(row["dx_um"]),
                    dz_um=float(row["dz_um"]),
                )
                result = _quantify(
                    maps.sv_raw,
                    localization,
                    config,
                    background_excluded_side=background_excluded_side,
                )
            except (OSError, ValueError) as error:
                invalid_record = _frame_processing_error_record(
                    row,
                    frame_id=frame_id,
                    source_path=source_path,
                    config=config,
                    error=error,
                )
                frame_records.append(invalid_record)
                localization_records.append(
                    {
                        "scan_id": scan_id,
                        "frame_id": frame_id,
                        "resolved_source_file": str(source_path),
                        "localization_valid": False,
                        "invalid_reason": invalid_record["invalid_reason"],
                    }
                )
                _write_json(
                    logs_dir / f"{frame_stem}_error.json",
                    {
                        "scan_id": scan_id,
                        "frame_id": frame_id,
                        "status": "invalid_frame",
                        "error_type": type(error).__name__,
                        "error": str(error),
                    },
                )
                if config.detection.enabled:
                    detection_records.append(
                        {
                            "scan_id": scan_id,
                            "frame_id": frame_id,
                            "status": "not_run_frame_invalid",
                            "detected": False,
                            "detectable_length_um": np.nan,
                            "right_censored": False,
                            "evaluable": False,
                            "not_evaluable_reason": invalid_record["invalid_reason"],
                            "blank_profile_key": None,
                        }
                    )
                continue
            common = _scalar_manifest_values(row)
            frame_records.append(
                {
                    **common,
                    "frame_id": frame_id,
                    "resolved_source_file": str(source_path),
                    **map_metadata_record,
                    **_manifest_completeness_record(row),
                    **_localization_record(localization),
                    "geometry_qc_valid": (
                        localization.lateral.valid and localization.top_edge.valid
                    ),
                    "background_qc_valid": result.background_complete,
                    "background_modes": ";".join(
                        sorted(set(result.background.row_mode.tolist()))
                    ),
                    "window_qc_valid": result.tail_window_complete,
                    **result.summary(),
                }
            )
            location_record = {
                "scan_id": scan_id,
                "frame_id": frame_id,
                "resolved_source_file": str(source_path),
                **_localization_record(localization),
                **asdict(localization.geometry),
                "z_bottom_edge_px": localization.geometry.z_bottom_edge_px,
            }
            localization_records.append(location_record)
            current_profiles = result.profile_records(scan_id)
            for profile_record in current_profiles:
                profile_record["frame_id"] = frame_id
                profile_record["bscan_index"] = int(row["bscan_index"])
                profile_record["flow_speed_mm_s"] = float(row["flow_speed_mm_s"])
                profile_records.append(profile_record)
            pd.DataFrame.from_records(current_profiles).to_csv(
                profiles_dir / f"{frame_stem}.csv", index=False
            )
            sensitivity_records.extend(
                _sensitivity_records(
                    scan_id,
                    frame_id,
                    maps.sv_raw,
                    localization,
                    config,
                    background_excluded_side,
                )
            )
            _save_frame_archive(
                arrays_dir / f"{frame_stem}.mat",
                maps=maps,
                result=result,
                localization=localization,
                row=row,
                config=config,
            )
            save_qc_figure(
                qc_dir / f"{frame_stem}_QC01-03.png",
                scan_id=frame_id,
                maps=maps,
                localization=localization,
                result=result,
            )
            detection_upstream_reasons: list[str] = []
            if not localization.source_qc_valid:
                detection_upstream_reasons.append("source_qc_failed")
            if not result.background_complete:
                detection_upstream_reasons.append("background_qc_failed")
            if not result.tail_window_complete:
                detection_upstream_reasons.append("tail_window_incomplete")
            if config.detection.enabled and detection_upstream_reasons:
                detection_records.append(
                    {
                        "scan_id": scan_id,
                        "frame_id": frame_id,
                        "status": "not_run_upstream_qc_failed",
                        "detected": False,
                        "detectable_length_um": np.nan,
                        "right_censored": False,
                        "evaluable": False,
                        "not_evaluable_reason": ";".join(detection_upstream_reasons),
                        "blank_profile_key": None,
                    }
                )
            elif config.detection.enabled and blank_data is not None:
                blank_key = next(
                    (key for key in (frame_id, scan_id, "all") if key in blank_data),
                    None,
                )
                if blank_key is not None:
                    signal_profile = _tail_centre_profile(result)
                    selected_blank_profiles = blank_data[blank_key]
                    detection = detect_tail_extent(
                        signal_profile,
                        selected_blank_profiles,
                        dz_um=config.calibration.dz_um,
                        bin_rows=config.detection.bin_rows,
                        threshold_mad_multiplier=config.detection.threshold_mad_multiplier,
                        minimum_consecutive_bins=config.detection.minimum_consecutive_bins,
                        minimum_blank_samples_per_bin=(
                            config.detection.minimum_blank_samples_per_bin
                        ),
                    )
                    detection_records.append(
                        {
                            "scan_id": scan_id,
                            "frame_id": frame_id,
                            "status": detection.invalid_reason,
                            "detected": detection.detected,
                            "detectable_length_um": detection.detectable_length_um,
                            "right_censored": detection.right_censored,
                            "evaluable": detection.invalid_reason
                            in ("ok", "right_censored", "not_detected"),
                            "not_evaluable_reason": (
                                None
                                if detection.invalid_reason
                                in ("ok", "right_censored", "not_detected")
                                else detection.invalid_reason
                            ),
                            "blank_profile_key": blank_key,
                        }
                    )
                    save_detection_qc(
                        qc_dir / f"{frame_stem}_QC03_detection.png",
                        scan_id=frame_id,
                        detection=detection,
                    )
                    if blank_source is None or blank_source_sha256 is None:
                        raise AssertionError("loaded blank source metadata is missing")
                    _save_detection_archive(
                        arrays_dir / f"{frame_stem}_detection.mat",
                        signal_profile=signal_profile,
                        blank_profiles=selected_blank_profiles,
                        detection=detection,
                        blank_profile_key=blank_key,
                        blank_source_path=blank_source,
                        blank_source_sha256=blank_source_sha256,
                        config=config,
                    )
                    bins = detection.bin_records(scan_id)
                    for bin_record in bins:
                        bin_record["frame_id"] = frame_id
                    detection_bin_records.extend(bins)
                else:
                    detection_records.append(
                        {
                            "scan_id": scan_id,
                            "frame_id": frame_id,
                            "status": "not_run_missing_scan_or_all_blank_key",
                            "detected": False,
                            "detectable_length_um": np.nan,
                            "right_censored": False,
                            "evaluable": False,
                            "not_evaluable_reason": (
                                "no frame_id, scan_id, or all key in blank NPZ"
                            ),
                            "blank_profile_key": None,
                        }
                    )
            elif config.detection.enabled:
                detection_records.append(
                    {
                        "scan_id": scan_id,
                        "frame_id": frame_id,
                        "status": "not_run_missing_matched_blank_profiles",
                        "detected": False,
                        "detectable_length_um": np.nan,
                        "right_censored": False,
                        "evaluable": False,
                        "not_evaluable_reason": "matched blank profiles not supplied",
                        "blank_profile_key": None,
                    }
                )
    except Exception as error:
        _write_json(
            output / "run_failed.json",
            {
                "status": "failed",
                "error_type": type(error).__name__,
                "error": str(error),
                "completed_frame_count": len(frame_records),
            },
        )
        raise
    frame_table = pd.DataFrame.from_records(frame_records)
    profile_table = pd.DataFrame.from_records(profile_records)
    localization_table = pd.DataFrame.from_records(localization_records)
    sensitivity_table = pd.DataFrame.from_records(sensitivity_records)
    detection_table = pd.DataFrame.from_records(
        detection_records,
        columns=[
            "scan_id",
            "frame_id",
            "status",
            "detected",
            "detectable_length_um",
            "right_censored",
            "evaluable",
            "not_evaluable_reason",
            "blank_profile_key",
        ],
    )
    detection_bins_table = pd.DataFrame.from_records(
        detection_bin_records,
        columns=[
            "scan_id",
            "frame_id",
            "bin_index_0based",
            "bin_start_um",
            "bin_stop_um",
            "signal_bin_mean",
            "blank_median",
            "blank_scaled_mad",
            "threshold",
            "blank_count",
            "exceeds_threshold",
            "sustained_detection",
        ],
    )
    manual_adjustments_table = pd.DataFrame.from_records(
        manual_adjustment_records,
        columns=[
            "scan_id",
            "frame_id",
            "adjustment_type",
            "value",
            "reason",
            "source",
        ],
    )
    if not detection_table.empty:
        detection_status = detection_table[
            [
                "frame_id",
                "status",
                "detected",
                "detectable_length_um",
                "right_censored",
                "evaluable",
                "not_evaluable_reason",
                "blank_profile_key",
            ]
        ].rename(
            columns={
                "status": "detection_qc_status",
                "detected": "tail_detected",
                "right_censored": "detection_right_censored",
                "evaluable": "detection_qc_valid",
                "not_evaluable_reason": "detection_not_evaluable_reason",
            }
        )
        frame_table = frame_table.merge(detection_status, on="frame_id", how="left")
    else:
        frame_table["detection_qc_status"] = "disabled"
        frame_table["detection_qc_valid"] = False
        frame_table["tail_detected"] = False
        frame_table["detectable_length_um"] = np.nan
        frame_table["detection_right_censored"] = False
        frame_table["detection_not_evaluable_reason"] = "detection disabled"
        frame_table["blank_profile_key"] = None
    scan_summary_table = _build_scan_summary(frame_table, detection_table)
    frame_path = output / "frame_results.csv"
    profiles_path = output / "profiles.csv"
    frame_table.to_csv(frame_path, index=False)
    profile_table.to_csv(profiles_path, index=False)
    localization_table.to_csv(output / "localization.csv", index=False)
    sensitivity_table.to_csv(output / "sensitivity_results.csv", index=False)
    detection_table.to_csv(output / "detection_results.csv", index=False)
    detection_bins_table.to_csv(output / "detection_bins.csv", index=False)
    scan_summary_table.to_csv(output / "scan_summary.csv", index=False)
    manual_adjustments_table.to_csv(
        logs_dir / "manual_adjustments.csv", index=False
    )

    valid_count = int(frame_table["valid"].sum()) if len(frame_table) else 0
    metadata = {
        "status": "complete",
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "schema_version": config.schema_version,
        "parameter_status": config.parameter_status,
        "formal_signal": "var(abs(E), 1, 3)",
        "code_git": _git_state(),
        "config_sha256": _sha256(config_source),
        "manifest_sha256": _sha256(manifest_source),
        "frame_count": len(frame_table),
        "valid_frame_count": valid_count,
        "invalid_frame_count": len(frame_table) - valid_count,
        "frame_processing_error_count": int(
            frame_table["invalid_reason"]
            .astype(str)
            .str.startswith("frame_processing_error:")
            .sum()
        ),
        "identity_incomplete_frame_count": int(
            (~frame_table["manifest_identity_complete"].astype(bool)).sum()
        ),
        "acquisition_metadata_incomplete_frame_count": int(
            (~frame_table["manifest_acquisition_metadata_complete"].astype(bool)).sum()
        ),
        "manual_adjustment_count": len(manual_adjustments_table),
        "detection_blank_profiles": (
            None if blank_profiles_path is None else str(Path(blank_profiles_path).resolve())
        ),
        "detection_blank_profiles_sha256": (
            blank_source_sha256
        ),
        "environment": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "h5py": h5py.__version__,
            "matplotlib": matplotlib.__version__,
        },
        "products": [
            "frame_results.csv",
            "scan_summary.csv",
            "profiles.csv",
            "localization.csv",
            "sensitivity_results.csv",
            "detection_results.csv",
            "detection_bins.csv",
            "profiles/",
            "arrays/",
            "qc/",
            "logs/",
            "logs/manual_adjustments.csv",
        ],
    }
    _write_json(output / "run_complete.json", metadata)
    _write_json(logs_dir / "run_complete.json", metadata)
    return BatchRunResult(
        output_dir=output,
        frame_count=len(frame_table),
        valid_frame_count=valid_count,
        frame_results_path=frame_path,
        profiles_path=profiles_path,
    )
