"""Adapter between the mentor's full-volume tracker and our SV geometry.

The tracking image is used only to locate the vessel. Formal tail
quantification continues to use this project's linear sv_raw maps.
"""

from __future__ import annotations

import hashlib
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import pandas as pd

from .geometry import VesselGeometry
from .localization import (
    LocalizationResult,
    LocalBodyResult,
    MentorTrackingEvidence,
    TopEdgeResult,
)
from .mentor import tracking_core
from .mentor.xroi import (
    add_assessability_score,
    isolated_jump_correct,
    local_body_features,
)


PRIMARY_ALPHA = 0.15
FINAL_X_METHOD = "X4_centroid_isolated_jump_corrected"
LATERAL_WIDTH_METHOD = "X1_continuous_run_recentered_on_X4"
TOP_METHOD = "primary_alpha_z_upper"
REQUIRED_ASSESSABILITY = "assessable"

REQUIRED_TRACKING_COLUMNS = {
    "scan_id",
    "frame_index",
    "alpha",
    "x_center_px",
    "z_upper_px",
    "new_tracking_class",
    "valid_local_body",
    "x1_local_geometry_px",
    "x2_robust_centroid_px",
    "x4_centroid_isolated_jump_corrected_px",
    "x4_jump_corrected",
    "x1_fallback",
    "local_body_run_width_px",
    "expected_lateral_width_px",
    "local_body_background",
    "local_body_sigma",
    "local_body_peak_cnr",
    "local_body_axial_completeness",
    "assessability_score",
    "vessel_presence_prediction",
}


@dataclass(frozen=True)
class MentorTrackingBundle:
    primary_table_path: Path
    metadata_path: Path
    alpha_table_paths: tuple[Path, ...]
    frame_count: int


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _safe_name(value: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("._")
    return safe or "scan"


def _numeric(record: pd.Series, name: str) -> float:
    try:
        result = float(record.get(name))
    except (TypeError, ValueError) as error:
        raise ValueError(f"mentor tracking column {name} must be numeric") from error
    if not np.isfinite(result):
        raise ValueError(f"mentor tracking column {name} must be finite")
    return result


def _boolean(value: Any, *, name: str) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    normalized = str(value).strip().lower()
    if normalized in {"true", "1"}:
        return True
    if normalized in {"false", "0"}:
        return False
    raise ValueError(f"mentor tracking column {name} must be boolean")


def load_flow_dicom(path: str | Path) -> np.ndarray:
    """Load an uncompressed Flow DICOM as float32 [frame, z, x]."""

    import pydicom

    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(source)
    dataset = pydicom.dcmread(source)
    volume = np.asarray(dataset.pixel_array, dtype=np.float32)
    if volume.ndim == 2:
        volume = volume[np.newaxis, :, :]
    if volume.ndim != 3:
        raise ValueError(
            f"Flow DICOM must decode to [frame,z,x], got {volume.shape}"
        )
    expected_rows = int(dataset.Rows)
    expected_columns = int(dataset.Columns)
    if volume.shape[1:] != (expected_rows, expected_columns):
        raise ValueError(
            "Flow DICOM pixel order does not match Rows and Columns: "
            f"{volume.shape[1:]} versus {(expected_rows, expected_columns)}"
        )
    declared_frames = int(getattr(dataset, "NumberOfFrames", 1))
    if volume.shape[0] != declared_frames:
        raise ValueError(
            f"Flow DICOM frame count mismatch: {volume.shape[0]} versus "
            f"{declared_frames}"
        )
    if not np.isfinite(volume).all():
        raise ValueError("Flow DICOM contains non-finite pixels")
    return volume


def build_mentor_tracking_tables(
    volume: np.ndarray,
    *,
    scan_id: str,
    diameter_um: float,
    tracking_config: Mapping[str, Any] | None = None,
) -> tuple[pd.DataFrame, dict[float, pd.DataFrame], dict[str, Any]]:
    """Run full-volume Viterbi, X1/X2/X4 and assessability."""

    effective_config = tracking_core.merge_tracking_config(tracking_config)
    primary_alpha = float(effective_config["primary_alpha"])
    outputs, _, _ = tracking_core.track_volume(
        np.asarray(volume, dtype=np.float32),
        scan_id=scan_id,
        diameter_um=float(diameter_um),
        config=effective_config,
    )
    primary_key = min(outputs, key=lambda value: abs(value - primary_alpha))
    if not np.isclose(primary_key, primary_alpha, rtol=0.0, atol=1e-12):
        raise ValueError(f"primary alpha {primary_alpha} was not generated")
    tracking = outputs[primary_key].sort_values("frame_index").reset_index(drop=True)
    if tracking["frame_index"].duplicated().any():
        raise ValueError("mentor tracker generated duplicate frame indices")
    expected = np.arange(len(tracking), dtype=int)
    if not np.array_equal(tracking["frame_index"].to_numpy(dtype=int), expected):
        raise ValueError("mentor tracker frame indices must be contiguous from zero")

    feature_rows: list[dict[str, Any]] = []
    for frame_index, row in tracking.iterrows():
        feature = local_body_features(volume[frame_index], row)
        feature["frame_index"] = int(frame_index)
        feature_rows.append(feature)
    features = pd.DataFrame.from_records(feature_rows)
    neighbor = features["local_body_peak_cnr"].shift(1).combine(
        features["local_body_peak_cnr"].shift(-1),
        lambda left, right: np.nanmedian([left, right]),
    )
    features["neighbor_peak_cnr"] = neighbor.fillna(
        features["local_body_peak_cnr"]
    )
    features = add_assessability_score(features)
    corrected, changed = isolated_jump_correct(
        features["x2_robust_centroid_px"].to_numpy(float)
    )
    features["x4_centroid_isolated_jump_corrected_px"] = corrected
    features["x4_jump_corrected"] = changed
    feature_columns = [
        column
        for column in features.columns
        if column != "frame_index" and column not in tracking.columns
    ]
    combined = pd.concat(
        [
            tracking.reset_index(drop=True),
            features[feature_columns].reset_index(drop=True),
        ],
        axis=1,
    )
    metadata = {
        "primary_alpha": primary_alpha,
        "final_x_method": FINAL_X_METHOD,
        "lateral_width_method": LATERAL_WIDTH_METHOD,
        "top_method": TOP_METHOD,
        "upper_edge_method": str(effective_config["upper_edge_method"]),
        "required_assessability": REQUIRED_ASSESSABILITY,
        "effective_tracking_config": effective_config,
    }
    return combined, outputs, metadata


def write_mentor_tracking_bundle(
    flow_dicom_path: str | Path,
    *,
    scan_id: str,
    diameter_um: float,
    output_dir: str | Path,
    tracking_config: Mapping[str, Any] | None = None,
) -> MentorTrackingBundle:
    """Run and save one auditable localization-only tracking bundle."""

    source = Path(flow_dicom_path).resolve()
    output = Path(output_dir).resolve()
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"output directory is not empty: {output}")
    volume = load_flow_dicom(source)
    primary, alpha_tables, metadata = build_mentor_tracking_tables(
        volume,
        scan_id=scan_id,
        diameter_um=diameter_um,
        tracking_config=tracking_config,
    )
    output.mkdir(parents=True, exist_ok=True)
    stem = _safe_name(scan_id)
    primary_path = output / f"{stem}_mentor_tracking.csv"
    primary.to_csv(primary_path, index=False)
    alpha_paths: list[Path] = []
    for alpha, table in sorted(alpha_tables.items()):
        alpha_tag = tracking_core.alpha_tag(float(alpha))
        path = output / f"{stem}_alpha{alpha_tag}_tracking.csv"
        table.to_csv(path, index=False)
        alpha_paths.append(path)
    core_path = Path(tracking_core.__file__).resolve()
    metadata.update(
        {
            "status": "complete",
            "scan_id": scan_id,
            "diameter_um": float(diameter_um),
            "flow_dicom_path": str(source),
            "flow_dicom_sha256": _sha256(source),
            "volume_shape_frame_z_x": list(volume.shape),
            "primary_table": primary_path.name,
            "primary_table_sha256": _sha256(primary_path),
            "alpha_tables": [path.name for path in alpha_paths],
            "mentor_tracking_core_path": str(core_path),
            "mentor_tracking_core_sha256": _sha256(core_path),
            "metric_boundary": (
                "localization only; no mentor tail AUC metric is computed"
            ),
        }
    )
    metadata_path = output / f"{stem}_mentor_tracking.metadata.json"
    with metadata_path.open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, ensure_ascii=False, indent=2, allow_nan=False)
        handle.write(chr(10))
    return MentorTrackingBundle(
        primary_table_path=primary_path,
        metadata_path=metadata_path,
        alpha_table_paths=tuple(alpha_paths),
        frame_count=len(primary),
    )


def select_tracking_record(
    path: str | Path,
    *,
    scan_id: str,
    frame_index: int,
) -> pd.Series:
    """Select exactly one auditable tracking row."""

    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(source)
    table = pd.read_csv(source)
    missing = sorted(REQUIRED_TRACKING_COLUMNS - set(table.columns))
    if missing:
        raise ValueError(
            "mentor tracking table missing columns: " + ", ".join(missing)
        )
    frame_values = pd.to_numeric(table["frame_index"], errors="coerce")
    selected = table.loc[
        table["scan_id"].astype(str).eq(str(scan_id))
        & frame_values.eq(int(frame_index))
    ]
    if len(selected) != 1:
        raise ValueError(
            "mentor tracking table must contain exactly one row for "
            f"scan_id={scan_id}, frame_index={frame_index}; found {len(selected)}"
        )
    return selected.iloc[0]


def build_localization_from_tracking(
    record: pd.Series,
    *,
    diameter_um: float,
    dx_um: float,
    dz_um: float,
    required_assessability: str = REQUIRED_ASSESSABILITY,
) -> LocalizationResult:
    """Map mentor X4/X1/z_upper outputs into our physical source geometry."""

    x_viterbi = _numeric(record, "x_center_px")
    z_upper = _numeric(record, "z_upper_px")
    x1_center = _numeric(record, "x1_local_geometry_px")
    x2_center = _numeric(record, "x2_robust_centroid_px")
    x4_center = _numeric(
        record, "x4_centroid_isolated_jump_corrected_px"
    )
    run_width_value = _numeric(record, "local_body_run_width_px")
    if run_width_value != round(run_width_value) or run_width_value < 1:
        raise ValueError("local_body_run_width_px must be a positive integer")
    run_width = int(round(run_width_value))
    expected_width = int(round(_numeric(record, "expected_lateral_width_px")))
    if expected_width < 1:
        raise ValueError("expected_lateral_width_px must be positive")
    valid_local_body = _boolean(
        record.get("valid_local_body"), name="valid_local_body"
    )
    x1_fallback = _boolean(record.get("x1_fallback"), name="x1_fallback")
    x4_changed = _boolean(
        record.get("x4_jump_corrected"), name="x4_jump_corrected"
    )
    tracking_class = str(
        record.get("new_tracking_class", record.get("tracking_class", ""))
    ).strip()
    assessment = str(record.get("vessel_presence_prediction", "")).strip()
    reasons: list[str] = []
    if tracking_class not in {"high_confidence", "model_assisted"}:
        reasons.append("mentor_tracking_failed")
    if not valid_local_body:
        reasons.append("mentor_local_body_invalid")
    if x1_fallback:
        reasons.append("mentor_x1_fallback")
    if assessment != required_assessability:
        reasons.append(f"mentor_assessability_{assessment or 'missing'}")
    qc_valid = not reasons

    original_left = int(round(x1_center - (run_width - 1) / 2.0))
    original_right = original_left + run_width - 1
    lateral = LocalBodyResult(
        valid=valid_local_body,
        invalid_reason="ok" if valid_local_body else "mentor_local_body_invalid",
        x_anchor_center_px=x_viterbi,
        x_left_center_px=original_left,
        x_right_center_px=original_right,
        x1_local_geometry_center_px=x1_center,
        fallback=x1_fallback,
        run_width_px=run_width,
        expected_lateral_width_px=expected_width,
        local_background=_numeric(record, "local_body_background"),
        local_sigma=_numeric(record, "local_body_sigma"),
        peak_cnr=_numeric(record, "local_body_peak_cnr"),
        axial_completeness=_numeric(
            record, "local_body_axial_completeness"
        ),
    )
    top_peak = record.get(
        "peak_snr", record.get("z_edge_local_robust_snr", np.nan)
    )
    try:
        top_peak_cnr = float(top_peak)
    except (TypeError, ValueError):
        top_peak_cnr = float("nan")
    top = TopEdgeResult(
        valid=True,
        invalid_reason="ok",
        z_anchor_center_px=z_upper,
        z_top_center_px=z_upper,
        fallback=False,
        local_background=float("nan"),
        local_sigma=float("nan"),
        peak_cnr=top_peak_cnr,
    )
    geometry = VesselGeometry(
        x_left_edge_px=x4_center - run_width / 2.0,
        x_right_edge_px=x4_center + run_width / 2.0,
        z_top_edge_px=z_upper - 0.5,
        diameter_um=float(diameter_um),
        dx_um=float(dx_um),
        dz_um=float(dz_um),
    )
    evidence = MentorTrackingEvidence(
        scan_id=str(record.get("scan_id", "")),
        frame_index=int(round(_numeric(record, "frame_index"))),
        primary_alpha=_numeric(record, "alpha"),
        viterbi_x_center_px=x_viterbi,
        z_upper_px=z_upper,
        tracking_class=tracking_class,
        x_path_confidence_class=str(
            record.get("x_path_confidence_class", "")
        ),
        z_edge_confidence_class=str(
            record.get("z_edge_confidence_class", "")
        ),
        x2_robust_centroid_px=x2_center,
        x4_centroid_isolated_jump_corrected_px=x4_center,
        x4_jump_corrected=x4_changed,
        assessability_score=_numeric(record, "assessability_score"),
        vessel_presence_prediction=assessment,
        valid_local_body=valid_local_body,
        x1_fallback=x1_fallback,
        qc_valid=qc_valid,
        invalid_reason="ok" if qc_valid else ";".join(reasons),
    )
    return LocalizationResult(
        geometry=geometry,
        lateral=lateral,
        top_edge=top,
        mentor_tracking=evidence,
    )
