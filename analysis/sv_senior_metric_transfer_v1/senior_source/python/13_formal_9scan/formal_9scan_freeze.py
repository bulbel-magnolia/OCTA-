#!/usr/bin/env python
"""Run the predeclared formal nine-scan OCTA tail freeze Gate.

This is deliberately a *freeze* runner, not a method-development program.
It consumes the already-frozen X4/assessability full-frame table and the
alpha=0.15 tracking tables, and never recomputes or tunes X, ROI, z, guard,
window, assessability thresholds, or normalization parameters.  The script
loads one exported Flow volume at a time, keeps all 500 frames, and writes
auditable frame-, scan-, spatial-, sensitivity-, denominator-, and Gate-level
results.  It does not perform velocity/diameter inference and it never
extends beyond the nine explicitly listed scans.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import importlib.util
import json
import math
import shutil
import sys
import uuid
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from tailq.io import (  # noqa: E402
    export_path,
    load_config,
    load_flow,
    load_manifest,
    output_root,
    tracking_csv_path,
)
from tailq.metrics import _interval_values, _normalizer, _window_pixels, extract_profiles  # noqa: E402


# This literal list is an intentional identity lock.  It is checked against
# the historical package list and the manifest before any Flow is read.
FROZEN_SCAN_IDS = (
    "d185_s1_scan32",
    "d185_s5_scan38",
    "d185_s10_scan5",
    "d235_s1_scan17",
    "d235_s5_scan23",
    "d235_s10_scan29",
    "d285_s1_scan2",
    "d285_s5_scan8",
    "d285_s10_scan14",
)
EXPECTED_CONDITIONS = {
    "d185_s1_scan32": (185.0, 1.0),
    "d185_s5_scan38": (185.0, 5.0),
    "d185_s10_scan5": (185.0, 10.0),
    "d235_s1_scan17": (235.0, 1.0),
    "d235_s5_scan23": (235.0, 5.0),
    "d235_s10_scan29": (235.0, 10.0),
    "d285_s1_scan2": (285.0, 1.0),
    "d285_s5_scan8": (285.0, 5.0),
    "d285_s10_scan14": (285.0, 10.0),
}
SLOW_AXIS_REGIONS = (
    ("front", 0, 167),
    ("middle", 167, 333),
    ("rear", 333, 500),
)
REGION_NAMES = tuple(item[0] for item in SLOW_AXIS_REGIONS)

PRIMARY_ALPHA = 0.15
PRIMARY_ROI_FRACTION = 0.40
PRIMARY_GUARD = 2
GUARDS = (0, 2, 4)
PRIMARY_WINDOW_UM = (0, 200)
SECONDARY_WINDOWS_UM = ((0, 100), (100, 200), (200, 300))
AXIAL_UM_PER_PX = 6.7
LATERAL_UM_PER_PX = 12.7
FRAMES_PER_SCAN = 500
N_Z = 351
N_X = 500
EPSILON = 1e-12

# These are Gate operationalizations, declared in code before the formal run.
# The 10/15% limits come directly from the user's predeclared Gate.  The
# 50%-valid rule reuses the frozen aggregation QC floor.  A 25% regional
# departure is only a severe spatial-dominance diagnostic, not a biological
# effect threshold.  P95 QC is intentionally conservative and only decides
# whether P95 remains secondary (Gate B); it cannot change the raw Gate.
VALID_FRACTION_FLOOR = 0.50
SENS_MEDIAN_LIMIT_PCT = 10.0
SENS_MATERIAL_LIMIT_PCT = 15.0
GUARD_MATERIAL_SCAN_COUNT_MAX = 1
SPATIAL_SEVERE_DEVIATION_PCT = 25.0
P95_VALID_RATE_FLOOR = 0.90
P95_ROBUST_CV_LIMIT = 0.50

DEFAULT_CONFIG = PROJECT_ROOT / "configs" / "pilot_20260414.json"
DEFAULT_FEATURES = (
    PROJECT_ROOT
    / "work"
    / "pilot_20260414"
    / "xroi_assessability_development_20260827"
    / "assessability_all_frames_features.csv"
)
DEFAULT_OUTPUT = (
    PROJECT_ROOT
    / "work"
    / "pilot_20260414"
    / "formal_9scan_freeze_20260902"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(DEFAULT_CONFIG))
    parser.add_argument("--features", default=str(DEFAULT_FEATURES))
    parser.add_argument("--tracking-root", default=None)
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT))
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _numeric(value: Any) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return np.nan
    return result if np.isfinite(result) else np.nan


def _finite_values(values: Iterable[Any]) -> np.ndarray:
    array = pd.to_numeric(pd.Series(list(values)), errors="coerce").to_numpy(float)
    return array[np.isfinite(array)]


def _summary(values: Iterable[Any]) -> dict[str, Any]:
    array = _finite_values(values)
    if not len(array):
        return {"n": 0, "median": np.nan, "iqr": np.nan, "p25": np.nan, "p75": np.nan}
    p25, median, p75 = np.quantile(array, [0.25, 0.50, 0.75])
    return {
        "n": int(len(array)),
        "median": float(median),
        "iqr": float(p75 - p25),
        "p25": float(p25),
        "p75": float(p75),
    }


def _robust_cv(values: Iterable[Any]) -> float:
    array = _finite_values(values)
    if not len(array):
        return np.nan
    median = float(np.median(array))
    if not np.isfinite(median) or abs(median) <= EPSILON:
        return np.nan
    mad = float(np.median(np.abs(array - median)))
    return float(1.4826 * mad / abs(median))


def _relative_change(reference: Any, comparison: Any) -> float:
    ref = _numeric(reference)
    comp = _numeric(comparison)
    if not np.isfinite([ref, comp]).all() or abs(ref) <= EPSILON:
        return np.nan
    return float(100.0 * abs(comp - ref) / abs(ref))


def _signed_change(reference: Any, comparison: Any) -> float:
    ref = _numeric(reference)
    comp = _numeric(comparison)
    if not np.isfinite([ref, comp]).all() or abs(ref) <= EPSILON:
        return np.nan
    return float(100.0 * (comp - ref) / abs(ref))


def _region_for_frame(frame_index: int) -> str:
    for name, start, stop in SLOW_AXIS_REGIONS:
        if start <= int(frame_index) < stop:
            return name
    raise ValueError(f"frame_index outside 0..499: {frame_index}")


def _portable_config(cfg: dict[str, Any]) -> tuple[dict[str, Any], Path, str]:
    """Use the relocated project output root without editing the JSON config."""

    configured = output_root(cfg)
    portable = PROJECT_ROOT / "work" / "pilot_20260414"
    configured_exports = configured / str(cfg["output"]["export_subdir"])
    portable_exports = portable / str(cfg["output"]["export_subdir"])
    if configured_exports.exists():
        return cfg, configured, "configured_root"
    if portable_exports.exists():
        effective = copy.deepcopy(cfg)
        effective["output"]["root"] = str(portable)
        return effective, portable, "relocated_project_root"
    return cfg, configured, "configured_root_missing_exports"


def _load_frozen_features(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Frozen assessability feature table not found: {path}")
    table = pd.read_csv(path)
    required = {
        "dataset_id",
        "frame_index_0based",
        "diameter_um",
        "velocity_mm_s",
        "assessability_score",
        "vessel_presence_prediction",
        "x4_centroid_isolated_jump_corrected_px",
        "z_upper_algorithm_px",
    }
    missing = required.difference(table.columns)
    if missing:
        raise ValueError(f"Frozen features missing required columns: {sorted(missing)}")
    if len(table) != len(FROZEN_SCAN_IDS) * FRAMES_PER_SCAN:
        raise ValueError(
            f"Frozen feature row count is {len(table)}, expected "
            f"{len(FROZEN_SCAN_IDS) * FRAMES_PER_SCAN} (9 scans x 500 frames)"
        )
    if set(table["dataset_id"].astype(str)) != set(FROZEN_SCAN_IDS):
        raise ValueError("Frozen feature table does not contain exactly the locked nine scans")
    table = table.copy()
    table["dataset_id"] = table["dataset_id"].astype(str)
    table["frame_index_0based"] = pd.to_numeric(table["frame_index_0based"], errors="raise").astype(int)
    for scan_id in FROZEN_SCAN_IDS:
        subset = table.loc[table["dataset_id"].eq(scan_id)]
        if len(subset) != FRAMES_PER_SCAN or set(subset["frame_index_0based"]) != set(range(FRAMES_PER_SCAN)):
            raise ValueError(f"Frozen features have incomplete/duplicated frames for {scan_id}")
    table["diameter_um"] = pd.to_numeric(table["diameter_um"], errors="coerce")
    table["velocity_mm_s"] = pd.to_numeric(table["velocity_mm_s"], errors="coerce")
    table["assessability_score"] = pd.to_numeric(table["assessability_score"], errors="coerce")
    table["x4_centroid_isolated_jump_corrected_px"] = pd.to_numeric(
        table["x4_centroid_isolated_jump_corrected_px"], errors="coerce"
    )
    table["z_upper_algorithm_px"] = pd.to_numeric(table["z_upper_algorithm_px"], errors="coerce")
    allowed = {"assessable", "uncertain", "not_assessable"}
    if not table["vessel_presence_prediction"].astype(str).isin(allowed).all():
        raise ValueError("Frozen features contain an unknown assessability class")
    table["vessel_presence_prediction"] = table["vessel_presence_prediction"].astype(str)
    table["slow_axis_region"] = table["frame_index_0based"].map(_region_for_frame)
    return table.sort_values(["dataset_id", "frame_index_0based"]).reset_index(drop=True)


def _validate_manifest_and_inputs(cfg: dict[str, Any], config_path: Path, manifest: pd.DataFrame, tracking_root: Path) -> None:
    if tuple(FROZEN_SCAN_IDS) != tuple(EXPECTED_CONDITIONS):  # pragma: no cover - identity guard
        raise AssertionError("Internal frozen scan identity/condition map changed")
    if set(manifest.index if manifest.index.name == "dataset_id" else manifest["dataset_id"]) < set(FROZEN_SCAN_IDS):
        raise ValueError("Manifest does not contain all frozen scan IDs")
    indexed = manifest.set_index("dataset_id", drop=False)
    for scan_id in FROZEN_SCAN_IDS:
        if scan_id not in indexed.index:
            raise ValueError(f"Frozen scan identity is not in manifest: {scan_id}")
        row = indexed.loc[scan_id]
        expected_d, expected_v = EXPECTED_CONDITIONS[scan_id]
        if not np.isclose(float(row["diameter_um"]), expected_d) or not np.isclose(float(row["velocity_mm_s"]), expected_v):
            raise ValueError(
                f"Manifest identity mismatch for {scan_id}: "
                f"got diameter={row['diameter_um']}, velocity={row['velocity_mm_s']}"
            )
        flow_path = export_path(cfg, scan_id)
        if not flow_path.exists():
            raise FileNotFoundError(f"Missing frozen Flow export for {scan_id}: {flow_path}")
        tracking_path = tracking_csv_path(tracking_root, scan_id, PRIMARY_ALPHA)
        if not tracking_path.exists():
            raise FileNotFoundError(f"Missing alpha=0.15 tracking for {scan_id}: {tracking_path}")


def _normalization_fields(profiles: Any, z_upper: int, diameter_px: int, metrics_cfg: dict[str, Any]) -> dict[str, Any]:
    p95, p95_valid, core_start, core_stop = _normalizer(
        profiles.excess,
        z_upper,
        diameter_px,
        metrics_cfg["core_fraction"],
        "p95",
        EPSILON,
    )
    # The public helper intentionally supports p95/median/none.  P90 is only
    # a QC audit control, so calculate it explicitly on the identical core.
    if 0 <= core_start < core_stop <= len(profiles.excess):
        core = profiles.excess[core_start:core_stop]
        p90 = float(np.quantile(core, 0.90)) if np.isfinite(core).all() else np.nan
        p90_valid = bool(np.isfinite(p90) and p90 > EPSILON)
    else:
        p90, p90_valid = np.nan, False
    return {
        "p95_denominator": float(p95) if np.isfinite(p95) else np.nan,
        "p95_denominator_valid": bool(p95_valid),
        "p90_denominator": float(p90) if np.isfinite(p90) else np.nan,
        "p90_denominator_valid": bool(p90_valid),
        "core_start_px": int(core_start),
        "core_stop_exclusive_px": int(core_stop),
    }


def _interval_result(profiles: Any, z_start: int, start_um: int, stop_um: int) -> dict[str, Any]:
    start_px = _window_pixels(start_um, AXIAL_UM_PER_PX)
    stop_px = _window_pixels(stop_um, AXIAL_UM_PER_PX)
    values = _interval_values(
        profiles,
        z_start=z_start,
        offset_start_px=start_px,
        offset_stop_px=stop_px,
        axial_um_per_px=AXIAL_UM_PER_PX,
    )
    name = f"{start_um}_{stop_um}"
    return {
        "name": name,
        "start_px_offset": int(start_px),
        "stop_px_offset": int(stop_px),
        "complete": bool(values["complete"]),
        "raw_auc": float(values["tail_auc_raw"]) if np.isfinite(values["tail_auc_raw"]) else np.nan,
        "central_profile_auc": float(values["central_profile_auc"]) if np.isfinite(values["central_profile_auc"]) else np.nan,
        "background_auc": float(values["background_auc"]) if np.isfinite(values["background_auc"]) else np.nan,
        "central_profile_median": float(values["central_profile_median"]) if np.isfinite(values["central_profile_median"]) else np.nan,
        "background_median": float(values["background_median"]) if np.isfinite(values["background_median"]) else np.nan,
        "tail_excess_median": float(values["tail_excess_median"]) if np.isfinite(values["tail_excess_median"]) else np.nan,
        "invalid_reason": "ok" if bool(values["complete"]) else "incomplete_window_or_nonfinite_profile",
    }


def _frame_metrics(
    frame: np.ndarray,
    *,
    x_x4: float,
    z_upper: float,
    diameter_um: float,
    metrics_cfg: dict[str, Any],
) -> dict[str, Any]:
    """Compute frozen geometry metrics for one frame, including guard QC."""

    result: dict[str, Any] = {
        "x_x4": x_x4,
        "z_upper": z_upper,
        "diameter_axial_px": np.nan,
        "tail_start_guard0_px": np.nan,
        "tail_start_px": np.nan,
        "tail_start": np.nan,
        "tail_start_um": np.nan,
        "tail_start_guard4_px": np.nan,
        "roi_x0": np.nan,
        "roi_x1_exclusive": np.nan,
        "side_left_x0": np.nan,
        "side_left_x1_exclusive": np.nan,
        "side_right_x0": np.nan,
        "side_right_x1_exclusive": np.nan,
        "roi_complete": False,
        "p95_denominator": np.nan,
        "p95_denominator_valid": False,
        "p90_denominator": np.nan,
        "p90_denominator_valid": False,
        "core_start_px": np.nan,
        "core_stop_exclusive_px": np.nan,
    }
    interval_names = ("0_100", "100_200", "200_300")
    for guard in GUARDS:
        result[f"raw_auc_0_200_guard{guard}"] = np.nan
        result[f"raw_auc_0_200_guard{guard}_valid"] = False
        result[f"invalid_reason_guard{guard}"] = "not_computed"
    # Canonical primary aliases are populated from guard=2 below.  Keeping
    # these names separate makes the CSV self-explanatory while preserving all
    # three guard sensitivity columns.
    result["raw_auc_0_200"] = np.nan
    result["raw_auc_0_200_valid"] = False
    result["invalid_reason_raw_auc_0_200"] = "not_computed"
    for name in interval_names:
        result[f"raw_auc_{name}"] = np.nan
        result[f"raw_auc_{name}_valid"] = False
        result[f"invalid_reason_raw_auc_{name}"] = "not_computed"
    result["p95_normalized_auc_0_200"] = np.nan
    result["p95_normalized_auc_0_200_valid"] = False
    result["invalid_reason_p95_normalized_auc_0_200"] = "not_computed"
    result["background_estimate"] = np.nan
    result["background_mad_estimate"] = np.nan

    reasons: list[str] = []
    if not np.isfinite([x_x4, z_upper, diameter_um]).all():
        result["invalid_reason"] = "missing_frozen_localization_or_condition"
        for guard in GUARDS:
            result[f"invalid_reason_guard{guard}"] = "missing_frozen_localization_or_condition"
        for name in interval_names:
            result[f"invalid_reason_raw_auc_{name}"] = "missing_frozen_localization_or_condition"
        result["invalid_reason_p95_normalized_auc_0_200"] = "missing_frozen_localization_or_condition"
        return {**result, "invalid_reason": "missing_frozen_localization_or_condition"}

    diameter_px = max(1, int(round(float(diameter_um) / AXIAL_UM_PER_PX)))
    z_int = int(round(float(z_upper)))
    result["diameter_axial_px"] = diameter_px
    result["z_upper"] = z_int
    result["tail_start_guard0_px"] = z_int + diameter_px + 0
    result["tail_start_px"] = z_int + diameter_px + PRIMARY_GUARD
    result["tail_start"] = result["tail_start_px"]
    result["tail_start_um"] = float(result["tail_start_px"] * AXIAL_UM_PER_PX)
    result["tail_start_guard4_px"] = z_int + diameter_px + 4

    local_metrics = dict(metrics_cfg)
    local_metrics["central_fraction"] = PRIMARY_ROI_FRACTION
    profiles = extract_profiles(
        frame,
        x_center=float(x_x4),
        diameter_um=float(diameter_um),
        lateral_um_per_px=LATERAL_UM_PER_PX,
        metrics_cfg=local_metrics,
    )
    result["roi_x0"], result["roi_x1_exclusive"] = profiles.central_bounds
    result["side_left_x0"], result["side_left_x1_exclusive"] = profiles.left_bounds
    result["side_right_x0"], result["side_right_x1_exclusive"] = profiles.right_bounds
    result["roi_complete"] = bool(profiles.roi_complete)
    if not profiles.roi_complete:
        reasons.append("roi_incomplete")
        for guard in GUARDS:
            result[f"invalid_reason_guard{guard}"] = "roi_incomplete"
        for name in interval_names:
            result[f"invalid_reason_raw_auc_{name}"] = "roi_incomplete"
        result["invalid_reason_p95_normalized_auc_0_200"] = "roi_incomplete"
        return {**result, "invalid_reason": ";".join(reasons)}

    normalization = _normalization_fields(profiles, z_int, diameter_px, local_metrics)
    result.update(normalization)
    core_slice = slice(int(normalization["core_start_px"]), int(normalization["core_stop_exclusive_px"]))
    if 0 <= core_slice.start < core_slice.stop <= len(profiles.background) and np.isfinite(profiles.background[core_slice]).all():
        result["background_estimate_core"] = float(np.median(profiles.background[core_slice]))
    else:
        result["background_estimate_core"] = np.nan

    for guard in GUARDS:
        z_start = z_int + diameter_px + guard
        item = _interval_result(profiles, z_start, 0, 200)
        result[f"raw_auc_0_200_guard{guard}"] = item["raw_auc"]
        result[f"raw_auc_0_200_guard{guard}_valid"] = bool(item["complete"])
        result[f"invalid_reason_guard{guard}"] = item["invalid_reason"]
        if guard == PRIMARY_GUARD:
            result["raw_auc_0_200"] = item["raw_auc"]
            result["raw_auc_0_200_valid"] = bool(item["complete"])
            result["invalid_reason_raw_auc_0_200"] = item["invalid_reason"]
            result["background_estimate"] = item["background_median"]
            if 0 <= z_start < len(profiles.background):
                stop = min(len(profiles.background), z_start + _window_pixels(200, AXIAL_UM_PER_PX))
                bg = profiles.background[z_start:stop]
                bg_mad = profiles.background_mad[z_start:stop]
                if len(bg) and np.isfinite(bg).all():
                    result["background_estimate"] = float(np.median(bg))
                if len(bg_mad) and np.isfinite(bg_mad).all():
                    result["background_mad_estimate"] = float(np.median(bg_mad))

    z_primary = z_int + diameter_px + PRIMARY_GUARD
    for start_um, stop_um in SECONDARY_WINDOWS_UM:
        item = _interval_result(profiles, z_primary, start_um, stop_um)
        name = f"{start_um}_{stop_um}"
        result[f"raw_auc_{name}"] = item["raw_auc"]
        result[f"raw_auc_{name}_valid"] = bool(item["complete"])
        result[f"invalid_reason_raw_auc_{name}"] = item["invalid_reason"]

    primary_item = _interval_result(profiles, z_primary, 0, 200)
    if not primary_item["complete"]:
        reasons.append("raw_auc_0_200_incomplete")
    if not normalization["p95_denominator_valid"]:
        reasons.append("p95_denominator_invalid")
    p95_norm_valid = bool(primary_item["complete"] and normalization["p95_denominator_valid"])
    result["p95_normalized_auc_0_200_valid"] = p95_norm_valid
    result["p95_normalized_auc_0_200"] = (
        float(primary_item["raw_auc"] / normalization["p95_denominator"])
        if p95_norm_valid
        else np.nan
    )
    result["invalid_reason_p95_normalized_auc_0_200"] = "ok" if p95_norm_valid else ";".join(
        [reason for reason in ("raw_auc_0_200_incomplete" if not primary_item["complete"] else "", "p95_denominator_invalid" if not normalization["p95_denominator_valid"] else "") if reason]
    )
    result["invalid_reason"] = ";".join(dict.fromkeys(reasons)) if reasons else "ok"
    return result


def _scan_row(frame_table: pd.DataFrame, scan_id: str, diameter_um: float, velocity_mm_s: float, technical_replicate: int) -> dict[str, Any]:
    assessable = frame_table["assessability_class"].eq("assessable")
    uncertain = frame_table["assessability_class"].eq("uncertain")
    not_assessable = frame_table["assessability_class"].eq("not_assessable")
    valid = frame_table["raw_auc_0_200_valid"].astype(bool)
    primary_mask = assessable & valid
    p95_mask = primary_mask & frame_table["p95_normalized_auc_0_200_valid"].astype(bool)
    raw = _summary(frame_table.loc[primary_mask, "raw_auc_0_200"])
    p95 = _summary(frame_table.loc[p95_mask, "p95_normalized_auc_0_200"])
    all_raw = _summary(frame_table.loc[valid, "raw_auc_0_200"])
    row: dict[str, Any] = {
        "scan_id": scan_id,
        "dataset_id": scan_id,
        "diameter_um": diameter_um,
        "velocity_mm_s": velocity_mm_s,
        "technical_replicate": technical_replicate,
        "n_frames_total": int(len(frame_table)),
        "assessable_count": int(assessable.sum()),
        "uncertain_count": int(uncertain.sum()),
        "not_assessable_count": int(not_assessable.sum()),
        "assessable_rate": float(assessable.mean()),
        "uncertain_rate": float(uncertain.mean()),
        "not_assessable_rate": float(not_assessable.mean()),
        "primary_valid_frame_count": int(primary_mask.sum()),
        "valid_frame_count": int(primary_mask.sum()),
        "assessable_frame_count": int(assessable.sum()),
        "uncertain_frame_count": int(uncertain.sum()),
        "not_assessable_frame_count": int(not_assessable.sum()),
        "primary_valid_fraction_of_assessable": float(primary_mask.sum() / assessable.sum()) if assessable.any() else np.nan,
        "all_class_valid_frame_count": int(valid.sum()),
        "all_class_valid_rate": float(valid.mean()),
        "raw_auc_0_200_primary_median": raw["median"],
        "raw_auc_0_200_primary_iqr": raw["iqr"],
        "raw_auc_0_200_primary_p25": raw["p25"],
        "raw_auc_0_200_primary_p75": raw["p75"],
        "p95_normalized_auc_0_200_primary_median": p95["median"],
        "p95_normalized_auc_0_200_primary_iqr": p95["iqr"],
        "p95_normalized_valid_frame_count": int(p95_mask.sum()),
        "p95_normalized_valid_fraction_of_primary": float(p95_mask.sum() / primary_mask.sum()) if primary_mask.any() else np.nan,
        "raw_auc_0_200_all_class_median": all_raw["median"],
        "raw_auc_0_200_all_class_iqr": all_raw["iqr"],
        "aggregation_rule": "pooled_median_of_all_valid_assessable_frames; no block medians",
    }
    return row


def _spatial_rows(frame_table: pd.DataFrame, scan_id: str, diameter_um: float, velocity_mm_s: float) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for region, start, stop in SLOW_AXIS_REGIONS:
        group = frame_table.loc[frame_table["slow_axis_region"].eq(region)]
        assessable = group["assessability_class"].eq("assessable")
        uncertain = group["assessability_class"].eq("uncertain")
        not_assessable = group["assessability_class"].eq("not_assessable")
        valid = assessable & group["raw_auc_0_200_valid"].astype(bool)
        p95_valid = valid & group["p95_normalized_auc_0_200_valid"].astype(bool)
        raw = _summary(group.loc[valid, "raw_auc_0_200"])
        p95 = _summary(group.loc[p95_valid, "p95_normalized_auc_0_200"])
        rows.append(
            {
                "scan_id": scan_id,
                "dataset_id": scan_id,
                "diameter_um": diameter_um,
                "velocity_mm_s": velocity_mm_s,
                "slow_axis_region": region,
                "frame_start_0based": start,
                "frame_stop_exclusive_0based": stop,
                "n_frames": int(len(group)),
                "assessable_count": int(assessable.sum()),
                "uncertain_count": int(uncertain.sum()),
                "not_assessable_count": int(not_assessable.sum()),
                "assessable_rate": float(assessable.mean()) if len(group) else np.nan,
                "uncertain_rate": float(uncertain.mean()) if len(group) else np.nan,
                "not_assessable_rate": float(not_assessable.mean()) if len(group) else np.nan,
                "primary_valid_count": int(valid.sum()),
                "raw_auc_0_200_median": raw["median"],
                "raw_auc_0_200_iqr": raw["iqr"],
                "p95_normalized_auc_0_200_median": p95["median"],
                "p95_normalized_auc_0_200_iqr": p95["iqr"],
            }
        )
    return rows


def _assessability_sensitivity(scan_table: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for scan_id, group in scan_table.groupby("scan_id", sort=False):
        primary = group.loc[(group["assessability_class"] == "assessable") & group["raw_auc_0_200_valid"].astype(bool), "raw_auc_0_200"]
        sensitivity = group.loc[group["assessability_class"].isin(["assessable", "uncertain"]) & group["raw_auc_0_200_valid"].astype(bool), "raw_auc_0_200"]
        p = _summary(primary)
        s = _summary(sensitivity)
        rows.append(
            {
                "scan_id": scan_id,
                "diameter_um": float(group["diameter_um"].iloc[0]),
                "velocity_mm_s": float(group["velocity_mm_s"].iloc[0]),
                "primary_definition": "assessable only",
                "sensitivity_definition": "assessable + uncertain",
                "primary_valid_count": p["n"],
                "sensitivity_valid_count": s["n"],
                "primary_raw_auc_median": p["median"],
                "sensitivity_raw_auc_median": s["median"],
                "absolute_relative_change_pct": _relative_change(p["median"], s["median"]),
                "signed_change_pct": _signed_change(p["median"], s["median"]),
                "status": "ok" if np.isfinite([p["median"], s["median"]]).all() else "invalid_scan_median",
            }
        )
    return pd.DataFrame(rows)


def _guard_sensitivity(scan_table: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for scan_id, group in scan_table.groupby("scan_id", sort=False):
        mask = group["assessability_class"].eq("assessable")
        medians: dict[int, float] = {}
        counts: dict[int, int] = {}
        for guard in GUARDS:
            valid_col = f"raw_auc_0_200_guard{guard}_valid"
            value_col = f"raw_auc_0_200_guard{guard}"
            values = group.loc[mask & group[valid_col].astype(bool), value_col]
            summary = _summary(values)
            medians[guard] = summary["median"]
            counts[guard] = summary["n"]
        rows.append(
            {
                "scan_id": scan_id,
                "diameter_um": float(group["diameter_um"].iloc[0]),
                "velocity_mm_s": float(group["velocity_mm_s"].iloc[0]),
                "guard0_primary_median": medians[0],
                "guard2_primary_median": medians[2],
                "guard4_primary_median": medians[4],
                "guard0_valid_count": counts[0],
                "guard2_valid_count": counts[2],
                "guard4_valid_count": counts[4],
                "guard0_vs_guard2_abs_relative_change_pct": _relative_change(medians[2], medians[0]),
                "guard4_vs_guard2_abs_relative_change_pct": _relative_change(medians[2], medians[4]),
                "guard0_vs_guard2_signed_change_pct": _signed_change(medians[2], medians[0]),
                "guard4_vs_guard2_signed_change_pct": _signed_change(medians[2], medians[4]),
            }
        )
    return pd.DataFrame(rows)


def _p95_qc(frame_table: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for scan_id, group in frame_table.groupby("scan_id", sort=False):
        for scope, mask in (
            ("all_frames", np.ones(len(group), dtype=bool)),
            ("primary_assessable", group["assessability_class"].eq("assessable").to_numpy(bool)),
        ):
            values = group.loc[mask & group["p95_denominator_valid"].astype(bool), "p95_denominator"]
            summary = _summary(values)
            all_n = int(mask.sum())
            rows.append(
                {
                    "scan_id": scan_id,
                    "diameter_um": float(group["diameter_um"].iloc[0]),
                    "velocity_mm_s": float(group["velocity_mm_s"].iloc[0]),
                    "scope": scope,
                    "n_frames_scope": all_n,
                    "valid_count": summary["n"],
                    "valid_rate": float(summary["n"] / all_n) if all_n else np.nan,
                    "median": summary["median"],
                    "iqr": summary["iqr"],
                    "robust_cv": _robust_cv(values),
                    "p5": float(np.quantile(_finite_values(values), 0.05)) if summary["n"] else np.nan,
                    "p90": float(np.quantile(_finite_values(values), 0.90)) if summary["n"] else np.nan,
                    "p95": float(np.quantile(_finite_values(values), 0.95)) if summary["n"] else np.nan,
                    "denominator_definition": "p95 of excess in frozen vessel core fraction [0.10, 0.45]",
                }
            )
    return pd.DataFrame(rows)


def _make_figures(scan_primary: pd.DataFrame, assess_sens: pd.DataFrame, guard_sens: pd.DataFrame, spatial: pd.DataFrame, p95_qc: pd.DataFrame, figure_dir: Path) -> None:
    figure_dir.mkdir(parents=True, exist_ok=True)
    labels = scan_primary["scan_id"].tolist()
    x = np.arange(len(labels))

    # 1. Stacked predicted-class rates.
    fig, ax = plt.subplots(figsize=(12, 5))
    bottom = np.zeros(len(scan_primary))
    colors = {"assessable_rate": "#2c7fb8", "uncertain_rate": "#fdae61", "not_assessable_rate": "#d7191c"}
    names = {"assessable_rate": "assessable", "uncertain_rate": "uncertain", "not_assessable_rate": "not_assessable"}
    for col in colors:
        values = scan_primary[col].to_numpy(float)
        ax.bar(x, values, bottom=bottom, label=names[col], color=colors[col])
        bottom += values
    ax.set_ylim(0, 1)
    ax.set_ylabel("fraction of 500 frames")
    ax.set_title("Frozen assessability composition (all 500 frames)")
    ax.set_xticks(x, labels, rotation=45, ha="right", fontsize=8)
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(figure_dir / "01_assessable_rate.png", dpi=220)
    plt.close(fig)

    # 2. Scan-level primary raw AUC.
    fig, ax = plt.subplots(figsize=(12, 5))
    values = scan_primary["raw_auc_0_200_primary_median"].to_numpy(float)
    ax.scatter(x, values, s=70, color="#2166ac")
    for xi, yi in zip(x, values):
        if np.isfinite(yi):
            ax.annotate(f"{yi:.0f}", (xi, yi), xytext=(0, 7), textcoords="offset points", ha="center", fontsize=8)
    ax.set_ylabel("raw AUC 0–200 μm (Flow·μm)")
    ax.set_title("Frozen scan-level primary raw tail burden")
    ax.set_xticks(x, labels, rotation=45, ha="right", fontsize=8)
    fig.tight_layout()
    fig.savefig(figure_dir / "02_scan_raw_auc_0_200.png", dpi=220)
    plt.close(fig)

    # 3. Assessability sensitivity paired plot.
    merged = assess_sens.set_index("scan_id").reindex(labels)
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.plot(x, merged["primary_raw_auc_median"], "o-", label="assessable only", color="#2166ac")
    ax.plot(x, merged["sensitivity_raw_auc_median"], "s--", label="assessable + uncertain", color="#e08214")
    ax.set_ylabel("raw AUC 0–200 μm median")
    ax.set_title("Assessability sensitivity (paired scan medians)")
    ax.set_xticks(x, labels, rotation=45, ha="right", fontsize=8)
    ax.legend()
    fig.tight_layout()
    fig.savefig(figure_dir / "03_assessability_sensitivity.png", dpi=220)
    plt.close(fig)

    # 4. Guard sensitivity.
    merged = guard_sens.set_index("scan_id").reindex(labels)
    fig, ax = plt.subplots(figsize=(12, 5))
    for guard, color in ((0, "#762a83"), (2, "#1b7837"), (4, "#af8dc3")):
        ax.plot(x, merged[f"guard{guard}_primary_median"], "o-", label=f"guard {guard} px", color=color)
    ax.set_ylabel("raw AUC 0–200 μm median")
    ax.set_title("Guard sensitivity (primary assessable frames)")
    ax.set_xticks(x, labels, rotation=45, ha="right", fontsize=8)
    ax.legend()
    fig.tight_layout()
    fig.savefig(figure_dir / "04_guard_sensitivity.png", dpi=220)
    plt.close(fig)

    # 5. Front/middle/rear spatial heterogeneity.
    fig, ax = plt.subplots(figsize=(12, 5))
    width = 0.25
    for index, region in enumerate(REGION_NAMES):
        subset = spatial.set_index(["scan_id", "slow_axis_region"]).reindex(pd.MultiIndex.from_product([labels, [region]], names=["scan_id", "slow_axis_region"]))
        vals = subset["raw_auc_0_200_median"].to_numpy(float)
        ax.bar(x + (index - 1) * width, vals, width=width, label=region)
    ax.set_ylabel("raw AUC 0–200 μm median")
    ax.set_title("Spatial front/middle/rear heterogeneity")
    ax.set_xticks(x, labels, rotation=45, ha="right", fontsize=8)
    ax.legend()
    fig.tight_layout()
    fig.savefig(figure_dir / "05_spatial_front_middle_rear.png", dpi=220)
    plt.close(fig)

    # 6. P95 denominator QC.  Keep P90 visible only as an audit curve.
    qc = p95_qc.loc[p95_qc["scope"].eq("all_frames")].set_index("scan_id").reindex(labels)
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.errorbar(x, qc["median"], yerr=qc["iqr"] / 2, fmt="o", capsize=3, label="P95 denominator median ± IQR/2")
    ax.plot(x, qc["p90"], "s--", label="P90 audit quantile", alpha=0.75)
    ax.set_ylabel("denominator (Flow units)")
    ax.set_title("P95 denominator QC (all 500 frames; P90 audit only)")
    ax.set_xticks(x, labels, rotation=45, ha="right", fontsize=8)
    ax.legend()
    fig.tight_layout()
    fig.savefig(figure_dir / "06_p95_denominator_qc.png", dpi=220)
    plt.close(fig)


def _gate_tables(scan_primary: pd.DataFrame, assess_sens: pd.DataFrame, guard_sens: pd.DataFrame, spatial: pd.DataFrame, p95_qc: pd.DataFrame) -> tuple[pd.DataFrame, str, dict[str, Any]]:
    criteria: list[dict[str, Any]] = []

    valid_fraction = scan_primary["primary_valid_fraction_of_assessable"].to_numpy(float)
    c1 = bool(np.isfinite(valid_fraction).all() and (valid_fraction >= VALID_FRACTION_FLOOR).all())
    criteria.append({"criterion_id": "C1_primary_valid_rate", "pass": c1, "observed": float(np.nanmin(valid_fraction)) if np.isfinite(valid_fraction).any() else np.nan, "threshold": f">={VALID_FRACTION_FLOOR:.2f} in every scan", "notes": "assessable-only raw AUC 0-200 valid fraction"})

    assess_changes = assess_sens["absolute_relative_change_pct"].to_numpy(float)
    assess_median = float(np.nanmedian(assess_changes)) if np.isfinite(assess_changes).any() else np.nan
    assess_max = float(np.nanmax(assess_changes)) if np.isfinite(assess_changes).any() else np.nan
    c2 = bool(np.isfinite(assess_changes).all() and assess_median <= SENS_MEDIAN_LIMIT_PCT and assess_max <= SENS_MATERIAL_LIMIT_PCT)
    criteria.append({"criterion_id": "C2_assessability_sensitivity", "pass": c2, "observed": f"median={assess_median:.3f}%; max={assess_max:.3f}%" if np.isfinite([assess_median, assess_max]).all() else "nonfinite", "threshold": f"median<={SENS_MEDIAN_LIMIT_PCT:.1f}%; max<={SENS_MATERIAL_LIMIT_PCT:.1f}%", "notes": "assessable-only vs assessable+uncertain"})

    guard0 = guard_sens["guard0_vs_guard2_abs_relative_change_pct"].to_numpy(float)
    guard4 = guard_sens["guard4_vs_guard2_abs_relative_change_pct"].to_numpy(float)
    g0_med = float(np.nanmedian(guard0)) if np.isfinite(guard0).any() else np.nan
    g4_med = float(np.nanmedian(guard4)) if np.isfinite(guard4).any() else np.nan
    g0_over = int(np.sum(guard0 > SENS_MATERIAL_LIMIT_PCT))
    g4_over = int(np.sum(guard4 > SENS_MATERIAL_LIMIT_PCT))
    c3 = bool(np.isfinite(guard0).all() and np.isfinite(guard4).all() and g0_med <= SENS_MEDIAN_LIMIT_PCT and g4_med <= SENS_MEDIAN_LIMIT_PCT and g0_over <= GUARD_MATERIAL_SCAN_COUNT_MAX and g4_over <= GUARD_MATERIAL_SCAN_COUNT_MAX)
    criteria.append({"criterion_id": "C3_guard_sensitivity", "pass": c3, "observed": f"guard0 median={g0_med:.3f}%, >15% n={g0_over}; guard4 median={g4_med:.3f}%, >15% n={g4_over}" if np.isfinite([g0_med, g4_med]).all() else "nonfinite", "threshold": f"each median<={SENS_MEDIAN_LIMIT_PCT:.1f}%; each >{SENS_MATERIAL_LIMIT_PCT:.1f}% count<={GUARD_MATERIAL_SCAN_COUNT_MAX}", "notes": "guard 2 px remains primary"})

    d235 = assess_sens.loc[assess_sens["diameter_um"].eq(235.0), "absolute_relative_change_pct"].to_numpy(float)
    d235_med = float(np.nanmedian(d235)) if np.isfinite(d235).any() else np.nan
    d235_max = float(np.nanmax(d235)) if np.isfinite(d235).any() else np.nan
    c4 = bool(np.isfinite(d235).all() and d235_med <= SENS_MEDIAN_LIMIT_PCT and d235_max <= SENS_MATERIAL_LIMIT_PCT)
    criteria.append({"criterion_id": "C4_D235_QC_selection", "pass": c4, "observed": f"median={d235_med:.3f}%; max={d235_max:.3f}%" if np.isfinite([d235_med, d235_max]).all() else "nonfinite", "threshold": f"D235 median<={SENS_MEDIAN_LIMIT_PCT:.1f}%; max<={SENS_MATERIAL_LIMIT_PCT:.1f}%", "notes": "QC sensitivity only; no biological interpretation"})

    spatial_flags: list[str] = []
    for scan_id, group in spatial.groupby("scan_id", sort=False):
        overall = float(scan_primary.loc[scan_primary["scan_id"].eq(scan_id), "raw_auc_0_200_primary_median"].iloc[0])
        region_vals = group["raw_auc_0_200_median"].to_numpy(float)
        if not np.isfinite(overall) or not np.isfinite(region_vals).all():
            spatial_flags.append(scan_id + ":missing_region")
            continue
        deviations = np.abs(region_vals - overall) / max(abs(overall), EPSILON) * 100.0
        if float(np.nanmax(deviations)) > SPATIAL_SEVERE_DEVIATION_PCT:
            spatial_flags.append(scan_id + f":max_dev={np.nanmax(deviations):.1f}%")
    c5 = len(spatial_flags) == 0
    criteria.append({"criterion_id": "C5_spatial_dominance", "pass": c5, "observed": "none" if c5 else "; ".join(spatial_flags), "threshold": f"no region deviation >{SPATIAL_SEVERE_DEVIATION_PCT:.1f}% from scan median and no missing region median", "notes": "severe spatial QC diagnostic, not a biological effect test"})

    all_p95 = p95_qc.loc[p95_qc["scope"].eq("all_frames")]
    p95_bad = all_p95.loc[(all_p95["valid_rate"] < P95_VALID_RATE_FLOOR) | (all_p95["robust_cv"] > P95_ROBUST_CV_LIMIT) | ~np.isfinite(all_p95["median"])]
    p95_flag = not p95_bad.empty
    criteria.append({"criterion_id": "P95_QC_secondary", "pass": not p95_flag, "observed": "pass" if not p95_flag else ", ".join(p95_bad["scan_id"].astype(str).tolist()), "threshold": f"all-frame valid_rate>={P95_VALID_RATE_FLOOR:.2f}, robust CV<={P95_ROBUST_CV_LIMIT:.2f}, finite median", "notes": "failure makes P95 exploratory only; it cannot overturn raw AUC"})

    raw_pass = bool(c1 and c2 and c3 and c4 and c5)
    gate = "A" if raw_pass and not p95_flag else "B" if raw_pass and p95_flag else "C"
    gate_reason = {
        "A": "raw AUC primary passed all predeclared stability/QC criteria; P95 denominator QC also passed",
        "B": "raw AUC primary passed; P95 denominator QC is sensitive/low-quality, so P95 remains exploratory",
        "C": "raw AUC primary failed at least one predeclared stability/QC criterion; do not enter 27 scans",
    }[gate]
    summary = {"gate": gate, "raw_gate_pass": raw_pass, "p95_qc_flag": p95_flag, "gate_reason": gate_reason, "spatial_flags": spatial_flags}
    criteria.append({"criterion_id": "FINAL_GATE", "pass": gate in ("A", "B"), "observed": gate, "threshold": "A/B permits later 27-scan entry; C blocks it", "notes": gate_reason})
    table = pd.DataFrame(criteria)
    return table, gate, summary


def _write_report(path: Path, *, metadata: dict[str, Any], scan_primary: pd.DataFrame, qc: pd.DataFrame, assess_sens: pd.DataFrame, guard_sens: pd.DataFrame, p95_qc: pd.DataFrame, spatial: pd.DataFrame, gate_table: pd.DataFrame, gate: str, gate_summary: dict[str, Any]) -> None:
    def fmt(value: Any) -> str:
        if isinstance(value, (float, np.floating)):
            return "NA" if not np.isfinite(value) else f"{float(value):.4g}"
        return str(value)

    lines = [
        "# FORMAL 9-SCAN FREEZE REPORT (2026-09-02)",
        "",
        f"**预声明 Gate：{gate}** — {gate_summary['gate_reason']}",
        "",
        "## 1. 本轮边界与冻结声明",
        "",
        "本轮只处理预先锁定的 9 个代表扫描、每卷全部 500 帧。未扩展到 27/45 卷；未进行流速/直径显著性检验，也未根据结果修改 X4、central 40% ROI、当前 z upper、assessability 规则、guard、窗口或归一化方案。500 帧不是 500 个独立实验样本。",
        "",
        "冻结方法：X4 横向定位；central 40% vessel ROI；当前 tracking 的 z_upper；guard=2 px 为主结果（0/4 px 仅敏感性）；raw AUC 0–200 μm 为主指标；primary 只纳入 predicted assessable；sensitivity 纳入 assessable+uncertain；not_assessable 排除且不补零/插值。",
        "",
        "## 2. 输入审计",
        "",
        f"- 配置文件：`{metadata['config_path']}`；SHA256：`{metadata['config_sha256']}`。",
        f"- 有效输出根：`{metadata['effective_output_root']}`（{metadata['output_root_mode']}）。",
        f"- 冻结全帧特征：`{metadata['features_path']}`；SHA256：`{metadata['features_sha256']}`；行数={metadata['feature_rows']}。",
        f"- 9 卷身份严格核验：{', '.join(FROZEN_SCAN_IDS)}。每卷 Flow 形状均应为 `(500, 351, 500)`，tracking alpha={PRIMARY_ALPHA:.2f} 且帧号 0–499。",
        "",
        "## 3. 扫描级主结果（描述性，仅 n=9 卷）",
        "",
        "scan-level 主值是该卷所有 valid + assessable 帧的 pooled median；IQR 来自同一帧集合，不使用 block median。",
        "",
        "| scan_id | D (μm) | v (mm/s) | assessable | uncertain | not_assessable | primary valid n | raw AUC median | IQR | P95-normalized median |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for _, row in scan_primary.iterrows():
        lines.append(
            "| " + " | ".join(
                [
                    str(row["scan_id"]),
                    fmt(row["diameter_um"]),
                    fmt(row["velocity_mm_s"]),
                    f"{float(row['assessable_rate']):.3f}",
                    f"{float(row['uncertain_rate']):.3f}",
                    f"{float(row['not_assessable_rate']):.3f}",
                    str(int(row["primary_valid_frame_count"])),
                    fmt(row["raw_auc_0_200_primary_median"]),
                    fmt(row["raw_auc_0_200_primary_iqr"]),
                    fmt(row["p95_normalized_auc_0_200_primary_median"]),
                ]
            )
            + " |"
        )
    lines += [
        "",
        "## 4. D185/D235/D285 QC 分层",
        "",
        "下表只描述预测可评估性和 QC 选择，不把直径解释为不可见的生物学原因。重点是 D235 的 assessable-only 与 assessable+uncertain 是否改变同一卷的 raw AUC。",
        "",
        "| diameter | frame n | assessable rate | uncertain rate | not_assessable rate |",
        "|---:|---:|---:|---:|---:|",
    ]
    diameter_summary = qc.loc[qc["level"].eq("diameter")]
    for _, row in diameter_summary.iterrows():
        lines.append(f"| {fmt(row['diameter_um'])} | {int(row['n_frames_total'])} | {float(row['assessable_rate']):.3f} | {float(row['uncertain_rate']):.3f} | {float(row['not_assessable_rate']):.3f} |")
    lines += [
        "",
        "## 5. 必要敏感性结果",
        "",
        "### Assessability",
        "",
        "| scan_id | assessable-only median | assessable+uncertain median | absolute relative change |",
        "|---|---:|---:|---:|",
    ]
    for _, row in assess_sens.iterrows():
        lines.append(f"| {row['scan_id']} | {fmt(row['primary_raw_auc_median'])} | {fmt(row['sensitivity_raw_auc_median'])} | {fmt(row['absolute_relative_change_pct'])}% |")
    lines += [
        "",
        "### Guard",
        "",
        "| scan_id | guard 0 median | guard 2 median (primary) | guard 4 median | guard0 vs 2 | guard4 vs 2 |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for _, row in guard_sens.iterrows():
        lines.append(f"| {row['scan_id']} | {fmt(row['guard0_primary_median'])} | {fmt(row['guard2_primary_median'])} | {fmt(row['guard4_primary_median'])} | {fmt(row['guard0_vs_guard2_abs_relative_change_pct'])}% | {fmt(row['guard4_vs_guard2_abs_relative_change_pct'])}% |")
    lines += [
        "",
        "### P95 denominator QC",
        "",
        "P95 denominator 只作为 relative tail-to-vessel 辅助指标的 QC。raw AUC 表示 absolute tail burden；如果 P95-normalized 与 raw 方向不同，不能据此判定其中一个错误。P90 仅作稳定性审计，不参与分母选择。",
        "",
        "| scan_id | scope | valid rate | median | IQR | robust CV | P5 | P90 audit | P95 |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for _, row in p95_qc.iterrows():
        lines.append(f"| {row['scan_id']} | {row['scope']} | {fmt(row['valid_rate'])} | {fmt(row['median'])} | {fmt(row['iqr'])} | {fmt(row['robust_cv'])} | {fmt(row['p5'])} | {fmt(row['p90'])} | {fmt(row['p95'])} |")
    lines += [
        "",
        "## 6. Front/middle/rear",
        "",
        "每个区域保留原始帧数和 assessability 分类；区域 raw AUC 仍使用 assessable-only 的 pooled median。此处用于发现空间异质性和单一区域主导风险，不进行生物学显著性推断。详表见 `07_spatial_front_middle_rear.csv`。",
        "",
        "## 7. Gate 判定（阈值预先声明）",
        "",
        f"- C1：每卷 assessable-only 的 0–200 μm raw AUC 有效比例 ≥ {VALID_FRACTION_FLOOR:.2f}。",
        f"- C2：assessable-only vs assessable+uncertain 的 scan-level raw AUC 绝对相对变化，跨卷 median ≤ {SENS_MEDIAN_LIMIT_PCT:.1f}%，且任何卷 ≤ {SENS_MATERIAL_LIMIT_PCT:.1f}%。",
        f"- C3：guard 0/4 相对 guard 2 的跨卷 median ≤ {SENS_MEDIAN_LIMIT_PCT:.1f}%；> {SENS_MATERIAL_LIMIT_PCT:.1f}% 的卷数不超过 {GUARD_MATERIAL_SCAN_COUNT_MAX}（2 px 始终为主结果）。",
        f"- C4：D235 子集同样满足 C2；只做 QC 选择敏感性判断，不做直径生物学解释。",
        f"- C5：任一 slow-axis 区域相对该卷主 median 偏离不超过 {SPATIAL_SEVERE_DEVIATION_PCT:.1f}%，且三个区域均有有效区域 median。",
        f"- P95 辅助 QC：all-frame valid rate ≥ {P95_VALID_RATE_FLOOR:.2f}、robust CV ≤ {P95_ROBUST_CV_LIMIT:.2f}。失败只触发 Gate B，不改变 raw Gate。",
        "",
        "| criterion | pass | observed | threshold |",
        "|---|---|---|---|",
    ]
    for _, row in gate_table.iterrows():
        lines.append(f"| {row['criterion_id']} | {row['pass']} | {row['observed']} | {row['threshold']} |")
    lines += [
        "",
        f"**最终 Gate = {gate}。** 本报告到此停止；不会自动运行 27 卷或 45 卷。Gate 结果也不是流速或直径的生物学结论。",
        "",
        "## 8. 输出文件",
        "",
        "- `01_frame_level_metrics.csv`：4500 帧，保留原始分类、定位、ROI、背景、raw/P95、0/2/4 px guard 和每个 invalid_reason。",
        "- `02_scan_level_primary.csv`：9 卷 pooled median/IQR 主汇总。",
        "- `03_qc_summary.csv`：9 卷及 D185/D235/D285 QC 分层。",
        "- `04_assessability_sensitivity.csv`、`05_guard_sensitivity.csv`、`06_p95_denominator_qc.csv`、`07_spatial_front_middle_rear.csv`、`08_gate_summary.csv`。",
        "- `figures/01_assessable_rate.png` 至 `figures/06_p95_denominator_qc.png`。",
        "",
        "## 9. 实验单位纪律",
        "",
        "当前每个条件只有一个代表扫描；帧级值只用于该扫描的 pooled 描述和 QC，不能把 500 帧当作 n=500 的独立重复。若 Gate A/B 允许后续扩展，27 卷阶段仍应以独立扫描/独立实验作为统计单位。",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _safe_replace_directory(building: Path, target: Path) -> None:
    resolved_target = target.resolve()
    allowed_parent = (PROJECT_ROOT / "work" / "pilot_20260414").resolve()
    if resolved_target.parent != allowed_parent or resolved_target.name != "formal_9scan_freeze_20260902":
        raise ValueError(f"Refusing to replace an unexpected output directory: {resolved_target}")
    if target.exists():
        shutil.rmtree(target)
    building.replace(target)


def main() -> int:
    args = parse_args()
    config_path = Path(args.config).resolve()
    cfg, _ = load_config(config_path)
    effective_cfg, effective_root, root_mode = _portable_config(cfg)
    manifest = load_manifest(effective_cfg, config_path)
    features_path = Path(args.features).resolve()
    frozen_features = _load_frozen_features(features_path)
    tracking_root = Path(args.tracking_root).resolve() if args.tracking_root else effective_root / str(effective_cfg["output"]["tracking_subdir"] if "tracking_subdir" in effective_cfg["output"] else "tracking")
    _validate_manifest_and_inputs(effective_cfg, config_path, manifest, tracking_root)

    output_dir = Path(args.output_dir).resolve()
    if output_dir.exists() and not args.overwrite:
        raise FileExistsError(f"Output exists; rerun with --overwrite: {output_dir}")
    building = output_dir.parent / f".{output_dir.name}.building-{uuid.uuid4().hex}"
    building.mkdir(parents=True, exist_ok=False)
    figure_dir = building / "figures"
    all_rows: list[dict[str, Any]] = []
    scan_rows: list[dict[str, Any]] = []
    spatial_rows: list[dict[str, Any]] = []
    indexed_manifest = manifest.set_index("dataset_id", drop=False)
    metrics_cfg = dict(effective_cfg["metrics"])
    metrics_cfg["central_fraction"] = PRIMARY_ROI_FRACTION

    try:
        for scan_order, scan_id in enumerate(FROZEN_SCAN_IDS, start=1):
            manifest_row = indexed_manifest.loc[scan_id]
            diameter_um = float(manifest_row["diameter_um"])
            velocity_mm_s = float(manifest_row["velocity_mm_s"])
            technical_replicate = int(manifest_row["technical_replicate"])
            flow_path = export_path(effective_cfg, scan_id)
            flow, source_frames, flow_metadata = load_flow(flow_path)
            if flow.shape != (FRAMES_PER_SCAN, N_Z, N_X):
                raise ValueError(f"{scan_id}: expected Flow shape (500,351,500), got {flow.shape}")
            if set(np.asarray(source_frames, dtype=int).tolist()) != set(range(FRAMES_PER_SCAN)):
                raise ValueError(f"{scan_id}: source frame mapping is not exactly 0..499")
            tracking_path = tracking_csv_path(tracking_root, scan_id, PRIMARY_ALPHA)
            # This audited normalizer is intentionally loaded dynamically so
            # the formal runner uses the exact confidence/localization schema
            # already used in the held-out audit.
            legacy_path = PROJECT_ROOT / "06_blinded_qc" / "generate_blinded_qc.py"
            legacy_spec = importlib.util.spec_from_file_location(f"legacy_norm_{scan_order}", legacy_path)
            if legacy_spec is None or legacy_spec.loader is None:
                raise ImportError(f"Cannot load tracking normalizer: {legacy_path}")
            legacy = importlib.util.module_from_spec(legacy_spec)
            legacy_spec.loader.exec_module(legacy)
            tracking = legacy.normalize_tracking_table(tracking_path, require_new_confidence=True).set_index("frame_index", drop=False)
            frozen = frozen_features.loc[frozen_features["dataset_id"].eq(scan_id)].set_index("frame_index_0based", drop=False)
            if len(tracking) != FRAMES_PER_SCAN or len(frozen) != FRAMES_PER_SCAN:
                raise ValueError(f"{scan_id}: tracking/features do not contain 500 frames")

            scan_frame_rows: list[dict[str, Any]] = []
            for frame_index in range(FRAMES_PER_SCAN):
                feature = frozen.loc[frame_index]
                track = tracking.loc[frame_index]
                x_x4 = _numeric(feature["x4_centroid_isolated_jump_corrected_px"])
                # The current z upper is the frozen alpha=0.15 tracking z;
                # feature z is cross-checked, never substituted silently.
                z_track = _numeric(track["z_upper_px"])
                z_feature = _numeric(feature["z_upper_algorithm_px"])
                if np.isfinite([z_track, z_feature]).all() and int(round(z_track)) != int(round(z_feature)):
                    raise ValueError(f"{scan_id} frame {frame_index}: frozen z mismatch between tracking and features")
                z_upper = z_track
                metric = _frame_metrics(flow[frame_index], x_x4=x_x4, z_upper=z_upper, diameter_um=diameter_um, metrics_cfg=metrics_cfg)
                row: dict[str, Any] = {
                    "scan_id": scan_id,
                    "dataset_id": scan_id,
                    "scan_order": scan_order,
                    "diameter_um": diameter_um,
                    "velocity_mm_s": velocity_mm_s,
                    "technical_replicate": technical_replicate,
                    "frame_index": frame_index,
                    "frame_index_0based": frame_index,
                    "frame_index_1based": frame_index + 1,
                    "slow_axis_region": str(feature["slow_axis_region"]),
                    "assessability_class": str(feature["vessel_presence_prediction"]),
                    "assessability_score": _numeric(feature["assessability_score"]),
                    "frozen_feature_x4_px": x_x4,
                    "frozen_feature_z_upper_px": z_feature,
                    "tracking_z_upper_px": z_track,
                    "alpha": PRIMARY_ALPHA,
                    "roi_central_fraction": PRIMARY_ROI_FRACTION,
                    "primary_guard_px": PRIMARY_GUARD,
                    "raw_auc_window_start_um": PRIMARY_WINDOW_UM[0],
                    "raw_auc_window_stop_um": PRIMARY_WINDOW_UM[1],
                    "axial_um_per_px": AXIAL_UM_PER_PX,
                    "lateral_um_per_px": LATERAL_UM_PER_PX,
                    "source_frame_index": int(source_frames[frame_index]),
                }
                row.update(metric)
                row["primary_frame_included"] = bool(row["assessability_class"] == "assessable" and row["raw_auc_0_200_guard2_valid"])
                row["sensitivity_frame_included"] = bool(row["assessability_class"] in ("assessable", "uncertain") and row["raw_auc_0_200_guard2_valid"])
                row["primary_exclusion_reason"] = "included" if row["primary_frame_included"] else ("not_assessable_or_uncertain" if row["assessability_class"] != "assessable" else row["invalid_reason"])
                row["sensitivity_exclusion_reason"] = "included" if row["sensitivity_frame_included"] else ("not_assessable" if row["assessability_class"] == "not_assessable" else row["invalid_reason"])
                scan_frame_rows.append(row)

            frame_table = pd.DataFrame(scan_frame_rows)
            frame_table["scan_id"] = scan_id
            all_rows.extend(scan_frame_rows)
            scan_rows.append(_scan_row(frame_table, scan_id, diameter_um, velocity_mm_s, technical_replicate))
            spatial_rows.extend(_spatial_rows(frame_table, scan_id, diameter_um, velocity_mm_s))

        frame_level = pd.DataFrame(all_rows)
        if len(frame_level) != len(FROZEN_SCAN_IDS) * FRAMES_PER_SCAN:
            raise AssertionError("Formal frame output does not contain exactly 4500 rows")
        scan_primary = pd.DataFrame(scan_rows)
        spatial = pd.DataFrame(spatial_rows)
        assess_sens = _assessability_sensitivity(frame_level)
        guard_sens = _guard_sensitivity(frame_level)
        p95_qc = _p95_qc(frame_level)

        # Diameter-level QC is pooled descriptively over frames; it is not a
        # condition-level inferential test and does not replace scan medians.
        qc_rows: list[dict[str, Any]] = []
        for _, row in scan_primary.iterrows():
            qc_rows.append({"level": "scan", **row.to_dict()})
        for diameter, group in frame_level.groupby("diameter_um", sort=True):
            counts = group["assessability_class"].value_counts()
            valid = group["raw_auc_0_200_valid"].astype(bool)
            qc_rows.append(
                {
                    "level": "diameter",
                    "scan_id": f"D{int(diameter)}_pooled_QC",
                    "dataset_id": f"D{int(diameter)}_pooled_QC",
                    "diameter_um": float(diameter),
                    "velocity_mm_s": np.nan,
                    "technical_replicate": np.nan,
                    "n_frames_total": int(len(group)),
                    "assessable_count": int(counts.get("assessable", 0)),
                    "uncertain_count": int(counts.get("uncertain", 0)),
                    "not_assessable_count": int(counts.get("not_assessable", 0)),
                    "assessable_rate": float((group["assessability_class"] == "assessable").mean()),
                    "uncertain_rate": float((group["assessability_class"] == "uncertain").mean()),
                    "not_assessable_rate": float((group["assessability_class"] == "not_assessable").mean()),
                    "primary_valid_frame_count": int((group["assessability_class"].eq("assessable") & valid).sum()),
                    "valid_frame_count": int((group["assessability_class"].eq("assessable") & valid).sum()),
                    "assessable_frame_count": int((group["assessability_class"] == "assessable").sum()),
                    "uncertain_frame_count": int((group["assessability_class"] == "uncertain").sum()),
                    "not_assessable_frame_count": int((group["assessability_class"] == "not_assessable").sum()),
                    "primary_valid_fraction_of_assessable": float((group["assessability_class"].eq("assessable") & valid).sum() / max(int((group["assessability_class"] == "assessable").sum()), 1)),
                    "all_class_valid_frame_count": int(valid.sum()),
                    "all_class_valid_rate": float(valid.mean()),
                    "aggregation_rule": "pooled descriptive QC across the three frozen scans; no inferential test",
                }
            )
        qc = pd.DataFrame(qc_rows)
        gate_table, gate, gate_summary = _gate_tables(scan_primary, assess_sens, guard_sens, spatial, p95_qc)

        metadata = {
            "schema_version": "1.0.0",
            "created_local": datetime.now().astimezone().isoformat(),
            "config_path": str(config_path),
            "config_sha256": _sha256(config_path),
            "effective_output_root": str(effective_root),
            "output_root_mode": root_mode,
            "features_path": str(features_path),
            "features_sha256": _sha256(features_path),
            "feature_rows": int(len(frozen_features)),
            "tracking_root": str(tracking_root),
            "frozen_scan_ids": list(FROZEN_SCAN_IDS),
            "frames_per_scan": FRAMES_PER_SCAN,
            "flow_shape_expected": [FRAMES_PER_SCAN, N_Z, N_X],
            "frozen_method": {
                "x": "X4_centroid_isolated_jump_corrected_px from frozen development table",
                "roi_fraction": PRIMARY_ROI_FRACTION,
                "z": "alpha=0.15 current tracking z_upper_px; feature z cross-checked",
                "primary_guard_px": PRIMARY_GUARD,
                "sensitivity_guards_px": list(GUARDS),
                "primary_raw_auc_window_um": list(PRIMARY_WINDOW_UM),
                "assessability": "primary=assessable only; sensitivity=assessable+uncertain; not_assessable excluded",
                "aggregation": "pooled median of valid assessable frames; no block medians",
            },
            "gate_operationalization": {
                "valid_fraction_floor": VALID_FRACTION_FLOOR,
                "sensitivity_median_limit_pct": SENS_MEDIAN_LIMIT_PCT,
                "sensitivity_material_limit_pct": SENS_MATERIAL_LIMIT_PCT,
                "guard_material_scan_count_max": GUARD_MATERIAL_SCAN_COUNT_MAX,
                "spatial_severe_deviation_pct": SPATIAL_SEVERE_DEVIATION_PCT,
                "p95_valid_rate_floor": P95_VALID_RATE_FLOOR,
                "p95_robust_cv_limit": P95_ROBUST_CV_LIMIT,
            },
            "gate": gate,
            "gate_summary": gate_summary,
            "no_biological_inference": True,
            "status": "formal_nine_scan_freeze_complete_stopped_before_27_scan_extension",
        }

        frame_level.to_csv(building / "01_frame_level_metrics.csv", index=False, encoding="utf-8-sig")
        scan_primary.to_csv(building / "02_scan_level_primary.csv", index=False, encoding="utf-8-sig")
        qc.to_csv(building / "03_qc_summary.csv", index=False, encoding="utf-8-sig")
        assess_sens.to_csv(building / "04_assessability_sensitivity.csv", index=False, encoding="utf-8-sig")
        guard_sens.to_csv(building / "05_guard_sensitivity.csv", index=False, encoding="utf-8-sig")
        p95_qc.to_csv(building / "06_p95_denominator_qc.csv", index=False, encoding="utf-8-sig")
        spatial.to_csv(building / "07_spatial_front_middle_rear.csv", index=False, encoding="utf-8-sig")
        gate_table.to_csv(building / "08_gate_summary.csv", index=False, encoding="utf-8-sig")
        (building / "run_metadata.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
        _make_figures(scan_primary, assess_sens, guard_sens, spatial, p95_qc, figure_dir)
        _write_report(building / "FORMAL_9SCAN_FREEZE_REPORT_20260902.md", metadata=metadata, scan_primary=scan_primary, qc=qc, assess_sens=assess_sens, guard_sens=guard_sens, p95_qc=p95_qc, spatial=spatial, gate_table=gate_table, gate=gate, gate_summary=gate_summary)
        _safe_replace_directory(building, output_dir)
    except Exception:
        shutil.rmtree(building, ignore_errors=True)
        raise
    print(json.dumps({"output_dir": str(output_dir), "gate": gate, "scan_count": len(FROZEN_SCAN_IDS), "frame_count": len(FROZEN_SCAN_IDS) * FRAMES_PER_SCAN}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
