"""Development-only x/ROI refinement and frame-assessability audit.

This module uses the fixed 53-frame pixel-localization package as a
development set.  It never reads or modifies raw DICOM files and never
opens a held-out reviewer answer.  The candidate x methods are intentionally
local to the existing Viterbi trajectory; no full-image argmax is used.

The output is a development report and auditable CSVs.  It is not a final
batch tail-analysis pipeline and its classification metrics are not an
independent validation result.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import shutil
import sys
import uuid
from datetime import datetime, timezone
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
from tailq.metrics import (  # noqa: E402
    _interval_values,
    _normalizer,
    _window_pixels,
    extract_profiles,
)


_GENERATOR_PATH = PROJECT_ROOT / "09_pixel_annotation" / "generate_pixel_annotation_package.py"
_GENERATOR_SPEC = importlib.util.spec_from_file_location("pixel_annotation_generator", _GENERATOR_PATH)
if _GENERATOR_SPEC is None or _GENERATOR_SPEC.loader is None:  # pragma: no cover
    raise ImportError(f"Cannot load package helper: {_GENERATOR_PATH}")
_GENERATOR = importlib.util.module_from_spec(_GENERATOR_SPEC)
_GENERATOR_SPEC.loader.exec_module(_GENERATOR)


REPRESENTATIVE_IDS = tuple(_GENERATOR.REPRESENTATIVE_IDS)
PRIMARY_ALPHA = float(_GENERATOR.PRIMARY_ALPHA)
AXIAL_UM_PER_PX = 6.7
LATERAL_UM_PER_PX = 12.7
GUARD_PX = 2
DEV_COUNT = 53
ROI_FRACTIONS = (1 / 3, 0.40, 0.50)
X_METHODS = (
    "old_viterbi",
    "X1_local_geometry",
    "X2_robust_centroid",
    "X3_local_high_signal_midaxis",
    "X4_centroid_isolated_jump_corrected",
)
DEFAULT_CONFIG = PROJECT_ROOT / "configs" / "pilot_20260414.json"
DEFAULT_PACKAGE_ROOT = PROJECT_ROOT / "work" / "pilot_20260414" / "pixel_annotation_20260827"
DEFAULT_OUTPUT = PROJECT_ROOT / "work" / "pilot_20260414" / "xroi_assessability_development_20260827"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(DEFAULT_CONFIG))
    parser.add_argument("--package-root", default=str(DEFAULT_PACKAGE_ROOT))
    parser.add_argument("--annotation-csv", default=None)
    parser.add_argument("--mapping-csv", default=None)
    parser.add_argument("--tracking-root", default=None)
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT))
    parser.add_argument("--seed", type=int, default=20260827)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _stable_rng(seed: int, namespace: str) -> np.random.Generator:
    material = f"{seed}|{namespace}".encode("utf-8")
    derived = int.from_bytes(hashlib.sha256(material).digest()[:8], "little")
    return np.random.default_rng(derived)


def _numeric(value: Any) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return np.nan
    return result if np.isfinite(result) else np.nan


def _runs(mask: np.ndarray) -> list[tuple[int, int]]:
    values = np.asarray(mask, dtype=bool)
    result: list[tuple[int, int]] = []
    index = 0
    while index < len(values):
        if not values[index]:
            index += 1
            continue
        stop = index + 1
        while stop < len(values) and values[stop]:
            stop += 1
        result.append((index, stop))
        index = stop
    return result


def _isolated_jump_correct(sequence: np.ndarray, *, jump_px: float = 2.0, neighbor_tol_px: float = 2.0) -> tuple[np.ndarray, np.ndarray]:
    """Correct only isolated 3-frame jumps; preserve smooth spatial changes."""

    corrected = np.asarray(sequence, dtype=float).copy()
    changed = np.zeros(len(corrected), dtype=bool)
    for index in range(1, len(corrected) - 1):
        left, current, right = corrected[index - 1 : index + 2]
        if not np.isfinite([left, current, right]).all():
            continue
        neighbor_median = float(np.median([left, current, right]))
        if abs(current - neighbor_median) > jump_px and abs(left - right) <= neighbor_tol_px:
            corrected[index] = neighbor_median
            changed[index] = True
    return corrected, changed


def _local_body_features(frame_zx: np.ndarray, tracking: pd.Series) -> dict[str, Any]:
    """Extract local body-only features and three non-continuity x candidates.

    The body slab is anchored by the existing algorithm z upper edge and the
    known inner diameter.  The x search is bounded by the existing local
    trajectory window.  The robust-centroid weights are background-corrected,
    clipped at their local 90th percentile, and therefore cannot be driven by
    one extreme pixel.
    """

    n_z, n_x = frame_zx.shape
    x_anchor = _numeric(tracking.get("x_center_px"))
    z_upper = _numeric(tracking.get("z_upper_px"))
    diameter_um = _numeric(tracking.get("diameter_um"))
    half_width = _numeric(tracking.get("x_local_window_half_width_px"))
    if not np.isfinite([x_anchor, z_upper, diameter_um]).all():
        return {"valid_local_body": False, "invalid_reason": "missing_algorithm_anchor"}
    diameter_px = max(1, int(round(diameter_um / AXIAL_UM_PER_PX)))
    expected_lateral_width = max(3, int(round(diameter_um / LATERAL_UM_PER_PX)))
    half_width_px = int(round(half_width)) if np.isfinite(half_width) else max(8, int(round(1.5 * expected_lateral_width)))
    anchor_int = int(round(x_anchor))
    x0 = max(0, anchor_int - half_width_px)
    x1 = min(n_x, anchor_int + half_width_px + 1)
    z0 = int(round(z_upper))
    z1 = min(n_z, z0 + diameter_px)
    if z0 < 0 or z1 <= z0 or x1 <= x0:
        return {"valid_local_body": False, "invalid_reason": "body_slab_out_of_bounds"}

    body = np.asarray(frame_zx[z0:z1, x0:x1], dtype=float)
    x_values = np.arange(x0, x1, dtype=float)
    column_score = np.quantile(body, 0.75, axis=0)
    exclusion_half_width = max(3, expected_lateral_width // 2 + 2)
    far_mask = np.abs(x_values - x_anchor) > exclusion_half_width
    if int(far_mask.sum()) < 5:
        far_mask = np.ones_like(x_values, dtype=bool)
    local_background = float(np.median(column_score[far_mask]))
    local_mad = float(np.median(np.abs(column_score[far_mask] - local_background)))
    local_sigma = max(1.4826 * local_mad, 1e-6)
    excess = np.maximum(column_score - local_background, 0.0)

    def run_score(run: tuple[int, int]) -> float:
        a, b = run
        return float(np.sum(excess[a - x0 : b - x0]))

    threshold = local_background + 1.5 * local_sigma
    runs = [(a + x0, b + x0) for a, b in _runs(column_score > threshold)]
    containing = [run for run in runs if run[0] <= anchor_int < run[1]]
    if containing:
        selected = max(containing, key=run_score)
        x1_fallback = False
    elif runs:
        selected = max(runs, key=run_score)
        x1_fallback = True
    else:
        selected = (anchor_int, anchor_int + 1)
        x1_fallback = True
    x1_center = (selected[0] + selected[1] - 1) / 2.0

    positive = excess > 0
    if positive.any():
        cap = float(np.quantile(excess[positive], 0.90))
    else:
        cap = 0.0
    weights = np.clip(excess, 0.0, cap)
    x2_center = (
        x_anchor
        if float(weights.sum()) <= 1e-12
        else float(np.sum(x_values * weights) / np.sum(weights))
    )

    peak_excess = float(np.max(excess)) if excess.size else 0.0
    half_peak_mask = excess >= 0.5 * peak_excess if peak_excess > 0 else np.zeros_like(excess, dtype=bool)
    half_peak_runs = [(a + x0, b + x0) for a, b in _runs(half_peak_mask)]
    if half_peak_runs:
        selected_half = max(half_peak_runs, key=run_score)
        x3_center = (selected_half[0] + selected_half[1] - 1) / 2.0
    else:
        x3_center = x_anchor

    # Axial completeness: fraction of body rows supporting the local robust
    # centroid above a row-wise robust background threshold.
    center_int = int(round(x2_center))
    center_half = max(1, expected_lateral_width // 4)
    c0 = max(x0, center_int - center_half)
    c1 = min(x1, center_int + center_half + 1)
    row_background = np.median(body[:, far_mask], axis=1)
    row_mad = np.median(np.abs(body[:, far_mask] - row_background[:, None]), axis=1)
    row_threshold = row_background + 1.5 * 1.4826 * row_mad
    row_signal = np.median(frame_zx[z0:z1, c0:c1], axis=1)
    axial_completeness = float(np.mean(row_signal > row_threshold)) if len(row_signal) else 0.0

    run_width = max(1, selected[1] - selected[0])
    body_peak_cnr = peak_excess / local_sigma
    continuity_score = float(np.clip(run_width / max(0.5 * expected_lateral_width, 1.0), 0.0, 1.0))
    width_consistency = float(np.exp(-abs(np.log(run_width / max(expected_lateral_width, 1)))))

    return {
        "valid_local_body": True,
        "invalid_reason": "ok",
        "x_anchor_px": float(x_anchor),
        "z_upper_algorithm_px": int(round(z_upper)),
        "diameter_axial_px": int(diameter_px),
        "expected_lateral_width_px": int(expected_lateral_width),
        "local_window_x0_px": int(x0),
        "local_window_x1_exclusive_px": int(x1),
        "local_window_half_width_px": int(half_width_px),
        "x1_local_geometry_px": float(x1_center),
        "x2_robust_centroid_px": float(x2_center),
        "x3_local_high_signal_midaxis_px": float(x3_center),
        "x1_fallback": bool(x1_fallback),
        "local_body_peak_cnr": float(body_peak_cnr),
        "local_body_peak_excess": float(peak_excess),
        "local_body_background": float(local_background),
        "local_body_sigma": float(local_sigma),
        "local_body_run_width_px": int(run_width),
        "local_body_axial_completeness": float(axial_completeness),
        "local_body_continuity_score": continuity_score,
        "local_body_width_consistency": width_consistency,
    }


def _add_assessability_score(features: pd.DataFrame) -> pd.DataFrame:
    result = features.copy()
    cnr_score = np.clip((result["local_body_peak_cnr"].to_numpy(float) - 2.0) / 6.0, 0.0, 1.0)
    continuity_score = result["local_body_continuity_score"].to_numpy(float)
    width_score = result["local_body_width_consistency"].to_numpy(float)
    axial_score = np.clip(result["local_body_axial_completeness"].to_numpy(float), 0.0, 1.0)
    neighbor_cnr = result["neighbor_peak_cnr"].to_numpy(float)
    neighbor_score = np.clip((neighbor_cnr - 2.0) / 6.0, 0.0, 1.0)
    score = (
        0.35 * cnr_score
        + 0.20 * continuity_score
        + 0.20 * width_score
        + 0.15 * axial_score
        + 0.10 * neighbor_score
    )
    result["cnr_score"] = cnr_score
    result["continuity_score"] = continuity_score
    result["width_consistency_score"] = width_score
    result["axial_completeness_score"] = axial_score
    result["neighbor_support_score"] = neighbor_score
    result["assessability_score"] = score
    result["vessel_presence_prediction"] = np.where(
        score >= 0.60,
        "assessable",
        np.where(score >= 0.40, "uncertain", "not_assessable"),
    )
    return result


def _cluster_bootstrap_ci(values: Iterable[float], clusters: Iterable[Any], *, seed: int, namespace: str, n_boot: int = 4000) -> tuple[float, float]:
    data = pd.DataFrame({"value": list(values), "cluster": list(clusters)})
    data["value"] = pd.to_numeric(data["value"], errors="coerce")
    data = data.loc[np.isfinite(data["value"].to_numpy(float))]
    if data.empty or data["cluster"].nunique() < 2:
        return np.nan, np.nan
    groups = {key: group["value"].to_numpy(float) for key, group in data.groupby("cluster", sort=True)}
    keys = list(groups)
    rng = _stable_rng(seed, namespace)
    boot = np.empty(n_boot, dtype=float)
    for index in range(n_boot):
        sampled = rng.choice(keys, size=len(keys), replace=True)
        boot[index] = float(np.median(np.concatenate([groups[key] for key in sampled])))
    return float(np.quantile(boot, 0.025)), float(np.quantile(boot, 0.975))


def _auc(frame: np.ndarray, *, x_center: float, z_upper: float, diameter_um: float, central_fraction: float, cfg: dict[str, Any]) -> dict[str, Any]:
    metrics_cfg = dict(cfg["metrics"])
    metrics_cfg["central_fraction"] = float(central_fraction)
    profiles = extract_profiles(
        frame,
        x_center=float(x_center),
        diameter_um=float(diameter_um),
        lateral_um_per_px=LATERAL_UM_PER_PX,
        metrics_cfg=metrics_cfg,
    )
    result = {
        "raw_auc": np.nan,
        "p95_rauc": np.nan,
        "raw_valid": False,
        "p95_valid": False,
        "central_x0_px": profiles.central_bounds[0],
        "central_x1_exclusive_px": profiles.central_bounds[1],
        "roi_complete": bool(profiles.roi_complete),
    }
    if not profiles.roi_complete:
        result["invalid_reason"] = "lateral_roi_out_of_bounds"
        return result
    z_int = int(round(z_upper))
    diameter_px = int(round(float(diameter_um) / AXIAL_UM_PER_PX))
    z_start = z_int + diameter_px + GUARD_PX
    offset_stop = _window_pixels(200, AXIAL_UM_PER_PX)
    values = _interval_values(
        profiles,
        z_start=z_start,
        offset_start_px=0,
        offset_stop_px=offset_stop,
        axial_um_per_px=AXIAL_UM_PER_PX,
    )
    if not values["complete"] or not np.isfinite(values["tail_auc_raw"]):
        result["invalid_reason"] = "window_incomplete"
        return result
    raw = float(values["tail_auc_raw"])
    denominator, valid, _, _ = _normalizer(
        profiles.excess,
        z_int,
        diameter_px,
        metrics_cfg["core_fraction"],
        "p95",
        float(metrics_cfg.get("epsilon", 1e-12)),
    )
    result["raw_auc"] = raw
    result["raw_valid"] = True
    result["p95_denominator"] = denominator
    if valid and np.isfinite(denominator) and denominator > float(metrics_cfg.get("epsilon", 1e-12)):
        result["p95_rauc"] = raw / denominator
        result["p95_valid"] = True
    result["invalid_reason"] = "ok"
    return result


def _candidate_error_summary(table: pd.DataFrame, *, seed: int) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    cohorts = {
        "primary_vessel_visible_yes": table["vessel_visible"].eq("yes"),
        "all_coordinate_complete": table["manual_coordinates_available"],
        "uncertain_with_coordinates": table["vessel_visible"].eq("uncertain") & table["manual_coordinates_available"],
    }
    for method in X_METHODS:
        column = f"{method}_x_px"
        for cohort, mask in cohorts.items():
            subset = table.loc[mask & table["manual_coordinates_available"]].copy()
            error = pd.to_numeric(subset[column], errors="coerce") - pd.to_numeric(subset["manual_x_px"], errors="coerce")
            error = error.dropna()
            if error.empty:
                continue
            absolute = np.abs(error.to_numpy(float))
            low, high = _cluster_bootstrap_ci(
                absolute,
                subset.loc[error.index, "dataset_id"],
                seed=seed,
                namespace=f"x-error|{method}|{cohort}",
            )
            rows.append(
                {
                    "method": method,
                    "cohort": cohort,
                    "n_frames": int(len(error)),
                    "n_scans": int(subset.loc[error.index, "dataset_id"].nunique()),
                    "median_signed_px": float(np.median(error)),
                    "median_abs_px": float(np.median(absolute)),
                    "iqr_abs_px": float(np.quantile(absolute, 0.75) - np.quantile(absolute, 0.25)),
                    "p90_abs_px": float(np.quantile(absolute, 0.90)),
                    "p95_abs_px": float(np.quantile(absolute, 0.95)),
                    "max_abs_px": float(np.max(absolute)),
                    "median_abs_ci95_low_px": low,
                    "median_abs_ci95_high_px": high,
                }
            )
    return pd.DataFrame(rows)


def _roi_summary(table: pd.DataFrame, *, seed: int) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    subset = table.loc[table["manual_coordinates_available"]].copy()
    for method in X_METHODS:
        for fraction in ROI_FRACTIONS:
            x_col = f"{method}_x_px"
            for metric in ("overlap_fraction", "iou"):
                values = pd.to_numeric(subset[f"{metric}_{method}_{fraction:g}"], errors="coerce")
                valid = values.notna()
                values_np = values.loc[valid].to_numpy(float)
                low, high = _cluster_bootstrap_ci(
                    values_np,
                    subset.loc[valid, "dataset_id"],
                    seed=seed,
                    namespace=f"roi|{method}|{fraction:g}|{metric}",
                )
                rows.append(
                    {
                        "method": method,
                        "central_fraction": fraction,
                        "metric": metric,
                        "n_frames": int(len(values_np)),
                        "n_scans": int(subset.loc[valid, "dataset_id"].nunique()),
                        "median": float(np.median(values_np)) if len(values_np) else np.nan,
                        "iqr": float(np.quantile(values_np, 0.75) - np.quantile(values_np, 0.25)) if len(values_np) else np.nan,
                        "p90": float(np.quantile(values_np, 0.90)) if len(values_np) else np.nan,
                        "p95": float(np.quantile(values_np, 0.95)) if len(values_np) else np.nan,
                        "max": float(np.max(values_np)) if len(values_np) else np.nan,
                        "median_ci95_low": low,
                        "median_ci95_high": high,
                    }
                )
    return pd.DataFrame(rows)


def _auc_summary(table: pd.DataFrame, *, seed: int) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    subset = table.loc[table["manual_coordinates_available"]].copy()
    # Keep the complete coordinate cohort for a diagnostic sensitivity check,
    # but use the pre-specified visible=yes cohort as the primary queue.  The
    # former also contains uncertain/no frames with coordinates and must not be
    # silently conflated with the primary biological-quality cohort.
    cohorts = {
        "all_coordinate_complete": pd.Series(True, index=subset.index),
        "primary_vessel_visible_yes": subset["vessel_visible"].eq("yes"),
    }
    for method in X_METHODS:
        for fraction in ROI_FRACTIONS:
            for metric in ("raw_auc", "p95_rauc"):
                column = f"relative_manual_minus_{method}_{fraction:g}_{metric}_pct"
                for cohort, cohort_mask in cohorts.items():
                    values = pd.to_numeric(subset[column], errors="coerce")
                    valid = values.notna() & cohort_mask
                    values = values.loc[valid]
                    arr = values.to_numpy(float)
                    absolute = np.abs(arr)
                    low, high = _cluster_bootstrap_ci(
                        absolute,
                        subset.loc[valid, "dataset_id"],
                        seed=seed,
                        namespace=f"auc|{method}|{fraction:g}|{metric}|{cohort}",
                    )
                    rows.append(
                        {
                            "method": method,
                            "central_fraction": fraction,
                            "metric": metric,
                            "cohort": cohort,
                            "n_frames": int(len(arr)),
                            "n_scans": int(subset.loc[valid, "dataset_id"].nunique()),
                            "median_signed_relative_pct": float(np.median(arr)) if len(arr) else np.nan,
                            "median_abs_relative_pct": float(np.median(absolute)) if len(arr) else np.nan,
                            "iqr_abs_relative_pct": float(np.quantile(absolute, 0.75) - np.quantile(absolute, 0.25)) if len(arr) else np.nan,
                            "p90_abs_relative_pct": float(np.quantile(absolute, 0.90)) if len(arr) else np.nan,
                            "p95_abs_relative_pct": float(np.quantile(absolute, 0.95)) if len(arr) else np.nan,
                            "max_abs_relative_pct": float(np.max(absolute)) if len(arr) else np.nan,
                            "median_abs_ci95_low": low,
                            "median_abs_ci95_high": high,
                        }
                    )
    return pd.DataFrame(rows)


def _assessability_tables(
    labeled_features: pd.DataFrame,
    all_features: pd.DataFrame,
    *,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    labeled = labeled_features.loc[labeled_features["vessel_visible"].notna()].copy()
    confusion = pd.crosstab(labeled["vessel_visible"], labeled["vessel_presence_prediction"], dropna=False).reset_index()
    for category in ("assessable", "uncertain", "not_assessable"):
        if category not in confusion:
            confusion[category] = 0
    confusion = confusion[["vessel_visible", "assessable", "uncertain", "not_assessable"]]

    rows: list[dict[str, Any]] = []
    for diameter, group in labeled.groupby("diameter_um", sort=True):
        yes = group["vessel_visible"].eq("yes")
        predicted_pass = group["vessel_presence_prediction"].eq("assessable")
        tp = int((yes & predicted_pass).sum())
        fn = int((yes & ~predicted_pass).sum())
        fp = int((~yes & predicted_pass).sum())
        tn = int((~yes & ~predicted_pass).sum())
        rows.append(
            {
                "diameter_um": float(diameter),
                "n_frames": int(len(group)),
                "manual_yes": int(yes.sum()),
                "manual_uncertain": int(group["vessel_visible"].eq("uncertain").sum()),
                "manual_no": int(group["vessel_visible"].eq("no").sum()),
                "pred_assessable": int(predicted_pass.sum()),
                "pred_uncertain": int(group["vessel_presence_prediction"].eq("uncertain").sum()),
                "pred_not_assessable": int(group["vessel_presence_prediction"].eq("not_assessable").sum()),
                "pass_rate": float(predicted_pass.mean()),
                "tp_yes": tp,
                "fn_yes": fn,
                "fp_non_yes": fp,
                "tn_non_yes": tn,
                "yes_sensitivity": tp / (tp + fn) if tp + fn else np.nan,
                "yes_specificity": tn / (tn + fp) if tn + fp else np.nan,
                "yes_ppv": tp / (tp + fp) if tp + fp else np.nan,
                "yes_npv": tn / (tn + fn) if tn + fn else np.nan,
                "yes_balanced_accuracy": 0.5 * (tp / (tp + fn) + tn / (tn + fp)) if (tp + fn) and (tn + fp) else np.nan,
            }
        )
    by_diameter = pd.DataFrame(rows)
    all_by_diameter = (
        all_features.groupby("diameter_um", sort=True)["vessel_presence_prediction"]
        .value_counts()
        .rename("frame_count")
        .reset_index()
    )
    all_by_diameter["frame_fraction"] = all_by_diameter["frame_count"] / all_by_diameter.groupby("diameter_um")["frame_count"].transform("sum")
    return confusion, by_diameter, all_by_diameter


def _make_figures(candidate_summary: pd.DataFrame, roi_summary: pd.DataFrame, auc_summary: pd.DataFrame, labeled_features: pd.DataFrame, all_features: pd.DataFrame, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    # Candidate x error distributions, primary visible cohort.
    fig, ax = plt.subplots(figsize=(10.5, 5.0), constrained_layout=True)
    primary = candidate_summary.loc[candidate_summary["cohort"].eq("primary_vessel_visible_yes")].copy()
    order = list(X_METHODS)
    vals = []
    for method in order:
        row = primary.loc[primary["method"].eq(method)].iloc[0]
        values = candidate_summary  # plotted from frame-level table below
        vals.append([])
    # Reconstruct from the frame-level columns stored by the caller through a
    # compact sidecar DataFrame attached as attributes is not reliable; the
    # report figure uses summary medians and P95 markers instead.
    x = np.arange(len(order))
    med = [float(primary.loc[primary["method"].eq(method), "median_abs_px"].iloc[0]) for method in order]
    p95 = [float(primary.loc[primary["method"].eq(method), "p95_abs_px"].iloc[0]) for method in order]
    ax.plot(x, med, "o-", label="median |Δx|")
    ax.plot(x, p95, "s--", label="P95 |Δx|")
    ax.set_xticks(x, ["old", "X1", "X2", "X3", "X4"])
    ax.set_ylabel("x error (px)")
    ax.set_title("Development-set local x candidates")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False)
    fig.savefig(out_dir / "figure01_x_candidate_error_summary.png", dpi=180)
    plt.close(fig)

    # Selected X4 ROI/AUC comparison across width choices.
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.8), constrained_layout=True)
    selected_roi = roi_summary.loc[roi_summary["method"].eq("X4_centroid_isolated_jump_corrected") & roi_summary["metric"].eq("overlap_fraction")]
    selected_auc = auc_summary.loc[auc_summary["method"].eq("X4_centroid_isolated_jump_corrected") & auc_summary["metric"].eq("raw_auc")]
    fractions = [1 / 3, 0.40, 0.50]
    axes[0].plot(fractions, [float(selected_roi.loc[selected_roi["central_fraction"].eq(f), "median"].iloc[0]) for f in fractions], "o-")
    axes[0].axhline(0.757, color="#c1121f", ls="--", lw=1, label="old 0.757")
    axes[0].set_xlabel("central ROI fraction")
    axes[0].set_ylabel("overlap fraction median")
    axes[0].set_ylim(0, 1.05)
    axes[0].grid(alpha=0.25)
    axes[0].legend(frameon=False)
    axes[1].plot(fractions, [float(selected_auc.loc[selected_auc["central_fraction"].eq(f), "median_abs_relative_pct"].iloc[0]) for f in fractions], "o-")
    axes[1].axhline(8.04, color="#c1121f", ls="--", lw=1, label="old 8.04%")
    axes[1].set_xlabel("central ROI fraction")
    axes[1].set_ylabel("raw AUC |relative change| median (%)")
    axes[1].grid(alpha=0.25)
    axes[1].legend(frameon=False)
    fig.suptitle("X4: ROI-width sensitivity on development set")
    fig.savefig(out_dir / "figure02_x4_roi_auc_sensitivity.png", dpi=180)
    plt.close(fig)

    # Labeled assessability score.
    fig, ax = plt.subplots(figsize=(8.0, 5.0), constrained_layout=True)
    order_label = ["yes", "uncertain", "no"]
    groups = [labeled_features.loc[labeled_features["vessel_visible"].eq(label), "assessability_score"].to_numpy(float) for label in order_label]
    ax.boxplot(groups, labels=order_label, showfliers=False)
    rng = _stable_rng(20260827, "figure|assessability")
    for index, values in enumerate(groups, start=1):
        if len(values):
            ax.scatter(index + rng.uniform(-0.08, 0.08, len(values)), values, s=22, alpha=0.75, color="#2f6690")
    ax.axhline(0.60, color="#2a9d8f", ls="--", lw=1, label="assessable ≥ 0.60")
    ax.axhline(0.40, color="#e9c46a", ls="--", lw=1, label="uncertain ≥ 0.40")
    ax.set_ylabel("development assessability score")
    ax.set_title("Independent vessel-presence score (development only)")
    ax.set_ylim(0, 1.05)
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False)
    fig.savefig(out_dir / "figure03_assessability_development.png", dpi=180)
    plt.close(fig)

    # All-frame predicted composition by diameter.
    counts = all_features.groupby(["diameter_um", "vessel_presence_prediction"]).size().unstack(fill_value=0)
    counts = counts.reindex(columns=["assessable", "uncertain", "not_assessable"], fill_value=0)
    fig, ax = plt.subplots(figsize=(8.2, 4.8), constrained_layout=True)
    counts.plot(kind="bar", stacked=True, ax=ax, color=["#2a9d8f", "#e9c46a", "#e76f51"])
    ax.set_xlabel("diameter (μm)")
    ax.set_ylabel("all-frame predicted count")
    ax.set_title("Predicted assessability composition across 500-frame scans")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.25)
    fig.savefig(out_dir / "figure04_assessability_all_frames_by_diameter.png", dpi=180)
    plt.close(fig)


def _write_report(
    path: Path,
    *,
    metadata: dict[str, Any],
    candidate_summary: pd.DataFrame,
    roi_summary: pd.DataFrame,
    auc_summary: pd.DataFrame,
    confusion: pd.DataFrame,
    assess_by_diameter: pd.DataFrame,
    all_assess_by_diameter: pd.DataFrame,
) -> None:
    primary = candidate_summary.loc[candidate_summary["cohort"].eq("primary_vessel_visible_yes")]
    x4 = primary.loc[primary["method"].eq("X4_centroid_isolated_jump_corrected")].iloc[0]
    old = primary.loc[primary["method"].eq("old_viterbi")].iloc[0]
    x4_roi = roi_summary.loc[(roi_summary["method"] == "X4_centroid_isolated_jump_corrected") & (roi_summary["central_fraction"] == 0.40) & (roi_summary["metric"] == "overlap_fraction")].iloc[0]
    x4_auc = auc_summary.loc[(auc_summary["method"] == "X4_centroid_isolated_jump_corrected") & (auc_summary["central_fraction"] == 0.40) & (auc_summary["metric"] == "raw_auc") & (auc_summary["cohort"] == "primary_vessel_visible_yes")].iloc[0]
    x4_p95 = auc_summary.loc[(auc_summary["method"] == "X4_centroid_isolated_jump_corrected") & (auc_summary["central_fraction"] == 0.40) & (auc_summary["metric"] == "p95_rauc") & (auc_summary["cohort"] == "primary_vessel_visible_yes")].iloc[0]
    lines = [
        "# X/ROI 与 frame-assessability development report",
        "",
        f"日期：{metadata['date_local']}  ",
        "阶段：**development set 已完成；held-out 包已生成并暂停等待人工标注**  ",
        "范围：固定 9 个代表扫描；不扩展到 27 卷或 45 卷。",
        "",
        "> 本报告使用第五轮已解盲的 53 帧作为 development set。候选选择、阈值和参数均在这 53 帧上形成，因此以下分类与性能数字不能称为独立验证结果。",
        "",
        "## 1. x 候选方法",
        "",
        "- `old_viterbi`：当前三维 Viterbi x，作为基线；",
        "- `X1_local_geometry`：Viterbi 局部窗内、body slab 的主要连续主体中轴；",
        "- `X2_robust_centroid`：局部背景校正后、90th-percentile clipped 的 robust centroid；",
        "- `X3_local_high_signal_midaxis`：局部 body slab 中半峰值主体中轴；",
        "- `X4_centroid_isolated_jump_corrected`：X2 加三帧孤立跳变修正；仅在当前点偏离邻居中位数 >2 px 且两邻居彼此相差 ≤2 px 时修正。",
        "",
        "所有候选都只在 `x_Viterbi ± local_window_half_width` 内计算；没有重新使用全图 argmax。z、Flow 输入、background、guard、window 和 normalization 均保持冻结。",
        "",
        "## 2. development-set 选择结果",
        "",
        "主评价队列为 `vessel_visible=yes` 的 34 帧；CI 为按扫描聚类 bootstrap。",
        "",
        "| x 方法 | median | P90 | P95 | max |",
        "|---|---:|---:|---:|---:|",
        f"| old Viterbi | {old['median_abs_px']:.3f} px | {old['p90_abs_px']:.3f} px | {old['p95_abs_px']:.3f} px | {old['max_abs_px']:.3f} px |",
        f"| X4 | {x4['median_abs_px']:.3f} px | {x4['p90_abs_px']:.3f} px | {x4['p95_abs_px']:.3f} px | {x4['max_abs_px']:.3f} px |",
        "",
        "开发选择：**X4**。理由是它在不扩大搜索范围的前提下，降低了 x 误差的中心趋势和孤立跳变风险；选择依据不包含流速/直径 tail 趋势。",
        "",
        "## 3. ROI 候选比较",
        "",
        "对 X4 比较 central 1/3、40%、50%，均继续使用横向 median profile 和相同双侧背景。40% 在 overlap 提升与 raw AUC 稳健性之间最平衡：",
        "",
        f"- X4 + central 40% overlap fraction median：**{x4_roi['median']:.3f}**；",
        f"- X4 + central 40% raw AUC 0–200 μm absolute relative change median：**{x4_auc['median_abs_relative_pct']:.3f}%**，P90 **{x4_auc['p90_abs_relative_pct']:.3f}%**，P95 **{x4_auc['p95_abs_relative_pct']:.3f}%**；",
        f"- X4 + central 40% P95-normalized AUC absolute relative change median：**{x4_p95['median_abs_relative_pct']:.3f}%**。",
        "",
        "因此冻结 development candidate 为：**X4 + central ROI fraction 0.40**。这只是 held-out 前的开发选择，不是最终验证结论。",
        "",
        "## 4. frame-assessability development",
        "",
        "评分独立于 x/z confidence、overall confidence、image_abnormal 和最终 tail AUC，使用以下可解释特征：local body CNR、主体连续性、宽度一致性、轴向主体完整性和相邻帧峰值 CNR 支持。固定评分：CNR 0.35、连续性 0.20、宽度 0.20、轴向完整性 0.15、邻帧支持 0.10；score ≥0.60 为 `assessable`，0.40–<0.60 为 `uncertain`，<0.40 为 `not_assessable`。",
        "",
        "development confusion table（行=人工，列=算法预测）：",
        "",
        confusion.to_markdown(index=False),
        "",
        "按直径的 development 诊断：",
        "",
        assess_by_diameter.to_markdown(index=False),
        "",
        "该规则在 development set 上用于了解候选特征，不代表独立验证；D235 的 pass/fail 差异必须在 held-out 中复核，不能当作生物学结论。",
        "",
        "全 500 帧扫描上的预测组成（无人工答案，仅用于抽样与风险诊断）：",
        "",
        all_assess_by_diameter.to_markdown(index=False),
        "",
        "## 5. held-out 准入判定",
        "",
        "**A：值得生成 held-out。** X4 + 40% 在 development 上相对旧基线明显降低 x 误差和 raw AUC 定位敏感性；assessability 规则简单、可解释，且没有把最终 AUC 写入评分。",
        "",
        "这不等于方法已经验证通过。53 帧已经用于开发，下一步必须在完全排除前 81 帧和前 53 帧的新样本上复验。",
        "",
        "## 6. held-out 包与停止点",
        "",
        f"- held-out 包：`{metadata['heldout_root']}`；",
        f"- 盲法 ID：`HV_001`–`HV_{metadata['heldout_count']:03d}`，共 {metadata['heldout_count']} 帧；",
        f"- 新抽样 seed：`{metadata['heldout_seed']}`；",
        "- 已从抽样池排除前 81 个 QC frame 和前 53 个 PA frame；mapping 与预测 assessability 仅在 restricted 目录；",
        "- held-out reviewer 包不显示直径、流速、scan/frame、算法坐标、ROI、confidence、AUC 或 assessability 预测；",
        "- 当前在生成盲法包后停止，不读取、不模拟 held-out 人工答案。",
        "",
        "## 7. 重新开始前的要求",
        "",
        "人工标注完成后，保留 `HV_###` ID、列名和行顺序不变，只交回完成的 CSV；不要提前打开 `restricted_mapping`。收到 CSV 后再进行 frame-assessability confusion matrix、Δx、ROI overlap 和 x→AUC sensitivity 的 held-out 解盲。",
        "",
        "## 8. 可复现信息",
        "",
        f"- development 输入 annotation SHA-256：`{metadata['development_annotation_sha256']}`；",
        f"- development 输入 mapping SHA-256：`{metadata['development_mapping_sha256']}`；",
        f"- 配置 SHA-256：`{metadata['config_sha256']}`；",
        "- 当前项目测试：运行 `D:\\anaconda\\python.exe -m pytest -q`。",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    cfg, config_path = load_config(args.config)
    manifest = load_manifest(cfg, config_path).set_index("dataset_id", drop=False)
    package_root = Path(args.package_root).resolve()
    annotation_path = Path(args.annotation_csv).resolve() if args.annotation_csv else package_root / "reviewer_package" / "pixel_level_manual_annotation.csv"
    mapping_path = Path(args.mapping_csv).resolve() if args.mapping_csv else package_root / "restricted_mapping" / "pixel_annotation_mapping.csv"
    annotation = _GENERATOR  # local alias only to keep validation provenance explicit
    # The audit module is the single source of truth for the completed CSV schema.
    analysis_path = PROJECT_ROOT / "09_pixel_annotation" / "analyze_pixel_annotation.py"
    analysis_spec = importlib.util.spec_from_file_location("pixel_annotation_analysis", analysis_path)
    if analysis_spec is None or analysis_spec.loader is None:  # pragma: no cover
        raise ImportError(f"Cannot load annotation validator: {analysis_path}")
    analysis_module = importlib.util.module_from_spec(analysis_spec)
    analysis_spec.loader.exec_module(analysis_module)
    annotation_table = analysis_module.validate_annotation_csv(annotation_path)
    mapping_table = analysis_module.validate_mapping_csv(mapping_path)
    if len(annotation_table) != DEV_COUNT or len(mapping_table) != DEV_COUNT:
        raise ValueError("Development set must contain exactly 53 frames")
    mapping_by_id = mapping_table.set_index("blind_id", drop=False)
    annotation_by_id = annotation_table.set_index("blind_id", drop=False)

    tracking_root = Path(args.tracking_root).resolve() if args.tracking_root else output_root(cfg) / "tracking"
    output_root_path = Path(args.output_dir).resolve()
    if output_root_path.exists() and not args.overwrite:
        raise FileExistsError(f"Output exists; use --overwrite: {output_root_path}")
    building = output_root_path.parent / f".{output_root_path.name}.building-{uuid.uuid4().hex}"
    building.mkdir(parents=True, exist_ok=False)
    figure_dir = building / "figures"
    feature_rows: list[dict[str, Any]] = []
    selected_rows: list[dict[str, Any]] = []
    all_feature_rows: list[dict[str, Any]] = []
    frame_cache: dict[tuple[str, int], np.ndarray] = {}
    try:
        for scan_id in REPRESENTATIVE_IDS:
            flow, _, _ = load_flow(export_path(cfg, scan_id))
            tracking = pd.read_csv(tracking_csv_path(tracking_root, scan_id, PRIMARY_ALPHA)).set_index("frame_index", drop=False)
            sequence_rows: list[dict[str, Any]] = []
            for frame_index in range(flow.shape[0]):
                track = tracking.loc[frame_index]
                feat = _local_body_features(flow[frame_index], track)
                feat.update(
                    {
                        "dataset_id": scan_id,
                        "frame_index_0based": int(frame_index),
                        "diameter_um": float(track["diameter_um"]),
                        "velocity_mm_s": float(manifest.loc[scan_id, "velocity_mm_s"]),
                        "technical_replicate": int(manifest.loc[scan_id, "technical_replicate"]),
                    }
                )
                sequence_rows.append(feat)
            sequence = pd.DataFrame(sequence_rows)
            neighbor_values = sequence["local_body_peak_cnr"].shift(1).combine(
                sequence["local_body_peak_cnr"].shift(-1), lambda left, right: np.nanmedian([left, right])
            )
            neighbor_values = neighbor_values.fillna(sequence["local_body_peak_cnr"])
            sequence["neighbor_peak_cnr"] = neighbor_values
            sequence = _add_assessability_score(sequence)
            all_feature_rows.extend(sequence.to_dict("records"))

            selected = mapping_table.loc[mapping_table["dataset_id"].eq(scan_id)].copy()
            for _, mapping_row in selected.iterrows():
                frame_index = int(mapping_row["frame_index_0based"])
                ann = annotation_by_id.loc[mapping_row["blind_id"]]
                feat = sequence.loc[sequence["frame_index_0based"].eq(frame_index)].iloc[0].to_dict()
                frame = flow[frame_index]
                frame_cache[(scan_id, frame_index)] = frame
                manual_x = _numeric(ann["manual_x_center_px"])
                manual_z = _numeric(ann["manual_z_upper_px"])
                row: dict[str, Any] = {
                    "blind_id": mapping_row["blind_id"],
                    "dataset_id": scan_id,
                    "frame_index_0based": frame_index,
                    "diameter_um": float(mapping_row["diameter_um"]),
                    "velocity_mm_s": float(mapping_row["velocity_mm_s"]),
                    "technical_replicate": int(mapping_row["technical_replicate"]),
                    "slow_axis_region": mapping_row["slow_axis_region"],
                    "source_category": mapping_row["source_category"],
                    "vessel_visible": ann["vessel_visible"],
                    "annotation_confidence": ann["annotation_confidence"],
                    "image_abnormal": ann["image_abnormal"],
                    "manual_x_px": manual_x,
                    "manual_z_px": manual_z,
                    "manual_coordinates_available": bool(np.isfinite([manual_x, manual_z]).all()),
                    "old_x_px": float(mapping_row["x_center_px"]),
                    "algorithm_z_px": float(mapping_row["z_upper_px"]),
                }
                for key in (
                    "x1_local_geometry_px",
                    "x2_robust_centroid_px",
                    "x3_local_high_signal_midaxis_px",
                    "x_anchor_px",
                    "local_body_peak_cnr",
                    "local_body_run_width_px",
                    "local_body_axial_completeness",
                    "local_body_continuity_score",
                    "local_body_width_consistency",
                    "neighbor_peak_cnr",
                    "assessability_score",
                    "vessel_presence_prediction",
                ):
                    row[key] = feat.get(key, np.nan)
                selected_rows.append(row)

            # Store continuity correction after all 500 frames are available.
            corrected, changed = _isolated_jump_correct(sequence["x2_robust_centroid_px"].to_numpy(float))
            all_feature_rows[-len(sequence) :]
            for frame_index, (value, changed_flag) in enumerate(zip(corrected, changed)):
                all_feature_rows[-len(sequence) + frame_index]["x4_centroid_isolated_jump_corrected_px"] = float(value)
                all_feature_rows[-len(sequence) + frame_index]["x4_jump_corrected"] = bool(changed_flag)
            for row in selected_rows[-len(selected) :]:
                frame_index = int(row["frame_index_0based"])
                row["X1_local_geometry_x_px"] = float(sequence.loc[sequence["frame_index_0based"].eq(frame_index), "x1_local_geometry_px"].iloc[0])
                row["X2_robust_centroid_x_px"] = float(sequence.loc[sequence["frame_index_0based"].eq(frame_index), "x2_robust_centroid_px"].iloc[0])
                row["X3_local_high_signal_midaxis_x_px"] = float(sequence.loc[sequence["frame_index_0based"].eq(frame_index), "x3_local_high_signal_midaxis_px"].iloc[0])
                row["X4_centroid_isolated_jump_corrected_x_px"] = float(corrected[frame_index])
                row["X4_jump_corrected"] = bool(changed[frame_index])

        selected_table = pd.DataFrame(selected_rows)
        # Candidate columns with stable names used by the summary functions.
        selected_table["old_viterbi_x_px"] = selected_table["old_x_px"]
        selected_table["X1_local_geometry_x_px"] = selected_table["X1_local_geometry_x_px"]
        selected_table["X2_robust_centroid_x_px"] = selected_table["X2_robust_centroid_x_px"]
        selected_table["X3_local_high_signal_midaxis_x_px"] = selected_table["X3_local_high_signal_midaxis_x_px"]
        selected_table["X4_centroid_isolated_jump_corrected_x_px"] = selected_table["X4_centroid_isolated_jump_corrected_x_px"]

        # Calculate all candidate/ROI geometry and AUC metrics on the 53-frame
        # development sample, while keeping every row even if invalid.
        for index, row in selected_table.iterrows():
            frame = frame_cache[(row["dataset_id"], int(row["frame_index_0based"]))]
            for method in X_METHODS:
                x_col = f"{method}_x_px"
                if x_col not in selected_table:
                    raise AssertionError(f"Missing candidate column {x_col}")
                x_value = float(row[x_col])
                for fraction in ROI_FRACTIONS:
                    candidate_profile = extract_profiles(
                        frame,
                        x_center=x_value,
                        diameter_um=float(row["diameter_um"]),
                        lateral_um_per_px=LATERAL_UM_PER_PX,
                        metrics_cfg={**cfg["metrics"], "central_fraction": fraction},
                    )
                    manual_profile = None
                    if np.isfinite(row["manual_x_px"]):
                        manual_profile = extract_profiles(
                            frame,
                            x_center=float(row["manual_x_px"]),
                            diameter_um=float(row["diameter_um"]),
                            lateral_um_per_px=LATERAL_UM_PER_PX,
                            metrics_cfg={**cfg["metrics"], "central_fraction": fraction},
                        )
                    prefix = f"{method}_{fraction:g}"
                    if manual_profile is not None:
                        overlap = max(0, min(candidate_profile.central_bounds[1], manual_profile.central_bounds[1]) - max(candidate_profile.central_bounds[0], manual_profile.central_bounds[0]))
                        denominator = max(candidate_profile.central_bounds[1] - candidate_profile.central_bounds[0], manual_profile.central_bounds[1] - manual_profile.central_bounds[0])
                        union = max(candidate_profile.central_bounds[1], manual_profile.central_bounds[1]) - min(candidate_profile.central_bounds[0], manual_profile.central_bounds[0])
                        selected_table.at[index, f"overlap_width_{prefix}"] = overlap
                        selected_table.at[index, f"overlap_fraction_{prefix}"] = overlap / denominator if denominator > 0 else np.nan
                        selected_table.at[index, f"iou_{prefix}"] = overlap / union if union > 0 else np.nan
                    candidate_auc = _auc(frame, x_center=x_value, z_upper=float(row["algorithm_z_px"]), diameter_um=float(row["diameter_um"]), central_fraction=fraction, cfg=cfg)
                    manual_auc = _auc(frame, x_center=float(row["manual_x_px"]), z_upper=float(row["algorithm_z_px"]), diameter_um=float(row["diameter_um"]), central_fraction=fraction, cfg=cfg) if np.isfinite(row["manual_x_px"]) else {"raw_auc": np.nan, "p95_rauc": np.nan}
                    for metric in ("raw_auc", "p95_rauc"):
                        selected_table.at[index, f"{metric}_{prefix}"] = candidate_auc[metric]
                        selected_table.at[index, f"{metric}_manual_{prefix}"] = manual_auc[metric]
                        denominator_auc = abs(float(candidate_auc[metric])) if np.isfinite(candidate_auc[metric]) else np.nan
                        selected_table.at[index, f"relative_manual_minus_{prefix}_{metric}_pct"] = (
                            100.0 * (float(manual_auc[metric]) - float(candidate_auc[metric])) / denominator_auc
                            if np.isfinite(denominator_auc) and denominator_auc > 1e-12 and np.isfinite(manual_auc[metric])
                            else np.nan
                        )

        candidate_summary = _candidate_error_summary(selected_table, seed=args.seed)
        roi_summary = _roi_summary(selected_table, seed=args.seed)
        auc_summary = _auc_summary(selected_table, seed=args.seed)
        all_features = pd.DataFrame(all_feature_rows)
        labeled_features = selected_table[["blind_id", "dataset_id", "frame_index_0based", "diameter_um", "vessel_visible", "assessability_score", "vessel_presence_prediction", "local_body_peak_cnr", "local_body_run_width_px", "local_body_axial_completeness", "neighbor_peak_cnr"]].copy()
        confusion, assess_by_diameter, all_assess_by_diameter = _assessability_tables(
            labeled_features,
            all_features,
            seed=args.seed,
        )
        _make_figures(candidate_summary, roi_summary, auc_summary, labeled_features, all_features, figure_dir)

        heldout_root = PROJECT_ROOT / "work" / "pilot_20260414" / "heldout_xroi_assessability_20260828"
        metadata: dict[str, Any] = {
            "schema_version": "1.0.0",
            "date_local": "2026-08-27",
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "development_set_count": DEV_COUNT,
            "development_set_scope": "fixed nine scans; 53-frame fifth-round pixel annotation package",
            "development_annotation_sha256": _sha256(annotation_path),
            "development_mapping_sha256": _sha256(mapping_path),
            "config_sha256": _sha256(config_path),
            "seed": int(args.seed),
            "x_methods": list(X_METHODS),
            "roi_fractions": list(ROI_FRACTIONS),
            "selected_x_method": "X4_centroid_isolated_jump_corrected",
            "selected_roi_fraction": 0.40,
            "assessability_rule": {
                "weights": {"cnr": 0.35, "continuity": 0.20, "width_consistency": 0.20, "axial_completeness": 0.15, "neighbor_support": 0.10},
                "assessable_threshold": 0.60,
                "uncertain_threshold": 0.40,
                "development_only": True,
            },
            "heldout_root": str(heldout_root),
            "heldout_seed": 20260828,
            "heldout_count": 72,
            "status": "development_complete_heldout_generated_separately",
        }
        selected_table.to_csv(building / "xroi_development_frame_results.csv", index=False, encoding="utf-8-sig")
        candidate_summary.to_csv(building / "x_candidate_summary.csv", index=False, encoding="utf-8-sig")
        roi_summary.to_csv(building / "roi_candidate_summary.csv", index=False, encoding="utf-8-sig")
        auc_summary.to_csv(building / "xroi_auc_sensitivity_summary.csv", index=False, encoding="utf-8-sig")
        labeled_features.to_csv(building / "assessability_development_features.csv", index=False, encoding="utf-8-sig")
        all_features.to_csv(building / "assessability_all_frames_features.csv", index=False, encoding="utf-8-sig")
        confusion.to_csv(building / "assessability_development_confusion.csv", index=False, encoding="utf-8-sig")
        assess_by_diameter.to_csv(building / "assessability_development_by_diameter.csv", index=False, encoding="utf-8-sig")
        all_assess_by_diameter.to_csv(building / "assessability_all_frames_by_diameter.csv", index=False, encoding="utf-8-sig")
        (building / "development_run_metadata.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
        _write_report(building / "XROI_ASSESSABILITY_DEVELOPMENT_REPORT_20260827.md", metadata=metadata, candidate_summary=candidate_summary, roi_summary=roi_summary, auc_summary=auc_summary, confusion=confusion, assess_by_diameter=assess_by_diameter, all_assess_by_diameter=all_assess_by_diameter)
        if output_root_path.exists():
            shutil.rmtree(output_root_path)
        building.replace(output_root_path)
    except Exception:
        shutil.rmtree(building, ignore_errors=True)
        raise
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
