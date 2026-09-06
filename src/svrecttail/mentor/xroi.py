"""Frozen X1/X2/X4 and assessability subset from the mentor delivery.

Source: python/10_xroi_refinement/develop_xroi_assessability.py in
Bscan_tail_auto_quantification_for_junior_20260903(1).

Only localization and frame-assessability functions are included. The
mentor package's tail AUC metrics are intentionally outside this module.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd


AXIAL_UM_PER_PX = 6.7
LATERAL_UM_PER_PX = 12.7


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


def isolated_jump_correct(
    sequence: np.ndarray,
    *,
    jump_px: float = 2.0,
    neighbor_tol_px: float = 2.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Correct only isolated 3-frame jumps; preserve smooth spatial changes."""

    corrected = np.asarray(sequence, dtype=float).copy()
    changed = np.zeros(len(corrected), dtype=bool)
    for index in range(1, len(corrected) - 1):
        left, current, right = corrected[index - 1 : index + 2]
        if not np.isfinite([left, current, right]).all():
            continue
        neighbor_median = float(np.median([left, current, right]))
        if (
            abs(current - neighbor_median) > jump_px
            and abs(left - right) <= neighbor_tol_px
        ):
            corrected[index] = neighbor_median
            changed[index] = True
    return corrected, changed


def local_body_features(
    frame_zx: np.ndarray, tracking: pd.Series
) -> dict[str, Any]:
    """Extract the mentor X1, X2 and X3 local body candidates."""

    n_z, n_x = frame_zx.shape
    x_anchor = _numeric(tracking.get("x_center_px"))
    z_upper = _numeric(tracking.get("z_upper_px"))
    diameter_um = _numeric(tracking.get("diameter_um"))
    half_width = _numeric(tracking.get("x_local_window_half_width_px"))
    if not np.isfinite([x_anchor, z_upper, diameter_um]).all():
        return {
            "valid_local_body": False,
            "invalid_reason": "missing_algorithm_anchor",
        }
    diameter_px = max(1, int(round(diameter_um / AXIAL_UM_PER_PX)))
    expected_lateral_width = max(
        3, int(round(diameter_um / LATERAL_UM_PER_PX))
    )
    half_width_px = (
        int(round(half_width))
        if np.isfinite(half_width)
        else max(8, int(round(1.5 * expected_lateral_width)))
    )
    anchor_int = int(round(x_anchor))
    x0 = max(0, anchor_int - half_width_px)
    x1 = min(n_x, anchor_int + half_width_px + 1)
    z0 = int(round(z_upper))
    z1 = min(n_z, z0 + diameter_px)
    if z0 < 0 or z1 <= z0 or x1 <= x0:
        return {
            "valid_local_body": False,
            "invalid_reason": "body_slab_out_of_bounds",
        }

    body = np.asarray(frame_zx[z0:z1, x0:x1], dtype=float)
    x_values = np.arange(x0, x1, dtype=float)
    column_score = np.quantile(body, 0.75, axis=0)
    exclusion_half_width = max(3, expected_lateral_width // 2 + 2)
    far_mask = np.abs(x_values - x_anchor) > exclusion_half_width
    if int(far_mask.sum()) < 5:
        far_mask = np.ones_like(x_values, dtype=bool)
    local_background = float(np.median(column_score[far_mask]))
    local_mad = float(
        np.median(np.abs(column_score[far_mask] - local_background))
    )
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
    cap = float(np.quantile(excess[positive], 0.90)) if positive.any() else 0.0
    weights = np.clip(excess, 0.0, cap)
    x2_center = (
        x_anchor
        if float(weights.sum()) <= 1e-12
        else float(np.sum(x_values * weights) / np.sum(weights))
    )

    peak_excess = float(np.max(excess)) if excess.size else 0.0
    half_peak_mask = (
        excess >= 0.5 * peak_excess
        if peak_excess > 0
        else np.zeros_like(excess, dtype=bool)
    )
    half_peak_runs = [
        (a + x0, b + x0) for a, b in _runs(half_peak_mask)
    ]
    if half_peak_runs:
        selected_half = max(half_peak_runs, key=run_score)
        x3_center = (selected_half[0] + selected_half[1] - 1) / 2.0
    else:
        x3_center = x_anchor

    center_int = int(round(x2_center))
    center_half = max(1, expected_lateral_width // 4)
    c0 = max(x0, center_int - center_half)
    c1 = min(x1, center_int + center_half + 1)
    row_background = np.median(body[:, far_mask], axis=1)
    row_mad = np.median(
        np.abs(body[:, far_mask] - row_background[:, None]), axis=1
    )
    row_threshold = row_background + 1.5 * 1.4826 * row_mad
    row_signal = np.median(frame_zx[z0:z1, c0:c1], axis=1)
    axial_completeness = (
        float(np.mean(row_signal > row_threshold)) if len(row_signal) else 0.0
    )

    run_width = max(1, selected[1] - selected[0])
    body_peak_cnr = peak_excess / local_sigma
    continuity_score = float(
        np.clip(
            run_width / max(0.5 * expected_lateral_width, 1.0),
            0.0,
            1.0,
        )
    )
    width_consistency = float(
        np.exp(-abs(np.log(run_width / max(expected_lateral_width, 1))))
    )

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


def add_assessability_score(features: pd.DataFrame) -> pd.DataFrame:
    """Apply the mentor's frozen five-component assessability score."""

    result = features.copy()
    cnr_score = np.clip(
        (result["local_body_peak_cnr"].to_numpy(float) - 2.0) / 6.0,
        0.0,
        1.0,
    )
    continuity_score = result["local_body_continuity_score"].to_numpy(float)
    width_score = result["local_body_width_consistency"].to_numpy(float)
    axial_score = np.clip(
        result["local_body_axial_completeness"].to_numpy(float), 0.0, 1.0
    )
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
