"""Per-frame OCTA tail metrics from floating-point pre-filter Flow.

The implementation intentionally keeps geometry, background correction and
normalisation explicit.  It never sums laterally: the central vessel third is
collapsed with a median depth profile before axial integration.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable

import numpy as np
import pandas as pd


MAD_SCALE = 1.4826


# Nominal physical intervals used for the robustness diagnostic.  Pixel
# endpoints deliberately use the same round(um / axial_scale) convention as
# the frozen rNT-AUC implementation, so the decomposition adds no new geometry
# rule and the cumulative values remain directly comparable with legacy rows.
DECOMPOSITION_INTERVALS_UM: tuple[tuple[int, int], ...] = (
    (0, 100),
    (100, 200),
    (200, 300),
    (0, 200),
    (0, 300),
)


@dataclass(frozen=True)
class FrameProfiles:
    profile: np.ndarray
    background: np.ndarray
    background_mad: np.ndarray
    excess: np.ndarray
    central_bounds: tuple[int, int]
    left_bounds: tuple[int, int]
    right_bounds: tuple[int, int]
    roi_complete: bool


def _centered_bounds(center: int, width: int) -> tuple[int, int]:
    """Return a deterministic half-open interval of exactly ``width`` pixels."""

    start = int(center) - int(width) // 2
    return start, start + int(width)


def _profile_geometry(
    x_center: float,
    diameter_um: float,
    lateral_um_per_px: float,
    n_x: int,
    metrics_cfg: dict[str, Any],
) -> tuple[tuple[int, int], tuple[int, int], tuple[int, int], bool]:
    vessel_width = max(1, int(round(diameter_um / lateral_um_per_px)))
    central_width = max(
        int(metrics_cfg.get("minimum_profile_width_px", 3)),
        int(round(vessel_width * float(metrics_cfg["central_fraction"]))),
    )
    background_width = max(
        1,
        int(round(vessel_width * float(
            metrics_cfg.get("side_background_width_fraction_of_diameter", 1 / 3)
        ))),
    )
    gap = max(
        1,
        int(round(vessel_width * float(
            metrics_cfg.get("side_background_gap_fraction_of_diameter", 1 / 4)
        ))),
    )
    center = int(round(float(x_center)))
    vessel_left, vessel_right = _centered_bounds(center, vessel_width)
    central = _centered_bounds(center, central_width)
    left = (vessel_left - gap - background_width, vessel_left - gap)
    right = (vessel_right + gap, vessel_right + gap + background_width)
    all_bounds = (central, left, right)
    complete = all(0 <= lo < hi <= n_x for lo, hi in all_bounds)
    return central, left, right, complete


def extract_profiles(
    frame_zx: np.ndarray,
    *,
    x_center: float,
    diameter_um: float,
    lateral_um_per_px: float,
    metrics_cfg: dict[str, Any],
) -> FrameProfiles:
    if frame_zx.ndim != 2:
        raise ValueError(f"Expected frame[z,x], got {frame_zx.shape}")
    n_z, n_x = frame_zx.shape
    central, left, right, complete = _profile_geometry(
        x_center, diameter_um, lateral_um_per_px, n_x, metrics_cfg
    )
    if not complete:
        nan = np.full(n_z, np.nan, dtype=float)
        return FrameProfiles(nan, nan.copy(), nan.copy(), nan.copy(), central, left, right, False)

    c0, c1 = central
    l0, l1 = left
    r0, r1 = right
    central_values = np.asarray(frame_zx[:, c0:c1], dtype=float)
    background_values = np.concatenate(
        (np.asarray(frame_zx[:, l0:l1], dtype=float), np.asarray(frame_zx[:, r0:r1], dtype=float)),
        axis=1,
    )
    minimum_bg = int(metrics_cfg.get("minimum_background_pixels_per_depth", 6))
    if background_values.shape[1] < minimum_bg:
        nan = np.full(n_z, np.nan, dtype=float)
        return FrameProfiles(nan, nan.copy(), nan.copy(), nan.copy(), central, left, right, False)

    profile = np.median(central_values, axis=1)
    background = np.median(background_values, axis=1)
    background_mad = np.median(np.abs(background_values - background[:, None]), axis=1)
    excess = np.maximum(profile - background, 0.0)
    return FrameProfiles(
        profile,
        background,
        background_mad,
        excess,
        central,
        left,
        right,
        True,
    )


def _normalizer(
    excess: np.ndarray,
    z_upper: int,
    diameter_px: int,
    core_fraction: Iterable[float],
    mode: str,
    epsilon: float,
) -> tuple[float, bool, int, int]:
    fractions = list(core_fraction)
    if len(fractions) != 2 or not 0 <= fractions[0] < fractions[1] <= 1:
        raise ValueError(f"Invalid core_fraction={fractions}")
    core_start = z_upper + int(np.ceil(fractions[0] * diameter_px))
    core_stop = z_upper + int(np.floor(fractions[1] * diameter_px)) + 1
    if core_start < 0 or core_stop > len(excess) or core_stop <= core_start:
        return np.nan, False, core_start, core_stop
    values = excess[core_start:core_stop]
    if mode == "p95":
        value = float(np.quantile(values, 0.95))
    elif mode == "median":
        value = float(np.median(values))
    elif mode == "none":
        return 1.0, True, core_start, core_stop
    else:
        raise ValueError(f"Unknown normalisation mode: {mode}")
    return value, bool(np.isfinite(value) and value > epsilon), core_start, core_stop


def _window_pixels(window_um: float, axial_um_per_px: float) -> int:
    return max(1, int(round(float(window_um) / float(axial_um_per_px))))


def _tail_length_pixels(
    above: np.ndarray,
    *,
    minimum_true: int,
    allowed_gap: int,
    stop_below: int,
) -> tuple[float, int | None]:
    """Conservative extent of the threshold-connected tail from its start.

    The component must start at the operational tail boundary (up to the
    configured one-pixel gap), contain at least ``minimum_true`` positive
    pixels, and cannot cross a gap longer than ``allowed_gap``.  A run of
    ``stop_below`` negatives independently confirms termination.
    """

    above = np.asarray(above, dtype=bool)
    if above.size == 0:
        return np.nan, None
    last_true: int | None = None
    consecutive_true = 0
    gap = 0
    confirmed = False
    disconnected = False
    for index, value in enumerate(above):
        if value:
            # Once a gap is longer than the permitted interruption, a later
            # supra-threshold island is not part of the tail connected to
            # z_s.  Do not let that island move the endpoint deeper.
            if disconnected:
                break
            consecutive_true += 1
            last_true = index
            gap = 0
            confirmed = confirmed or consecutive_true >= minimum_true
        else:
            consecutive_true = 0
            gap += 1
            disconnected = disconnected or gap > allowed_gap
            if gap >= stop_below:
                break
    if not confirmed or last_true is None:
        return 0.0, None
    return float(last_true + 1), last_true


def _tail_cnr(
    frame_zx: np.ndarray,
    profiles: FrameProfiles,
    z_start: int,
    n_depth: int,
    epsilon: float,
) -> tuple[float, bool]:
    z_stop = z_start + n_depth
    if z_start < 0 or z_stop > frame_zx.shape[0]:
        return np.nan, False
    c0, c1 = profiles.central_bounds
    l0, l1 = profiles.left_bounds
    r0, r1 = profiles.right_bounds
    tail = np.asarray(frame_zx[z_start:z_stop, c0:c1], dtype=float).ravel()
    background = np.concatenate(
        (
            np.asarray(frame_zx[z_start:z_stop, l0:l1], dtype=float).ravel(),
            np.asarray(frame_zx[z_start:z_stop, r0:r1], dtype=float).ravel(),
        )
    )
    median_background = float(np.median(background))
    mad_background = float(np.median(np.abs(background - median_background)))
    denominator = MAD_SCALE * mad_background
    if not np.isfinite(denominator) or denominator <= epsilon:
        return np.nan, False
    cnr = (float(np.median(tail)) - median_background) / denominator
    return float(cnr), bool(np.isfinite(cnr))


def one_factor_parameter_sets(metrics_cfg: dict[str, Any], alphas: Iterable[float]) -> list[dict[str, Any]]:
    """Pre-specified one-factor-at-a-time matrix; no outcome-driven tuning."""

    primary_alpha = float(metrics_cfg.get("primary_alpha", 0.15))
    primary_guard = int(metrics_cfg["primary_guard_px"])
    primary_window = float(metrics_cfg["primary_window_um"])
    primary_norm = str(metrics_cfg["primary_normalization"])
    sets: list[dict[str, Any]] = [
        {
            "parameter_set_id": "primary",
            "alpha": primary_alpha,
            "guard_px": primary_guard,
            "window_um": primary_window,
            "normalization": primary_norm,
            "changed_factor": "none",
        }
    ]
    for alpha in alphas:
        alpha = float(alpha)
        if not np.isclose(alpha, primary_alpha):
            sets.append({**sets[0], "parameter_set_id": f"alpha_{alpha:.2f}".replace(".", "p"), "alpha": alpha, "changed_factor": "alpha"})
    for guard in metrics_cfg["guards_px"]:
        guard = int(guard)
        if guard != primary_guard:
            sets.append({**sets[0], "parameter_set_id": f"guard_{guard}px", "guard_px": guard, "changed_factor": "guard"})
    for window in metrics_cfg["windows_um"]:
        window = float(window)
        if not np.isclose(window, primary_window):
            sets.append({**sets[0], "parameter_set_id": f"window_{window:g}um", "window_um": window, "changed_factor": "window"})
    for mode in metrics_cfg["normalizations"]:
        mode = str(mode)
        if mode != primary_norm:
            sets.append({**sets[0], "parameter_set_id": f"normalization_{mode}", "normalization": mode, "changed_factor": "normalization"})
    identifiers = [item["parameter_set_id"] for item in sets]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("Duplicate sensitivity parameter_set_id")
    return sets


def decomposition_formula_metadata(
    scales_cfg: dict[str, Any], metrics_cfg: dict[str, Any]
) -> dict[str, Any]:
    """Machine-readable definitions for the frame decomposition CSV.

    This is written beside every scan output.  In particular, it records that
    an incomplete interval is represented as missing rather than zero-filled,
    interpolated or extrapolated.
    """

    axial = float(scales_cfg["axial_um_per_px"])
    guard = int(metrics_cfg["primary_guard_px"])
    endpoints = {
        str(depth): _window_pixels(float(depth), axial) for depth in (100, 200, 300)
    }
    return {
        "schema_version": "1.0.0",
        "coordinate_convention": "zero_based_half_open_intervals",
        "z_tail_start_formula": "z_upper_px + round(diameter_um / axial_um_per_px) + primary_guard_px",
        "primary_guard_px": guard,
        "z_max_definition": (
            "z_last_usable_px is the deepest row with finite central/background profiles; "
            "z_max_usable_boundary_px is its exclusive lower pixel boundary"
        ),
        "margin_px_formula": "z_last_usable_px - z_tail_start_px",
        "margin_um_formula": "margin_px * axial_um_per_px",
        "available_tail_samples_formula": "z_max_usable_boundary_px - z_tail_start_px",
        "axial_um_per_px": axial,
        "nominal_endpoint_pixels": endpoints,
        "intervals_um": [f"{start}_{stop}" for start, stop in DECOMPOSITION_INTERVALS_UM],
        "complete_interval_rule": (
            "both half-open bounds lie inside the image and every central/background "
            "profile sample in the interval is finite"
        ),
        "incomplete_interval_policy": "NaN; never zero-fill, interpolate, extrapolate, or move z_tail_start_px",
        "central_profile": "lateral median within the central vessel-third ROI",
        "background_profile": "per-depth median of the bilateral local background ROIs",
        "tail_excess_formula": "maximum(central_profile - background_profile, 0)",
        "tail_auc_raw_formula": "axial_um_per_px * sum(tail_excess over the complete interval)",
        "central_profile_auc_formula": "axial_um_per_px * sum(central_profile over the complete interval)",
        "background_auc_formula": "axial_um_per_px * sum(background_profile over the complete interval)",
        "vessel_core_interval": list(metrics_cfg["core_fraction"]),
        "vessel_denominators": ["p95", "p90", "median"],
        "normalized_rauc_formula": "tail_auc_raw / vessel_core_denominator",
        "primary_rauc_alias": "rAUC_<interval> equals rAUC_p95_<interval>",
    }


def _finite_profile_boundary(profiles: FrameProfiles) -> tuple[int | None, int | None]:
    """Return deepest finite row and its exclusive boundary.

    Raw pilot volumes are finite throughout.  Keeping this explicit prevents a
    future malformed/partially masked row from silently being counted as
    available depth.
    """

    finite = (
        np.isfinite(profiles.profile)
        & np.isfinite(profiles.background)
        & np.isfinite(profiles.background_mad)
        & np.isfinite(profiles.excess)
    )
    indices = np.flatnonzero(finite)
    if not len(indices):
        return None, None
    last = int(indices[-1])
    return last, last + 1


def _interval_values(
    profiles: FrameProfiles,
    *,
    z_start: int,
    offset_start_px: int,
    offset_stop_px: int,
    axial_um_per_px: float,
) -> dict[str, float | bool]:
    """Calculate one exact interval or return NaNs when it is incomplete."""

    start = int(z_start + offset_start_px)
    stop = int(z_start + offset_stop_px)
    n_z = len(profiles.profile)
    complete = bool(
        0 <= start < stop <= n_z
        and np.isfinite(profiles.profile[start:stop]).all()
        and np.isfinite(profiles.background[start:stop]).all()
        and np.isfinite(profiles.excess[start:stop]).all()
    )
    if not complete:
        return {
            "complete": False,
            "tail_auc_raw": np.nan,
            "central_profile_auc": np.nan,
            "background_auc": np.nan,
            "central_profile_median": np.nan,
            "background_median": np.nan,
            "tail_excess_median": np.nan,
        }
    central = profiles.profile[start:stop]
    background = profiles.background[start:stop]
    excess = profiles.excess[start:stop]
    return {
        "complete": True,
        "tail_auc_raw": float(axial_um_per_px * np.sum(excess)),
        "central_profile_auc": float(axial_um_per_px * np.sum(central)),
        "background_auc": float(axial_um_per_px * np.sum(background)),
        "central_profile_median": float(np.median(central)),
        "background_median": float(np.median(background)),
        "tail_excess_median": float(np.median(excess)),
    }


def compute_scan_auc_decomposition(
    flow_fzx: np.ndarray,
    source_frames: np.ndarray,
    primary_tracking: pd.DataFrame,
    *,
    dataset_id: str,
    diameter_um: float,
    velocity_mm_s: float,
    technical_replicate: int,
    scales_cfg: dict[str, Any],
    metrics_cfg: dict[str, Any],
    tracking_cfg: dict[str, Any],
) -> pd.DataFrame:
    """Return one diagnostic row per frame at the frozen primary geometry.

    Unlike the long one-factor-at-a-time table, this table stores numerator,
    denominator, incremental shells, cumulative intervals and depth margin in
    one row.  No value is calculated from a partial depth interval.
    """

    if flow_fzx.ndim != 3:
        raise ValueError(f"Expected Flow[frame,z,x], got {flow_fzx.shape}")
    source_frames = np.asarray(source_frames, dtype=int)
    if flow_fzx.shape[0] != len(source_frames):
        raise ValueError("source_frames length differs from Flow frame dimension")
    axial = float(scales_cfg["axial_um_per_px"])
    lateral = float(scales_cfg["lateral_um_per_px"])
    diameter_px = int(round(float(diameter_um) / axial))
    guard_px = int(metrics_cfg["primary_guard_px"])
    epsilon = float(metrics_cfg.get("epsilon", 1e-12))
    primary_alpha = float(tracking_cfg["primary_alpha"])
    core_fraction = metrics_cfg["core_fraction"]
    tracking_map = {
        int(row.frame_index): row for _, row in primary_tracking.iterrows()
    }
    endpoint_px = {
        depth: _window_pixels(float(depth), axial) for depth in (0, 100, 200, 300)
    }
    endpoint_px[0] = 0
    rows: list[dict[str, Any]] = []

    confidence_columns = (
        "old_tracking_class",
        "x_path_confidence",
        "x_path_confidence_class",
        "z_edge_confidence",
        "z_edge_confidence_class",
        "overall_tracking_confidence",
        "overall_tracking_confidence_class",
        "qc_valid",
    )

    for local_index, source_frame in enumerate(source_frames):
        frame = flow_fzx[local_index]
        track = tracking_map.get(int(source_frame))
        row: dict[str, Any] = {
            "dataset_id": dataset_id,
            "diameter_um": float(diameter_um),
            "velocity_mm_s": float(velocity_mm_s),
            "technical_replicate": int(technical_replicate),
            "frame_index": int(source_frame),
            "alpha": primary_alpha,
            "primary_guard_px": guard_px,
            "tracking_class": None if track is None else str(track.get("tracking_class", "")),
            "x_center_px": np.nan if track is None else track.get("x_center_px", np.nan),
            "z_upper_px": np.nan if track is None else track.get("z_upper_px", np.nan),
            "diameter_axial_px": diameter_px,
            "z_lower_op_px": np.nan,
            "z_tail_start_px": np.nan,
            "z_last_usable_px": np.nan,
            "deepest_valid_index_px": np.nan,
            "z_max_usable_boundary_px": np.nan,
            "available_tail_samples": np.nan,
            "margin_px": np.nan,
            "margin_um": np.nan,
            "core_start_px": np.nan,
            "core_stop_exclusive_px": np.nan,
            "vessel_core_p95": np.nan,
            "vessel_core_p90": np.nan,
            "vessel_core_median": np.nan,
            "normalizer_p95_valid": False,
            "normalizer_p90_valid": False,
            "normalizer_median_valid": False,
            "decomposition_valid": False,
            "invalid_reason": "tracking_missing_or_failed",
        }
        for name in confidence_columns:
            row[name] = np.nan if track is None or name not in track.index else track[name]
        for depth in (100, 200, 300):
            row[f"complete_{depth}um"] = False
        for interval_start, interval_stop in DECOMPOSITION_INTERVALS_UM:
            label = f"{interval_start}_{interval_stop}"
            row[f"shell_complete_{label}"] = False
            row[f"tail_auc_raw_{label}"] = np.nan
            row[f"central_profile_auc_{label}"] = np.nan
            row[f"background_auc_{label}"] = np.nan
            row[f"central_profile_median_{label}"] = np.nan
            row[f"background_median_{label}"] = np.nan
            row[f"tail_excess_median_{label}"] = np.nan
            for mode in ("p95", "p90", "median"):
                row[f"rAUC_{mode}_{label}"] = np.nan
                row[f"rAUC_{mode}_{label}_valid"] = False
            # Backward-friendly primary-denominator alias requested for the
            # diagnostic package, distinct from legacy rnt_auc_um.
            row[f"rAUC_{label}"] = np.nan
            row[f"rAUC_{label}_valid"] = False

        tracking_valid = bool(
            track is not None
            and str(track.get("tracking_class", "failed")) != "failed"
            and np.isfinite(track.get("x_center_px", np.nan))
            and np.isfinite(track.get("z_upper_px", np.nan))
            and ("qc_valid" not in track.index or bool(track["qc_valid"]))
        )
        if not tracking_valid:
            rows.append(row)
            continue

        z_upper = int(round(float(track["z_upper_px"])))
        profiles = extract_profiles(
            frame,
            x_center=float(track["x_center_px"]),
            diameter_um=float(diameter_um),
            lateral_um_per_px=lateral,
            metrics_cfg=metrics_cfg,
        )
        if not profiles.roi_complete:
            row["invalid_reason"] = "lateral_roi_out_of_bounds"
            rows.append(row)
            continue

        z_lower = z_upper + diameter_px
        z_start = z_lower + guard_px
        last_usable, boundary = _finite_profile_boundary(profiles)
        row.update(
            {
                "z_upper_px": z_upper,
                "z_lower_op_px": z_lower,
                "z_tail_start_px": z_start,
                "z_last_usable_px": np.nan if last_usable is None else last_usable,
                "deepest_valid_index_px": np.nan if last_usable is None else last_usable,
                "z_max_usable_boundary_px": np.nan if boundary is None else boundary,
                "available_tail_samples": np.nan if boundary is None else int(boundary - z_start),
                "margin_px": np.nan if last_usable is None else int(last_usable - z_start),
                "margin_um": np.nan if last_usable is None else float((last_usable - z_start) * axial),
            }
        )
        if not (0 <= z_start < frame.shape[0]) or boundary is None:
            row["invalid_reason"] = "tail_start_out_of_bounds_or_no_finite_depth"
            rows.append(row)
            continue

        # Extract all denominator summaries from the identical background-
        # corrected upper-core samples used by the legacy normalizer.
        p95, p95_valid, core_start, core_stop = _normalizer(
            profiles.excess, z_upper, diameter_px, core_fraction, "p95", epsilon
        )
        median, median_valid, _, _ = _normalizer(
            profiles.excess, z_upper, diameter_px, core_fraction, "median", epsilon
        )
        if 0 <= core_start < core_stop <= len(profiles.excess):
            p90 = float(np.quantile(profiles.excess[core_start:core_stop], 0.90))
            p90_valid = bool(np.isfinite(p90) and p90 > epsilon)
        else:
            p90 = np.nan
            p90_valid = False
        row.update(
            {
                "core_start_px": core_start,
                "core_stop_exclusive_px": core_stop,
                "vessel_core_p95": p95,
                "vessel_core_p90": p90,
                "vessel_core_median": median,
                "normalizer_p95_valid": p95_valid,
                "normalizer_p90_valid": p90_valid,
                "normalizer_median_valid": median_valid,
            }
        )

        denominators = {"p95": p95, "p90": p90, "median": median}
        denominator_valid = {
            "p95": p95_valid, "p90": p90_valid, "median": median_valid
        }
        for depth in (100, 200, 300):
            stop = z_start + endpoint_px[depth]
            row[f"complete_{depth}um"] = bool(
                stop <= frame.shape[0]
                and np.isfinite(profiles.profile[z_start:stop]).all()
                and np.isfinite(profiles.background[z_start:stop]).all()
                and np.isfinite(profiles.excess[z_start:stop]).all()
            )

        for interval_start, interval_stop in DECOMPOSITION_INTERVALS_UM:
            label = f"{interval_start}_{interval_stop}"
            values = _interval_values(
                profiles,
                z_start=z_start,
                offset_start_px=endpoint_px[interval_start],
                offset_stop_px=endpoint_px[interval_stop],
                axial_um_per_px=axial,
            )
            row[f"shell_complete_{label}"] = bool(values["complete"])
            for field in (
                "tail_auc_raw",
                "central_profile_auc",
                "background_auc",
                "central_profile_median",
                "background_median",
                "tail_excess_median",
            ):
                row[f"{field}_{label}"] = values[field]
            for mode, denominator in denominators.items():
                valid = bool(values["complete"] and denominator_valid[mode])
                value = (
                    float(values["tail_auc_raw"]) / (float(denominator) + epsilon)
                    if valid
                    else np.nan
                )
                row[f"rAUC_{mode}_{label}"] = value
                row[f"rAUC_{mode}_{label}_valid"] = valid
            row[f"rAUC_{label}"] = row[f"rAUC_p95_{label}"]
            row[f"rAUC_{label}_valid"] = row[f"rAUC_p95_{label}_valid"]

        row["decomposition_valid"] = True
        row["invalid_reason"] = "ok"
        rows.append(row)

    return pd.DataFrame(rows)


def compute_scan_metrics(
    flow_fzx: np.ndarray,
    source_frames: np.ndarray,
    tracking_by_alpha: dict[float, pd.DataFrame],
    *,
    dataset_id: str,
    diameter_um: float,
    velocity_mm_s: float,
    technical_replicate: int,
    scales_cfg: dict[str, Any],
    metrics_cfg: dict[str, Any],
    tracking_cfg: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return long rNT-AUC and auxiliary frame tables for one scan."""

    if flow_fzx.ndim != 3:
        raise ValueError(f"Expected Flow[frame,z,x], got {flow_fzx.shape}")
    if flow_fzx.shape[0] != len(source_frames):
        raise ValueError("source_frames length differs from Flow frame dimension")
    axial = float(scales_cfg["axial_um_per_px"])
    lateral = float(scales_cfg["lateral_um_per_px"])
    diameter_px = int(round(float(diameter_um) / axial))
    epsilon = float(metrics_cfg.get("epsilon", 1e-12))
    primary_alpha = float(tracking_cfg["primary_alpha"])
    metrics_cfg = {**metrics_cfg, "primary_alpha": primary_alpha}
    parameter_sets = one_factor_parameter_sets(metrics_cfg, sorted(tracking_by_alpha))

    tracking_maps: dict[float, dict[int, pd.Series]] = {}
    for alpha, table in tracking_by_alpha.items():
        tracking_maps[float(alpha)] = {
            int(row.frame_index): row for _, row in table.iterrows()
        }

    auc_rows: list[dict[str, Any]] = []
    auxiliary_rows: list[dict[str, Any]] = []
    tail_ks = [float(value) for value in metrics_cfg["tail_length_k"]]
    primary_guard = int(metrics_cfg["primary_guard_px"])
    aux_guards = [int(value) for value in metrics_cfg["guards_px"]]
    cnr_pixels = _window_pixels(metrics_cfg["tail_cnr_window_um"], axial)

    for local_index, source_frame in enumerate(np.asarray(source_frames, dtype=int)):
        frame = flow_fzx[local_index]
        frame_profiles: dict[float, tuple[pd.Series | None, FrameProfiles | None, int | None]] = {}
        for alpha, mapping in tracking_maps.items():
            track = mapping.get(int(source_frame))
            if track is None or str(track["tracking_class"]) == "failed":
                frame_profiles[alpha] = (track, None, None)
                continue
            if not np.isfinite(track["x_center_px"]) or not np.isfinite(track["z_upper_px"]):
                frame_profiles[alpha] = (track, None, None)
                continue
            if "qc_valid" in track and not bool(track["qc_valid"]):
                frame_profiles[alpha] = (track, None, None)
                continue
            z_upper = int(round(float(track["z_upper_px"])))
            profiles = extract_profiles(
                frame,
                x_center=float(track["x_center_px"]),
                diameter_um=float(diameter_um),
                lateral_um_per_px=lateral,
                metrics_cfg=metrics_cfg,
            )
            frame_profiles[alpha] = (track, profiles, z_upper)

        for params in parameter_sets:
            alpha = float(params["alpha"])
            track, profiles, z_upper = frame_profiles.get(alpha, (None, None, None))
            base = {
                "dataset_id": dataset_id,
                "diameter_um": float(diameter_um),
                "velocity_mm_s": float(velocity_mm_s),
                "technical_replicate": int(technical_replicate),
                "frame_index": int(source_frame),
                **params,
            }
            valid = track is not None and profiles is not None and profiles.roi_complete and z_upper is not None
            reason = "ok" if valid else (
                "tracking_missing_or_failed" if track is None or profiles is None else "lateral_roi_out_of_bounds"
            )
            result = {
                **base,
                "metric_units": (
                    "um_relative_to_vessel_core"
                    if str(params["normalization"]) != "none"
                    else "flow_intensity_times_um"
                ),
                "tracking_class": None if track is None else str(track["tracking_class"]),
                "x_center_px": np.nan if track is None else track["x_center_px"],
                "z_upper_px": np.nan if z_upper is None else z_upper,
                "diameter_axial_px": diameter_px,
                "z_lower_op_px": np.nan,
                "z_tail_start_px": np.nan,
                "normalizer": np.nan,
                "normalizer_valid": False,
                "window_px": _window_pixels(params["window_um"], axial),
                "window_complete": False,
                "rnt_auc_um": np.nan,
                "metric_valid": False,
                "invalid_reason": reason,
            }
            if valid:
                z_lower = z_upper + diameter_px
                z_start = z_lower + int(params["guard_px"])
                n_window = int(result["window_px"])
                window_complete = 0 <= z_start and z_start + n_window <= frame.shape[0]
                normalizer, normalizer_valid, core_start, core_stop = _normalizer(
                    profiles.excess,
                    z_upper,
                    diameter_px,
                    metrics_cfg["core_fraction"],
                    str(params["normalization"]),
                    epsilon,
                )
                metric_valid = bool(window_complete and normalizer_valid)
                auc = (
                    axial * float(np.sum(profiles.excess[z_start:z_start + n_window])) / (normalizer + epsilon)
                    if metric_valid else np.nan
                )
                result.update(
                    {
                        "z_lower_op_px": z_lower,
                        "z_tail_start_px": z_start,
                        "core_start_px": core_start,
                        "core_stop_exclusive_px": core_stop,
                        "normalizer": normalizer,
                        "normalizer_valid": normalizer_valid,
                        "window_complete": window_complete,
                        "rnt_auc_um": auc,
                        "metric_valid": metric_valid,
                        "invalid_reason": "ok" if metric_valid else (
                            "depth_window_out_of_bounds" if not window_complete else "normalizer_nonpositive"
                        ),
                    }
                )
            auc_rows.append(result)

        # Auxiliary metrics use the primary alpha; guard and k sensitivities
        # are explicit but are not crossed with the AUC OFAT matrix.
        track, profiles, z_upper = frame_profiles.get(primary_alpha, (None, None, None))
        for guard in aux_guards:
            for tail_k in tail_ks:
                base = {
                    "dataset_id": dataset_id,
                    "diameter_um": float(diameter_um),
                    "velocity_mm_s": float(velocity_mm_s),
                    "technical_replicate": int(technical_replicate),
                    "frame_index": int(source_frame),
                    "alpha": primary_alpha,
                    "guard_px": guard,
                    "tail_length_k": tail_k,
                    "tracking_class": None if track is None else str(track["tracking_class"]),
                    "z_tail_start_px": np.nan,
                    "tail_length_px": np.nan,
                    "tail_length_um": np.nan,
                    "tail_endpoint_px": np.nan,
                    "tail_cnr": np.nan,
                    "tail_cnr_valid": False,
                    "aux_valid": False,
                    "invalid_reason": "tracking_missing_or_failed",
                }
                valid = track is not None and profiles is not None and profiles.roi_complete and z_upper is not None
                if valid:
                    z_start = z_upper + diameter_px + guard
                    if 0 <= z_start < frame.shape[0]:
                        threshold = profiles.background + tail_k * MAD_SCALE * profiles.background_mad
                        length_px, endpoint_rel = _tail_length_pixels(
                            profiles.profile[z_start:] > threshold[z_start:],
                            minimum_true=int(metrics_cfg["tail_length_min_run_px"]),
                            allowed_gap=int(metrics_cfg["tail_length_allowed_gap_px"]),
                            stop_below=int(metrics_cfg["tail_length_stop_below_px"]),
                        )
                        cnr, cnr_valid = _tail_cnr(frame, profiles, z_start, cnr_pixels, epsilon)
                        base.update(
                            {
                                "z_tail_start_px": z_start,
                                "tail_length_px": length_px,
                                "tail_length_um": length_px * axial if np.isfinite(length_px) else np.nan,
                                "tail_endpoint_px": (
                                    z_start + endpoint_rel if endpoint_rel is not None else np.nan
                                ),
                                "tail_cnr": cnr,
                                "tail_cnr_valid": cnr_valid,
                                "aux_valid": True,
                                "invalid_reason": "ok",
                            }
                        )
                    else:
                        base["invalid_reason"] = "tail_start_out_of_bounds"
                auxiliary_rows.append(base)

    return pd.DataFrame(auc_rows), pd.DataFrame(auxiliary_rows)
