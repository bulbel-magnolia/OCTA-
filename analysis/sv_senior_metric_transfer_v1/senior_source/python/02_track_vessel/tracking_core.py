"""Auditable slow-axis-constrained tracking of a phantom vessel in OMAG Flow.

The public coordinate convention is zero-based ``[frame, z, x]``.  Local
measurements provide candidate centres/upper edges, but the reported track is
regularised using the complete slow-axis sequence.  A weak frame is therefore
never silently discarded: it is labelled ``model_assisted`` or ``failed``.

This module deliberately does not use structural OCT to define vessel borders.
The operational lower border is ``z_upper + round(diameter / axial_pitch)``.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence
import hashlib
import json
import math

import h5py
import numpy as np
import pandas as pd
from scipy.ndimage import (
    distance_transform_edt,
    gaussian_filter,
    gaussian_filter1d,
    median_filter,
    uniform_filter1d,
)


CLASS_HIGH = "high_confidence"
CLASS_ASSISTED = "model_assisted"
CLASS_FAILED = "failed"
EVIDENCE_STRONG = "strong_evidence"
EVIDENCE_WEAK = "weak_evidence"
EVIDENCE_UNAVAILABLE = "unavailable"


# Defaults live in one place.  The formal pilot JSON should repeat/freeze every
# value so that a run does not depend on a future code-default change.
DEFAULT_TRACKING_CONFIG: dict[str, Any] = {
    "upper_edge_alphas": [0.10, 0.15, 0.20],
    "primary_alpha": 0.15,
    "axial_um_per_px": 6.7,
    "lateral_um_per_px": 12.7,
    "primary_guard_px": 2,
    "primary_window_um": 300.0,
    "x_search_fraction": [0.10, 0.90],
    "z_search_px": [150, 330],
    "localization_depth_sigma_px": 1.2,
    "localization_slow_sigma_frames": 2.0,
    "x_core_depth_fraction_of_diameter": 0.40,
    "x_viterbi_max_jump_px": 8,
    "x_viterbi_jump_penalty": 0.10,
    "z_viterbi_max_jump_px": 8,
    "z_viterbi_jump_penalty": 0.12,
    "x_seed_min_robust_snr": 3.0,
    "x_seed_max_path_residual_px": 10.0,
    "confidence_revision_version": "local_x_and_z_components_v1",
    # Revised x-confidence is evaluated only near the already-established
    # Viterbi path.  The central search band is one known vessel radius; the
    # bilateral comparison bands and trajectory-change padding define the
    # local (not full-field) evidence window.
    "x_confidence_core_radius_fraction_of_diameter": 0.50,
    "x_confidence_flank_width_fraction_of_diameter": 0.75,
    "x_confidence_trajectory_change_multiplier": 2.0,
    "x_confidence_min_trajectory_pad_px": 2,
    "x_confidence_min_flank_points": 6,
    "x_confidence_min_robust_snr": 3.0,
    "x_confidence_alignment_tolerance_fraction_of_diameter": 0.35,
    "x_confidence_min_alignment_tolerance_px": 3.0,
    # A score reaches 0.5 when all component evidence is exactly at its
    # categorical threshold and approaches 1 by this multiple of threshold.
    "confidence_full_score_multiple": 2.0,
    "upper_edge_noise_multiplier": 3.0,
    "upper_edge_min_component_px": 3,
    "peak_min_snr": 3.0,
    "central_width_fraction_of_lateral_diameter": 1.0 / 3.0,
    "central_min_width_px": 3,
    "background_side_gap_px": 6,
    "background_side_width_px": 16,
    "background_min_total_width_px": 16,
    "background_depth_trend_sigma_px": 8.0,
    "trajectory_median_window_frames": 15,
    "trajectory_gaussian_sigma_frames": 3.0,
    "trajectory_outlier_min_px": 3.0,
    "trajectory_outlier_mad_multiplier": 3.0,
    "trajectory_min_seed_frames": 8,
    "trajectory_min_seed_fraction": 0.05,
    "max_model_assisted_gap_frames": 75,
    "qc_overlay_frames": 12,
    "input_dataset_candidates": ["p_bld_ed", "p_bld_ed_sparse"],
    "expected_z_px": 351,
    "expected_x_px": 500,
    "expected_frames": None,
}


class TrackingError(RuntimeError):
    """Raised when a scan cannot be interpreted or tracked defensibly."""


@dataclass(frozen=True)
class FlowVolume:
    data: np.ndarray
    dataset_name: str
    source_axis_order: tuple[str, str, str]
    axis_order_source: str
    source_shape: tuple[int, int, int]
    source_path: str


@dataclass(frozen=True)
class ProfileBundle:
    central: np.ndarray
    background: np.ndarray
    excess: np.ndarray
    excess_localization: np.ndarray
    background_sigma: np.ndarray
    central_left: np.ndarray
    central_right: np.ndarray
    side_left_start: np.ndarray
    side_left_stop: np.ndarray
    side_right_start: np.ndarray
    side_right_stop: np.ndarray
    roi_valid: np.ndarray


def robust_sigma(values: np.ndarray, axis: int | tuple[int, ...] | None = None) -> np.ndarray:
    """Scaled MAD, with NaN support and no hidden positive floor."""
    med = np.nanmedian(values, axis=axis, keepdims=True)
    return 1.4826 * np.nanmedian(np.abs(values - med), axis=axis)


def _as_plain_axis_order(value: Any) -> str | None:
    if value is None:
        return None
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="strict")
    if isinstance(value, str):
        return value
    arr = np.asarray(value)
    if arr.dtype.kind in "SU":
        flat = arr.ravel().tolist()
        return "".join(
            item.decode("utf-8") if isinstance(item, bytes) else str(item)
            for item in flat
        )
    return None


def parse_axis_order(value: str | Sequence[str]) -> tuple[str, str, str]:
    if isinstance(value, str):
        cleaned = value.lower().replace("[", "").replace("]", "")
        bits = tuple(x.strip() for x in cleaned.replace(";", ",").split(",") if x.strip())
    else:
        bits = tuple(str(x).lower().strip() for x in value)
    aliases = {"y": "frame", "frames": "frame", "depth": "z", "lateral": "x"}
    bits = tuple(aliases.get(x, x) for x in bits)
    if len(bits) != 3 or set(bits) != {"frame", "z", "x"}:
        raise TrackingError(
            f"axis_order must contain frame,z,x exactly once; received {value!r}"
        )
    return bits  # type: ignore[return-value]


def _sidecar_axis_order(path: Path) -> str | None:
    candidates = [
        path.with_suffix(path.suffix + ".json"),
        Path(str(path) + ".metadata.json"),
        path.with_suffix(".json"),
        path.parent / f"{path.stem}_metadata.json",
    ]
    for candidate in candidates:
        if not candidate.exists():
            continue
        try:
            obj = json.loads(candidate.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        # Exporter records how h5py actually sees the MATLAB v7.3 dataset.
        # Prefer that over any logical/MATLAB-order field in the same sidecar.
        for key in (
            "axis_order_h5py",
            "h5_axis_order",
            "dataset_axis_order",
            "axis_order",
        ):
            if key in obj:
                value = obj[key]
                return ",".join(value) if isinstance(value, list) else str(value)
        metadata = obj.get("metadata", {})
        for key in (
            "axis_order_h5py",
            "h5_axis_order",
            "dataset_axis_order",
            "axis_order",
        ):
            if key in metadata:
                value = metadata[key]
                return ",".join(value) if isinstance(value, list) else str(value)
    return None


def _unique_shape_axis_order(
    shape: tuple[int, int, int],
    expected_z: int | None,
    expected_x: int | None,
    expected_frames: int | None,
) -> tuple[str, str, str] | None:
    """Return an order only when dimensions make it mathematically unique."""
    import itertools

    candidates: list[tuple[str, str, str]] = []
    for order in itertools.permutations(("frame", "z", "x")):
        sizes = dict(zip(order, shape))
        if expected_z is not None and sizes["z"] != expected_z:
            continue
        if expected_x is not None and sizes["x"] != expected_x:
            continue
        if expected_frames is not None and sizes["frame"] != expected_frames:
            continue
        candidates.append(order)
    return candidates[0] if len(candidates) == 1 else None


def load_flow_volume(
    path: str | Path,
    *,
    dataset_candidates: Sequence[str] = ("p_bld_ed", "p_bld_ed_sparse"),
    axis_order: str | Sequence[str] | None = None,
    expected_z: int | None = 351,
    expected_x: int | None = 500,
    expected_frames: int | None = None,
) -> FlowVolume:
    """Load floating Flow and strictly establish its axis order.

    Attribute/sidecar/CLI metadata wins.  Shape inference is accepted only if
    exactly one permutation matches expected dimensions.  A full 500x500x351
    MATLAB v7.3 array is intentionally rejected without metadata because frame
    and x are ambiguous.
    """
    source = Path(path).resolve()
    if not source.exists():
        raise TrackingError(f"Flow input does not exist: {source}")
    if source.suffix.lower() == ".npy":
        raw = np.load(source, mmap_mode=None)
        dataset_name = "npy"
        attr_order = None
    else:
        try:
            h5 = h5py.File(source, "r")
        except OSError as exc:
            raise TrackingError(f"Cannot open Flow HDF5/MAT file {source}: {exc}") from exc
        with h5:
            selected = next((key for key in dataset_candidates if key in h5), None)
            if selected is None:
                available = sorted(k for k in h5.keys() if not k.startswith("#"))
                raise TrackingError(
                    f"No Flow dataset among {list(dataset_candidates)} in {source}; "
                    f"available={available}"
                )
            ds = h5[selected]
            if not isinstance(ds, h5py.Dataset) or ds.ndim != 3:
                raise TrackingError(f"Dataset {selected!r} must be a 3-D numeric array")
            attr_order = None
            for owner in (ds.attrs, h5.attrs):
                for key in ("axis_order", "h5_axis_order", "dataset_axis_order"):
                    if key in owner:
                        attr_order = _as_plain_axis_order(owner[key])
                        if attr_order:
                            break
                if attr_order:
                    break
            raw = np.asarray(ds)
            dataset_name = selected
    if raw.ndim != 3:
        raise TrackingError(f"Flow input must be 3-D; got shape {raw.shape}")
    if raw.dtype.kind not in "fiu":
        raise TrackingError(f"Flow input must be numeric; got dtype {raw.dtype}")

    order_source: str
    if axis_order is not None:
        order = parse_axis_order(axis_order)
        order_source = "cli_or_caller"
    elif attr_order:
        order = parse_axis_order(attr_order)
        order_source = "dataset_attribute"
    else:
        sidecar_order = _sidecar_axis_order(source)
        if sidecar_order:
            order = parse_axis_order(sidecar_order)
            order_source = "metadata_sidecar"
        else:
            unique = _unique_shape_axis_order(
                tuple(int(x) for x in raw.shape), expected_z, expected_x, expected_frames
            )
            if unique is None:
                raise TrackingError(
                    "Axis order is ambiguous. Add dataset attribute axis_order, a metadata "
                    "sidecar, or --axis-order. No axis permutation was guessed. "
                    f"shape={raw.shape}, expected(frame,z,x)="
                    f"({expected_frames},{expected_z},{expected_x})"
                )
            order = unique
            order_source = "unique_shape_match"

    permutation = tuple(order.index(name) for name in ("frame", "z", "x"))
    data = np.transpose(raw, permutation).astype(np.float32, copy=False)
    if expected_z is not None and data.shape[1] != expected_z:
        raise TrackingError(f"Expected z={expected_z}, got normalized shape {data.shape}")
    if expected_x is not None and data.shape[2] != expected_x:
        raise TrackingError(f"Expected x={expected_x}, got normalized shape {data.shape}")
    if expected_frames is not None and data.shape[0] != expected_frames:
        raise TrackingError(
            f"Expected {expected_frames} frames, got normalized shape {data.shape}"
        )
    finite_fraction = float(np.isfinite(data).mean())
    if finite_fraction < 0.999:
        raise TrackingError(
            f"Flow has too many non-finite pixels: finite_fraction={finite_fraction:.6f}"
        )
    return FlowVolume(
        data=data,
        dataset_name=dataset_name,
        source_axis_order=order,
        axis_order_source=order_source,
        source_shape=tuple(int(x) for x in raw.shape),
        source_path=str(source),
    )


def merge_tracking_config(config: Mapping[str, Any] | None) -> dict[str, Any]:
    merged = dict(DEFAULT_TRACKING_CONFIG)
    if config:
        tracking = config.get("tracking", config)
        if not isinstance(tracking, Mapping):
            raise TrackingError("configuration field 'tracking' must be an object")
        raw_tracking = dict(tracking)
        # Accept both the concise study-config vocabulary and the explicit
        # internal names written into each run summary.
        aliases = {
            "alphas": "upper_edge_alphas",
            "minimum_axial_component_px": "upper_edge_min_component_px",
            "background_mad_multiplier": "upper_edge_noise_multiplier",
        }
        for source, target in aliases.items():
            if source in raw_tracking:
                raw_tracking[target] = raw_tracking[source]
        merged.update({k: v for k, v in raw_tracking.items() if k in DEFAULT_TRACKING_CONFIG})
        calibration = config.get("spatial_calibration", config.get("scales", {}))
        if isinstance(calibration, Mapping):
            aliases = {
                "axial_um_per_px": "axial_um_per_px",
                "lateral_um_per_px": "lateral_um_per_px",
            }
            for src, dst in aliases.items():
                if src in calibration:
                    merged[dst] = calibration[src]
    validate_tracking_config(merged)
    return merged


def validate_tracking_config(cfg: Mapping[str, Any]) -> None:
    alphas = [float(x) for x in cfg["upper_edge_alphas"]]
    if not alphas or any(not 0 < x < 1 for x in alphas):
        raise TrackingError("upper_edge_alphas must be non-empty and in (0,1)")
    if float(cfg["primary_alpha"]) not in alphas:
        raise TrackingError("primary_alpha must be included in upper_edge_alphas")
    zlo, zhi = (int(x) for x in cfg["z_search_px"])
    if zlo < 0 or zhi <= zlo:
        raise TrackingError("z_search_px must be an increasing non-negative pair")
    xlo, xhi = (float(x) for x in cfg["x_search_fraction"])
    if not 0 <= xlo < xhi <= 1:
        raise TrackingError("x_search_fraction must lie within [0,1]")
    for key in (
        "axial_um_per_px",
        "lateral_um_per_px",
        "upper_edge_noise_multiplier",
        "peak_min_snr",
        "x_confidence_min_robust_snr",
        "x_confidence_core_radius_fraction_of_diameter",
        "x_confidence_flank_width_fraction_of_diameter",
        "x_confidence_alignment_tolerance_fraction_of_diameter",
        "confidence_full_score_multiple",
    ):
        if float(cfg[key]) <= 0:
            raise TrackingError(f"{key} must be positive")


def _row_robust_zscore(score: np.ndarray) -> np.ndarray:
    med = np.nanmedian(score, axis=1, keepdims=True)
    sig = robust_sigma(score, axis=1)[:, None]
    fallback = np.nanstd(score, axis=1, keepdims=True)
    sig = np.where(sig > np.finfo(float).eps, sig, fallback)
    sig = np.where(sig > np.finfo(float).eps, sig, 1.0)
    return np.clip((score - med) / sig, -8.0, 30.0)


def viterbi_continuous_path(
    score: np.ndarray, max_jump: int, jump_penalty: float
) -> np.ndarray:
    """Maximum-score path with a finite slow-axis jump neighbourhood."""
    if score.ndim != 2 or score.shape[0] < 1 or score.shape[1] < 2:
        raise TrackingError(f"Viterbi score must be [frame,state], got {score.shape}")
    n, states = score.shape
    clean = np.nan_to_num(score.astype(np.float64), nan=-1e6, neginf=-1e6, posinf=1e6)
    back = np.zeros((n, states), dtype=np.int16 if states < 32767 else np.int32)
    previous = clean[0].copy()
    state_idx = np.arange(states)
    for frame in range(1, n):
        best = np.full(states, -np.inf)
        best_prev = np.zeros(states, dtype=back.dtype)
        for delta in range(-int(max_jump), int(max_jump) + 1):
            target_start = max(0, delta)
            target_stop = min(states, states + delta)
            prev_start = target_start - delta
            prev_stop = target_stop - delta
            candidate = previous[prev_start:prev_stop] - jump_penalty * (delta**2)
            replace = candidate > best[target_start:target_stop]
            if np.any(replace):
                sl = best[target_start:target_stop]
                sl[replace] = candidate[replace]
                best[target_start:target_stop] = sl
                bp = best_prev[target_start:target_stop]
                bp[replace] = state_idx[prev_start:prev_stop][replace]
                best_prev[target_start:target_stop] = bp
        previous = clean[frame] + best
        back[frame] = best_prev
    path = np.empty(n, dtype=np.int32)
    path[-1] = int(np.argmax(previous))
    for frame in range(n - 1, 0, -1):
        path[frame - 1] = back[frame, path[frame]]
    return path


def _odd_window(value: int, maximum: int) -> int:
    value = max(1, min(int(value), int(maximum)))
    if value % 2 == 0:
        value = value - 1 if value == maximum else value + 1
    return max(1, value)


def locate_lateral_track(
    volume: np.ndarray, diameter_um: float, cfg: Mapping[str, Any]
) -> dict[str, np.ndarray]:
    nframe, nz, nx = volume.shape
    x0 = max(0, int(math.floor(float(cfg["x_search_fraction"][0]) * nx)))
    x1 = min(nx, int(math.ceil(float(cfg["x_search_fraction"][1]) * nx)))
    z0 = max(0, int(cfg["z_search_px"][0]))
    z1 = min(nz, int(cfg["z_search_px"][1]))
    if x1 - x0 < 20 or z1 - z0 < 20:
        raise TrackingError("Configured x/z tracking search region is too small")
    diameter_z = diameter_um / float(cfg["axial_um_per_px"])
    core_depth = max(
        3, int(round(diameter_z * float(cfg["x_core_depth_fraction_of_diameter"])))
    )
    scores = np.empty((nframe, x1 - x0), dtype=np.float32)
    local_seed = np.empty(nframe, dtype=np.int32)
    for frame in range(nframe):
        bscan = volume[frame].astype(np.float32, copy=False)
        smooth = gaussian_filter(
            bscan, sigma=(float(cfg["localization_depth_sigma_px"]), 1.0)
        )
        lateral_background = np.nanmedian(smooth[:, x0:x1], axis=1)
        excess = np.maximum(smooth[:, x0:x1] - lateral_background[:, None], 0.0)
        vertical_core = uniform_filter1d(excess, size=core_depth, axis=0, mode="nearest")
        scores[frame] = np.nanmax(vertical_core[z0:z1], axis=0)
        local_seed[frame] = int(np.nanargmax(scores[frame])) + x0
    standardized = _row_robust_zscore(scores)
    relative_path = viterbi_continuous_path(
        standardized,
        int(cfg["x_viterbi_max_jump_px"]),
        float(cfg["x_viterbi_jump_penalty"]),
    )
    path = relative_path.astype(float) + x0
    smoothed = median_filter(
        path,
        size=_odd_window(int(cfg["trajectory_median_window_frames"]), nframe),
        mode="nearest",
    )
    smoothed = gaussian_filter1d(
        smoothed, sigma=float(cfg["trajectory_gaussian_sigma_frames"]), mode="nearest"
    )
    best = np.nanmax(scores, axis=1)
    base = np.nanmedian(scores, axis=1)
    sig = robust_sigma(scores, axis=1)
    sig = np.where(sig > 0, sig, np.nanstd(scores, axis=1))
    local_snr = (best - base) / np.maximum(sig, np.finfo(np.float32).eps)
    residual = local_seed.astype(float) - smoothed
    seed_high = (
        (local_snr >= float(cfg["x_seed_min_robust_snr"]))
        & (np.abs(residual) <= float(cfg["x_seed_max_path_residual_px"]))
    )
    min_seed = max(
        int(cfg["trajectory_min_seed_frames"]),
        int(math.ceil(float(cfg["trajectory_min_seed_fraction"]) * nframe)),
    )
    if int(seed_high.sum()) < min_seed:
        raise TrackingError(
            "Insufficient high-confidence lateral anchors for a defensible 3-D track: "
            f"{int(seed_high.sum())}/{nframe}, required={min_seed}"
        )
    distance = distance_transform_edt(~seed_high)
    supported = distance <= int(cfg["max_model_assisted_gap_frames"])

    # New x evidence: the Viterbi path defines only the centre of a local
    # inspection window.  Evidence itself comes from the unsmoothed image
    # score within a diameter-defined central band versus bilateral local
    # flanks.  No full-field per-frame argmax participates in this confidence.
    lateral_diameter_px = diameter_um / float(cfg["lateral_um_per_px"])
    core_radius = max(
        1,
        int(
            math.ceil(
                lateral_diameter_px
                * float(cfg["x_confidence_core_radius_fraction_of_diameter"])
            )
        ),
    )
    flank_width = max(
        int(cfg["x_confidence_min_flank_points"]),
        int(
            math.ceil(
                lateral_diameter_px
                * float(cfg["x_confidence_flank_width_fraction_of_diameter"])
            )
        ),
    )
    path_change = np.zeros(nframe, dtype=float)
    if nframe > 1:
        steps = np.abs(np.diff(smoothed))
        path_change[0] = steps[0]
        path_change[-1] = steps[-1]
        if nframe > 2:
            path_change[1:-1] = np.maximum(steps[:-1], steps[1:])
    trajectory_pad = np.maximum(
        int(cfg["x_confidence_min_trajectory_pad_px"]),
        np.ceil(
            float(cfg["x_confidence_trajectory_change_multiplier"]) * path_change
        ).astype(int),
    )
    local_half_width = core_radius + flank_width + trajectory_pad
    alignment_tolerance = max(
        float(cfg["x_confidence_min_alignment_tolerance_px"]),
        lateral_diameter_px
        * float(cfg["x_confidence_alignment_tolerance_fraction_of_diameter"]),
    )
    local_peak = np.full(nframe, np.nan, dtype=float)
    local_peak_score = np.full(nframe, np.nan, dtype=float)
    local_baseline = np.full(nframe, np.nan, dtype=float)
    local_noise = np.full(nframe, np.nan, dtype=float)
    local_snr = np.full(nframe, np.nan, dtype=float)
    local_residual = np.full(nframe, np.nan, dtype=float)
    local_window_start = np.full(nframe, -1, dtype=np.int32)
    local_window_stop = np.full(nframe, -1, dtype=np.int32)
    local_available = np.zeros(nframe, dtype=bool)
    for frame in range(nframe):
        xc = float(smoothed[frame])
        half = int(local_half_width[frame])
        win0 = max(x0, int(math.floor(xc - half)))
        win1 = min(x1, int(math.ceil(xc + half)) + 1)
        core0 = max(win0, int(math.floor(xc - core_radius)))
        core1 = min(win1, int(math.ceil(xc + core_radius)) + 1)
        local_window_start[frame] = win0
        local_window_stop[frame] = win1
        if core1 <= core0:
            continue
        central_score = scores[frame, core0 - x0 : core1 - x0]
        if not np.isfinite(central_score).any():
            continue
        peak_relative = int(np.nanargmax(central_score))
        peak_x = core0 + peak_relative
        flank = np.concatenate(
            (
                scores[frame, win0 - x0 : core0 - x0],
                scores[frame, core1 - x0 : win1 - x0],
            )
        )
        flank = flank[np.isfinite(flank)]
        if flank.size < int(cfg["x_confidence_min_flank_points"]):
            continue
        baseline = float(np.nanmedian(flank))
        noise = float(robust_sigma(flank))
        if not np.isfinite(noise) or noise <= np.finfo(np.float32).eps:
            noise = float(np.nanstd(flank))
        if not np.isfinite(noise) or noise <= np.finfo(np.float32).eps:
            continue
        peak_value = float(scores[frame, peak_x - x0])
        local_peak[frame] = peak_x
        local_peak_score[frame] = peak_value
        local_baseline[frame] = baseline
        local_noise[frame] = noise
        local_snr[frame] = (peak_value - baseline) / noise
        local_residual[frame] = peak_x - xc
        local_available[frame] = True
    local_strong = (
        local_available
        & (local_snr >= float(cfg["x_confidence_min_robust_snr"]))
        & (np.abs(local_residual) <= alignment_tolerance)
    )
    full_multiple = float(cfg["confidence_full_score_multiple"])
    snr_quality = np.clip(
        local_snr
        / (
            full_multiple
            * max(float(cfg["x_confidence_min_robust_snr"]), np.finfo(float).eps)
        ),
        0.0,
        1.0,
    )
    alignment_quality = np.clip(
        1.0 - np.abs(local_residual) / (full_multiple * alignment_tolerance),
        0.0,
        1.0,
    )
    local_confidence = np.sqrt(snr_quality * alignment_quality)
    local_confidence[~local_available] = 0.0
    return {
        "x_center": smoothed,
        "x_local_seed": local_seed.astype(float),
        "x_seed_high": seed_high,
        "x_supported": supported,
        "x_robust_snr": local_snr,
        "x_residual": residual,
        "x_score": best,
        "x_local_evidence_available": local_available,
        "x_local_evidence_strong": local_strong,
        "x_path_confidence": local_confidence,
        "x_local_peak": local_peak,
        "x_local_peak_score": local_peak_score,
        "x_local_baseline": local_baseline,
        "x_local_noise_sigma": local_noise,
        "x_local_robust_snr": local_snr,
        "x_local_path_residual": local_residual,
        "x_alignment_tolerance_px": np.full(nframe, alignment_tolerance),
        "x_local_window_half_width_px": local_half_width.astype(float),
        "x_local_window_start_px": local_window_start,
        "x_local_window_stop_px": local_window_stop,
        "x_local_path_change_px": path_change,
    }


def extract_profiles(
    volume: np.ndarray,
    x_center: np.ndarray,
    diameter_um: float,
    cfg: Mapping[str, Any],
) -> ProfileBundle:
    nframe, nz, nx = volume.shape
    lateral_diameter = diameter_um / float(cfg["lateral_um_per_px"])
    width = max(
        int(cfg["central_min_width_px"]),
        int(round(lateral_diameter * float(cfg["central_width_fraction_of_lateral_diameter"]))),
    )
    if width % 2 == 0:
        width += 1
    half = width // 2
    side_gap = int(cfg["background_side_gap_px"])
    side_width = int(cfg["background_side_width_px"])
    vessel_radius = int(math.ceil(lateral_diameter / 2.0))

    central = np.full((nframe, nz), np.nan, np.float32)
    background = np.full_like(central, np.nan)
    bg_sigma = np.full_like(central, np.nan)
    bounds = [np.full(nframe, -1, np.int32) for _ in range(6)]
    roi_valid = np.zeros(nframe, dtype=bool)
    for frame in range(nframe):
        xc = int(round(float(x_center[frame])))
        c0, c1 = xc - half, xc + half + 1
        l1 = xc - vessel_radius - side_gap
        l0 = l1 - side_width
        r0 = xc + vessel_radius + side_gap + 1
        r1 = r0 + side_width
        bounds[0][frame], bounds[1][frame] = c0, c1
        bounds[2][frame], bounds[3][frame] = l0, l1
        bounds[4][frame], bounds[5][frame] = r0, r1
        if c0 < 0 or c1 > nx or l0 < 0 or r1 > nx:
            continue
        sides = np.concatenate((volume[frame, :, l0:l1], volume[frame, :, r0:r1]), axis=1)
        if sides.shape[1] < int(cfg["background_min_total_width_px"]):
            continue
        central[frame] = np.nanmedian(volume[frame, :, c0:c1], axis=1)
        background[frame] = np.nanmedian(sides, axis=1)
        bg_sigma[frame] = robust_sigma(sides, axis=1)
        roi_valid[frame] = np.isfinite(central[frame]).all() and np.isfinite(
            background[frame]
        ).all()
    excess = np.maximum(central - background, 0.0)
    loc = gaussian_filter(
        np.nan_to_num(excess, nan=0.0),
        sigma=(
            float(cfg["localization_slow_sigma_frames"]),
            float(cfg["localization_depth_sigma_px"]),
        ),
        mode="nearest",
    )
    return ProfileBundle(
        central=central,
        background=background,
        excess=excess,
        excess_localization=loc,
        background_sigma=bg_sigma,
        central_left=bounds[0],
        central_right=bounds[1],
        side_left_start=bounds[2],
        side_left_stop=bounds[3],
        side_right_start=bounds[4],
        side_right_stop=bounds[5],
        roi_valid=roi_valid,
    )


def locate_peak_path(profiles: ProfileBundle, cfg: Mapping[str, Any]) -> np.ndarray:
    nframe, nz = profiles.excess_localization.shape
    z0 = max(0, int(cfg["z_search_px"][0]))
    z1 = min(nz, int(cfg["z_search_px"][1]))
    noise = np.maximum(profiles.background_sigma[:, z0:z1], np.finfo(np.float32).eps)
    score = profiles.excess_localization[:, z0:z1] / noise
    score = _row_robust_zscore(score)
    relative = viterbi_continuous_path(
        score,
        int(cfg["z_viterbi_max_jump_px"]),
        float(cfg["z_viterbi_jump_penalty"]),
    )
    return relative + z0


def _component_containing(mask: np.ndarray, index: int) -> tuple[int, int] | None:
    if index < 0 or index >= mask.size or not bool(mask[index]):
        return None
    start = index
    stop = index + 1
    while start > 0 and mask[start - 1]:
        start -= 1
    while stop < mask.size and mask[stop]:
        stop += 1
    return start, stop


def _robust_trajectory(
    candidate: np.ndarray,
    initially_valid: np.ndarray,
    cfg: Mapping[str, Any],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    n = candidate.size
    accepted = initially_valid & np.isfinite(candidate)
    min_seed = max(
        int(cfg["trajectory_min_seed_frames"]),
        int(math.ceil(float(cfg["trajectory_min_seed_fraction"]) * n)),
    )
    if int(accepted.sum()) < min_seed:
        raise TrackingError(
            "Insufficient upper-edge anchors for slow-axis modelling: "
            f"{int(accepted.sum())}/{n}, required={min_seed}"
        )
    positions = np.arange(n)
    model = np.full(n, np.nan)
    for _ in range(3):
        good = np.flatnonzero(accepted)
        filled = np.interp(positions, good, candidate[good])
        window = _odd_window(int(cfg["trajectory_median_window_frames"]), n)
        model = median_filter(filled, size=window, mode="nearest")
        model = gaussian_filter1d(
            model, sigma=float(cfg["trajectory_gaussian_sigma_frames"]), mode="nearest"
        )
        residual = candidate - model
        sigma = float(robust_sigma(residual[accepted]))
        tolerance = max(
            float(cfg["trajectory_outlier_min_px"]),
            float(cfg["trajectory_outlier_mad_multiplier"]) * sigma,
        )
        revised = initially_valid & np.isfinite(candidate) & (np.abs(residual) <= tolerance)
        if int(revised.sum()) < min_seed:
            break
        if np.array_equal(revised, accepted):
            accepted = revised
            break
        accepted = revised
    residual = candidate - model
    support_distance = distance_transform_edt(~accepted)
    supported = support_distance <= int(cfg["max_model_assisted_gap_frames"])
    return model, accepted, residual, support_distance


def track_one_alpha(
    *,
    scan_id: str,
    diameter_um: float,
    alpha: float,
    xtrack: Mapping[str, np.ndarray],
    peak_path: np.ndarray,
    profiles: ProfileBundle,
    cfg: Mapping[str, Any],
) -> pd.DataFrame:
    nframe, nz = profiles.excess.shape
    seed_upper = np.full(nframe, np.nan)
    peak_z = np.full(nframe, np.nan)
    peak_excess = np.full(nframe, np.nan)
    noise_sigma = np.full(nframe, np.nan)
    threshold = np.full(nframe, np.nan)
    component_width = np.zeros(nframe, dtype=np.int32)
    seed_valid = np.zeros(nframe, dtype=bool)
    peak_snr = np.full(nframe, np.nan)
    min_component = int(cfg["upper_edge_min_component_px"])
    diameter_px_float = diameter_um / float(cfg["axial_um_per_px"])

    for frame in range(nframe):
        if not profiles.roi_valid[frame]:
            continue
        pk0 = int(round(float(peak_path[frame])))
        lo = max(0, pk0 - 2)
        hi = min(nz, pk0 + 3)
        if hi <= lo:
            continue
        pk = lo + int(np.argmax(profiles.excess_localization[frame, lo:hi]))
        peak = float(profiles.excess_localization[frame, pk])
        core0 = max(0, int(round(pk - 0.5 * diameter_px_float)))
        core1 = min(nz, int(round(pk + 0.5 * diameter_px_float)) + 1)
        # ``B(z)`` can have a large non-zero level even when its local noise is
        # modest.  The threshold formula calls for MAD of background *noise*,
        # not the lateral dispersion of raw positive OMAG values.  Estimate it
        # from the high-frequency residual of the same-frame, same-depth flank
        # profile after removing its slow axial trend.
        bg_profile = profiles.background[frame]
        bg_residual = bg_profile - gaussian_filter1d(
            bg_profile,
            sigma=float(cfg["background_depth_trend_sigma_px"]),
            mode="nearest",
        )
        sigma = float(robust_sigma(bg_residual[core0:core1]))
        if not np.isfinite(sigma) or sigma <= 0 or not np.isfinite(peak):
            continue
        thr = max(float(alpha) * peak, float(cfg["upper_edge_noise_multiplier"]) * sigma)
        mask = profiles.excess_localization[frame] >= thr
        component = _component_containing(mask, pk)
        peak_z[frame] = pk
        peak_excess[frame] = peak
        noise_sigma[frame] = sigma
        threshold[frame] = thr
        peak_snr[frame] = peak / sigma
        if component is None:
            continue
        start, stop = component
        component_width[frame] = stop - start
        if stop - start < min_component:
            continue
        seed_upper[frame] = start
        seed_valid[frame] = peak_snr[frame] >= float(cfg["peak_min_snr"])

    model, accepted, z_residual, support_distance = _robust_trajectory(
        seed_upper, seed_valid, cfg
    )
    z_upper = np.rint(model)
    diameter_px = int(round(diameter_px_float))
    z_lower = z_upper + diameter_px
    z_start = z_lower + int(cfg["primary_guard_px"])
    window_px = int(round(float(cfg["primary_window_um"]) / float(cfg["axial_um_per_px"])))
    window_complete = (z_start >= 0) & (z_start + window_px <= nz)
    x_supported = np.asarray(xtrack["x_supported"], bool)
    supported = (
        (support_distance <= int(cfg["max_model_assisted_gap_frames"]))
        & x_supported
        & profiles.roi_valid
        & np.isfinite(z_upper)
        & (z_upper >= 0)
        & (z_lower < nz)
    )
    high = accepted & np.asarray(xtrack["x_seed_high"], bool) & supported
    assisted = (~high) & supported
    classes = np.full(nframe, CLASS_FAILED, dtype=object)
    classes[assisted] = CLASS_ASSISTED
    classes[high] = CLASS_HIGH

    # A bounded audit score; the class remains the authoritative categorical QC.
    snr_part = np.clip(
        (peak_snr - float(cfg["peak_min_snr"])) / max(float(cfg["peak_min_snr"]), 1.0),
        0.0,
        1.0,
    )
    residual_scale = max(float(cfg["trajectory_outlier_min_px"]), 1.0)
    residual_part = np.exp(-np.abs(np.nan_to_num(z_residual, nan=999.0)) / residual_scale)
    x_part = np.clip(np.asarray(xtrack["x_robust_snr"]) / 10.0, 0.0, 1.0)
    confidence = snr_part * residual_part * x_part
    confidence[assisted] = np.minimum(confidence[assisted], 0.49)
    confidence[high] = np.maximum(confidence[high], 0.50)
    confidence[classes == CLASS_FAILED] = 0.0

    # Evidence-component confidence revision.  These quantities audit the
    # unchanged final trajectory; they do not move x_center or z_upper and do
    # not change qc_valid.  The old mixed class/score remain authoritative for
    # backward compatibility during this diagnostic round.
    x_available = (
        np.asarray(xtrack["x_local_evidence_available"], bool)
        & profiles.roi_valid
        & np.isfinite(np.asarray(xtrack["x_center"], float))
    )
    x_strong = (
        np.asarray(xtrack["x_local_evidence_strong"], bool) & x_available
    )
    x_evidence_class = np.full(nframe, EVIDENCE_UNAVAILABLE, dtype=object)
    x_evidence_class[x_available & ~x_strong] = EVIDENCE_WEAK
    x_evidence_class[x_strong] = EVIDENCE_STRONG
    x_path_confidence = np.asarray(xtrack["x_path_confidence"], float).copy()
    x_path_confidence[~x_available] = 0.0

    z_candidate_present = np.isfinite(seed_upper)
    z_available = (
        z_candidate_present
        & profiles.roi_valid
        & np.isfinite(peak_snr)
        & np.isfinite(z_residual)
    )
    accepted_residual = z_residual[accepted & np.isfinite(z_residual)]
    accepted_sigma = (
        float(robust_sigma(accepted_residual)) if accepted_residual.size else math.nan
    )
    if not np.isfinite(accepted_sigma):
        accepted_sigma = 0.0
    z_residual_tolerance = max(
        float(cfg["trajectory_outlier_min_px"]),
        float(cfg["trajectory_outlier_mad_multiplier"]) * accepted_sigma,
    )
    z_strong = accepted & z_available
    z_evidence_class = np.full(nframe, EVIDENCE_UNAVAILABLE, dtype=object)
    z_evidence_class[z_available & ~z_strong] = EVIDENCE_WEAK
    z_evidence_class[z_strong] = EVIDENCE_STRONG
    full_multiple = float(cfg["confidence_full_score_multiple"])
    z_snr_quality = np.clip(
        peak_snr
        / (full_multiple * max(float(cfg["peak_min_snr"]), np.finfo(float).eps)),
        0.0,
        1.0,
    )
    z_component_quality = np.clip(
        component_width
        / (full_multiple * max(float(min_component), np.finfo(float).eps)),
        0.0,
        1.0,
    )
    z_residual_quality = np.clip(
        1.0
        - np.abs(z_residual)
        / (full_multiple * max(z_residual_tolerance, np.finfo(float).eps)),
        0.0,
        1.0,
    )
    z_edge_confidence = np.cbrt(
        z_snr_quality * z_component_quality * z_residual_quality
    )
    z_edge_confidence[~z_available] = 0.0
    overall_tracking_confidence = np.sqrt(
        np.clip(x_path_confidence, 0.0, 1.0)
        * np.clip(z_edge_confidence, 0.0, 1.0)
    )
    new_geometry_valid = (
        profiles.roi_valid
        & np.isfinite(np.asarray(xtrack["x_center"], float))
        & np.isfinite(z_upper)
        & (z_upper >= 0)
        & (z_lower < nz)
    )
    overall_tracking_confidence[~new_geometry_valid] = 0.0
    # Deliberately exclude legacy x_seed_high/x_supported and z-anchor support
    # distances from the revised class.  Reusing them would reintroduce the
    # full-field-argmax bias this diagnostic is designed to isolate.
    new_high = x_strong & z_strong & new_geometry_valid
    new_assisted = (~new_high) & new_geometry_valid
    new_classes = np.full(nframe, CLASS_FAILED, dtype=object)
    new_classes[new_assisted] = CLASS_ASSISTED
    new_classes[new_high] = CLASS_HIGH

    flags: list[str] = []
    for frame in range(nframe):
        item: list[str] = []
        if not profiles.roi_valid[frame]:
            item.append("ROI_OUT_OF_BOUNDS")
        if not bool(xtrack["x_seed_high"][frame]):
            item.append("X_MODEL_ASSISTED")
        if not bool(x_available[frame]):
            item.append("X_LOCAL_EVIDENCE_UNAVAILABLE")
        elif not bool(x_strong[frame]):
            item.append("X_LOCAL_EVIDENCE_WEAK")
        if not np.isfinite(seed_upper[frame]):
            item.append("NO_VALID_UPPER_COMPONENT")
        elif not accepted[frame]:
            item.append("UPPER_SEED_TRAJECTORY_OUTLIER")
        if np.isfinite(peak_snr[frame]) and peak_snr[frame] < float(cfg["peak_min_snr"]):
            item.append("LOW_PEAK_SNR")
        if not bool(z_available[frame]):
            item.append("Z_EDGE_EVIDENCE_UNAVAILABLE")
        elif not bool(z_strong[frame]):
            item.append("Z_EDGE_EVIDENCE_WEAK")
        if support_distance[frame] > int(cfg["max_model_assisted_gap_frames"]):
            item.append("LONG_UNSUPPORTED_Z_GAP")
        if not bool(x_supported[frame]):
            item.append("LONG_UNSUPPORTED_X_GAP")
        if not bool(window_complete[frame]):
            item.append("PRIMARY_WINDOW_INCOMPLETE")
        if classes[frame] == CLASS_FAILED:
            item.append("TRACKING_FAILED")
        flags.append(";".join(item))

    frame_index = np.arange(nframe, dtype=np.int32)
    return pd.DataFrame(
        {
            "scan_id": scan_id,
            "frame_index": frame_index,
            "coordinate_base": 0,
            "alpha": float(alpha),
            "diameter_um": float(diameter_um),
            "diameter_axial_px": diameter_px,
            "x_center_px": np.asarray(xtrack["x_center"]),
            "z_peak_px": peak_z,
            "z_upper_px": z_upper,
            "z_lower_op_px": z_lower,
            "z_tail_start_px": z_start,
            "tracking_class": classes,
            "tracking_confidence": confidence,
            "old_tracking_class": classes,
            "old_tracking_confidence": confidence,
            "new_tracking_class": new_classes,
            "overall_tracking_confidence": overall_tracking_confidence,
            "x_path_confidence_class": x_evidence_class,
            "x_path_confidence": x_path_confidence,
            "x_local_evidence_available": x_available,
            "x_local_evidence_strong": x_strong,
            "x_local_peak_px": np.asarray(xtrack["x_local_peak"]),
            "x_local_peak_score": np.asarray(xtrack["x_local_peak_score"]),
            "x_local_baseline_score": np.asarray(xtrack["x_local_baseline"]),
            "x_local_noise_sigma": np.asarray(xtrack["x_local_noise_sigma"]),
            "x_local_peak_snr": np.asarray(xtrack["x_local_robust_snr"]),
            "x_local_path_residual_px": np.asarray(
                xtrack["x_local_path_residual"]
            ),
            "x_alignment_tolerance_px": np.asarray(
                xtrack["x_alignment_tolerance_px"]
            ),
            "x_local_window_half_width_px": np.asarray(
                xtrack["x_local_window_half_width_px"]
            ),
            "x_local_window_x0_px": np.asarray(xtrack["x_local_window_start_px"]),
            "x_local_window_x1_exclusive_px": np.asarray(
                xtrack["x_local_window_stop_px"]
            ),
            "x_local_path_change_px": np.asarray(
                xtrack["x_local_path_change_px"]
            ),
            "z_edge_confidence_class": z_evidence_class,
            "z_edge_confidence": z_edge_confidence,
            "z_candidate_present": z_candidate_present,
            "z_candidate_accepted": accepted,
            "z_residual_tolerance_px": z_residual_tolerance,
            "peak_excess": peak_excess,
            "background_sigma": noise_sigma,
            "peak_snr": peak_snr,
            "upper_threshold": threshold,
            "upper_component_width_px": component_width,
            "seed_x_px": np.asarray(xtrack["x_local_seed"]),
            "seed_z_upper_px": seed_upper,
            "model_z_upper_px": model,
            "residual_x_px": np.asarray(xtrack["x_residual"]),
            "residual_z_px": z_residual,
            "distance_to_z_anchor_frames": support_distance,
            "x_robust_snr": np.asarray(xtrack["x_robust_snr"]),
            "central_x0_px": profiles.central_left,
            "central_x1_exclusive_px": profiles.central_right,
            "side_left_x0_px": profiles.side_left_start,
            "side_left_x1_exclusive_px": profiles.side_left_stop,
            "side_right_x0_px": profiles.side_right_start,
            "side_right_x1_exclusive_px": profiles.side_right_stop,
            "window_300_fits": window_complete,
            "qc_valid": classes != CLASS_FAILED,
            "qc_flags": flags,
        }
    )


def track_volume(
    volume: np.ndarray,
    *,
    scan_id: str,
    diameter_um: float,
    config: Mapping[str, Any] | None = None,
) -> tuple[dict[float, pd.DataFrame], ProfileBundle, dict[str, np.ndarray]]:
    cfg = merge_tracking_config(config)
    if volume.ndim != 3:
        raise TrackingError(f"Expected [frame,z,x], got {volume.shape}")
    if not np.isfinite(volume).all():
        raise TrackingError("Tracking input contains non-finite values")
    if volume.shape[0] < int(cfg["trajectory_min_seed_frames"]):
        raise TrackingError("Too few frames to establish a slow-axis trajectory")
    z1 = int(cfg["z_search_px"][1])
    if z1 > volume.shape[1]:
        raise TrackingError(
            f"z_search_px ends at {z1}, beyond input depth {volume.shape[1]}"
        )
    xtrack = locate_lateral_track(volume, diameter_um, cfg)
    profiles = extract_profiles(volume, xtrack["x_center"], diameter_um, cfg)
    if int(profiles.roi_valid.sum()) < int(cfg["trajectory_min_seed_frames"]):
        raise TrackingError(
            f"Only {int(profiles.roi_valid.sum())}/{volume.shape[0]} frames have valid "
            "central and bilateral background ROIs"
        )
    peak_path = locate_peak_path(profiles, cfg)
    xtrack = dict(xtrack)
    xtrack["z_peak_path"] = peak_path.astype(float)
    outputs: dict[float, pd.DataFrame] = {}
    for alpha in [float(x) for x in cfg["upper_edge_alphas"]]:
        outputs[alpha] = track_one_alpha(
            scan_id=scan_id,
            diameter_um=diameter_um,
            alpha=alpha,
            xtrack=xtrack,
            peak_path=peak_path,
            profiles=profiles,
            cfg=cfg,
        )
    return outputs, profiles, xtrack


def alpha_tag(alpha: float) -> str:
    return f"{int(round(alpha * 100)):03d}"


def dataframe_summary(df: pd.DataFrame) -> dict[str, Any]:
    n = int(len(df))
    counts = df["tracking_class"].value_counts().to_dict()
    high = int(counts.get(CLASS_HIGH, 0))
    assisted = int(counts.get(CLASS_ASSISTED, 0))
    failed = int(counts.get(CLASS_FAILED, 0))
    new_counts = df["new_tracking_class"].value_counts().to_dict()
    x_counts = df["x_path_confidence_class"].value_counts().to_dict()
    z_counts = df["z_edge_confidence_class"].value_counts().to_dict()
    return {
        "scan_id": str(df["scan_id"].iloc[0]),
        "alpha": float(df["alpha"].iloc[0]),
        "frame_count": n,
        "high_confidence_frames": high,
        "model_assisted_frames": assisted,
        "failed_frames": failed,
        "high_confidence_fraction": high / n if n else math.nan,
        "model_assisted_fraction": assisted / n if n else math.nan,
        "failed_fraction": failed / n if n else math.nan,
        "old_high_confidence_fraction": high / n if n else math.nan,
        "old_model_assisted_fraction": assisted / n if n else math.nan,
        "old_failed_fraction": failed / n if n else math.nan,
        "new_high_confidence_fraction": int(new_counts.get(CLASS_HIGH, 0)) / n
        if n
        else math.nan,
        "new_model_assisted_fraction": int(new_counts.get(CLASS_ASSISTED, 0)) / n
        if n
        else math.nan,
        "new_failed_fraction": int(new_counts.get(CLASS_FAILED, 0)) / n
        if n
        else math.nan,
        "x_strong_evidence_fraction": int(x_counts.get(EVIDENCE_STRONG, 0)) / n
        if n
        else math.nan,
        "x_weak_evidence_fraction": int(x_counts.get(EVIDENCE_WEAK, 0)) / n
        if n
        else math.nan,
        "x_unavailable_fraction": int(x_counts.get(EVIDENCE_UNAVAILABLE, 0)) / n
        if n
        else math.nan,
        "z_strong_evidence_fraction": int(z_counts.get(EVIDENCE_STRONG, 0)) / n
        if n
        else math.nan,
        "z_weak_evidence_fraction": int(z_counts.get(EVIDENCE_WEAK, 0)) / n
        if n
        else math.nan,
        "z_unavailable_fraction": int(z_counts.get(EVIDENCE_UNAVAILABLE, 0)) / n
        if n
        else math.nan,
        "x_path_confidence_median": float(df["x_path_confidence"].median()),
        "z_edge_confidence_median": float(df["z_edge_confidence"].median()),
        "overall_tracking_confidence_median": float(
            df["overall_tracking_confidence"].median()
        ),
        "qc_valid_fraction": float(df["qc_valid"].mean()) if n else math.nan,
        "window_300_complete_frames": int(df["window_300_fits"].sum()),
        "window_300_complete_fraction": float(df["window_300_fits"].mean()) if n else math.nan,
        "z_upper_median_px": float(df.loc[df.qc_valid, "z_upper_px"].median()),
        "z_upper_iqr_px": float(
            np.subtract(*np.nanpercentile(df.loc[df.qc_valid, "z_upper_px"], [75, 25]))
        ),
        "x_center_median_px": float(df.loc[df.qc_valid, "x_center_px"].median()),
        "peak_snr_median": float(df.loc[df.qc_valid, "peak_snr"].median()),
    }


def code_sha256(paths: Iterable[str | Path]) -> str:
    digest = hashlib.sha256()
    for path in sorted(Path(x).resolve() for x in paths):
        digest.update(path.name.encode("utf-8"))
        digest.update(path.read_bytes())
    return digest.hexdigest()
