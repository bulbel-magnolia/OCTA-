"""OMAG-only localization feeding the formal SV geometry."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.ndimage import gaussian_filter1d

from .geometry import VesselGeometry


FloatArray = NDArray[np.float64]


def _runs(mask: NDArray[np.bool_]) -> list[tuple[int, int]]:
    """Return true runs as half-open index intervals."""

    padded = np.pad(np.asarray(mask, dtype=np.int8), (1, 1))
    changes = np.diff(padded)
    starts = np.flatnonzero(changes == 1)
    stops = np.flatnonzero(changes == -1)
    return list(zip(starts.tolist(), stops.tolist(), strict=True))


@dataclass(frozen=True)
class SurfaceGuidedAnchorResult:
    valid: bool
    invalid_reason: str
    surface_z_center_px: float
    surface_to_vessel_top_um: float
    effective_refractive_index: float
    predicted_z_top_center_px: float
    x_anchor_center_px: float
    search_z_start_px: int
    search_z_stop_px: int
    peak_score: float
    background_score: float
    background_sigma: float
    peak_cnr: float

    @property
    def qc_valid(self) -> bool:
        return self.valid


@dataclass(frozen=True)
class LocalBodyResult:
    valid: bool
    invalid_reason: str
    x_anchor_center_px: float
    x_left_center_px: int
    x_right_center_px: int
    x1_local_geometry_center_px: float
    fallback: bool
    run_width_px: int
    expected_lateral_width_px: int
    local_background: float
    local_sigma: float
    peak_cnr: float
    axial_completeness: float

    @property
    def qc_valid(self) -> bool:
        minimum_width = max(2, int(np.floor(0.3 * self.expected_lateral_width_px)))
        return bool(
            self.valid
            and not self.fallback
            and self.run_width_px >= minimum_width
            and self.peak_cnr >= 1.5
            and self.axial_completeness >= 0.25
        )


@dataclass(frozen=True)
class TopEdgeResult:
    valid: bool
    invalid_reason: str
    z_anchor_center_px: float
    z_top_center_px: float
    fallback: bool
    local_background: float
    local_sigma: float
    peak_cnr: float

    @property
    def qc_valid(self) -> bool:
        return bool(self.valid and not self.fallback and self.peak_cnr >= 1.5)


@dataclass(frozen=True)
class MentorTrackingEvidence:
    """Full-volume mentor tracking evidence attached to one selected frame."""

    scan_id: str
    frame_index: int
    primary_alpha: float
    viterbi_x_center_px: float
    z_upper_px: float
    tracking_class: str
    x_path_confidence_class: str
    z_edge_confidence_class: str
    x2_robust_centroid_px: float
    x4_centroid_isolated_jump_corrected_px: float
    x4_jump_corrected: bool
    assessability_score: float
    vessel_presence_prediction: str
    valid_local_body: bool
    x1_fallback: bool
    qc_valid: bool
    invalid_reason: str


@dataclass(frozen=True)
class LocalizationResult:
    geometry: VesselGeometry
    lateral: LocalBodyResult
    top_edge: TopEdgeResult
    coarse_anchor: SurfaceGuidedAnchorResult | None = None
    mentor_tracking: MentorTrackingEvidence | None = None

    @property
    def source_qc_valid(self) -> bool:
        if self.mentor_tracking is not None:
            return self.mentor_tracking.qc_valid
        coarse_valid = (
            True if self.coarse_anchor is None else self.coarse_anchor.qc_valid
        )
        return coarse_valid and self.lateral.qc_valid and self.top_edge.qc_valid


def find_surface_guided_anchor(
    omag_frame: NDArray[np.floating],
    *,
    surface_z_center_px: float,
    surface_to_vessel_top_um: float,
    effective_refractive_index: float,
    diameter_um: float,
    dx_um: float,
    dz_um: float,
    axial_margin_above_px: int = 6,
    axial_margin_below_px: int = 10,
    lateral_smoothing_sigma_px: float = 2.5,
    edge_exclusion_px: int = 12,
    minimum_peak_cnr: float = 5.0,
) -> SurfaceGuidedAnchorResult:
    """Find a global lateral anchor in the known surface-relative depth band.

    The search uses a log copy of the co-registered OMAG frame. Each depth row
    is standardized by its lateral median and MAD so that the continuous vessel
    response outranks depth-dependent background and isolated bright pixels.
    Formal SV integration remains unchanged and uses linear sv_raw.
    """

    frame = np.asarray(omag_frame, dtype=np.float64)
    if frame.ndim != 2:
        raise ValueError("omag_frame must be a 2-D depth-by-A-line array")
    parameters = np.asarray(
        [
            surface_z_center_px,
            surface_to_vessel_top_um,
            effective_refractive_index,
            diameter_um,
            dx_um,
            dz_um,
            lateral_smoothing_sigma_px,
            minimum_peak_cnr,
        ],
        dtype=np.float64,
    )
    if not np.isfinite(parameters).all():
        raise ValueError("surface-guided localization parameters must be finite")
    if (
        surface_z_center_px < 0
        or surface_to_vessel_top_um <= 0
        or effective_refractive_index <= 0
        or diameter_um <= 0
        or dx_um <= 0
        or dz_um <= 0
        or axial_margin_above_px < 0
        or axial_margin_below_px < 0
        or lateral_smoothing_sigma_px <= 0
        or edge_exclusion_px < 0
        or minimum_peak_cnr < 0
    ):
        raise ValueError("surface-guided localization parameters are out of range")

    n_z, n_x = frame.shape
    predicted_top = (
        float(surface_z_center_px)
        + float(surface_to_vessel_top_um)
        * float(effective_refractive_index)
        / float(dz_um)
    )
    diameter_rows = max(3, int(round(float(diameter_um) / float(dz_um))))
    search_start = max(0, int(np.floor(predicted_top - axial_margin_above_px)))
    search_stop = min(
        n_z,
        int(np.ceil(predicted_top + diameter_rows + axial_margin_below_px)),
    )
    if search_stop - search_start < 3:
        raise ValueError("surface-guided depth band is outside the image")
    if 2 * edge_exclusion_px >= n_x - 2:
        raise ValueError("edge exclusion leaves no lateral search region")

    band = frame[search_start:search_stop]
    finite = np.isfinite(band)
    usable_rows = finite.any(axis=1)
    if not usable_rows.any():
        raise ValueError("surface-guided depth band has no finite OMAG values")
    band = band[usable_rows]
    finite = finite[usable_rows]
    transformed = np.full_like(band, np.nan, dtype=np.float64)
    transformed[finite] = np.log1p(np.maximum(band[finite], 0.0))
    row_background = np.nanmedian(transformed, axis=1, keepdims=True)
    row_mad = np.nanmedian(
        np.abs(transformed - row_background), axis=1, keepdims=True
    )
    row_sigma = np.maximum(1.4826 * row_mad, np.finfo(float).eps)
    standardized = (transformed - row_background) / row_sigma
    positive_support = np.clip(
        np.nan_to_num(standardized, nan=0.0, posinf=8.0, neginf=0.0),
        0.0,
        8.0,
    )
    column_score = np.mean(positive_support, axis=0)
    smoothed_score = gaussian_filter1d(
        column_score,
        sigma=float(lateral_smoothing_sigma_px),
        mode="nearest",
    )

    allowed = np.ones(n_x, dtype=bool)
    if edge_exclusion_px:
        allowed[:edge_exclusion_px] = False
        allowed[n_x - edge_exclusion_px :] = False
    ranked_score = np.where(allowed, smoothed_score, -np.inf)
    x_anchor = int(np.argmax(ranked_score))

    expected_width = max(3, int(round(float(diameter_um) / float(dx_um))))
    background_mask = allowed.copy()
    exclusion = max(3, 2 * expected_width)
    background_mask[
        max(0, x_anchor - exclusion) : min(n_x, x_anchor + exclusion + 1)
    ] = False
    background_values = smoothed_score[background_mask]
    if background_values.size < 5:
        raise ValueError("too few lateral background columns for coarse localization")
    background = float(np.median(background_values))
    background_mad = float(np.median(np.abs(background_values - background)))
    background_sigma = max(1.4826 * background_mad, np.finfo(float).eps)
    peak_score = float(smoothed_score[x_anchor])
    peak_cnr = (peak_score - background) / background_sigma
    valid = bool(np.isfinite(peak_cnr) and peak_cnr >= minimum_peak_cnr)
    return SurfaceGuidedAnchorResult(
        valid=valid,
        invalid_reason="ok" if valid else "coarse_peak_cnr_below_threshold",
        surface_z_center_px=float(surface_z_center_px),
        surface_to_vessel_top_um=float(surface_to_vessel_top_um),
        effective_refractive_index=float(effective_refractive_index),
        predicted_z_top_center_px=float(predicted_top),
        x_anchor_center_px=float(x_anchor),
        search_z_start_px=search_start,
        search_z_stop_px=search_stop,
        peak_score=peak_score,
        background_score=background,
        background_sigma=background_sigma,
        peak_cnr=float(peak_cnr),
    )


def _invalid_lateral(reason: str, x_anchor: float, expected_width: int) -> LocalBodyResult:
    return LocalBodyResult(
        valid=False,
        invalid_reason=reason,
        x_anchor_center_px=float(x_anchor),
        x_left_center_px=-1,
        x_right_center_px=-1,
        x1_local_geometry_center_px=float("nan"),
        fallback=True,
        run_width_px=0,
        expected_lateral_width_px=expected_width,
        local_background=float("nan"),
        local_sigma=float("nan"),
        peak_cnr=float("nan"),
        axial_completeness=0.0,
    )


def localize_lateral_body(
    omag_frame: NDArray[np.floating],
    *,
    x_anchor_center_px: float,
    z_top_anchor_center_px: float,
    diameter_um: float,
    dx_um: float,
    dz_um: float,
    local_half_width_px: int | None = None,
    threshold_sigma: float = 1.5,
) -> LocalBodyResult:
    """Port the mentor X1 local-body rule and expose its selected boundaries."""

    frame = np.asarray(omag_frame, dtype=np.float64)
    if frame.ndim != 2:
        raise ValueError("omag_frame must be a 2-D depth-by-A-line array")
    expected_width = max(3, int(round(diameter_um / dx_um)))
    if not np.isfinite([x_anchor_center_px, z_top_anchor_center_px]).all():
        return _invalid_lateral("missing_algorithm_anchor", x_anchor_center_px, expected_width)
    n_z, n_x = frame.shape
    diameter_rows = max(1, int(round(diameter_um / dz_um)))
    half_width = (
        int(local_half_width_px)
        if local_half_width_px is not None
        else max(8, int(round(1.5 * expected_width)))
    )
    anchor = int(round(x_anchor_center_px))
    x0 = max(0, anchor - half_width)
    x1 = min(n_x, anchor + half_width + 1)
    z0 = int(round(z_top_anchor_center_px))
    z1 = min(n_z, z0 + diameter_rows)
    if z0 < 0 or z1 <= z0 or x1 <= x0:
        return _invalid_lateral("body_slab_out_of_bounds", x_anchor_center_px, expected_width)
    body = frame[z0:z1, x0:x1]
    if not np.isfinite(body).any():
        return _invalid_lateral("body_slab_nonfinite", x_anchor_center_px, expected_width)

    column_score = np.nanquantile(body, 0.75, axis=0)
    x_values = np.arange(x0, x1, dtype=np.float64)
    exclusion_half_width = max(3, expected_width // 2 + 2)
    far_mask = np.abs(x_values - x_anchor_center_px) > exclusion_half_width
    far_mask &= np.isfinite(column_score)
    if int(far_mask.sum()) < 5:
        far_mask = np.isfinite(column_score)
    if int(far_mask.sum()) < 1:
        return _invalid_lateral("no_finite_local_background", x_anchor_center_px, expected_width)
    background = float(np.median(column_score[far_mask]))
    mad = float(np.median(np.abs(column_score[far_mask] - background)))
    sigma = max(1.4826 * mad, np.finfo(float).eps)
    excess = np.maximum(np.nan_to_num(column_score - background, nan=0.0), 0.0)
    threshold = background + threshold_sigma * sigma
    runs = [(a + x0, b + x0) for a, b in _runs(column_score > threshold)]

    def run_score(run: tuple[int, int]) -> float:
        start, stop = run
        return float(np.sum(excess[start - x0 : stop - x0]))

    containing = [run for run in runs if run[0] <= anchor < run[1]]
    if containing:
        selected = max(containing, key=run_score)
        fallback = False
    elif runs:
        selected = max(runs, key=run_score)
        fallback = True
    else:
        selected = (max(0, min(n_x - 1, anchor)), max(1, min(n_x, anchor + 1)))
        fallback = True

    run_width = max(1, selected[1] - selected[0])
    x_center = (selected[0] + selected[1] - 1) / 2.0
    peak_excess = float(np.max(excess)) if excess.size else 0.0
    peak_cnr = peak_excess / sigma

    center_half = max(1, expected_width // 4)
    center_int = int(round(x_center))
    c0 = max(x0, center_int - center_half)
    c1 = min(x1, center_int + center_half + 1)
    far_body = body[:, far_mask]
    row_background = np.nanmedian(far_body, axis=1)
    row_mad = np.nanmedian(np.abs(far_body - row_background[:, None]), axis=1)
    row_threshold = row_background + threshold_sigma * 1.4826 * row_mad
    row_signal = np.nanmedian(frame[z0:z1, c0:c1], axis=1)
    finite_rows = np.isfinite(row_signal) & np.isfinite(row_threshold)
    axial_completeness = (
        float(np.mean(row_signal[finite_rows] > row_threshold[finite_rows]))
        if finite_rows.any()
        else 0.0
    )
    return LocalBodyResult(
        valid=True,
        invalid_reason="ok",
        x_anchor_center_px=float(x_anchor_center_px),
        x_left_center_px=int(selected[0]),
        x_right_center_px=int(selected[1] - 1),
        x1_local_geometry_center_px=float(x_center),
        fallback=fallback,
        run_width_px=run_width,
        expected_lateral_width_px=expected_width,
        local_background=background,
        local_sigma=sigma,
        peak_cnr=peak_cnr,
        axial_completeness=axial_completeness,
    )


def reestablish_z_top(
    omag_frame: NDArray[np.floating],
    *,
    x_center_px: float,
    z_anchor_center_px: float,
    diameter_um: float,
    dz_um: float,
    centre_half_width_px: int = 1,
    threshold_sigma: float = 1.5,
) -> TopEdgeResult:
    """Re-establish the top row on the local X1 centre line.

    The supplied tracking position limits the search. A robust upper-flank
    baseline defines the threshold, and the nearest sustained rising run gives
    the new top edge.
    """

    frame = np.asarray(omag_frame, dtype=np.float64)
    if frame.ndim != 2:
        raise ValueError("omag_frame must be a 2-D depth-by-A-line array")
    if not np.isfinite([x_center_px, z_anchor_center_px]).all():
        return TopEdgeResult(
            False, "missing_algorithm_anchor", z_anchor_center_px, float("nan"),
            True, float("nan"), float("nan"), float("nan")
        )
    n_z, n_x = frame.shape
    x_int = int(round(x_center_px))
    x0 = max(0, x_int - centre_half_width_px)
    x1 = min(n_x, x_int + centre_half_width_px + 1)
    anchor = int(round(z_anchor_center_px))
    diameter_rows = max(3, int(round(diameter_um / dz_um)))
    search_radius = max(3, int(np.ceil(0.45 * diameter_rows)))
    search0 = max(0, anchor - search_radius)
    search1 = min(n_z, anchor + search_radius + 1)
    baseline0 = max(0, anchor - diameter_rows)
    baseline1 = max(baseline0, search0 - 1)
    if x1 <= x0 or search1 <= search0:
        return TopEdgeResult(
            False, "central_line_out_of_bounds", z_anchor_center_px, float("nan"),
            True, float("nan"), float("nan"), float("nan")
        )
    profile = np.nanmedian(frame[:, x0:x1], axis=1)
    baseline = profile[baseline0:baseline1]
    baseline = baseline[np.isfinite(baseline)]
    if baseline.size < 3:
        fallback_stop = min(search1, search0 + max(3, search_radius // 2))
        baseline = profile[search0:fallback_stop]
        baseline = baseline[np.isfinite(baseline)]
    if baseline.size == 0:
        return TopEdgeResult(
            False, "no_finite_top_background", z_anchor_center_px, float("nan"),
            True, float("nan"), float("nan"), float("nan")
        )
    background = float(np.median(baseline))
    mad = float(np.median(np.abs(baseline - background)))
    sigma = max(1.4826 * mad, np.finfo(float).eps)
    search_profile = profile[search0:search1]
    threshold = background + threshold_sigma * sigma
    all_runs = _runs(np.isfinite(search_profile) & (search_profile > threshold))
    sustained = [(a + search0, b + search0) for a, b in all_runs if b - a >= 2]
    if sustained:
        selected = min(sustained, key=lambda run: abs(run[0] - z_anchor_center_px))
        z_top = float(selected[0])
        fallback = False
        reason = "ok"
    else:
        z_top = float(anchor)
        fallback = True
        reason = "no_sustained_top_edge"
    peak = float(np.nanmax(search_profile) - background) if np.isfinite(search_profile).any() else 0.0
    return TopEdgeResult(
        valid=True,
        invalid_reason=reason,
        z_anchor_center_px=float(z_anchor_center_px),
        z_top_center_px=z_top,
        fallback=fallback,
        local_background=background,
        local_sigma=sigma,
        peak_cnr=peak / sigma,
    )


def localize_geometry(
    omag_frame: NDArray[np.floating],
    *,
    x_anchor_center_px: float,
    z_anchor_center_px: float,
    diameter_um: float,
    dx_um: float,
    dz_um: float,
    local_half_width_px: int | None = None,
) -> LocalizationResult:
    """Create the physical source geometry from a co-registered OMAG frame."""

    lateral = localize_lateral_body(
        omag_frame,
        x_anchor_center_px=x_anchor_center_px,
        z_top_anchor_center_px=z_anchor_center_px,
        diameter_um=diameter_um,
        dx_um=dx_um,
        dz_um=dz_um,
        local_half_width_px=local_half_width_px,
    )
    if not lateral.valid:
        raise ValueError(f"lateral localization failed: {lateral.invalid_reason}")
    top = reestablish_z_top(
        omag_frame,
        x_center_px=lateral.x1_local_geometry_center_px,
        z_anchor_center_px=z_anchor_center_px,
        diameter_um=diameter_um,
        dz_um=dz_um,
    )
    if not top.valid:
        raise ValueError(f"top-edge localization failed: {top.invalid_reason}")
    geometry = VesselGeometry.from_inclusive_centres(
        x_left_center_px=lateral.x_left_center_px,
        x_right_center_px=lateral.x_right_center_px,
        z_top_center_px=top.z_top_center_px,
        diameter_um=diameter_um,
        dx_um=dx_um,
        dz_um=dz_um,
    )
    return LocalizationResult(geometry=geometry, lateral=lateral, top_edge=top)


def localize_geometry_from_surface(
    omag_frame: NDArray[np.floating],
    *,
    surface_z_center_px: float,
    surface_to_vessel_top_um: float,
    effective_refractive_index: float,
    diameter_um: float,
    dx_um: float,
    dz_um: float,
    axial_margin_above_px: int = 6,
    axial_margin_below_px: int = 10,
    lateral_smoothing_sigma_px: float = 2.5,
    edge_exclusion_px: int = 12,
    minimum_peak_cnr: float = 5.0,
    local_half_width_px: int | None = None,
) -> LocalizationResult:
    """Run global surface-guided anchoring followed by the local X1 rules."""

    coarse = find_surface_guided_anchor(
        omag_frame,
        surface_z_center_px=surface_z_center_px,
        surface_to_vessel_top_um=surface_to_vessel_top_um,
        effective_refractive_index=effective_refractive_index,
        diameter_um=diameter_um,
        dx_um=dx_um,
        dz_um=dz_um,
        axial_margin_above_px=axial_margin_above_px,
        axial_margin_below_px=axial_margin_below_px,
        lateral_smoothing_sigma_px=lateral_smoothing_sigma_px,
        edge_exclusion_px=edge_exclusion_px,
        minimum_peak_cnr=minimum_peak_cnr,
    )
    if not coarse.valid:
        raise ValueError(
            f"surface-guided localization failed: {coarse.invalid_reason}"
        )
    localized = localize_geometry(
        omag_frame,
        x_anchor_center_px=coarse.x_anchor_center_px,
        z_anchor_center_px=coarse.predicted_z_top_center_px,
        diameter_um=diameter_um,
        dx_um=dx_um,
        dz_um=dz_um,
        local_half_width_px=local_half_width_px,
    )
    return LocalizationResult(
        geometry=localized.geometry,
        lateral=localized.lateral,
        top_edge=localized.top_edge,
        coarse_anchor=coarse,
    )


def shifted_geometry(
    geometry: VesselGeometry, *, x_shift_px: float = 0.0, z_shift_px: float = 0.0
) -> VesselGeometry:
    """Apply a registered sensitivity perturbation to the frozen geometry."""

    return VesselGeometry(
        x_left_edge_px=geometry.x_left_edge_px + x_shift_px,
        x_right_edge_px=geometry.x_right_edge_px + x_shift_px,
        z_top_edge_px=geometry.z_top_edge_px + z_shift_px,
        diameter_um=geometry.diameter_um,
        dx_um=geometry.dx_um,
        dz_um=geometry.dz_um,
    )
