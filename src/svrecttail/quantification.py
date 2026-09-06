"""Formal SV source and rectangular-tail integrals."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray

from .background import BackgroundProfiles, estimate_background
from .geometry import (
    VesselGeometry,
    ellipse_is_complete,
    ellipse_weights,
    interval_is_complete,
    interval_overlap_weights,
)


FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class QuantificationResult:
    """Formal scalars, profiles, weights, and QC state for one frame."""

    geometry: VesselGeometry
    tail_gap_um: float
    requested_tail_length_um: float
    q_vessel: float
    q_tail: float
    ratio_tail_to_vessel: float
    q_tail_direct_check: float
    source_area_um2: float
    tail_area_um2: float
    requested_tail_area_um2: float
    source_mean: float
    tail_window_complete: bool
    source_window_complete: bool
    source_qc_valid: bool
    background_complete: bool
    valid: bool
    invalid_reason: str
    vessel_profile: FloatArray
    background: BackgroundProfiles
    tail_contrast_profile: FloatArray
    tail_linear_density: FloatArray
    corrected_sv: FloatArray
    vessel_x_weights: FloatArray
    vessel_ellipse_weights: FloatArray
    tail_z_weights: FloatArray

    def summary(self) -> dict[str, Any]:
        """Flatten the scalar result for a frame-level table."""

        data: dict[str, Any] = {
            "q_vessel": self.q_vessel,
            "q_tail": self.q_tail,
            "ratio_tail_to_vessel": self.ratio_tail_to_vessel,
            "q_tail_direct_check": self.q_tail_direct_check,
            "source_area_um2": self.source_area_um2,
            "tail_area_um2": self.tail_area_um2,
            "requested_tail_area_um2": self.requested_tail_area_um2,
            "source_mean": self.source_mean,
            "tail_window_complete": self.tail_window_complete,
            "source_window_complete": self.source_window_complete,
            "source_qc_valid": self.source_qc_valid,
            "background_complete": self.background_complete,
            "background_excluded_side": self.background.excluded_side or "",
            "valid": self.valid,
            "invalid_reason": self.invalid_reason,
            "tail_gap_um": self.tail_gap_um,
            "requested_tail_length_um": self.requested_tail_length_um,
            "lateral_width_um": self.geometry.lateral_width_um,
        }
        data.update(asdict(self.geometry))
        data["z_bottom_edge_px"] = self.geometry.z_bottom_edge_px
        return data

    def profile_records(self, scan_id: str) -> list[dict[str, Any]]:
        """Return one auditable record per axial pixel centre."""

        distance_um = (
            np.arange(self.vessel_profile.size, dtype=np.float64)
            - self.geometry.z_bottom_edge_px
        ) * self.geometry.dz_um
        records: list[dict[str, Any]] = []
        for z_index in range(self.vessel_profile.size):
            in_source = bool(self.vessel_ellipse_weights[z_index].sum() > 0)
            in_tail = bool(self.tail_z_weights[z_index] > 0)
            if in_source and in_tail:
                region = "source_tail_boundary"
            elif in_source:
                region = "source"
            elif in_tail:
                region = "tail"
            else:
                region = "other"
            values = (
                self.vessel_profile[z_index],
                self.background.combined[z_index],
                self.tail_contrast_profile[z_index],
                self.tail_linear_density[z_index],
            )
            records.append(
                {
                    "scan_id": scan_id,
                    "z_index_0based": z_index,
                    "z_um": z_index * self.geometry.dz_um,
                    "r_um": distance_um[z_index],
                    "distance_from_vessel_bottom_um": distance_um[z_index],
                    "V": self.vessel_profile[z_index],
                    "B_left": self.background.left[z_index],
                    "B_right": self.background.right[z_index],
                    "B": self.background.combined[z_index],
                    "B_standard_deviation": self.background.standard_deviation[z_index],
                    "B_median_diagnostic": self.background.median[z_index],
                    "B_scaled_mad_diagnostic": self.background.scaled_mad[z_index],
                    "background_mode": self.background.row_mode[z_index],
                    "T": self.tail_contrast_profile[z_index],
                    "P": self.tail_linear_density[z_index],
                    "tail_z_fraction": self.tail_z_weights[z_index],
                    "background_left_valid_count": self.background.left_valid_count[z_index],
                    "background_right_valid_count": self.background.right_valid_count[z_index],
                    "region": region,
                    "validity": "valid" if np.isfinite(values).all() else "nonfinite_profile",
                }
            )
        return records


def _weighted_row_mean(image: FloatArray, weights: FloatArray) -> FloatArray:
    finite = np.isfinite(image)
    weighted = np.where(finite, image, 0.0) * weights[None, :]
    denominator = (finite * weights[None, :]).sum(axis=1, dtype=np.float64)
    result = np.full(image.shape[0], np.nan, dtype=np.float64)
    np.divide(
        weighted.sum(axis=1, dtype=np.float64),
        denominator,
        out=result,
        where=denominator > 0,
    )
    return result


def _weighted_integral(
    image: FloatArray, weights: FloatArray, pixel_area_um2: float
) -> float:
    support = weights > 0
    if not np.any(support) or not np.isfinite(image[support]).all():
        return float("nan")
    return float(np.sum(image[support] * weights[support]) * pixel_area_um2)


def quantify_frame(
    sv_raw: NDArray[np.floating],
    geometry: VesselGeometry,
    *,
    tail_gap_um: float = 0.0,
    tail_length_um: float = 500.0,
    background_skip_columns: int = 3,
    background_strip_width_columns: int = 5,
    background_excluded_side: str | None = None,
    ellipse_supersample: int = 16,
    source_qc_valid: bool = True,
) -> QuantificationResult:
    """Quantify one linear SV frame using the frozen rectangular protocol.

    ``sv_raw`` must be the unnormalised signal ``var(abs(E), 1, 3)``. The
    row background is subtracted without clipping; negative residuals remain.
    """

    image = np.asarray(sv_raw, dtype=np.float64)
    if image.ndim != 2 or min(image.shape) < 1:
        raise ValueError("sv_raw must be a non-empty 2-D array")
    if tail_gap_um < 0 or tail_length_um <= 0:
        raise ValueError("tail gap must be non-negative and length positive")
    if np.isinf(image).any():
        raise ValueError("sv_raw cannot contain infinity")

    nz, nx = image.shape
    background = estimate_background(
        image,
        geometry,
        skip_columns=background_skip_columns,
        strip_width_columns=background_strip_width_columns,
        excluded_side=background_excluded_side,
    )
    corrected = image - background.combined[:, None]
    x_weights = interval_overlap_weights(
        nx, geometry.x_left_edge_px, geometry.x_right_edge_px
    )
    vessel_profile = _weighted_row_mean(image, x_weights)
    tail_contrast = vessel_profile - background.combined
    finite_corrected = np.isfinite(corrected)
    tail_linear_density = np.sum(
        np.where(finite_corrected, corrected, 0.0) * x_weights[None, :],
        axis=1,
        dtype=np.float64,
    ) * geometry.dx_um
    if np.any(x_weights > 0):
        row_has_full_signal = np.all(np.isfinite(image[:, x_weights > 0]), axis=1)
    else:
        row_has_full_signal = np.zeros(nz, dtype=bool)
    tail_linear_density[~row_has_full_signal | ~np.isfinite(background.combined)] = np.nan
    effective_width_um = float(x_weights.sum() * geometry.dx_um)
    profile_identity_support = np.isfinite(tail_linear_density) & np.isfinite(tail_contrast)
    if profile_identity_support.any() and not np.allclose(
        tail_linear_density[profile_identity_support],
        effective_width_um * tail_contrast[profile_identity_support],
        rtol=2e-12,
        atol=1e-12,
    ):
        raise AssertionError("P(z) and effective_width*T(z) disagree")

    gap_px = tail_gap_um / geometry.dz_um
    length_px = tail_length_um / geometry.dz_um
    tail_top = geometry.z_bottom_edge_px + gap_px
    tail_bottom = tail_top + length_px
    z_weights = interval_overlap_weights(nz, tail_top, tail_bottom)
    source_weights = ellipse_weights(
        image.shape, geometry, supersample=ellipse_supersample
    )
    source_complete = ellipse_is_complete(image.shape, geometry)
    tail_complete = interval_is_complete(nz, tail_top, tail_bottom) and interval_is_complete(
        nx, geometry.x_left_edge_px, geometry.x_right_edge_px
    )
    required_rows = (source_weights.sum(axis=1) > 0) | (z_weights > 0)
    if background.excluded_side == "left":
        strip_complete = background.right_columns.size == background_strip_width_columns
        rows_complete = (background.right_valid_count[required_rows] > 0).all()
    elif background.excluded_side == "right":
        strip_complete = background.left_columns.size == background_strip_width_columns
        rows_complete = (background.left_valid_count[required_rows] > 0).all()
    else:
        strip_complete = bool(
            background.left_columns.size == background_strip_width_columns
            and background.right_columns.size == background_strip_width_columns
        )
        rows_complete = bool(
            (background.left_valid_count[required_rows] > 0).all()
            and (background.right_valid_count[required_rows] > 0).all()
        )
    background_complete = bool(
        strip_complete
        and rows_complete
        and np.isfinite(background.combined[required_rows]).all()
    )

    pixel_area_um2 = geometry.dx_um * geometry.dz_um
    source_area_um2 = float(source_weights.sum() * pixel_area_um2)
    tail_area_weights = np.multiply.outer(z_weights, x_weights)
    tail_area_um2 = float(tail_area_weights.sum() * pixel_area_um2)
    requested_tail_area_um2 = geometry.lateral_width_um * tail_length_um

    q_vessel = _weighted_integral(
        corrected,
        source_weights,
        pixel_area_um2,
    )
    profile_support = z_weights > 0
    if profile_support.any() and np.isfinite(tail_linear_density[profile_support]).all():
        q_tail_profile = float(
            np.sum(tail_linear_density[profile_support] * z_weights[profile_support])
            * geometry.dz_um
        )
    else:
        q_tail_profile = float("nan")
    q_tail_direct = _weighted_integral(
        corrected,
        tail_area_weights,
        pixel_area_um2,
    )
    if np.isfinite(q_tail_profile) and np.isfinite(q_tail_direct):
        if not np.isclose(q_tail_profile, q_tail_direct, rtol=2e-12, atol=1e-12):
            raise AssertionError("profile and direct tail integrals disagree")

    reasons: list[str] = []
    if not source_complete:
        reasons.append("source_window_out_of_frame")
    if not source_qc_valid:
        reasons.append("source_qc_failed")
    if not tail_complete:
        reasons.append("tail_window_out_of_frame")
    if not background_complete:
        reasons.append("background_incomplete")
    if not np.isfinite(q_vessel):
        reasons.append("q_vessel_nonfinite")
    elif q_vessel <= 0:
        reasons.append("q_vessel_nonpositive")
    if not np.isfinite(q_tail_profile):
        reasons.append("q_tail_nonfinite")

    valid = len(reasons) == 0
    q_tail = q_tail_profile if tail_complete and background_complete else float("nan")
    ratio = q_tail / q_vessel if valid else float("nan")
    source_mean = (
        q_vessel / source_area_um2
        if np.isfinite(q_vessel) and source_area_um2 > 0
        else float("nan")
    )
    return QuantificationResult(
        geometry=geometry,
        tail_gap_um=float(tail_gap_um),
        requested_tail_length_um=float(tail_length_um),
        q_vessel=q_vessel,
        q_tail=q_tail,
        ratio_tail_to_vessel=float(ratio),
        q_tail_direct_check=(
            q_tail_direct if tail_complete and background_complete else float("nan")
        ),
        source_area_um2=source_area_um2,
        tail_area_um2=tail_area_um2,
        requested_tail_area_um2=float(requested_tail_area_um2),
        source_mean=float(source_mean),
        tail_window_complete=tail_complete,
        source_window_complete=source_complete,
        source_qc_valid=bool(source_qc_valid),
        background_complete=background_complete,
        valid=valid,
        invalid_reason="ok" if valid else ";".join(reasons),
        vessel_profile=vessel_profile,
        background=background,
        tail_contrast_profile=tail_contrast,
        tail_linear_density=tail_linear_density,
        corrected_sv=corrected,
        vessel_x_weights=x_weights,
        vessel_ellipse_weights=source_weights,
        tail_z_weights=z_weights,
    )
