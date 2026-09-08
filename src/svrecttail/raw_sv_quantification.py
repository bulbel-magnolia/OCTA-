"""No-background linear raw-SV integration; independent of historical quantification."""
from __future__ import annotations

import numpy as np

from .geometry import (
    VesselGeometry, ellipse_weights, ellipse_is_complete,
    interval_overlap_weights, interval_is_complete,
)

WINDOWS_UM = (100, 200, 300, 500)
TARGETS_UM = np.arange(0, 501, 10, dtype=float)
METRICS = ("source_area_um2", "source_q_raw", "source_mean_raw") + tuple(
    f"{name}_{length}um" for length in WINDOWS_UM
    for name in ("tail_area_um2", "tail_q_raw", "tail_mean_raw", "ri_tail")
) + ("ri_tail",)


def relative_error(actual: float, expected: float) -> float:
    return abs(actual - expected) / max(abs(expected), np.finfo(float).tiny)


def quantify_raw_sv(image, geometry: VesselGeometry | None, *, geometry_qc_valid=True):
    """Return scalar metrics and nearest native-row anchors without signal transforms.

    Invalid geometry has NA metrics and no anchors. Identity validation belongs to
    the caller and also runs for invalid frames. Incomplete geometry marked valid
    is a hard error, as are nonfinite images or a nonpositive source denominator.
    """
    if not geometry_qc_valid:
        return {name: np.nan for name in METRICS}, []
    g = geometry
    if g is None:
        raise ValueError("Valid entry requires geometry")
    sv = np.asarray(image, dtype=np.float64)
    if sv.ndim != 2 or not np.isfinite(sv).all():
        raise ValueError("Expected finite two-dimensional raw SV")
    if not ellipse_is_complete(sv.shape, g) or not interval_is_complete(
        sv.shape[0], g.z_bottom_edge_px, g.z_bottom_edge_px + 500 / g.dz_um
    ):
        raise ValueError("Valid geometry must contain full source and tail")
    pixel_area = g.dx_um * g.dz_um
    sw = ellipse_weights(sv.shape, g, supersample=16)
    mask = sw > 0
    source_area = float(sw.sum() * pixel_area)
    source_q = float(np.sum(sv[mask] * sw[mask], dtype=np.float64) * pixel_area)
    source_mean = source_q / source_area
    if not np.isfinite(source_mean) or source_mean <= 0:
        raise ValueError("Raw source mean must be finite and positive")
    result = dict(source_area_um2=source_area, source_q_raw=source_q,
                  source_mean_raw=source_mean)
    xw = interval_overlap_weights(sv.shape[1], g.x_left_edge_px, g.x_right_edge_px)
    row_mean = np.sum(sv * xw[None, :], axis=1, dtype=np.float64) / xw.sum()
    area_errors = []
    profile_errors = []
    for length in WINDOWS_UM:
        zw = interval_overlap_weights(sv.shape[0], g.z_bottom_edge_px,
                                      g.z_bottom_edge_px + length / g.dz_um)
        tw = np.multiply.outer(zw, xw)
        mask = tw > 0
        area = float(tw.sum() * pixel_area)
        q = float(np.sum(sv[mask] * tw[mask], dtype=np.float64) * pixel_area)
        mean = q / area
        profile_q = float(np.sum(row_mean * zw) * xw.sum() * pixel_area)
        err = relative_error(q, profile_q)
        ae = relative_error(area, g.lateral_width_um * length)
        if err > 1e-12 or ae > 1e-12:
            raise ValueError("Direct/profile integration or rectangle area mismatch")
        profile_errors.append(err)
        area_errors.append(ae)
        result.update({f"tail_area_um2_{length}um": area, f"tail_q_raw_{length}um": q,
                       f"tail_mean_raw_{length}um": mean, f"ri_tail_{length}um": mean / source_mean})
    result["ri_tail"] = result["ri_tail_500um"]
    result["direct_profile_max_relative_error"] = max(profile_errors)
    result["tail_area_max_relative_error"] = max(area_errors)
    # Native pixels with fractional overlap, including boundary pixels whose
    # centres can lie just outside the interval: identical to frozen D128.
    candidates = np.flatnonzero(zw > 0)
    rr = (np.arange(sv.shape[0]) - g.z_bottom_edge_px) * g.dz_um
    selected = candidates[np.argmin(np.abs(rr[candidates, None] - TARGETS_UM), axis=0)]
    errors = np.abs(rr[selected] - TARGETS_UM)
    if errors.max() > g.dz_um / 2 + 1e-9:
        raise ValueError("Native depth sampling exceeds half-pixel tolerance")
    anchors = [dict(target_r_um=float(t), selected_z_index_0based=int(z),
                    selected_r_um=float(rr[z]), abs_sampling_error_um=float(e),
                    tail_z_fraction=float(zw[z]), S_tail_raw_r=float(row_mean[z]),
                    S_vessel_raw=source_mean, RI_r=float(row_mean[z] / source_mean))
               for t, z, e in zip(TARGETS_UM, selected, errors)]
    if not all(np.isfinite(result[name]) for name in METRICS):
        raise ValueError("Nonfinite formal metric")
    return result, anchors
