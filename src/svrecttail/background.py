"""Fixed bilateral background sampling for rectangular-tail quantification."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from .geometry import VesselGeometry, vessel_centre_bounds


FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class BackgroundProfiles:
    combined: FloatArray
    left: FloatArray
    right: FloatArray
    standard_deviation: FloatArray
    median: FloatArray
    scaled_mad: FloatArray
    left_valid_count: NDArray[np.int64]
    right_valid_count: NDArray[np.int64]
    left_columns: NDArray[np.int64]
    right_columns: NDArray[np.int64]
    row_mode: NDArray[np.str_]
    excluded_side: str | None

    @property
    def complete_rows(self) -> NDArray[np.bool_]:
        return np.isfinite(self.combined)


def background_columns(
    nx: int,
    geometry: VesselGeometry,
    *,
    skip_columns: int = 3,
    strip_width_columns: int = 5,
) -> tuple[NDArray[np.int64], NDArray[np.int64]]:
    """Return clipped left/right columns from the fixed protocol.

    Defaults yield ``iL-8:iL-4`` and ``iR+4:iR+8`` (inclusive), where
    ``iL`` and ``iR`` are the body's first and last centre indices.
    """

    if nx < 1:
        raise ValueError("nx must be positive")
    if skip_columns < 0 or strip_width_columns < 1:
        raise ValueError("skip must be non-negative and strip width positive")
    i_left, i_right = vessel_centre_bounds(geometry)
    left = np.arange(
        i_left - skip_columns - strip_width_columns,
        i_left - skip_columns,
        dtype=np.int64,
    )
    right = np.arange(
        i_right + skip_columns + 1,
        i_right + skip_columns + strip_width_columns + 1,
        dtype=np.int64,
    )
    return left[(left >= 0) & (left < nx)], right[(right >= 0) & (right < nx)]


def _finite_row_mean(values: FloatArray) -> tuple[FloatArray, NDArray[np.int64]]:
    finite = np.isfinite(values)
    count = finite.sum(axis=1, dtype=np.int64)
    total = np.where(finite, values, 0.0).sum(axis=1, dtype=np.float64)
    mean = np.full(values.shape[0], np.nan, dtype=np.float64)
    np.divide(total, count, out=mean, where=count > 0)
    return mean, count


def estimate_background(
    sv_raw: NDArray[np.floating],
    geometry: VesselGeometry,
    *,
    skip_columns: int = 3,
    strip_width_columns: int = 5,
    excluded_side: str | None = None,
) -> BackgroundProfiles:
    """Compute bilateral mean plus non-formal row diagnostics.

    Both sides contribute according to finite pixel count. An explicitly
    excluded side is permitted only through the audited caller input; the
    routine never chooses the darker side. Median, standard deviation and
    scaled MAD remain diagnostic outputs only.
    """

    image = np.asarray(sv_raw, dtype=np.float64)
    if image.ndim != 2:
        raise ValueError("sv_raw must be a 2-D depth-by-A-line array")
    if excluded_side not in (None, "left", "right"):
        raise ValueError("excluded_side must be None, 'left', or 'right'")
    left_columns, right_columns = background_columns(
        image.shape[1],
        geometry,
        skip_columns=skip_columns,
        strip_width_columns=strip_width_columns,
    )
    left_values = image[:, left_columns]
    right_values = image[:, right_columns]
    left_mean, left_count = _finite_row_mean(left_values)
    right_mean, right_count = _finite_row_mean(right_values)

    if excluded_side == "left":
        values = right_values
    elif excluded_side == "right":
        values = left_values
    else:
        values = np.concatenate((left_values, right_values), axis=1)
    combined, _ = _finite_row_mean(values)
    finite = np.isfinite(values)
    no_values = ~finite.any(axis=1)
    if values.shape[1] == 0:
        standard_deviation = np.full(image.shape[0], np.nan, dtype=np.float64)
        median = standard_deviation.copy()
        scaled_mad = standard_deviation.copy()
    else:
        diagnostics = np.where(finite, values, np.nan)
        with np.errstate(invalid="ignore"):
            standard_deviation = np.nanstd(diagnostics, axis=1, ddof=0)
            median = np.nanmedian(diagnostics, axis=1)
            scaled_mad = 1.4826 * np.nanmedian(
                np.abs(diagnostics - median[:, None]), axis=1
            )
        standard_deviation[no_values] = np.nan
        median[no_values] = np.nan
        scaled_mad[no_values] = np.nan
    row_mode = np.full(image.shape[0], "unavailable", dtype="<U24")
    has_left = left_count > 0
    has_right = right_count > 0
    if excluded_side == "left":
        row_mode[has_right] = "right_only_explicit"
    elif excluded_side == "right":
        row_mode[has_left] = "left_only_explicit"
    else:
        row_mode[has_left & ~has_right] = "left_only_unapproved"
        row_mode[~has_left & has_right] = "right_only_unapproved"
        row_mode[has_left & has_right] = "bilateral"
    return BackgroundProfiles(
        combined=combined,
        left=left_mean,
        right=right_mean,
        standard_deviation=standard_deviation,
        median=median,
        scaled_mad=scaled_mad,
        left_valid_count=left_count,
        right_valid_count=right_count,
        left_columns=left_columns,
        right_columns=right_columns,
        row_mode=row_mode,
        excluded_side=excluded_side,
    )
