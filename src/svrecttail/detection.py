"""Exploratory tail-extent detection against matched blank profiles."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray


FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class DetectionResult:
    detected: bool
    detectable_length_um: float
    right_censored: bool
    invalid_reason: str
    signal_bin_mean: FloatArray
    blank_median: FloatArray
    blank_scaled_mad: FloatArray
    threshold: FloatArray
    blank_count: NDArray[np.int64]
    exceeds_threshold: NDArray[np.bool_]
    sustained_detection: NDArray[np.bool_]
    bin_rows: int
    dz_um: float

    def bin_records(self, scan_id: str) -> list[dict[str, object]]:
        records: list[dict[str, object]] = []
        for index in range(self.signal_bin_mean.size):
            records.append(
                {
                    "scan_id": scan_id,
                    "bin_index_0based": index,
                    "bin_start_um": index * self.bin_rows * self.dz_um,
                    "bin_stop_um": (index + 1) * self.bin_rows * self.dz_um,
                    "signal_bin_mean": self.signal_bin_mean[index],
                    "blank_median": self.blank_median[index],
                    "blank_scaled_mad": self.blank_scaled_mad[index],
                    "threshold": self.threshold[index],
                    "blank_count": self.blank_count[index],
                    "exceeds_threshold": self.exceeds_threshold[index],
                    "sustained_detection": self.sustained_detection[index],
                }
            )
        return records


def _bin_profile(profile: FloatArray, bin_rows: int, n_bins: int) -> FloatArray:
    trimmed = profile[: n_bins * bin_rows].reshape(n_bins, bin_rows)
    valid = np.isfinite(trimmed).all(axis=1)
    result = np.full(n_bins, np.nan, dtype=np.float64)
    result[valid] = trimmed[valid].mean(axis=1)
    return result


def _sustained_mask(mask: NDArray[np.bool_], minimum_length: int) -> NDArray[np.bool_]:
    output = np.zeros(mask.shape, dtype=bool)
    start = 0
    while start < mask.size:
        if not mask[start]:
            start += 1
            continue
        stop = start + 1
        while stop < mask.size and mask[stop]:
            stop += 1
        if stop - start >= minimum_length:
            output[start:stop] = True
        start = stop
    return output


def detect_tail_extent(
    signal_profile: NDArray[np.floating],
    blank_profiles: NDArray[np.floating],
    *,
    dz_um: float,
    bin_rows: int = 5,
    threshold_mad_multiplier: float = 3.0,
    minimum_consecutive_bins: int = 2,
    minimum_blank_samples_per_bin: int = 3,
) -> DetectionResult:
    """Detect the furthest sustained 5-row signal using blank nuisance data.

    ``signal_profile`` and every row of ``blank_profiles`` must start at their
    matched vessel-bottom/virtual-strip origin. Only complete non-overlapping
    bins are evaluated. Reaching the final evaluable bin is right-censored.
    """

    signal = np.asarray(signal_profile, dtype=np.float64).reshape(-1)
    blanks = np.asarray(blank_profiles, dtype=np.float64)
    if blanks.ndim == 1:
        blanks = blanks[None, :]
    if blanks.ndim != 2:
        raise ValueError("blank_profiles must be a 2-D sample-by-row array")
    if dz_um <= 0 or bin_rows < 1 or minimum_consecutive_bins < 1:
        raise ValueError("dz_um and bin counts must be positive")
    if threshold_mad_multiplier <= 0 or minimum_blank_samples_per_bin < 1:
        raise ValueError("threshold multiplier and minimum blank count must be positive")
    available_rows = min(signal.size, blanks.shape[1])
    n_bins = available_rows // bin_rows
    if n_bins < minimum_consecutive_bins:
        empty = np.empty(0, dtype=np.float64)
        return DetectionResult(
            detected=False,
            detectable_length_um=float("nan"),
            right_censored=False,
            invalid_reason="insufficient_complete_bins",
            signal_bin_mean=empty,
            blank_median=empty.copy(),
            blank_scaled_mad=empty.copy(),
            threshold=empty.copy(),
            blank_count=np.empty(0, dtype=np.int64),
            exceeds_threshold=np.empty(0, dtype=bool),
            sustained_detection=np.empty(0, dtype=bool),
            bin_rows=bin_rows,
            dz_um=float(dz_um),
        )

    signal_bins = _bin_profile(signal, bin_rows, n_bins)
    blank_bins = np.vstack([_bin_profile(row, bin_rows, n_bins) for row in blanks])
    blank_count = np.isfinite(blank_bins).sum(axis=0, dtype=np.int64)
    blank_median = np.full(n_bins, np.nan, dtype=np.float64)
    blank_mad = np.full(n_bins, np.nan, dtype=np.float64)
    for index in range(n_bins):
        values = blank_bins[:, index]
        values = values[np.isfinite(values)]
        if values.size < minimum_blank_samples_per_bin:
            continue
        blank_median[index] = np.median(values)
        blank_mad[index] = 1.4826 * np.median(np.abs(values - blank_median[index]))
    threshold = blank_median + threshold_mad_multiplier * blank_mad
    insufficient_blank = blank_count < minimum_blank_samples_per_bin
    degenerate_scale = (~insufficient_blank) & (
        ~np.isfinite(blank_mad) | (blank_mad <= np.finfo(float).eps)
    )
    nonfinite_signal = ~np.isfinite(signal_bins)
    if insufficient_blank.any():
        forced_reason = "insufficient_blank_samples"
    elif degenerate_scale.any():
        forced_reason = "degenerate_blank_scale"
    elif nonfinite_signal.any():
        forced_reason = "nonfinite_signal_bins"
    else:
        forced_reason = None
    evaluable = (
        np.isfinite(signal_bins)
        & np.isfinite(threshold)
        & ~insufficient_blank
        & ~degenerate_scale
    )
    exceeds = evaluable & (signal_bins > threshold)
    sustained = _sustained_mask(exceeds, minimum_consecutive_bins)
    detected = bool(sustained.any()) and forced_reason is None
    if forced_reason is not None:
        length_um = float("nan")
        right_censored = False
        reason = forced_reason
        sustained[:] = False
    elif detected:
        final_index = int(np.flatnonzero(sustained)[-1])
        length_um = (final_index + 1) * bin_rows * dz_um
        right_censored = final_index == n_bins - 1
        reason = "right_censored" if right_censored else "ok"
    else:
        length_um = float("nan")
        right_censored = False
        reason = "not_detected"
    return DetectionResult(
        detected=detected,
        detectable_length_um=float(length_um),
        right_censored=right_censored,
        invalid_reason=reason,
        signal_bin_mean=signal_bins,
        blank_median=blank_median,
        blank_scaled_mad=blank_mad,
        threshold=threshold,
        blank_count=blank_count,
        exceeds_threshold=exceeds,
        sustained_detection=sustained,
        bin_rows=bin_rows,
        dz_um=float(dz_um),
    )
