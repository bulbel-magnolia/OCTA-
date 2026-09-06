"""Formal signal construction from reconstructed complex repeats."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


def compute_sv_maps(
    reconstructed_complex: NDArray[np.complexfloating],
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Return ``sv_raw`` and ``stru_amp`` for a depth-x-A-line-x-repeat array.

    The variance uses denominator ``N`` (NumPy ``ddof=0``), matching MATLAB
    ``var(abs(E), 1, 3)`` exactly.
    """

    field = np.asarray(reconstructed_complex)
    if field.ndim != 3 or field.shape[2] < 1:
        raise ValueError("reconstructed_complex must be depth-by-A-line-by-repeat")
    amplitude = np.abs(field).astype(np.float64, copy=False)
    sv_raw = np.var(amplitude, axis=2, ddof=0, dtype=np.float64)
    stru_amp = np.mean(amplitude, axis=2, dtype=np.float64)
    return sv_raw, stru_amp
