import numpy as np

from svrecttail.signal import compute_sv_maps


def test_sv_uses_population_variance_denominator_n() -> None:
    amplitudes = np.array([1.0, 2.0, 3.0, 4.0])
    field = amplitudes.reshape(1, 1, 4).astype(np.complex128)
    sv_raw, stru_amp = compute_sv_maps(field)
    assert sv_raw[0, 0] == 1.25
    assert stru_amp[0, 0] == 2.5
    assert sv_raw[0, 0] != np.var(amplitudes, ddof=1)
