"""Physical geometry and fractional pixel-area weights.

Coordinates are zero-based pixel-centre coordinates. Pixel ``i`` occupies
``[i - 0.5, i + 0.5]`` and an image of length ``n`` occupies
``[-0.5, n - 0.5]``. Interval bounds use image-edge coordinates.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray


FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class VesselGeometry:
    """Local vessel and tail geometry in pixel-edge coordinates."""

    x_left_edge_px: float
    x_right_edge_px: float
    z_top_edge_px: float
    diameter_um: float
    dx_um: float
    dz_um: float

    def __post_init__(self) -> None:
        values = (
            self.x_left_edge_px,
            self.x_right_edge_px,
            self.z_top_edge_px,
            self.diameter_um,
            self.dx_um,
            self.dz_um,
        )
        if not np.all(np.isfinite(values)):
            raise ValueError("geometry values must be finite")
        if self.x_right_edge_px <= self.x_left_edge_px:
            raise ValueError("x_right_edge_px must exceed x_left_edge_px")
        if self.diameter_um <= 0 or self.dx_um <= 0 or self.dz_um <= 0:
            raise ValueError("diameter and pixel pitches must be positive")

    @property
    def z_bottom_edge_px(self) -> float:
        """Physical lower vessel edge; deliberately not rounded."""

        return self.z_top_edge_px + self.diameter_um / self.dz_um

    @property
    def x_center_px(self) -> float:
        return (self.x_left_edge_px + self.x_right_edge_px) / 2.0

    @property
    def z_center_px(self) -> float:
        return (self.z_top_edge_px + self.z_bottom_edge_px) / 2.0

    @property
    def lateral_width_um(self) -> float:
        return (self.x_right_edge_px - self.x_left_edge_px) * self.dx_um

    @classmethod
    def from_inclusive_centres(
        cls,
        *,
        x_left_center_px: int,
        x_right_center_px: int,
        z_top_center_px: float,
        diameter_um: float,
        dx_um: float,
        dz_um: float,
    ) -> "VesselGeometry":
        """Build geometry from inclusive lateral centre indices."""

        if x_right_center_px < x_left_center_px:
            raise ValueError("right centre index must not precede left")
        return cls(
            x_left_edge_px=float(x_left_center_px) - 0.5,
            x_right_edge_px=float(x_right_center_px) + 0.5,
            z_top_edge_px=float(z_top_center_px) - 0.5,
            diameter_um=float(diameter_um),
            dx_um=float(dx_um),
            dz_um=float(dz_um),
        )


def interval_overlap_weights(
    length: int, lower_edge_px: float, upper_edge_px: float
) -> FloatArray:
    """Return fractional overlap of every pixel with an edge interval."""

    if length < 0:
        raise ValueError("length must be non-negative")
    if not np.isfinite([lower_edge_px, upper_edge_px]).all():
        raise ValueError("interval bounds must be finite")
    if upper_edge_px < lower_edge_px:
        raise ValueError("upper bound must not be below lower bound")
    centres = np.arange(length, dtype=np.float64)
    lower = np.maximum(centres - 0.5, lower_edge_px)
    upper = np.minimum(centres + 0.5, upper_edge_px)
    return np.clip(upper - lower, 0.0, 1.0)


def interval_is_complete(length: int, lower_edge_px: float, upper_edge_px: float) -> bool:
    """Whether an interval lies entirely inside an image axis."""

    tolerance = 1e-12
    return bool(
        lower_edge_px >= -0.5 - tolerance
        and upper_edge_px <= length - 0.5 + tolerance
        and upper_edge_px > lower_edge_px
    )


def rectangle_weights(
    shape: tuple[int, int],
    *,
    x_left_edge_px: float,
    x_right_edge_px: float,
    z_top_edge_px: float,
    z_bottom_edge_px: float,
) -> FloatArray:
    """Fractional area weights for an axis-aligned rectangle."""

    nz, nx = shape
    x_weight = interval_overlap_weights(nx, x_left_edge_px, x_right_edge_px)
    z_weight = interval_overlap_weights(nz, z_top_edge_px, z_bottom_edge_px)
    return np.multiply.outer(z_weight, x_weight)


def ellipse_weights(
    shape: tuple[int, int],
    geometry: VesselGeometry,
    *,
    supersample: int = 16,
) -> FloatArray:
    """Return deterministic sub-pixel area fractions for the source ellipse."""

    if supersample < 1:
        raise ValueError("supersample must be at least 1")
    nz, nx = shape
    weights = np.zeros((nz, nx), dtype=np.float64)
    x_overlap = interval_overlap_weights(
        nx, geometry.x_left_edge_px, geometry.x_right_edge_px
    )
    z_overlap = interval_overlap_weights(
        nz, geometry.z_top_edge_px, geometry.z_bottom_edge_px
    )
    x_indices = np.flatnonzero(x_overlap > 0)
    z_indices = np.flatnonzero(z_overlap > 0)
    if x_indices.size == 0 or z_indices.size == 0:
        return weights

    offsets = (np.arange(supersample, dtype=np.float64) + 0.5) / supersample - 0.5
    x_samples = x_indices[:, None] + offsets[None, :]
    z_samples = z_indices[:, None] + offsets[None, :]
    x_radius_um = geometry.lateral_width_um / 2.0
    z_radius_um = geometry.diameter_um / 2.0
    x_term = ((x_samples - geometry.x_center_px) * geometry.dx_um / x_radius_um) ** 2
    z_term = ((z_samples - geometry.z_center_px) * geometry.dz_um / z_radius_um) ** 2
    inside = z_term[:, None, :, None] + x_term[None, :, None, :] <= 1.0
    local_weights = inside.mean(axis=(2, 3), dtype=np.float64)
    weights[np.ix_(z_indices, x_indices)] = local_weights
    return weights


def ellipse_is_complete(shape: tuple[int, int], geometry: VesselGeometry) -> bool:
    nz, nx = shape
    return interval_is_complete(
        nx, geometry.x_left_edge_px, geometry.x_right_edge_px
    ) and interval_is_complete(nz, geometry.z_top_edge_px, geometry.z_bottom_edge_px)


def vessel_centre_bounds(geometry: VesselGeometry) -> tuple[int, int]:
    """Return first and last centre indices assigned to the detected body."""

    first = int(np.ceil(geometry.x_left_edge_px))
    last = int(np.ceil(geometry.x_right_edge_px)) - 1
    return first, last
