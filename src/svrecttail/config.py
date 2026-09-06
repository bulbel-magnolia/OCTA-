"""Validated JSON configuration for frozen analysis runs."""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping


@dataclass(frozen=True)
class CalibrationConfig:
    diameter_um: float
    dx_um: float
    dz_um: float


@dataclass(frozen=True)
class GeometryConfig:
    tail_gap_um: float
    tail_length_um: float
    ellipse_supersample: int


@dataclass(frozen=True)
class BackgroundConfig:
    skip_columns: int
    strip_width_columns: int


@dataclass(frozen=True)
class DetectionConfig:
    enabled: bool
    bin_rows: int
    threshold_mad_multiplier: float
    minimum_consecutive_bins: int
    minimum_blank_samples_per_bin: int


@dataclass(frozen=True)
class RunConfig:
    schema_version: str
    parameter_status: str
    calibration: CalibrationConfig
    geometry: GeometryConfig
    background: BackgroundConfig
    detection: DetectionConfig
    raw: Mapping[str, Any]


def _positive(name: str, value: Any) -> float:
    result = float(value)
    if not math.isfinite(result) or result <= 0:
        raise ValueError(f"{name} must be finite and positive")
    return result


def _nonnegative(name: str, value: Any) -> float:
    result = float(value)
    if not math.isfinite(result) or result < 0:
        raise ValueError(f"{name} must be finite and non-negative")
    return result


def _integer(name: str, value: Any, *, minimum: int) -> int:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be an integer")
    numeric = float(value)
    if not math.isfinite(numeric) or numeric != round(numeric):
        raise ValueError(f"{name} must be an integer")
    result = int(numeric)
    if result < minimum:
        raise ValueError(f"{name} must be at least {minimum}")
    return result


def parse_config(data: Mapping[str, Any]) -> RunConfig:
    """Validate the protocol-critical fields of a config mapping."""

    signal = data.get("signal", {})
    if signal.get("formal_input") != "sv_raw":
        raise ValueError("signal.formal_input must be 'sv_raw'")
    if signal.get("definition") != "var(abs(E), 1, 3)":
        raise ValueError("formal SV definition must be var(abs(E), 1, 3)")
    if signal.get("variance_denominator") != "N":
        raise ValueError("formal variance denominator must be N")
    if signal.get("negative_background_residuals") != "retain":
        raise ValueError("negative background residuals must be retained")

    calibration = data["calibration"]
    geometry = data["geometry"]
    background = data["background"]
    detection = data["detection"]
    coordinate_base = _integer(
        "coordinate_base", geometry.get("coordinate_base", -1), minimum=0
    )
    if coordinate_base != 0:
        raise ValueError("coordinate_base must be zero")
    formal_estimator = background.get("formal_estimator")
    if formal_estimator != "mean":
        raise ValueError("formal background estimator must be mean")

    schema_version = str(data["schema_version"]).strip()
    parameter_status = str(data["parameter_status"]).strip()
    if not schema_version:
        raise ValueError("schema_version must be non-empty")
    if parameter_status not in ("pilot", "frozen"):
        raise ValueError("parameter_status must be pilot or frozen")
    if not isinstance(detection["enabled"], bool):
        raise ValueError("detection.enabled must be a JSON boolean")

    result = RunConfig(
        schema_version=schema_version,
        parameter_status=parameter_status,
        calibration=CalibrationConfig(
            diameter_um=_positive("diameter_um", calibration["diameter_um"]),
            dx_um=_positive("dx_um", calibration["dx_um"]),
            dz_um=_positive("dz_um", calibration["dz_um"]),
        ),
        geometry=GeometryConfig(
            tail_gap_um=_nonnegative("tail_gap_um", geometry["tail_gap_um"]),
            tail_length_um=_positive("tail_length_um", geometry["tail_length_um"]),
            ellipse_supersample=_integer(
                "ellipse_supersample", geometry["ellipse_supersample"], minimum=1
            ),
        ),
        background=BackgroundConfig(
            skip_columns=_integer(
                "skip_columns", background["skip_columns"], minimum=0
            ),
            strip_width_columns=_integer(
                "strip_width_columns",
                background["strip_width_columns"],
                minimum=1,
            ),
        ),
        detection=DetectionConfig(
            enabled=detection["enabled"],
            bin_rows=_integer("bin_rows", detection["bin_rows"], minimum=1),
            threshold_mad_multiplier=_positive(
                "threshold_mad_multiplier", detection["threshold_mad_multiplier"]
            ),
            minimum_consecutive_bins=_integer(
                "minimum_consecutive_bins",
                detection["minimum_consecutive_bins"],
                minimum=1,
            ),
            minimum_blank_samples_per_bin=_integer(
                "minimum_blank_samples_per_bin",
                detection["minimum_blank_samples_per_bin"],
                minimum=1,
            ),
        ),
        raw=data,
    )
    return result


def load_config(path: str | Path) -> RunConfig:
    config_path = Path(path)
    with config_path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    return parse_config(data)
