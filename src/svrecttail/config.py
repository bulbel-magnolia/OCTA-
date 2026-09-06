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
class LocalizationConfig:
    mode: str
    fixed_surface_z_center_px: float | None
    surface_to_vessel_top_um: float | None
    mentor_primary_alpha: float | None
    mentor_final_x_method: str | None
    mentor_lateral_width_method: str | None
    mentor_top_method: str | None
    mentor_required_assessability: str | None
    effective_refractive_index: float
    axial_margin_above_px: int
    axial_margin_below_px: int
    lateral_smoothing_sigma_px: float
    edge_exclusion_px: int
    minimum_peak_cnr: float


@dataclass(frozen=True)
class RunConfig:
    schema_version: str
    parameter_status: str
    calibration: CalibrationConfig
    geometry: GeometryConfig
    background: BackgroundConfig
    detection: DetectionConfig
    localization: LocalizationConfig
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
    localization = data.get("localization", {"mode": "manifest_anchor"})
    localization_mode = str(localization.get("mode", "manifest_anchor")).strip()
    if localization_mode not in (
        "manifest_anchor",
        "fixed_surface_global_x",
        "mentor_tracking",
    ):
        raise ValueError(
            "localization.mode must be manifest_anchor, "
            "fixed_surface_global_x, or mentor_tracking"
        )
    if localization_mode == "fixed_surface_global_x":
        fixed_surface_z_center_px = _nonnegative(
            "fixed_surface_z_center_px",
            localization["fixed_surface_z_center_px"],
        )
        surface_to_vessel_top_um = _positive(
            "surface_to_vessel_top_um",
            localization["surface_to_vessel_top_um"],
        )
    else:
        fixed_surface_z_center_px = None
        surface_to_vessel_top_um = None
    if localization_mode == "mentor_tracking":
        mentor_primary_alpha = _positive(
            "mentor_primary_alpha",
            localization.get("mentor_primary_alpha", 0.15),
        )
        mentor_final_x_method = str(
            localization.get(
                "mentor_final_x_method",
                "X4_centroid_isolated_jump_corrected",
            )
        ).strip()
        mentor_lateral_width_method = str(
            localization.get(
                "mentor_lateral_width_method",
                "X1_continuous_run_recentered_on_X4",
            )
        ).strip()
        mentor_top_method = str(
            localization.get("mentor_top_method", "primary_alpha_z_upper")
        ).strip()
        mentor_required_assessability = str(
            localization.get("mentor_required_assessability", "assessable")
        ).strip()
        expected_methods = {
            "mentor_final_x_method": (
                mentor_final_x_method,
                "X4_centroid_isolated_jump_corrected",
            ),
            "mentor_lateral_width_method": (
                mentor_lateral_width_method,
                "X1_continuous_run_recentered_on_X4",
            ),
            "mentor_top_method": (
                mentor_top_method,
                "primary_alpha_z_upper",
            ),
            "mentor_required_assessability": (
                mentor_required_assessability,
                "assessable",
            ),
        }
        mismatched_methods = [
            name
            for name, (actual, expected) in expected_methods.items()
            if actual != expected
        ]
        if mismatched_methods:
            raise ValueError(
                "unsupported mentor localization method: "
                + ", ".join(mismatched_methods)
            )
    else:
        mentor_primary_alpha = None
        mentor_final_x_method = None
        mentor_lateral_width_method = None
        mentor_top_method = None
        mentor_required_assessability = None
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
        localization=LocalizationConfig(
            mode=localization_mode,
            fixed_surface_z_center_px=fixed_surface_z_center_px,
            surface_to_vessel_top_um=surface_to_vessel_top_um,
            mentor_primary_alpha=mentor_primary_alpha,
            mentor_final_x_method=mentor_final_x_method,
            mentor_lateral_width_method=mentor_lateral_width_method,
            mentor_top_method=mentor_top_method,
            mentor_required_assessability=mentor_required_assessability,
            effective_refractive_index=_positive(
                "effective_refractive_index",
                localization.get("effective_refractive_index", 1.0),
            ),
            axial_margin_above_px=_integer(
                "axial_margin_above_px",
                localization.get("axial_margin_above_px", 6),
                minimum=0,
            ),
            axial_margin_below_px=_integer(
                "axial_margin_below_px",
                localization.get("axial_margin_below_px", 10),
                minimum=0,
            ),
            lateral_smoothing_sigma_px=_positive(
                "lateral_smoothing_sigma_px",
                localization.get("lateral_smoothing_sigma_px", 2.5),
            ),
            edge_exclusion_px=_integer(
                "edge_exclusion_px",
                localization.get("edge_exclusion_px", 12),
                minimum=0,
            ),
            minimum_peak_cnr=_nonnegative(
                "minimum_peak_cnr",
                localization.get("minimum_peak_cnr", 5.0),
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
