"""Per-frame visual quality-control figures."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse, Rectangle

from .detection import DetectionResult
from .io import FrameMaps
from .localization import LocalizationResult
from .quantification import QuantificationResult


def _display_log(image: np.ndarray) -> tuple[np.ndarray, float, float]:
    finite = image[np.isfinite(image)]
    if finite.size == 0:
        return np.zeros_like(image), 0.0, 1.0
    floor = max(0.0, float(np.percentile(finite, 0.5)))
    shown = np.log1p(np.maximum(image - floor, 0.0))
    finite_shown = shown[np.isfinite(shown)]
    return shown, float(np.percentile(finite_shown, 1)), float(np.percentile(finite_shown, 99.5))


def _overlay_geometry(ax: plt.Axes, result: QuantificationResult) -> None:
    geometry = result.geometry
    width = geometry.x_right_edge_px - geometry.x_left_edge_px
    height = geometry.diameter_um / geometry.dz_um
    ellipse = Ellipse(
        (geometry.x_center_px, geometry.z_center_px),
        width=width,
        height=height,
        fill=False,
        edgecolor="#00ffff",
        linewidth=1.3,
    )
    ax.add_patch(ellipse)
    ax.axvline(geometry.x_center_px, color="#00ff88", lw=0.9, alpha=0.9)
    ax.axhline(geometry.z_top_edge_px, color="#00ffff", lw=0.8, alpha=0.8)
    ax.axhline(geometry.z_bottom_edge_px, color="#00ffff", lw=0.8, alpha=0.8)
    tail_top = geometry.z_bottom_edge_px + result.tail_gap_um / geometry.dz_um
    tail_height = result.requested_tail_length_um / geometry.dz_um
    ax.add_patch(
        Rectangle(
            (geometry.x_left_edge_px, tail_top),
            width,
            tail_height,
            fill=False,
            edgecolor="#ffcc00",
            linewidth=1.3,
        )
    )
    for columns, color in (
        (result.background.left_columns, "#ff55ff"),
        (result.background.right_columns, "#ff55ff"),
    ):
        if columns.size:
            ax.axvspan(columns[0] - 0.5, columns[-1] + 0.5, color=color, alpha=0.15)


def save_qc_figure(
    path: str | Path,
    *,
    scan_id: str,
    maps: FrameMaps,
    localization: LocalizationResult,
    result: QuantificationResult,
) -> None:
    """Save structure/localization, profile, and tail QC in one figure."""

    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 3, figsize=(17, 9), constrained_layout=True)
    structure_source = maps.stru_amp if maps.stru_amp is not None else maps.omag_raw
    structure_show, structure_min, structure_max = _display_log(structure_source)
    omag_show, omag_min, omag_max = _display_log(maps.omag_raw)
    sv_show, sv_min, sv_max = _display_log(maps.sv_raw)
    physical_aspect = result.geometry.dz_um / result.geometry.dx_um
    axes[0, 0].imshow(
        structure_show,
        cmap="gray",
        aspect=physical_aspect,
        vmin=structure_min,
        vmax=structure_max,
    )
    axes[0, 0].set_title(
        "Structural amplitude (log1p)"
        if maps.stru_amp is not None
        else "Structural amplitude unavailable; OMAG shown"
    )
    axes[0, 1].imshow(
        omag_show,
        cmap="gray",
        aspect=physical_aspect,
        vmin=omag_min,
        vmax=omag_max,
    )
    axes[0, 1].set_title("OMAG localization display (log1p)")
    axes[0, 2].imshow(
        sv_show,
        cmap="magma",
        aspect=physical_aspect,
        vmin=sv_min,
        vmax=sv_max,
    )
    axes[0, 2].set_title("Formal input sv_raw (display log1p)")
    corrected = result.corrected_sv
    finite_corrected = np.abs(corrected[np.isfinite(corrected)])
    limit = float(np.percentile(finite_corrected, 99)) if finite_corrected.size else 1.0
    limit = max(limit, np.finfo(float).eps)
    axes[1, 0].imshow(
        corrected,
        cmap="coolwarm",
        aspect=physical_aspect,
        vmin=-limit,
        vmax=limit,
    )
    axes[1, 0].set_title("Background-corrected C (linear; negatives retained)")
    image_axes = [axes[0, 0], axes[0, 1], axes[0, 2], axes[1, 0]]
    geometry = result.geometry
    x_pad = max(12.0, 1.5 * (geometry.x_right_edge_px - geometry.x_left_edge_px))
    tail_stop = (
        geometry.z_bottom_edge_px
        + result.tail_gap_um / geometry.dz_um
        + result.requested_tail_length_um / geometry.dz_um
    )
    x_limits = (
        max(-0.5, geometry.x_left_edge_px - x_pad),
        min(maps.sv_raw.shape[1] - 0.5, geometry.x_right_edge_px + x_pad),
    )
    z_limits = (
        max(-0.5, geometry.z_top_edge_px - 10),
        min(maps.sv_raw.shape[0] - 0.5, tail_stop + 10),
    )
    for ax in image_axes:
        _overlay_geometry(ax, result)
        ax.set_xlim(*x_limits)
        ax.set_ylim(z_limits[1], z_limits[0])
        ax.set_xlabel("A-line x (0-based pixel centre)")
        ax.set_ylabel("Depth z (0-based pixel centre)")

    distance = (
        np.arange(result.tail_linear_density.size) - result.geometry.z_bottom_edge_px
    ) * result.geometry.dz_um
    profile_axis = axes[1, 1]
    profile_axis.plot(distance, result.vessel_profile, lw=1.0, label="V")
    profile_axis.plot(distance, result.background.left, lw=0.9, label="B_left")
    profile_axis.plot(distance, result.background.right, lw=0.9, label="B_right")
    profile_axis.plot(distance, result.background.combined, lw=1.2, label="B")
    profile_axis.plot(distance, result.tail_contrast_profile, lw=1.2, label="T")
    profile_axis.axvspan(
        -result.geometry.diameter_um,
        0.0,
        color="#00ffff",
        alpha=0.08,
        label="source axial extent",
    )
    profile_axis.axvline(0.0, color="black", lw=0.7)
    profile_axis.set_xlabel("r from vessel bottom (um)")
    profile_axis.set_ylabel("Linear SV")
    profile_axis.set_title("QC02 full-width profiles")
    profile_axis.legend(loc="best", fontsize=7, ncol=2)
    profile_axis.grid(alpha=0.2)

    tail_axis = axes[1, 2]
    tail_axis.plot(distance, result.tail_linear_density, color="#005f73", lw=1.4, label="P(z)")
    tail_axis.axvspan(
        result.tail_gap_um,
        result.tail_gap_um + result.requested_tail_length_um,
        color="#ffcc00",
        alpha=0.14,
        label="formal tail window",
    )
    tail_axis.axhline(0.0, color="black", lw=0.7)
    tail_axis.set_xlabel("r from vessel bottom (um)")
    tail_axis.set_ylabel("P(z), SV x um")
    tail_axis.set_title("QC03 signed tail profile")
    tail_axis.legend(loc="best", fontsize=8)
    tail_axis.grid(alpha=0.2)

    qc_text = (
        f"scan={scan_id} | valid={result.valid} | reason={result.invalid_reason}\n"
        f"X1={localization.lateral.x1_local_geometry_center_px:.2f}, "
        f"x=[{localization.lateral.x_left_center_px},{localization.lateral.x_right_center_px}], "
        f"z_top={localization.top_edge.z_top_center_px:.2f}, "
        f"source_qc={localization.source_qc_valid}\n"
        f"Qv={result.q_vessel:.6g}, Qtail={result.q_tail:.6g}, "
        f"R={result.ratio_tail_to_vessel:.6g}"
    )
    fig.suptitle(qc_text, fontsize=10)
    fig.savefig(output, dpi=160)
    plt.close(fig)


def save_detection_qc(
    path: str | Path, *, scan_id: str, detection: DetectionResult
) -> None:
    """Save the independent binned threshold and censoring diagnostic."""

    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    centres = (
        np.arange(detection.signal_bin_mean.size) + 0.5
    ) * detection.bin_rows * detection.dz_um
    fig, ax = plt.subplots(figsize=(9, 5), constrained_layout=True)
    ax.plot(centres, detection.signal_bin_mean, marker="o", label="T(r) bin mean")
    ax.plot(centres, detection.blank_median, lw=1.0, label="blank median")
    ax.plot(centres, detection.threshold, lw=1.4, label="median + 3 scaled MAD")
    for centre, sustained in zip(centres, detection.sustained_detection, strict=True):
        if sustained:
            half = 0.5 * detection.bin_rows * detection.dz_um
            ax.axvspan(centre - half, centre + half, color="#2a9d8f", alpha=0.16)
    ax.axhline(0.0, color="black", lw=0.7)
    ax.set_xlabel("r bin centre (um)")
    ax.set_ylabel("T(r), linear SV")
    ax.set_title(
        f"{scan_id} | {detection.invalid_reason} | "
        f"D_detect={detection.detectable_length_um:.6g} um"
    )
    ax.grid(alpha=0.2)
    ax.legend(loc="best")
    fig.savefig(output, dpi=160)
    plt.close(fig)
