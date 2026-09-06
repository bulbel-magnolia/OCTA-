"""Plot representative frames around an axial tracking excursion."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC = REPO_ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from svrecttail.mentor import tracking_core  # noqa: E402
from svrecttail.mentor_tracking import load_flow_dicom  # noqa: E402


OLD_COLOR = "#e6b800"
NEW_COLOR = "#d62728"
SEED_COLOR = "#00bcd4"
PEAK_COLOR = "#2667ff"
SUPPORT_COLOR = "#7b2cbf"


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--flow-dicom", required=True, type=Path)
    parser.add_argument("--new-tracking", required=True, type=Path)
    parser.add_argument("--old-tracking", required=True, type=Path)
    parser.add_argument("--tracking-config", required=True, type=Path)
    parser.add_argument("--frames", nargs=3, required=True, type=int)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--summary-output", required=True, type=Path)
    return parser.parse_args()


def _display_omag(image: np.ndarray) -> tuple[np.ndarray, float, float]:
    """Keep zero-valued DICOM background black while showing non-zero detail."""

    shown = image.astype(float, copy=False)
    positive = shown[np.isfinite(shown) & (shown > 0)]
    if positive.size == 0:
        return shown, 0.0, 1.0
    return (
        shown,
        float(np.percentile(positive, 1.0)),
        float(np.percentile(positive, 99.0)),
    )


def _row(table: pd.DataFrame, frame: int) -> pd.Series:
    selected = table.loc[table["frame_index"].astype(int).eq(frame)]
    if len(selected) != 1:
        raise ValueError(f"expected one tracking row for frame {frame}")
    return selected.iloc[0]


def _candidate_label(row: pd.Series) -> str:
    if bool(row["z_candidate_accepted"]):
        return "local candidate accepted"
    if bool(row["z_candidate_present"]):
        return "local candidate rejected"
    return "no local candidate; trajectory-assisted"


def main() -> int:
    args = arguments()
    with args.tracking_config.open("r", encoding="utf-8") as handle:
        raw_config = json.load(handle)
    config = tracking_core.merge_tracking_config(raw_config)
    volume = load_flow_dicom(args.flow_dicom)
    new = pd.read_csv(args.new_tracking)
    old = pd.read_csv(args.old_tracking)
    frames = [int(value) for value in args.frames]
    if len(set(frames)) != 3:
        raise ValueError("--frames must contain three distinct frame indices")

    dx_um = float(config["lateral_um_per_px"])
    dz_um = float(config["axial_um_per_px"])
    diameter_um = float(new["diameter_um"].iloc[0])
    diameter_z_px = diameter_um / dz_um
    z_view = (165.0, 285.0)

    fig = plt.figure(figsize=(15.5, 12.5), constrained_layout=True)
    grid = fig.add_gridspec(3, 3, height_ratios=(2.15, 1.45, 1.10))
    summary_rows: list[dict[str, object]] = []
    frame_colors = ("#009E73", "#D55E00", "#CC79A7")

    for column, (frame, frame_color) in enumerate(zip(frames, frame_colors)):
        new_row = _row(new, frame)
        old_row = _row(old, frame)
        x_center = float(new_row["x4_centroid_isolated_jump_corrected_px"])
        width = float(new_row["local_body_run_width_px"])
        x_left = x_center - width / 2.0
        x_right = x_center + width / 2.0
        old_top = float(old_row["z_upper_px"])
        new_top = float(new_row["z_upper_px"])
        seed = float(new_row["seed_z_upper_px"])
        peak = float(new_row["z_peak_px"])

        image_axis = fig.add_subplot(grid[0, column])
        shown, vmin, vmax = _display_omag(volume[frame])
        image_axis.imshow(
            shown,
            cmap="gray",
            aspect=dz_um / dx_um,
            vmin=vmin,
            vmax=vmax,
        )
        image_axis.axhline(
            old_top, color=OLD_COLOR, linewidth=1.7, linestyle="--"
        )
        image_axis.axhline(new_top, color=NEW_COLOR, linewidth=1.8)
        if np.isfinite(seed):
            image_axis.axhline(
                seed, color=SEED_COLOR, linewidth=1.3, linestyle=":"
            )
        if np.isfinite(peak):
            image_axis.scatter(
                [x_center],
                [peak],
                s=42,
                marker="v",
                color=PEAK_COLOR,
                edgecolor="white",
                linewidth=0.7,
                zorder=5,
            )
        image_axis.axvline(x_center, color="#00A676", linewidth=1.1)
        image_axis.axvline(x_left, color="#00A6D6", linewidth=0.9)
        image_axis.axvline(x_right, color="#00A6D6", linewidth=0.9)
        image_axis.add_patch(
            Ellipse(
                (x_center, new_top + diameter_z_px / 2.0),
                width=width,
                height=diameter_z_px,
                fill=False,
                color="#00A6D6",
                linewidth=1.5,
            )
        )
        image_axis.set_xlim(x_center - 34, x_center + 34)
        image_axis.set_ylim(z_view[1], z_view[0])
        image_axis.set_xlabel("A-line x (pixel centre)")
        if column == 0:
            image_axis.set_ylabel("Depth z (pixel centre)")
        image_axis.set_title(
            f"B-scan {frame}: {_candidate_label(new_row)}"
            + chr(10)
            + f"v1={old_top:.0f}, v2 model={new_top:.0f}, "
            + (
                f"local seed={seed:.0f}, peak={peak:.0f}"
                if np.isfinite(seed)
                else f"local seed=missing, peak={peak:.0f}"
            ),
            fontsize=10,
        )

        one_frame = volume[frame : frame + 1]
        profiles = tracking_core.extract_profiles(
            one_frame,
            np.asarray([x_center]),
            diameter_um,
            config,
        )
        row_support = tracking_core.extract_persistent_row_support(
            one_frame,
            np.asarray([x_center]),
            diameter_um,
            config,
        )[0]
        excess = profiles.excess[0].astype(float)
        local_slice = slice(int(z_view[0]), int(z_view[1]) + 1)
        scale = float(np.nanpercentile(np.maximum(excess[local_slice], 0.0), 99.0))
        if not np.isfinite(scale) or scale <= np.finfo(float).eps:
            scale = 1.0
        normalized = np.clip(excess / scale, -0.25, 1.25)
        depth = np.arange(excess.size)

        profile_axis = fig.add_subplot(grid[1, column])
        profile_axis.plot(
            normalized,
            depth,
            color="#333333",
            linewidth=1.1,
            label="central excess / p99",
        )
        profile_axis.plot(
            row_support,
            depth,
            color=SUPPORT_COLOR,
            linewidth=1.1,
            label="lateral support fraction",
        )
        profile_axis.axhline(
            old_top, color=OLD_COLOR, linewidth=1.4, linestyle="--"
        )
        profile_axis.axhline(new_top, color=NEW_COLOR, linewidth=1.5)
        if np.isfinite(seed):
            profile_axis.axhline(
                seed, color=SEED_COLOR, linewidth=1.2, linestyle=":"
            )
        profile_axis.scatter(
            [0.0],
            [peak],
            s=34,
            marker="v",
            color=PEAK_COLOR,
            zorder=4,
        )
        profile_axis.set_xlim(-0.25, 1.25)
        profile_axis.set_ylim(z_view[1], z_view[0])
        profile_axis.grid(alpha=0.18)
        profile_axis.set_xlabel("Normalized evidence / supported x fraction")
        if column == 0:
            profile_axis.set_ylabel("Depth z (pixel centre)")
        profile_axis.set_title(
            f"candidate accepted={bool(new_row['z_candidate_accepted'])}; "
            + f"assessment={new_row['vessel_presence_prediction']}",
            fontsize=9,
        )
        if column == 0:
            profile_axis.legend(loc="lower left", fontsize=8)

        summary_rows.append(
            {
                "frame_index": frame,
                "v1_z_upper_px": old_top,
                "v2_model_z_upper_px": new_top,
                "v2_local_seed_z_upper_px": seed,
                "v2_peak_z_px": peak,
                "z_candidate_present": bool(new_row["z_candidate_present"]),
                "z_candidate_accepted": bool(new_row["z_candidate_accepted"]),
                "new_tracking_class": new_row["new_tracking_class"],
                "assessability_score": float(new_row["assessability_score"]),
                "vessel_presence_prediction": new_row[
                    "vessel_presence_prediction"
                ],
                "upper_top_contrast_snr": new_row[
                    "upper_top_contrast_snr"
                ],
                "upper_bottom_contrast_snr": new_row[
                    "upper_bottom_contrast_snr"
                ],
                "upper_support_difference": new_row[
                    "upper_support_difference"
                ],
                "upper_balance_fraction": new_row[
                    "upper_balance_fraction"
                ],
            }
        )

    trajectory_axis = fig.add_subplot(grid[2, :])
    window = new["frame_index"].between(400, 485)
    new_window = new.loc[window]
    old_window = old.loc[old["frame_index"].between(400, 485)]
    trajectory_axis.plot(
        old_window["frame_index"],
        old_window["z_upper_px"],
        color=OLD_COLOR,
        linewidth=1.6,
        linestyle="--",
        label="v1 z_upper",
    )
    trajectory_axis.plot(
        new_window["frame_index"],
        new_window["z_upper_px"],
        color=NEW_COLOR,
        linewidth=1.8,
        label="v2 completed trajectory",
    )
    present = new_window["z_candidate_present"].astype(bool)
    accepted = new_window["z_candidate_accepted"].astype(bool)
    trajectory_axis.scatter(
        new_window.loc[accepted, "frame_index"],
        new_window.loc[accepted, "seed_z_upper_px"],
        s=18,
        color=SEED_COLOR,
        label="accepted local seed",
        zorder=4,
    )
    rejected = present & ~accepted
    trajectory_axis.scatter(
        new_window.loc[rejected, "frame_index"],
        new_window.loc[rejected, "seed_z_upper_px"],
        s=23,
        facecolor="none",
        edgecolor="#e66101",
        linewidth=0.9,
        label="rejected local seed",
        zorder=4,
    )
    trajectory_axis.axvspan(425, 465, color="#d62728", alpha=0.07)
    for frame, frame_color in zip(frames, frame_colors):
        trajectory_axis.axvline(frame, color=frame_color, linewidth=1.2)
        row = _row(new, frame)
        trajectory_axis.scatter(
            [frame],
            [float(row["z_upper_px"])],
            s=58,
            facecolor=frame_color,
            edgecolor="black",
            linewidth=0.6,
            zorder=5,
        )
    trajectory_axis.set_xlim(400, 485)
    trajectory_axis.invert_yaxis()
    trajectory_axis.set_xlabel("B-scan index")
    trajectory_axis.set_ylabel("z_upper (pixel centre; deeper is downward)")
    trajectory_axis.set_title(
        "Context: sparse/false deep candidates pull the completed v2 trajectory; "
        "the v1 path stays near the visible vessel"
    )
    trajectory_axis.grid(alpha=0.2)
    trajectory_axis.legend(loc="upper left", ncol=4, fontsize=8)

    handles = [
        Line2D([0], [0], color=OLD_COLOR, linestyle="--", label="v1 top"),
        Line2D([0], [0], color=NEW_COLOR, label="v2 completed top"),
        Line2D([0], [0], color=SEED_COLOR, linestyle=":", label="v2 local seed"),
        Line2D(
            [0],
            [0],
            color=PEAK_COLOR,
            marker="v",
            linestyle="none",
            label="v2 selected peak",
        ),
    ]
    fig.legend(handles=handles, loc="outside lower center", ncol=4, fontsize=9)
    fig.suptitle(
        "flow03 axial excursion diagnostic on raw OMAG | physical display aspect "
        f"dz/dx={dz_um:g}/{dx_um:g}",
        fontsize=15,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame.from_records(summary_rows).to_csv(args.summary_output, index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
