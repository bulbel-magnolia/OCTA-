"""Plot the 15-frame hybrid-localization pilot and all 500-frame tracks."""

from __future__ import annotations

import argparse
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

from svrecttail.io import load_frame_maps, resolve_source_path


def _arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--localization", required=True, type=Path)
    parser.add_argument("--old-localization", type=Path)
    parser.add_argument("--overview-output", required=True, type=Path)
    parser.add_argument("--trajectory-output", required=True, type=Path)
    parser.add_argument("--comparison-output", type=Path)
    parser.add_argument("--comparison-figure-output", type=Path)
    parser.add_argument(
        "--reference-label",
        default="reference",
        help="Short label for --old-localization in comparison outputs",
    )
    return parser.parse_args()


def _display_log(image: np.ndarray) -> tuple[np.ndarray, float, float]:
    finite = image[np.isfinite(image)]
    floor = max(0.0, float(np.percentile(finite, 0.5)))
    shown = np.log1p(np.maximum(image - floor, 0.0))
    finite_shown = shown[np.isfinite(shown)]
    return (
        shown,
        float(np.percentile(finite_shown, 1.0)),
        float(np.percentile(finite_shown, 99.7)),
    )


def _localization_row(
    localization: pd.DataFrame, scan_id: str, frame_index: int
) -> pd.Series:
    frame_values = pd.to_numeric(
        localization["mentor_frame_index"], errors="coerce"
    )
    selected = localization.loc[
        localization["scan_id"].astype(str).eq(scan_id)
        & frame_values.eq(frame_index)
    ]
    if len(selected) != 1:
        raise ValueError(
            f"expected one localization row for {scan_id} frame {frame_index}"
        )
    return selected.iloc[0]


def plot_overview(
    manifest_path: Path,
    manifest: pd.DataFrame,
    localization: pd.DataFrame,
    output_path: Path,
) -> None:
    speeds = sorted(
        manifest["flow_speed_mm_s"].astype(float).unique().tolist()
    )
    positions = ["front", "middle", "rear"]
    fig, axes = plt.subplots(
        len(speeds),
        len(positions),
        figsize=(15, 15),
        constrained_layout=True,
    )
    for row_index, speed in enumerate(speeds):
        for column_index, position in enumerate(positions):
            ax = axes[row_index, column_index]
            selected = manifest.loc[
                np.isclose(manifest["flow_speed_mm_s"].astype(float), speed)
                & manifest["position_label"].astype(str).eq(position)
            ]
            if len(selected) != 1:
                raise ValueError(
                    f"expected one frame for speed={speed}, position={position}"
                )
            manifest_row = selected.iloc[0]
            scan_id = str(manifest_row["scan_id"])
            frame_index = int(manifest_row["bscan_index"])
            maps = load_frame_maps(
                resolve_source_path(
                    manifest_path, str(manifest_row["source_file"])
                )
            )
            location = _localization_row(
                localization, scan_id, frame_index
            )
            shown, vmin, vmax = _display_log(maps.omag_raw)
            dx_um = float(location["dx_um"])
            dz_um = float(location["dz_um"])
            ax.imshow(
                shown,
                cmap="gray",
                aspect=dz_um / dx_um,
                vmin=vmin,
                vmax=vmax,
            )
            x_left = float(location["x_left_edge_px"])
            x_right = float(location["x_right_edge_px"])
            x_center = (x_left + x_right) / 2.0
            z_top = float(location["z_top_edge_px"])
            z_bottom = float(location["z_bottom_edge_px"])
            width = x_right - x_left
            height = z_bottom - z_top
            ax.add_patch(
                Ellipse(
                    (x_center, (z_top + z_bottom) / 2.0),
                    width=width,
                    height=height,
                    fill=False,
                    color="#00e5ff",
                    linewidth=1.8,
                )
            )
            ax.axvline(x_center, color="#00ff72", linewidth=1.3)
            ax.axvline(x_left, color="#00e5ff", linewidth=1.0)
            ax.axvline(x_right, color="#00e5ff", linewidth=1.0)
            ax.axhline(z_top, color="#ff4d4d", linewidth=1.2)
            ax.axhline(z_bottom, color="#ff9f1c", linewidth=1.2)
            x_pad = max(18.0, 1.7 * width)
            ax.set_xlim(
                max(-0.5, x_left - x_pad),
                min(maps.omag_raw.shape[1] - 0.5, x_right + x_pad),
            )
            ax.set_ylim(
                min(maps.omag_raw.shape[0] - 0.5, z_bottom + 25),
                max(-0.5, z_top - 25),
            )
            assessment = str(
                location["mentor_vessel_presence_prediction"]
            )
            score = float(location["mentor_assessability_score"])
            ax.set_title(
                f"{speed:g} mm/s | {position} | B-scan {frame_index}"
                + chr(10)
                + f"X4={x_center:.2f}, z_top={z_top + 0.5:.1f}, "
                f"width={width:.0f}px | {assessment} {score:.2f}",
                fontsize=9,
            )
            if row_index == len(speeds) - 1:
                ax.set_xlabel("A-line x (pixel centre)")
            if column_index == 0:
                ax.set_ylabel("Depth z (pixel centre)")
    handles = [
        Line2D([0], [0], color="#00ff72", label="X4 central A-line"),
        Line2D([0], [0], color="#00e5ff", label="X1-width left/right"),
        Line2D([0], [0], color="#ff4d4d", label="z_upper top"),
        Line2D([0], [0], color="#ff9f1c", label="128 um physical bottom"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=9)
    fig.suptitle(
        "Hybrid localization: mentor full-volume tracking, our SV geometry",
        fontsize=15,
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def plot_trajectories(
    manifest_path: Path,
    manifest: pd.DataFrame,
    output_path: Path,
) -> None:
    scan_rows = (
        manifest.sort_values("flow_speed_mm_s")
        .drop_duplicates("scan_id")
        .reset_index(drop=True)
    )
    fig, axes = plt.subplots(
        len(scan_rows), 2, figsize=(15, 13), sharex=True, constrained_layout=True
    )
    for row_index, manifest_row in scan_rows.iterrows():
        scan_id = str(manifest_row["scan_id"])
        speed = float(manifest_row["flow_speed_mm_s"])
        tracking_path = resolve_source_path(
            manifest_path, str(manifest_row["tracking_file"])
        )
        tracking = pd.read_csv(tracking_path)
        frame = tracking["frame_index"].to_numpy(float)
        x_viterbi = tracking["x_center_px"].to_numpy(float)
        x4 = tracking[
            "x4_centroid_isolated_jump_corrected_px"
        ].to_numpy(float)
        z_upper = tracking["z_upper_px"].to_numpy(float)
        assessment = tracking["vessel_presence_prediction"].astype(str)
        selected_frames = (
            manifest.loc[
                manifest["scan_id"].astype(str).eq(scan_id),
                "bscan_index",
            ]
            .astype(int)
            .to_numpy()
        )
        x_axis, z_axis = axes[row_index]
        x_axis.plot(
            frame,
            x_viterbi,
            color="#999999",
            linewidth=0.8,
            label="Viterbi",
        )
        x_axis.plot(frame, x4, color="#0066cc", linewidth=1.1, label="X4")
        weak = ~assessment.eq("assessable")
        if weak.any():
            x_axis.scatter(
                frame[weak],
                x4[weak],
                s=7,
                color="#e76f51",
                label="uncertain/not assessable",
                zorder=3,
            )
        x_axis.scatter(
            selected_frames,
            x4[selected_frames],
            s=28,
            marker="o",
            facecolor="#ffcc00",
            edgecolor="black",
            linewidth=0.5,
            label="15-frame pilot sample",
            zorder=4,
        )
        z_axis.plot(frame, z_upper, color="#7b2cbf", linewidth=1.1)
        z_axis.scatter(
            selected_frames,
            z_upper[selected_frames],
            s=28,
            marker="o",
            facecolor="#ffcc00",
            edgecolor="black",
            linewidth=0.5,
            zorder=4,
        )
        z_axis.invert_yaxis()
        counts = assessment.value_counts()
        x_axis.set_ylabel(
            f"{scan_id} ({speed:g} mm/s)" + chr(10) + "x (px)"
        )
        z_axis.set_ylabel("z_upper (px)")
        x_axis.grid(alpha=0.2)
        z_axis.grid(alpha=0.2)
        z_axis.set_title(
            " | ".join(
                f"{name}={int(counts.get(name, 0))}"
                for name in ("assessable", "uncertain", "not_assessable")
            ),
            fontsize=9,
        )
        if row_index == 0:
            x_axis.set_title("Lateral full-volume trajectory")
            z_axis.set_title(
                "Axial upper-edge trajectory | " + z_axis.get_title()
            )
        if row_index == len(scan_rows) - 1:
            x_axis.set_xlabel("B-scan index")
            z_axis.set_xlabel("B-scan index")
    axes[0, 0].legend(loc="best", fontsize=7)
    fig.suptitle(
        "Mentor slow-axis tracking across every 500-frame volume",
        fontsize=15,
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def plot_localization_comparison(
    manifest_path: Path,
    manifest: pd.DataFrame,
    old_localization: pd.DataFrame,
    new_localization: pd.DataFrame,
    output_path: Path,
    reference_label: str,
) -> None:
    """Overlay the reference and current axial borders on the same OMAG maps."""

    speeds = sorted(manifest["flow_speed_mm_s"].astype(float).unique())
    positions = ["front", "middle", "rear"]
    fig, axes = plt.subplots(
        len(speeds),
        len(positions),
        figsize=(15, 15),
        constrained_layout=True,
    )
    for row_index, speed in enumerate(speeds):
        for column_index, position in enumerate(positions):
            ax = axes[row_index, column_index]
            selected = manifest.loc[
                np.isclose(manifest["flow_speed_mm_s"].astype(float), speed)
                & manifest["position_label"].astype(str).eq(position)
            ]
            if len(selected) != 1:
                raise ValueError(
                    f"expected one frame for speed={speed}, position={position}"
                )
            manifest_row = selected.iloc[0]
            scan_id = str(manifest_row["scan_id"])
            frame_index = int(manifest_row["bscan_index"])
            maps = load_frame_maps(
                resolve_source_path(
                    manifest_path, str(manifest_row["source_file"])
                )
            )
            old = _localization_row(old_localization, scan_id, frame_index)
            new = _localization_row(new_localization, scan_id, frame_index)
            shown, vmin, vmax = _display_log(maps.omag_raw)
            dx_um = float(new["dx_um"])
            dz_um = float(new["dz_um"])
            ax.imshow(
                shown,
                cmap="gray",
                aspect=dz_um / dx_um,
                vmin=vmin,
                vmax=vmax,
            )

            old_top = float(old["z_top_edge_px"])
            old_bottom = float(old["z_bottom_edge_px"])
            new_top = float(new["z_top_edge_px"])
            new_bottom = float(new["z_bottom_edge_px"])
            x_left = float(new["x_left_edge_px"])
            x_right = float(new["x_right_edge_px"])
            x_center = 0.5 * (x_left + x_right)
            width = x_right - x_left
            ax.add_patch(
                Ellipse(
                    (x_center, 0.5 * (new_top + new_bottom)),
                    width=width,
                    height=new_bottom - new_top,
                    fill=False,
                    color="#00e5ff",
                    linewidth=1.5,
                )
            )
            ax.axvline(x_center, color="#00ff72", linewidth=1.1)
            ax.axvline(x_left, color="#00e5ff", linewidth=0.9)
            ax.axvline(x_right, color="#00e5ff", linewidth=0.9)
            ax.axhline(old_top, color="#ffe066", linewidth=1.5, linestyle="--")
            ax.axhline(old_bottom, color="#ffe066", linewidth=1.0, linestyle=":")
            ax.axhline(new_top, color="#ff4d4d", linewidth=1.5)
            ax.axhline(new_bottom, color="#ff9f1c", linewidth=1.2)

            x_pad = max(18.0, 1.7 * width)
            z_min = min(old_top, new_top)
            z_max = max(old_bottom, new_bottom)
            ax.set_xlim(
                max(-0.5, x_left - x_pad),
                min(maps.omag_raw.shape[1] - 0.5, x_right + x_pad),
            )
            ax.set_ylim(
                min(maps.omag_raw.shape[0] - 0.5, z_max + 25),
                max(-0.5, z_min - 25),
            )
            delta_px = new_top - old_top
            ax.set_title(
                f"{speed:g} mm/s | {position} | B-scan {frame_index}"
                + chr(10)
                + f"top: {reference_label}={old_top + 0.5:.1f}, "
                + f"v2={new_top + 0.5:.1f}, delta={delta_px:+.0f}px "
                + f"({delta_px * dz_um:+.1f} um)",
                fontsize=9,
            )
            if row_index == len(speeds) - 1:
                ax.set_xlabel("A-line x (pixel centre)")
            if column_index == 0:
                ax.set_ylabel("Depth z (pixel centre)")
    handles = [
        Line2D(
            [0],
            [0],
            color="#ffe066",
            linestyle="--",
            label=f"{reference_label} top",
        ),
        Line2D([0], [0], color="#ff4d4d", label="v2 top"),
        Line2D([0], [0], color="#ff9f1c", label="v2 128 um bottom"),
        Line2D([0], [0], color="#00e5ff", label="v2 lateral geometry"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=9)
    fig.suptitle(
        "Upper-edge localization comparison on the same OMAG B-scans",
        fontsize=15,
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def write_comparison(
    old_path: Path,
    new_localization: pd.DataFrame,
    output_path: Path,
    reference_label: str,
) -> None:
    old = pd.read_csv(old_path)
    old = old.set_index("frame_id", drop=False)
    records: list[dict[str, object]] = []
    for _, new in new_localization.iterrows():
        frame_id = str(new["frame_id"])
        previous = old.loc[frame_id]
        old_x = (
            float(previous["x_left_edge_px"])
            + float(previous["x_right_edge_px"])
        ) / 2.0
        new_x = (
            float(new["x_left_edge_px"]) + float(new["x_right_edge_px"])
        ) / 2.0
        old_z = float(previous["z_top_edge_px"]) + 0.5
        new_z = float(new["z_top_edge_px"]) + 0.5
        records.append(
            {
                "scan_id": new["scan_id"],
                "frame_id": frame_id,
                "bscan_index": int(new["mentor_frame_index"]),
                "reference_label": reference_label,
                "reference_x_center_px": old_x,
                "v2_X4_center_px": new_x,
                "delta_x_px_new_minus_old": new_x - old_x,
                "delta_x_um_new_minus_old": (
                    new_x - old_x
                ) * float(new["dx_um"]),
                "reference_z_top_center_px": old_z,
                "v2_z_upper_center_px": new_z,
                "delta_z_px_new_minus_old": new_z - old_z,
                "delta_z_um_new_minus_old": (
                    new_z - old_z
                ) * float(new["dz_um"]),
                "old_source_qc_valid": previous["source_qc_valid"],
                "new_source_qc_valid": new["source_qc_valid"],
                "new_assessability": new[
                    "mentor_vessel_presence_prediction"
                ],
            }
        )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame.from_records(records).to_csv(output_path, index=False)


def main() -> int:
    args = _arguments()
    manifest_path = args.manifest.resolve()
    manifest = pd.read_csv(manifest_path)
    localization = pd.read_csv(args.localization)
    plot_overview(
        manifest_path, manifest, localization, args.overview_output
    )
    plot_trajectories(manifest_path, manifest, args.trajectory_output)
    if args.old_localization is not None:
        if args.comparison_output is None:
            raise ValueError(
                "--comparison-output is required with --old-localization"
            )
        write_comparison(
            args.old_localization,
            localization,
            args.comparison_output,
            args.reference_label,
        )
        if args.comparison_figure_output is not None:
            plot_localization_comparison(
                manifest_path,
                manifest,
                pd.read_csv(args.old_localization),
                localization,
                args.comparison_figure_output,
                args.reference_label,
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
