#!/usr/bin/env python
"""Generate a condition-blinded manual localization-QC package.

Sampling is deliberately independent of all tail metrics and condition trends.
For each of the nine predeclared representative scans, one frame is sampled
from every crossing of three slow-axis regions and three old-confidence strata
(high, model-assisted, low-or-region-edge), giving 81 images in total.

The reviewer package and the restricted decoding table are written into
separate directories.  Images and reviewer-facing files contain only random
``QC_###`` identifiers.  New x/z confidence components are copied to the
restricted mapping when available, but never influence frame selection.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import shutil
import sys
import uuid
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from tailq.io import (  # noqa: E402
    export_path,
    load_config,
    load_manifest,
    output_root,
    read_export_metadata,
    tracking_csv_path,
)


REPRESENTATIVE_IDS = (
    "d185_s1_scan32",
    "d185_s5_scan38",
    "d185_s10_scan5",
    "d235_s1_scan17",
    "d235_s5_scan23",
    "d235_s10_scan29",
    "d285_s1_scan2",
    "d285_s5_scan8",
    "d285_s10_scan14",
)
DEFAULT_SEED = 20260818
PRIMARY_ALPHA = 0.15
SLOW_AXIS_REGIONS = (
    ("front", 0, 167),
    ("middle", 167, 333),
    ("rear", 333, 500),
)
SAMPLING_STRATA = ("old_high", "old_model_assisted", "old_low_or_region_edge")
REVIEW_COLUMNS = (
    "blind_id",
    "human_x_correct",
    "human_upper_edge_correct",
    "human_lower_boundary_reasonable",
    "human_roi_reasonable",
    "image_abnormal",
    "comments",
)
FORBIDDEN_REVIEWER_TOKENS = (
    "d185",
    "d235",
    "d285",
    "mm/s",
    "mm_per_s",
    "diameter",
    "velocity",
    "scan",
    "rnt",
    "auc",
    "trend",
    "high tail",
    "low tail",
)
BLIND_ID_PATTERN = re.compile(r"^QC_[0-9]{3}$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        default=str(PROJECT_ROOT / "configs" / "pilot_20260414.json"),
        help="Frozen pilot JSON configuration.",
    )
    parser.add_argument(
        "--tracking-root",
        default=None,
        help="Tracking directory; defaults to <configured output>/tracking.",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Package root; defaults to <configured output>/blinded_qc.",
    )
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument(
        "--require-new-confidence",
        action="store_true",
        help="Fail unless every tracking CSV contains the revised x/z confidence fields.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace an existing generated package. Never affects source data.",
    )
    return parser.parse_args()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _stable_rng(seed: int, namespace: str) -> np.random.Generator:
    material = f"{seed}|{namespace}".encode("utf-8")
    derived = int.from_bytes(hashlib.sha256(material).digest()[:8], "little")
    return np.random.default_rng(derived)


def _first_present(table: pd.DataFrame, names: Iterable[str], default: Any) -> pd.Series:
    for name in names:
        if name in table.columns:
            return table[name]
    return pd.Series([default] * len(table), index=table.index)


def normalize_tracking_table(path: Path, *, require_new_confidence: bool = False) -> pd.DataFrame:
    """Return one canonical row per zero-based frame.

    Old-only tracking files remain usable for package development.  A formal
    run can require revised confidence columns with ``--require-new-confidence``.
    """

    table = pd.read_csv(path)
    required = {
        "frame_index",
        "x_center_px",
        "z_upper_px",
        "z_lower_op_px",
        "z_tail_start_px",
        "central_x0_px",
        "central_x1_exclusive_px",
        "side_left_x0_px",
        "side_left_x1_exclusive_px",
        "side_right_x0_px",
        "side_right_x1_exclusive_px",
    }
    missing = required.difference(table.columns)
    if missing:
        raise ValueError(f"{path}: missing localization/ROI columns {sorted(missing)}")

    table = table.copy()
    table["frame_index"] = pd.to_numeric(table["frame_index"], errors="raise").astype(int)
    if table["frame_index"].duplicated().any():
        raise ValueError(f"{path}: duplicate frame indices")
    if set(table["frame_index"]) != set(range(500)):
        raise ValueError(f"{path}: expected exactly zero-based frames 0..499")

    table["old_tracking_class"] = _first_present(
        table, ("old_tracking_class", "tracking_class"), "unavailable"
    ).astype(str)
    table["old_tracking_confidence"] = pd.to_numeric(
        _first_present(table, ("old_tracking_confidence", "tracking_confidence"), np.nan),
        errors="coerce",
    )

    new_columns: Mapping[str, tuple[str, ...]] = {
        "x_path_confidence": ("x_path_confidence", "new_x_path_confidence"),
        "x_path_confidence_class": (
            "x_path_confidence_class",
            "new_x_path_confidence_class",
        ),
        "z_edge_confidence": ("z_edge_confidence", "new_z_edge_confidence"),
        "z_edge_confidence_class": (
            "z_edge_confidence_class",
            "new_z_edge_confidence_class",
        ),
        "overall_tracking_confidence": (
            "overall_tracking_confidence",
            "new_tracking_confidence",
        ),
        "new_tracking_class": ("new_tracking_class", "overall_tracking_confidence_class"),
    }
    for canonical, aliases in new_columns.items():
        default = np.nan if "class" not in canonical else ""
        table[canonical] = _first_present(table, aliases, default)
    for column in ("x_path_confidence", "z_edge_confidence", "overall_tracking_confidence"):
        table[column] = pd.to_numeric(table[column], errors="coerce")
    for column in ("x_path_confidence_class", "z_edge_confidence_class", "new_tracking_class"):
        table[column] = table[column].fillna("").astype(str)

    revised_complete = (
        table["x_path_confidence"].notna().all()
        and table["z_edge_confidence"].notna().all()
        and table["overall_tracking_confidence"].notna().all()
        and table["new_tracking_class"].ne("").all()
    )
    table.attrs["revised_confidence_complete"] = bool(revised_complete)
    if require_new_confidence and not revised_complete:
        raise ValueError(f"{path}: revised x/z/overall confidence is incomplete")

    numeric_required = [
        "x_center_px",
        "z_upper_px",
        "z_lower_op_px",
        "z_tail_start_px",
        "central_x0_px",
        "central_x1_exclusive_px",
        "side_left_x0_px",
        "side_left_x1_exclusive_px",
        "side_right_x0_px",
        "side_right_x1_exclusive_px",
    ]
    for column in numeric_required:
        table[column] = pd.to_numeric(table[column], errors="coerce")
    if not np.isfinite(table[numeric_required].to_numpy(float)).all():
        raise ValueError(f"{path}: non-finite localization or ROI coordinate")
    return table.sort_values("frame_index").reset_index(drop=True)


def _choose_one(
    candidates: pd.DataFrame,
    used: set[int],
    rng: np.random.Generator,
) -> tuple[pd.Series, int]:
    available = candidates.loc[~candidates["frame_index"].isin(used)]
    if available.empty:
        raise ValueError("Sampling cell has no unused candidates")
    position = int(rng.integers(0, len(available)))
    return available.iloc[position], len(available)


def sample_scan_frames(table: pd.DataFrame, scan_id: str, seed: int) -> pd.DataFrame:
    """Sample nine unique frames using old confidence only."""

    rng = _stable_rng(seed, f"sample|{scan_id}")
    used: set[int] = set()
    selected: list[dict[str, Any]] = []
    finite_conf = table["old_tracking_confidence"].replace([np.inf, -np.inf], np.nan)

    for region, start, stop in SLOW_AXIS_REGIONS:
        region_table = table.loc[table["frame_index"].between(start, stop - 1)].copy()
        if len(region_table) != stop - start:
            raise ValueError(f"{scan_id}: incomplete {region} slow-axis region")
        region_conf = finite_conf.loc[region_table.index]
        q20 = float(region_conf.quantile(0.20)) if region_conf.notna().any() else math.nan
        edge_width = max(6, int(round(0.05 * len(region_table))))
        at_region_edge = (
            (region_table["frame_index"] < start + edge_width)
            | (region_table["frame_index"] >= stop - edge_width)
        )
        low_evidence = region_table["old_tracking_confidence"] <= q20

        requested_pools: tuple[tuple[str, pd.DataFrame], ...] = (
            (
                "old_high",
                region_table.loc[region_table["old_tracking_class"] == "high_confidence"],
            ),
            (
                "old_model_assisted",
                region_table.loc[region_table["old_tracking_class"] == "model_assisted"],
            ),
            (
                "old_low_or_region_edge",
                region_table.loc[
                    (region_table["old_tracking_class"] == "failed")
                    | low_evidence
                    | at_region_edge
                ],
            ),
        )

        for stratum, requested_pool in requested_pools:
            fallback = "none"
            pool = requested_pool
            if pool.loc[~pool["frame_index"].isin(used)].empty:
                if stratum == "old_model_assisted":
                    pool = region_table.nsmallest(max(1, len(region_table) // 5), "old_tracking_confidence")
                    fallback = "no_old_model_assisted_in_region_used_lowest_confidence_quintile"
                elif stratum == "old_high":
                    pool = region_table
                    fallback = "no_old_high_in_region_used_any_region_frame"
                else:
                    pool = region_table
                    fallback = "no_low_or_edge_candidate_used_any_region_frame"
            row, pool_size = _choose_one(pool, used, rng)
            frame = int(row["frame_index"])
            used.add(frame)
            selected.append(
                {
                    "dataset_id": scan_id,
                    "frame_index": frame,
                    "slow_axis_region": region,
                    "sampling_stratum": stratum,
                    "sampling_fallback": fallback,
                    "sampling_pool_size": pool_size,
                    "sampling_seed": seed,
                }
            )

    result = pd.DataFrame(selected)
    if len(result) != 9 or result["frame_index"].nunique() != 9:
        raise AssertionError(f"{scan_id}: expected nine unique sampled frames")
    return result


def build_sampling_plan(
    tracking_tables: Mapping[str, pd.DataFrame], seed: int
) -> pd.DataFrame:
    rows = [sample_scan_frames(tracking_tables[scan_id], scan_id, seed) for scan_id in REPRESENTATIVE_IDS]
    plan = pd.concat(rows, ignore_index=True)
    order_rng = _stable_rng(seed, "global_blind_order")
    plan = plan.iloc[order_rng.permutation(len(plan))].reset_index(drop=True)
    plan.insert(0, "blind_id", [f"QC_{index:03d}" for index in range(1, len(plan) + 1)])
    return plan


def _dataset_name(handle: h5py.File) -> str:
    for candidate in ("p_bld_ed", "p_bld_ed_sparse"):
        if candidate in handle:
            return candidate
    raise KeyError("Flow MAT lacks p_bld_ed/p_bld_ed_sparse")


def read_selected_frames(path: Path, frames: Iterable[int]) -> dict[int, np.ndarray]:
    """Read only requested frames and return each as ``[z, x]`` float32."""

    metadata = read_export_metadata(path)
    axis_order = str(metadata.get("axis_order_h5py", "")).replace(" ", "").lower()
    if axis_order not in {"frame,x,z", "frame,z,x"}:
        raise ValueError(f"{path}: strict sidecar axis_order_h5py is required")
    source = metadata.get("source_frame_indices")
    with h5py.File(path, "r") as handle:
        dataset = handle[_dataset_name(handle)]
        if source is None:
            source_zero = np.arange(dataset.shape[0], dtype=int)
        else:
            source_zero = np.asarray(source, dtype=int)
            if source_zero.size and source_zero.min() >= 1:
                source_zero = source_zero - 1
        position_by_frame = {int(value): index for index, value in enumerate(source_zero)}
        output: dict[int, np.ndarray] = {}
        for frame in sorted(set(int(value) for value in frames)):
            if frame not in position_by_frame:
                raise IndexError(f"{path}: source frame {frame} is unavailable")
            raw = np.asarray(dataset[position_by_frame[frame]], dtype=np.float32)
            image = raw.T if axis_order == "frame,x,z" else raw
            if image.ndim != 2 or not np.isfinite(image).all():
                raise ValueError(f"{path}: invalid Flow frame {frame}")
            output[frame] = image
    return output


def _display_image(image: np.ndarray) -> tuple[np.ndarray, float, float]:
    finite = image[np.isfinite(image)]
    scale = max(float(np.percentile(finite, 75)), np.finfo(float).eps)
    transformed = np.arcsinh(np.maximum(image, 0) / scale)
    return (
        transformed,
        float(np.percentile(transformed, 2)),
        float(np.percentile(transformed, 99.7)),
    )


def _draw_overlay(ax: plt.Axes, row: pd.Series, z_size: int) -> None:
    xc = float(row["x_center_px"])
    z_upper = float(row["z_upper_px"])
    z_lower = float(row["z_lower_op_px"])
    z_start = float(row["z_tail_start_px"])
    ax.axvspan(
        float(row["central_x0_px"]),
        float(row["central_x1_exclusive_px"]),
        color="#00ffff",
        alpha=0.11,
        lw=0,
    )
    for left, right in (
        ("side_left_x0_px", "side_left_x1_exclusive_px"),
        ("side_right_x0_px", "side_right_x1_exclusive_px"),
    ):
        ax.axvspan(float(row[left]), float(row[right]), color="#ffe066", alpha=0.11, lw=0)
    ax.axvline(xc, color="#ffff00", lw=1.15, ls="--")
    ax.axhline(z_upper, color="#00ffff", lw=1.25)
    ax.axhline(z_lower, color="#ff9f1c", lw=1.25, ls="--")
    ax.axhline(z_start, color="#ff2e88", lw=1.35)

    shell_styles = (
        (100, "#7cff6b", (0, (2, 2))),
        (200, "#4cc9f0", (0, (5, 2))),
        (300, "#c77dff", (0, (1, 2))),
    )
    outside: list[int] = []
    for window_um, color, line_style in shell_styles:
        endpoint = z_start + round(window_um / 6.7)
        if endpoint <= z_size - 1:
            ax.axhline(endpoint, color=color, lw=1.05, ls=line_style)
        else:
            # The true endpoint is outside the acquired field; mark the image
            # boundary rather than inventing or clipping signal values.
            ax.plot(
                [xc],
                [z_size - 1],
                marker="v",
                ms=5,
                color=color,
                markeredgecolor="black",
                markeredgewidth=0.35,
                clip_on=True,
            )
            outside.append(window_um)
    if outside:
        label = "/".join(str(value) for value in outside) + " um endpoint outside FOV"
        ax.text(
            0.985,
            0.012,
            label,
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            color="white",
            fontsize=7.5,
            bbox={"facecolor": "black", "alpha": 0.58, "edgecolor": "none", "pad": 2},
        )


def render_blinded_image(image: np.ndarray, row: pd.Series, blind_id: str, output: Path) -> None:
    shown, vmin, vmax = _display_image(image)
    z_size, x_size = shown.shape
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 7.2))
    fig.subplots_adjust(left=0.06, right=0.99, bottom=0.22, top=0.88, wspace=0.15)
    for axis in axes:
        axis.imshow(shown, cmap="gray", vmin=vmin, vmax=vmax, origin="upper", aspect="auto")
        _draw_overlay(axis, row, z_size)
        axis.set_xlabel("lateral pixel")
        axis.set_ylabel("axial pixel")
    axes[0].set_title("Full B-scan")
    axes[0].set_xlim(0, x_size - 1)
    axes[0].set_ylim(z_size - 1, 0)

    x_values = [
        float(row["side_left_x0_px"]),
        float(row["side_right_x1_exclusive_px"]),
        float(row["x_center_px"]),
    ]
    xmin = max(0, min(x_values) - 12)
    xmax = min(x_size - 1, max(x_values) + 12)
    zmin = max(0, float(row["z_upper_px"]) - 25)
    axes[1].set_title("Local ROI view")
    axes[1].set_xlim(xmin, xmax)
    axes[1].set_ylim(z_size - 1, zmin)

    handles = [
        Line2D([0], [0], color="#ffff00", ls="--", lw=1.3, label="x center"),
        Line2D([0], [0], color="#00ffff", lw=1.3, label="Flow upper edge"),
        Line2D([0], [0], color="#ff9f1c", ls="--", lw=1.3, label="operational lower boundary"),
        Line2D([0], [0], color="#ff2e88", lw=1.3, label="tail start after 2-px guard"),
        Line2D([0], [0], color="#7cff6b", ls=(0, (2, 2)), lw=1.2, label="100 um endpoint"),
        Line2D([0], [0], color="#4cc9f0", ls=(0, (5, 2)), lw=1.2, label="200 um endpoint"),
        Line2D([0], [0], color="#c77dff", ls=(0, (1, 2)), lw=1.2, label="300 um endpoint"),
        Patch(facecolor="#00ffff", alpha=0.18, label="central ROI"),
        Patch(facecolor="#ffe066", alpha=0.18, label="left/right background ROIs"),
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.015),
        ncol=3,
        frameon=False,
        fontsize=8.5,
    )
    fig.suptitle(blind_id, y=0.965, fontsize=14, fontweight="bold")
    fig.savefig(
        output,
        dpi=180,
        bbox_inches="tight",
        metadata={"Title": blind_id, "Description": "Condition-blinded localization QC"},
    )
    plt.close(fig)


def reviewer_sheet(blind_ids: Iterable[str]) -> pd.DataFrame:
    rows = []
    for blind_id in blind_ids:
        row = {column: "" for column in REVIEW_COLUMNS}
        row["blind_id"] = blind_id
        rows.append(row)
    return pd.DataFrame(rows, columns=REVIEW_COLUMNS)


def _mapping_row(
    plan_row: pd.Series,
    tracking_row: pd.Series,
    manifest_row: pd.Series,
) -> dict[str, Any]:
    return {
        "blind_id": plan_row["blind_id"],
        "image_filename": f"{plan_row['blind_id']}.png",
        "dataset_id": plan_row["dataset_id"],
        "original_scan_number": int(manifest_row["scan_number"]),
        "diameter_um": float(manifest_row["diameter_um"]),
        "velocity_mm_s": float(manifest_row["velocity_mm_s"]),
        "technical_replicate": int(manifest_row["technical_replicate"]),
        "frame_index_0based": int(plan_row["frame_index"]),
        "frame_number_1based": int(plan_row["frame_index"]) + 1,
        "slow_axis_region": plan_row["slow_axis_region"],
        "sampling_stratum": plan_row["sampling_stratum"],
        "sampling_fallback": plan_row["sampling_fallback"],
        "sampling_pool_size": int(plan_row["sampling_pool_size"]),
        "sampling_seed": int(plan_row["sampling_seed"]),
        "old_tracking_class": tracking_row["old_tracking_class"],
        "old_tracking_confidence": tracking_row["old_tracking_confidence"],
        "x_path_confidence": tracking_row["x_path_confidence"],
        "x_path_confidence_class": tracking_row["x_path_confidence_class"],
        "z_edge_confidence": tracking_row["z_edge_confidence"],
        "z_edge_confidence_class": tracking_row["z_edge_confidence_class"],
        "overall_tracking_confidence": tracking_row["overall_tracking_confidence"],
        "new_tracking_class": tracking_row["new_tracking_class"],
        "x_center_px": tracking_row["x_center_px"],
        "z_upper_px": tracking_row["z_upper_px"],
        "z_lower_op_px": tracking_row["z_lower_op_px"],
        "z_tail_start_px": tracking_row["z_tail_start_px"],
    }


def _write_reviewer_instructions(path: Path) -> None:
    text = """# Blinded localization QC instructions

Review only the PNG files in `images/` and enter assessments in
`blinded_qc_sheet.csv`. Do not open any decoding/mapping table until all rows
have been scored and the completed CSV has been saved.

For each criterion use only: `yes`, `no`, or `uncertain`.

- `human_x_correct`: the marked lateral center follows the intended vessel.
- `human_upper_edge_correct`: the Flow upper-edge line is visually supported.
- `human_lower_boundary_reasonable`: the operational lower boundary is plausible.
- `human_roi_reasonable`: central and bilateral background ROIs avoid obvious contamination.
- `image_abnormal`: use `yes` for acquisition/processing abnormality unrelated to ordinary signal weakness.
- `comments`: optional short explanation; do not infer or enter experimental conditions.

Score every image independently. Do not delete, reorder, or selectively omit rows.
"""
    path.write_text(text, encoding="utf-8")


def _zip_reviewer_package(reviewer_dir: Path, archive: Path) -> None:
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as handle:
        for path in sorted(reviewer_dir.rglob("*")):
            if path.is_file():
                handle.write(path, path.relative_to(reviewer_dir.parent))


def self_check_package(package_root: Path, expected_count: int = 81) -> dict[str, Any]:
    reviewer = package_root / "reviewer_package"
    restricted = package_root / "restricted_mapping"
    images = sorted((reviewer / "images").glob("*.png"))
    sheet = pd.read_csv(reviewer / "blinded_qc_sheet.csv", keep_default_na=False)
    mapping = pd.read_csv(restricted / "blinded_qc_mapping.csv", keep_default_na=False)
    expected_ids = [f"QC_{index:03d}" for index in range(1, expected_count + 1)]

    if [path.stem for path in images] != expected_ids:
        raise AssertionError("Reviewer image filenames are incomplete or non-blinded")
    if sheet.columns.tolist() != list(REVIEW_COLUMNS):
        raise AssertionError("Reviewer sheet schema changed")
    if sheet["blind_id"].tolist() != expected_ids:
        raise AssertionError("Reviewer sheet blind IDs do not match images")
    if mapping["blind_id"].tolist() != expected_ids:
        raise AssertionError("Restricted mapping blind IDs do not match images")
    if not all(BLIND_ID_PATTERN.fullmatch(value) for value in sheet["blind_id"]):
        raise AssertionError("Invalid blind ID")
    response_columns = list(REVIEW_COLUMNS[1:])
    if sheet[response_columns].astype(str).ne("").any().any():
        raise AssertionError("Reviewer responses are not initially blank")

    if mapping["dataset_id"].nunique() != 9:
        raise AssertionError("Not all nine scans are represented")
    per_scan = mapping.groupby("dataset_id").size()
    if not (per_scan == 9).all():
        raise AssertionError("Every scan must contribute exactly nine images")
    crossed = mapping.groupby(["dataset_id", "slow_axis_region", "sampling_stratum"]).size()
    if len(crossed) != 9 * 3 * 3 or not (crossed == 1).all():
        raise AssertionError("Scan x region x sampling-stratum balance failed")
    if mapping.duplicated(["dataset_id", "frame_index_0based"]).any():
        raise AssertionError("A frame was selected twice within a scan")

    reviewer_names_and_text = []
    for path in reviewer.rglob("*"):
        reviewer_names_and_text.append(str(path.relative_to(reviewer)).lower())
        if path.suffix.lower() in {".csv", ".md", ".txt"}:
            reviewer_names_and_text.append(path.read_text(encoding="utf-8").lower())
    joined = "\n".join(reviewer_names_and_text)
    leaked = [token for token in FORBIDDEN_REVIEWER_TOKENS if token in joined]
    if leaked:
        raise AssertionError(f"Reviewer package contains forbidden token(s): {leaked}")

    with zipfile.ZipFile(package_root / "blinded_qc_reviewer_package.zip") as archive:
        archived_names = archive.namelist()
        if any("mapping" in name.lower() for name in archived_names):
            raise AssertionError("Restricted mapping leaked into reviewer ZIP")

    return {
        "status": "passed",
        "image_count": len(images),
        "sheet_rows": len(sheet),
        "mapping_rows": len(mapping),
        "scan_count": int(mapping["dataset_id"].nunique()),
        "per_scan_images": sorted(per_scan.astype(int).unique().tolist()),
        "sampling_fallback_count": int((mapping["sampling_fallback"] != "none").sum()),
        "reviewer_forbidden_token_hits": [],
        "reviewer_zip_excludes_mapping": True,
    }


def generate_package(args: argparse.Namespace) -> Path:
    cfg, config_path = load_config(args.config)
    manifest = load_manifest(cfg, config_path).set_index("dataset_id", drop=False)
    missing_manifest = set(REPRESENTATIVE_IDS).difference(manifest.index)
    if missing_manifest:
        raise ValueError(f"Manifest lacks representative scans: {sorted(missing_manifest)}")
    tracking_root = (
        Path(args.tracking_root).resolve()
        if args.tracking_root
        else output_root(cfg) / "tracking"
    )
    final_root = (
        Path(args.output_dir).resolve()
        if args.output_dir
        else output_root(cfg) / "blinded_qc"
    )
    if final_root.exists() and not args.overwrite:
        raise FileExistsError(f"Output already exists; use --overwrite to replace: {final_root}")

    tables: dict[str, pd.DataFrame] = {}
    tracking_paths: dict[str, Path] = {}
    for scan_id in REPRESENTATIVE_IDS:
        path = tracking_csv_path(tracking_root, scan_id, PRIMARY_ALPHA)
        tables[scan_id] = normalize_tracking_table(
            path, require_new_confidence=args.require_new_confidence
        )
        tracking_paths[scan_id] = path
    plan = build_sampling_plan(tables, args.seed)

    building = final_root.parent / f".{final_root.name}.building-{uuid.uuid4().hex}"
    reviewer_dir = building / "reviewer_package"
    image_dir = reviewer_dir / "images"
    restricted_dir = building / "restricted_mapping"
    image_dir.mkdir(parents=True, exist_ok=False)
    restricted_dir.mkdir(parents=True, exist_ok=False)

    mapping_rows: list[dict[str, Any]] = []
    try:
        for scan_id in REPRESENTATIVE_IDS:
            scan_plan = plan.loc[plan["dataset_id"] == scan_id]
            frames = read_selected_frames(export_path(cfg, scan_id), scan_plan["frame_index"])
            table = tables[scan_id].set_index("frame_index", drop=False)
            manifest_row = manifest.loc[scan_id]
            for _, plan_row in scan_plan.iterrows():
                frame = int(plan_row["frame_index"])
                tracking_row = table.loc[frame]
                blind_id = str(plan_row["blind_id"])
                render_blinded_image(
                    frames[frame], tracking_row, blind_id, image_dir / f"{blind_id}.png"
                )
                mapping_rows.append(_mapping_row(plan_row, tracking_row, manifest_row))

        mapping = pd.DataFrame(mapping_rows).sort_values("blind_id").reset_index(drop=True)
        mapping.to_csv(restricted_dir / "blinded_qc_mapping.csv", index=False, encoding="utf-8-sig")
        reviewer_sheet(mapping["blind_id"]).to_csv(
            reviewer_dir / "blinded_qc_sheet.csv",
            index=False,
            encoding="utf-8-sig",
            quoting=csv.QUOTE_MINIMAL,
        )
        _write_reviewer_instructions(reviewer_dir / "README_REVIEWER.md")

        run_metadata = {
            "schema_version": "1.0.0",
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "purpose": "condition-blinded manual localization QC",
            "sampling_seed": int(args.seed),
            "sampling_design": (
                "9 scans x 3 slow-axis regions x 3 old-confidence strata; "
                "new confidence does not affect sampling"
            ),
            "representative_scan_count": 9,
            "images_per_scan": 9,
            "total_images": 81,
            "primary_alpha": PRIMARY_ALPHA,
            "tracking_csv_sha256": {
                scan_id: _sha256(path) for scan_id, path in tracking_paths.items()
            },
            "config_path": str(config_path),
            "config_sha256": _sha256(config_path),
            "generator_path": str(Path(__file__).resolve()),
            "generator_sha256": _sha256(Path(__file__).resolve()),
            "new_confidence_complete_for_all_scans": all(
                table.attrs["revised_confidence_complete"] for table in tables.values()
            ),
            "reviewer_package_contains_mapping": False,
        }
        (restricted_dir / "blinded_qc_run_metadata.json").write_text(
            json.dumps(run_metadata, ensure_ascii=False, indent=2), encoding="utf-8"
        )
        _zip_reviewer_package(reviewer_dir, building / "blinded_qc_reviewer_package.zip")

        preliminary = self_check_package(building, expected_count=81)
        (restricted_dir / "blinded_qc_self_check.json").write_text(
            json.dumps(preliminary, ensure_ascii=False, indent=2), encoding="utf-8"
        )

        if final_root.exists():
            shutil.rmtree(final_root)
        building.replace(final_root)
    except Exception:
        shutil.rmtree(building, ignore_errors=True)
        raise
    return final_root


def main() -> int:
    args = parse_args()
    output = generate_package(args)
    check = self_check_package(output, expected_count=81)
    print(json.dumps({"output": str(output), "self_check": check}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
