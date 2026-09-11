#!/usr/bin/env python
"""Build the condition-blinded, pixel-level localization annotation package.

This module deliberately stops before any pixel-error or tail-metric calculation.
It samples the 35 previously disputed frames plus two deterministic, previous-
QC-yes controls per representative scan, exports anonymised raw Flow arrays and
display-only PNGs, and writes a blank reviewer CSV.  The decoding table is kept
in a separate restricted directory and is never included in the reviewer ZIP.

No source DICOM, raw export, tracking file, or existing pipeline file is edited.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import shutil
import sys
import uuid
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

# Reuse only the previously audited read/normalisation helpers.  The legacy
# module's overlay renderer is never called here.
_LEGACY_PATH = PROJECT_ROOT / "06_blinded_qc" / "generate_blinded_qc.py"
_LEGACY_SPEC = importlib.util.spec_from_file_location("tailq_legacy_blinded_qc", _LEGACY_PATH)
if _LEGACY_SPEC is None or _LEGACY_SPEC.loader is None:  # pragma: no cover - import guard
    raise ImportError(f"Cannot load audited helper module: {_LEGACY_PATH}")
_LEGACY = importlib.util.module_from_spec(_LEGACY_SPEC)
_LEGACY_SPEC.loader.exec_module(_LEGACY)

from tailq.io import export_path, load_config, load_manifest, output_root, tracking_csv_path  # noqa: E402


REPRESENTATIVE_IDS = tuple(_LEGACY.REPRESENTATIVE_IDS)
PRIMARY_ALPHA = float(_LEGACY.PRIMARY_ALPHA)
SLOW_AXIS_REGIONS = tuple(_LEGACY.SLOW_AXIS_REGIONS)
REGION_NAMES = tuple(region for region, _, _ in SLOW_AXIS_REGIONS)

DEFAULT_SEED = 20260827
DEFAULT_CONFIG = PROJECT_ROOT / "configs" / "pilot_20260414.json"
DEFAULT_PRIOR_EVALUATION = (
    PROJECT_ROOT
    / "work"
    / "pilot_20260414"
    / "unblinded_qc_evaluation_20260826"
    / "unblinded_frame_evaluation.csv"
)
DEFAULT_OUTPUT = PROJECT_ROOT / "work" / "pilot_20260414" / "pixel_annotation_20260827"

REVIEW_COLUMNS = (
    "blind_id",
    "vessel_visible",
    "manual_x_center_px",
    "manual_z_upper_px",
    "annotation_confidence",
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
    "viterbi",
    "algorithm",
    "guard",
    "rnt",
    "auc",
    "trend",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(DEFAULT_CONFIG))
    parser.add_argument("--prior-evaluation", default=str(DEFAULT_PRIOR_EVALUATION))
    parser.add_argument("--tracking-root", default=None)
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT))
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--controls-per-scan", type=int, default=2)
    parser.add_argument("--overwrite", action="store_true")
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


def _text(value: Any) -> str:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return ""
    return str(value).strip().lower()


def _validate_prior_table(prior: pd.DataFrame, path: Path) -> pd.DataFrame:
    required = {
        "blind_id",
        "dataset_id",
        "frame_index_0based",
        "slow_axis_region",
        "human_x_correct",
        "human_upper_edge_correct",
        "human_roi_reasonable",
        "image_abnormal",
    }
    missing = required.difference(prior.columns)
    if missing:
        raise ValueError(f"{path}: missing previous-QC columns {sorted(missing)}")
    table = prior.copy()
    if len(table) != 81:
        raise ValueError(f"{path}: expected the audited 81-row previous QC table, got {len(table)}")
    if table["blind_id"].duplicated().any():
        raise ValueError(f"{path}: duplicate previous blind IDs")
    table["dataset_id"] = table["dataset_id"].astype(str)
    table["frame_index_0based"] = pd.to_numeric(table["frame_index_0based"], errors="raise").astype(int)
    table["slow_axis_region"] = table["slow_axis_region"].astype(str).str.strip().str.lower()
    for column in (
        "human_x_correct",
        "human_upper_edge_correct",
        "human_roi_reasonable",
        "image_abnormal",
    ):
        table[column] = table[column].map(_text)
    keys = table["dataset_id"].astype(str) + "|" + table["frame_index_0based"].astype(str)
    if keys.duplicated().any():
        raise ValueError(f"{path}: duplicate dataset/frame keys")
    if set(table["dataset_id"]) != set(REPRESENTATIVE_IDS):
        raise ValueError(f"{path}: representative scan set does not match the frozen nine scans")
    if not table["slow_axis_region"].isin(REGION_NAMES).all():
        raise ValueError(f"{path}: unknown slow-axis region")
    return table.sort_values(["dataset_id", "frame_index_0based"]).reset_index(drop=True)


def _is_disputed(row: pd.Series) -> bool:
    return row["human_x_correct"] in {"no", "uncertain"} or row["human_upper_edge_correct"] in {
        "no",
        "uncertain",
    }


def _is_yes_control(row: pd.Series) -> bool:
    return (
        row["human_x_correct"] == "yes"
        and row["human_upper_edge_correct"] == "yes"
        and row["human_roi_reasonable"] == "yes"
        and row["image_abnormal"] == "no"
    )


def _choose_row(pool: pd.DataFrame, rng: np.random.Generator) -> pd.Series:
    if pool.empty:
        raise ValueError("Cannot sample from an empty pool")
    ordered = pool.sort_values(["frame_index_0based", "blind_id"]).reset_index(drop=True)
    return ordered.iloc[int(rng.integers(0, len(ordered)))]


def _select_controls(
    prior: pd.DataFrame,
    disputed: pd.DataFrame,
    *,
    seed: int,
    controls_per_scan: int,
) -> pd.DataFrame:
    """Select fixed yes controls with an explicit region-coverage rule.

    Two controls are selected per scan.  If the disputed set lacks a slow-axis
    region for that scan, one control is forced into the missing region; the
    remaining control(s) are randomised among the other available regions.  If
    all regions are already represented, distinct regions are chosen randomly.
    This rule is frozen before the new annotation and does not inspect any new
    annotation result.
    """

    if controls_per_scan < 1:
        raise ValueError("controls_per_scan must be positive")
    disputed_keys = set(zip(disputed["dataset_id"], disputed["frame_index_0based"]))
    pool = prior.loc[prior.apply(_is_yes_control, axis=1)].copy()
    pool = pool.loc[
        ~pool.apply(lambda row: (row["dataset_id"], int(row["frame_index_0based"])) in disputed_keys, axis=1)
    ]

    selected: list[dict[str, Any]] = []
    for scan_id in REPRESENTATIVE_IDS:
        scan_pool = pool.loc[pool["dataset_id"] == scan_id].copy()
        if len(scan_pool) < controls_per_scan:
            raise ValueError(f"{scan_id}: fewer than {controls_per_scan} eligible yes controls")
        disputed_regions = set(disputed.loc[disputed["dataset_id"] == scan_id, "slow_axis_region"])
        available_regions = [region for region in REGION_NAMES if not scan_pool.loc[scan_pool["slow_axis_region"] == region].empty]
        missing_regions = [region for region in REGION_NAMES if region not in disputed_regions and region in available_regions]
        rng = _stable_rng(seed, f"yes_control|{scan_id}")
        chosen_keys: set[tuple[str, int]] = set()
        targets: list[str] = []
        # First cover a missing region, if there is one.  This handles the one
        # known gap in the previous disputed union without hand-picking frames.
        if missing_regions:
            missing_order = list(np.asarray(missing_regions)[rng.permutation(len(missing_regions))])
            targets.append(str(missing_order[0]))
        remaining_regions = [region for region in available_regions if region not in targets]
        if controls_per_scan > len(targets):
            perm = list(np.asarray(remaining_regions)[rng.permutation(len(remaining_regions))]) if remaining_regions else []
            targets.extend(str(region) for region in perm[: controls_per_scan - len(targets)])
        # If a scan has fewer distinct available regions than controls, sample
        # a second frame from any unused region using the same fixed RNG.
        while len(targets) < controls_per_scan:
            targets.append(available_regions[int(rng.integers(0, len(available_regions)))])

        for target_region in targets:
            region_pool = scan_pool.loc[
                (scan_pool["slow_axis_region"] == target_region)
                & ~scan_pool.apply(
                    lambda row: (row["dataset_id"], int(row["frame_index_0based"])) in chosen_keys,
                    axis=1,
                )
            ]
            fallback = "none"
            if region_pool.empty:
                region_pool = scan_pool.loc[
                    ~scan_pool.apply(
                        lambda row: (row["dataset_id"], int(row["frame_index_0based"])) in chosen_keys,
                        axis=1,
                    )
                ]
                fallback = "target_region_exhausted_used_any_yes_control"
            chosen = _choose_row(region_pool, rng)
            key = (str(chosen["dataset_id"]), int(chosen["frame_index_0based"]))
            chosen_keys.add(key)
            selected.append(
                {
                    "dataset_id": key[0],
                    "frame_index_0based": key[1],
                    "slow_axis_region": str(chosen["slow_axis_region"]),
                    "source_category": "prior_yes_control",
                    "control_sampling_fallback": fallback,
                    "control_sampling_seed": int(seed),
                }
            )
    result = pd.DataFrame(selected)
    if len(result) != len(REPRESENTATIVE_IDS) * controls_per_scan:
        raise AssertionError("Unexpected number of controls")
    if result.duplicated(["dataset_id", "frame_index_0based"]).any():
        raise AssertionError("A control frame was selected more than once")
    return result


def build_sampling_plan(prior: pd.DataFrame, *, seed: int, controls_per_scan: int) -> pd.DataFrame:
    disputed = prior.loc[prior.apply(_is_disputed, axis=1)].copy()
    disputed = disputed.assign(
        source_category="prior_disputed",
        control_sampling_fallback="not_applicable",
        control_sampling_seed="not_applicable",
    )
    controls = _select_controls(
        prior,
        disputed,
        seed=seed,
        controls_per_scan=controls_per_scan,
    )
    selected = pd.concat(
        [
            disputed[
                [
                    "dataset_id",
                    "frame_index_0based",
                    "slow_axis_region",
                    "source_category",
                    "control_sampling_fallback",
                    "control_sampling_seed",
                ]
            ],
            controls,
        ],
        ignore_index=True,
    )
    if selected.duplicated(["dataset_id", "frame_index_0based"]).any():
        raise AssertionError("Disputed and control selections overlap")
    order_rng = _stable_rng(seed, "pixel_annotation_global_blind_order")
    selected = selected.iloc[order_rng.permutation(len(selected))].reset_index(drop=True)
    selected.insert(0, "blind_id", [f"PA_{index:03d}" for index in range(1, len(selected) + 1)])
    selected["sampling_seed"] = int(seed)
    selected["selection_order"] = np.arange(1, len(selected) + 1, dtype=int)
    if len(selected) != 53:
        raise AssertionError(f"Expected 53 frames (35 disputed + 18 controls), got {len(selected)}")
    return selected


def _display_image(image: np.ndarray) -> tuple[np.ndarray, float, float]:
    finite = image[np.isfinite(image)]
    if finite.size == 0:
        raise ValueError("Flow image has no finite values")
    scale = max(float(np.percentile(finite, 75)), np.finfo(float).eps)
    transformed = np.arcsinh(np.maximum(image, 0) / scale)
    return transformed, float(np.percentile(transformed, 2)), float(np.percentile(transformed, 99.7))


def render_raw_image(image: np.ndarray, blind_id: str, output: Path) -> None:
    """Render a raw-Flow display with no localization or condition overlays."""

    shown, vmin, vmax = _display_image(image)
    z_size, x_size = shown.shape
    fig, ax = plt.subplots(figsize=(11.4, 7.4))
    fig.subplots_adjust(left=0.085, right=0.985, bottom=0.105, top=0.87)
    ax.imshow(shown, cmap="gray", vmin=vmin, vmax=vmax, origin="upper", aspect="auto")
    ax.set_xlim(0, x_size - 1)
    ax.set_ylim(z_size - 1, 0)
    ax.set_xlabel("lateral pixel")
    ax.set_ylabel("axial pixel")
    ax.set_title("Raw Flow image")
    fig.suptitle(blind_id, y=0.965, fontsize=14, fontweight="bold")
    fig.savefig(
        output,
        dpi=180,
        bbox_inches="tight",
        metadata={"Title": blind_id, "Description": "Condition-blinded raw Flow image for pixel annotation"},
    )
    plt.close(fig)


def reviewer_sheet(blind_ids: Iterable[str]) -> pd.DataFrame:
    rows = []
    for blind_id in blind_ids:
        row = {column: "" for column in REVIEW_COLUMNS}
        row["blind_id"] = str(blind_id)
        rows.append(row)
    return pd.DataFrame(rows, columns=REVIEW_COLUMNS)


def _mapping_row(
    plan_row: pd.Series,
    prior_row: pd.Series,
    tracking_row: pd.Series,
    manifest_row: pd.Series,
) -> dict[str, Any]:
    return {
        "blind_id": str(plan_row["blind_id"]),
        "image_filename": f"{plan_row['blind_id']}.png",
        "frame_data_filename": f"{plan_row['blind_id']}.npy",
        "source_category": str(plan_row["source_category"]),
        "prior_qc_blind_id": str(prior_row["blind_id"]),
        "dataset_id": str(plan_row["dataset_id"]),
        "original_scan_number": int(manifest_row["scan_number"]),
        "diameter_um": float(manifest_row["diameter_um"]),
        "velocity_mm_s": float(manifest_row["velocity_mm_s"]),
        "technical_replicate": int(manifest_row["technical_replicate"]),
        "frame_index_0based": int(plan_row["frame_index_0based"]),
        "frame_number_1based": int(plan_row["frame_index_0based"]) + 1,
        "slow_axis_region": str(plan_row["slow_axis_region"]),
        "sampling_seed": int(plan_row["sampling_seed"]),
        "selection_order": int(plan_row["selection_order"]),
        "control_sampling_fallback": str(plan_row["control_sampling_fallback"]),
        "control_sampling_seed": str(plan_row["control_sampling_seed"]),
        "previous_human_x_correct": str(prior_row["human_x_correct"]),
        "previous_human_upper_edge_correct": str(prior_row["human_upper_edge_correct"]),
        "previous_human_roi_reasonable": str(prior_row["human_roi_reasonable"]),
        "previous_image_abnormal": str(prior_row["image_abnormal"]),
        "previous_comments": str(prior_row.get("comments", "")),
        "old_tracking_class": str(tracking_row["old_tracking_class"]),
        "old_tracking_confidence": float(tracking_row["old_tracking_confidence"]),
        "x_path_confidence": float(tracking_row["x_path_confidence"]),
        "x_path_confidence_class": str(tracking_row["x_path_confidence_class"]),
        "z_edge_confidence": float(tracking_row["z_edge_confidence"]),
        "z_edge_confidence_class": str(tracking_row["z_edge_confidence_class"]),
        "overall_tracking_confidence": float(tracking_row["overall_tracking_confidence"]),
        "new_tracking_class": str(tracking_row["new_tracking_class"]),
        "x_center_px": float(tracking_row["x_center_px"]),
        "z_upper_px": float(tracking_row["z_upper_px"]),
        "z_lower_op_px": float(tracking_row["z_lower_op_px"]),
        "z_tail_start_px": float(tracking_row["z_tail_start_px"]),
        "central_x0_px": float(tracking_row["central_x0_px"]),
        "central_x1_exclusive_px": float(tracking_row["central_x1_exclusive_px"]),
        "axial_um_per_px": 6.7,
        "lateral_um_per_px": 12.7,
    }


def _write_reviewer_instructions(path: Path) -> None:
    text = """# 像素级人工定位说明

请只查看 `images/` 中的 PNG，或使用项目目录中的交互工具读取
`frame_data/` 中的匿名数组；完成后只填写 `pixel_level_manual_annotation.csv`。
在全部行完成并保存副本前，不要打开任何解码或内部映射文件。

每一行对应一张图，不能删除、重排或选择性跳过。

## 推荐交互流程

1. 启动 `pixel_annotation_tool.py`。
2. 按 `x` 后点击血管主体的横向中心；程序记录最近像素。
3. 按 `z` 后点击 Flow 主体的近端上缘；程序记录最近像素。
4. 按 `1`、`2` 或 `3` 分别填写 `vessel_visible=yes/uncertain/no`。
5. 按 `h`、`m` 或 `l` 填写 `annotation_confidence=high/medium/low`。
6. 按 `i` 循环填写 `image_abnormal`（空白 → yes → no → 空白）。
7. 在窗口底部“备注”输入框中输入文字，按输入框内的 `Enter` 提交备注；主图中再按 `Enter` 保存当前行并进入下一张；`p`/`n` 可前后移动；`u` 撤销最近一次操作，`c` 清空当前行坐标、状态和备注；`s` 保存；`q` 保存并退出。

如果血管不可见或无法可靠判断，请选择 `no` 或 `uncertain`，不要强行填写坐标。
如果选择 `yes`，请填写两个坐标。坐标均为显示图像的像素坐标：横向为 x，轴向为 z，原点在左上角。
`comments` 只记录必要的可评价性说明，不填写任何实验条件或内部编号。

标注目标是“血管是否可评价、人工认为的 x 中心和 Flow 近端上缘”，不是评价信号强弱，也不是重新勾画拖尾。
"""
    path.write_text(text, encoding="utf-8")


def _zip_reviewer_package(reviewer_dir: Path, archive: Path) -> None:
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as handle:
        for path in sorted(reviewer_dir.rglob("*")):
            if path.is_file():
                handle.write(path, path.relative_to(reviewer_dir.parent))


def self_check_package(
    package_root: Path, *, expected_count: int = 53, require_blank_sheet: bool = True
) -> dict[str, Any]:
    reviewer = package_root / "reviewer_package"
    restricted = package_root / "restricted_mapping"
    images = sorted((reviewer / "images").glob("*.png"))
    arrays = sorted((reviewer / "frame_data").glob("*.npy"))
    sheet = pd.read_csv(reviewer / "pixel_level_manual_annotation.csv", keep_default_na=False)
    mapping = pd.read_csv(restricted / "pixel_annotation_mapping.csv", keep_default_na=False)
    expected_ids = [f"PA_{index:03d}" for index in range(1, expected_count + 1)]
    if [path.stem for path in images] != expected_ids:
        raise AssertionError("Reviewer PNG filenames are incomplete or not blinded")
    if [path.stem for path in arrays] != expected_ids:
        raise AssertionError("Reviewer frame arrays are incomplete or not blinded")
    if sheet.columns.tolist() != list(REVIEW_COLUMNS):
        raise AssertionError("Reviewer sheet schema changed")
    if sheet["blind_id"].tolist() != expected_ids:
        raise AssertionError("Reviewer sheet IDs do not match images")
    if mapping["blind_id"].tolist() != expected_ids:
        raise AssertionError("Restricted mapping IDs do not match reviewer order")
    if require_blank_sheet and sheet[list(REVIEW_COLUMNS[1:])].astype(str).ne("").any().any():
        raise AssertionError("Reviewer sheet is not initially blank")
    if mapping["dataset_id"].nunique() != len(REPRESENTATIVE_IDS):
        raise AssertionError("Not all nine representative scans are present")
    per_scan = mapping.groupby("dataset_id").size()
    if (per_scan < 3).any():
        raise AssertionError("A representative scan has too few selected frames")
    if (mapping.groupby("dataset_id")["source_category"].apply(lambda s: (s == "prior_yes_control").sum()) != 2).any():
        raise AssertionError("Each representative scan must contribute exactly two yes controls")
    if mapping.duplicated(["dataset_id", "frame_index_0based"]).any():
        raise AssertionError("A frame was selected more than once")
    if mapping["source_category"].value_counts().to_dict() != {
        "prior_disputed": 35,
        "prior_yes_control": 18,
    }:
        raise AssertionError("Unexpected disputed/control composition")
    selected_regions = mapping.groupby(["dataset_id", "slow_axis_region"]).size()
    for scan_id in REPRESENTATIVE_IDS:
        if any((scan_id, region) not in selected_regions.index for region in REGION_NAMES):
            raise AssertionError(f"{scan_id}: front/middle/rear coverage failed")

    reviewer_names_and_text: list[str] = []
    for path in reviewer.rglob("*"):
        reviewer_names_and_text.append(str(path.relative_to(reviewer)).lower())
        if path.suffix.lower() in {".csv", ".md", ".txt"}:
            reviewer_names_and_text.append(path.read_text(encoding="utf-8").lower())
    joined = "\n".join(reviewer_names_and_text)
    leaked = [token for token in FORBIDDEN_REVIEWER_TOKENS if token in joined]
    if leaked:
        raise AssertionError(f"Reviewer package contains forbidden token(s): {leaked}")

    archive_path = package_root / "pixel_annotation_reviewer_package.zip"
    with zipfile.ZipFile(archive_path) as archive:
        archived_names = archive.namelist()
        if any("mapping" in name.lower() or "restricted" in name.lower() for name in archived_names):
            raise AssertionError("Restricted mapping leaked into reviewer ZIP")

    return {
        "status": "passed",
        "image_count": len(images),
        "frame_array_count": len(arrays),
        "sheet_rows": len(sheet),
        "mapping_rows": len(mapping),
        "scan_count": int(mapping["dataset_id"].nunique()),
        "source_category_counts": mapping["source_category"].value_counts().to_dict(),
        "per_scan_counts": {str(k): int(v) for k, v in per_scan.items()},
        "front_middle_rear_coverage": True,
        "reviewer_forbidden_token_hits": [],
        "reviewer_zip_excludes_mapping": True,
    }


def generate_package(args: argparse.Namespace) -> Path:
    cfg, config_path = load_config(args.config)
    manifest = load_manifest(cfg, config_path).set_index("dataset_id", drop=False)
    if set(REPRESENTATIVE_IDS).difference(manifest.index):
        raise ValueError("Manifest lacks one or more frozen representative scans")
    prior_path = Path(args.prior_evaluation).resolve()
    prior = _validate_prior_table(pd.read_csv(prior_path, keep_default_na=False), prior_path)
    plan = build_sampling_plan(prior, seed=args.seed, controls_per_scan=args.controls_per_scan)

    tracking_root = Path(args.tracking_root).resolve() if args.tracking_root else output_root(cfg) / "tracking"
    final_root = Path(args.output_dir).resolve()
    if final_root.exists() and not args.overwrite:
        raise FileExistsError(f"Output already exists; use --overwrite to replace: {final_root}")

    tables: dict[str, pd.DataFrame] = {}
    tracking_paths: dict[str, Path] = {}
    for scan_id in REPRESENTATIVE_IDS:
        path = tracking_csv_path(tracking_root, scan_id, PRIMARY_ALPHA)
        tables[scan_id] = _LEGACY.normalize_tracking_table(path, require_new_confidence=True)
        tracking_paths[scan_id] = path

    prior_lookup = prior.set_index(["dataset_id", "frame_index_0based"], drop=False)
    building = final_root.parent / f".{final_root.name}.building-{uuid.uuid4().hex}"
    reviewer_dir = building / "reviewer_package"
    image_dir = reviewer_dir / "images"
    array_dir = reviewer_dir / "frame_data"
    restricted_dir = building / "restricted_mapping"
    image_dir.mkdir(parents=True, exist_ok=False)
    array_dir.mkdir(parents=True, exist_ok=False)
    restricted_dir.mkdir(parents=True, exist_ok=False)

    mapping_rows: list[dict[str, Any]] = []
    try:
        for scan_id in REPRESENTATIVE_IDS:
            scan_plan = plan.loc[plan["dataset_id"] == scan_id].copy()
            requested_frames = scan_plan["frame_index_0based"].astype(int).tolist()
            frames = _LEGACY.read_selected_frames(export_path(cfg, scan_id), requested_frames)
            table = tables[scan_id].set_index("frame_index", drop=False)
            manifest_row = manifest.loc[scan_id]
            for _, plan_row in scan_plan.iterrows():
                frame = int(plan_row["frame_index_0based"])
                prior_row = prior_lookup.loc[(scan_id, frame)]
                tracking_row = table.loc[frame]
                blind_id = str(plan_row["blind_id"])
                image = np.asarray(frames[frame], dtype=np.float32)
                np.save(array_dir / f"{blind_id}.npy", image, allow_pickle=False)
                render_raw_image(image, blind_id, image_dir / f"{blind_id}.png")
                mapping_rows.append(_mapping_row(plan_row, prior_row, tracking_row, manifest_row))

        mapping = pd.DataFrame(mapping_rows).sort_values("blind_id").reset_index(drop=True)
        mapping.to_csv(restricted_dir / "pixel_annotation_mapping.csv", index=False, encoding="utf-8-sig")
        reviewer_sheet(mapping["blind_id"]).to_csv(
            reviewer_dir / "pixel_level_manual_annotation.csv",
            index=False,
            encoding="utf-8-sig",
            quoting=csv.QUOTE_MINIMAL,
        )
        _write_reviewer_instructions(reviewer_dir / "README_REVIEWER.md")

        run_metadata: dict[str, Any] = {
            "schema_version": "1.0.0",
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "purpose": "condition-blinded pixel-level manual localization annotation",
            "sampling_seed": int(args.seed),
            "controls_per_scan": int(args.controls_per_scan),
            "selection_rule": (
                "all previous rows with x or Flow-upper human score no/uncertain, plus "
                "two previous x/upper/ROI-yes and non-abnormal controls per scan; "
                "missing slow-axis regions are covered by a deterministic region rule"
            ),
            "disputed_frame_count": int((mapping["source_category"] == "prior_disputed").sum()),
            "control_frame_count": int((mapping["source_category"] == "prior_yes_control").sum()),
            "total_frame_count": int(len(mapping)),
            "representative_scan_count": len(REPRESENTATIVE_IDS),
            "representative_scan_ids_restricted": list(REPRESENTATIVE_IDS),
            "display": {
                "transform": "arcsinh(max(raw_flow,0)/per_image_p75)",
                "clip_percentiles": [2.0, 99.7],
                "axes": "lateral pixel x axial pixel, origin upper-left",
                "algorithm_overlay": False,
                "condition_labels": False,
            },
            "coordinate_calibration": {"axial_um_per_px": 6.7, "lateral_um_per_px": 12.7},
            "prior_evaluation_csv": str(prior_path),
            "prior_evaluation_sha256": _sha256(prior_path),
            "tracking_csv_sha256": {scan_id: _sha256(path) for scan_id, path in tracking_paths.items()},
            "config_path": str(config_path),
            "config_sha256": _sha256(config_path),
            "generator_path": str(Path(__file__).resolve()),
            "generator_sha256": _sha256(Path(__file__).resolve()),
            "reviewer_package_contains_mapping": False,
            "next_state": "waiting for completed pixel_level_manual_annotation.csv",
        }
        (restricted_dir / "pixel_annotation_run_metadata.json").write_text(
            json.dumps(run_metadata, ensure_ascii=False, indent=2), encoding="utf-8"
        )
        _zip_reviewer_package(reviewer_dir, building / "pixel_annotation_reviewer_package.zip")
        preliminary = self_check_package(building, expected_count=len(mapping))
        (restricted_dir / "pixel_annotation_self_check.json").write_text(
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
    check = self_check_package(output)
    print(json.dumps({"output": str(output), "self_check": check}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
