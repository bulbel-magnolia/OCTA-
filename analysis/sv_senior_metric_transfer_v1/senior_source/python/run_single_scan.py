#!/usr/bin/env python3
"""Run one raw OCT .oct file through the frozen B-scan tail pipeline.

This is an orchestration wrapper.  The localization, feature, and metric
implementations are loaded from the copied audited source files; this wrapper
does not retune or reimplement those algorithms.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any

import numpy as np
import pandas as pd


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
PYTHON_ROOT = PACKAGE_ROOT / "python"
MATLAB_ROOT = PACKAGE_ROOT / "matlab"
DEFAULT_CONFIG = PACKAGE_ROOT / "config" / "example_config.json"
FRAMES = 500
N_Z = 351
N_X = 500
PRIMARY_ALPHA = 0.15
REGIONS = (("front", 0, 167), ("middle", 167, 333), ("rear", 333, 500))

if str(PYTHON_ROOT) not in sys.path:
    sys.path.insert(0, str(PYTHON_ROOT))

from tailq.io import load_flow, tracking_csv_path  # noqa: E402


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"无法加载交付包源码：{path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _region(frame_index: int) -> str:
    for name, start, stop in REGIONS:
        if start <= int(frame_index) < stop:
            return name
    raise ValueError(f"帧号不在 0..499：{frame_index}")


def _matlab_literal(value: Path | str) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def _resolve_matlab(explicit: str | None) -> Path:
    candidates: list[str] = []
    if explicit:
        candidates.append(explicit)
    env_value = os.environ.get("MATLAB_EXE")
    if env_value:
        candidates.append(env_value)
    candidates.append(r"C:\matlab\bin\matlab.exe")
    on_path = shutil.which("matlab")
    if on_path:
        candidates.append(on_path)
    seen: set[str] = set()
    for candidate in candidates:
        if candidate in seen:
            continue
        seen.add(candidate)
        path = Path(candidate)
        if path.is_file():
            return path.resolve()
    raise FileNotFoundError(
        "找不到 MATLAB 可执行文件。请使用 --matlab \"C:\\matlab\\bin\\matlab.exe\"，"
        "或设置环境变量 MATLAB_EXE。"
    )


def _read_config(path: Path) -> dict[str, Any]:
    config = json.loads(path.read_text(encoding="utf-8"))
    required = {"study_id", "source", "output", "scales", "metrics", "aggregation", "tracking"}
    missing = sorted(required.difference(config))
    if missing:
        raise ValueError(f"配置缺少 section：{missing}")
    return config


def _prepare_output(output: Path) -> None:
    if output.exists() and (output.is_file() or any(output.iterdir())):
        raise FileExistsError(f"输出目录已存在且非空，为避免覆盖请换一个目录：{output}")
    output.mkdir(parents=True, exist_ok=True)


def _run_matlab(raw: Path, flow_path: Path, log_path: Path, scan_id: str, matlab_exe: Path) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    expression = (
        "run_single_scan_export("
        f"{_matlab_literal(raw)},{_matlab_literal(flow_path)},{_matlab_literal(log_path)},"
        f"{_matlab_literal(scan_id)});"
    )
    command = [str(matlab_exe), "-batch", expression]
    with log_path.open("a", encoding="utf-8") as log_handle:
        log_handle.write(f"COMMAND: {command}\n")
        log_handle.flush()
        completed = subprocess.run(
            command,
            cwd=MATLAB_ROOT,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
    if completed.returncode != 0:
        raise RuntimeError(f"MATLAB 导出失败，返回码={completed.returncode}；请查看 {log_path}")


def _run_tracking(flow_path: Path, tracking_dir: Path, scan_id: str, diameter_um: float, config_path: Path) -> Path:
    tracker = PYTHON_ROOT / "02_track_vessel" / "track_vessel.py"
    scan_tracking_dir = tracking_dir / scan_id
    environment = os.environ.copy()
    environment["MPLBACKEND"] = "Agg"
    existing_pythonpath = environment.get("PYTHONPATH", "")
    environment["PYTHONPATH"] = str(PYTHON_ROOT) + (os.pathsep + existing_pythonpath if existing_pythonpath else "")
    command = [
        sys.executable,
        str(tracker),
        "--flow",
        str(flow_path),
        "--output-dir",
        str(scan_tracking_dir),
        "--scan-id",
        scan_id,
        "--diameter-um",
        str(float(diameter_um)),
        "--config",
        str(config_path),
        "--expected-frames",
        str(FRAMES),
    ]
    print("运行 Python 跟踪：", " ".join(f'\"{item}\"' if " " in item else item for item in command), flush=True)
    completed = subprocess.run(command, cwd=tracker.parent, env=environment, check=False)
    if completed.returncode != 0:
        raise RuntimeError(f"Python 跟踪失败，返回码={completed.returncode}；请查看 {scan_tracking_dir}")
    return tracking_csv_path(tracking_dir, scan_id, PRIMARY_ALPHA)


def _make_features(
    flow: np.ndarray,
    source_frames: np.ndarray,
    tracking_path: Path,
    scan_id: str,
    diameter_um: float,
    velocity_mm_s: float,
    dev: Any,
    output_path: Path,
) -> pd.DataFrame:
    tracking = pd.read_csv(tracking_path)
    if len(tracking) != FRAMES or tracking["frame_index"].duplicated().any():
        raise ValueError(f"{scan_id}: tracking 必须包含 500 个不重复帧")
    tracking = tracking.sort_values("frame_index").set_index("frame_index", drop=False)
    if set(tracking.index.astype(int)) != set(range(FRAMES)):
        raise ValueError(f"{scan_id}: tracking 帧号必须为 0..499")
    if flow.shape != (FRAMES, N_Z, N_X):
        raise ValueError(f"{scan_id}: 期望 flow 形状 {(FRAMES, N_Z, N_X)}，实际为 {flow.shape}")
    if set(np.asarray(source_frames, dtype=int).tolist()) != set(range(FRAMES)):
        raise ValueError(f"{scan_id}: source frame mapping 必须为 0..499")

    rows: list[dict[str, Any]] = []
    for frame_index in range(FRAMES):
        track = tracking.loc[frame_index]
        feature = dev._local_body_features(flow[frame_index], track)
        feature.update(
            {
                "dataset_id": scan_id,
                "scan_order": 1,
                "frame_index_0based": int(frame_index),
                "frame_index_1based": int(frame_index + 1),
                "diameter_um": float(diameter_um),
                "velocity_mm_s": float(velocity_mm_s),
                "alpha": PRIMARY_ALPHA,
                "slow_axis_region": _region(frame_index),
                "source_frame_index": int(source_frames[frame_index]),
            }
        )
        rows.append(feature)
    sequence = pd.DataFrame(rows)
    neighbor = sequence["local_body_peak_cnr"].shift(1).combine(
        sequence["local_body_peak_cnr"].shift(-1),
        lambda left, right: np.nanmedian([left, right]),
    )
    sequence["neighbor_peak_cnr"] = neighbor.fillna(sequence["local_body_peak_cnr"])
    sequence = dev._add_assessability_score(sequence)
    corrected, changed = dev._isolated_jump_correct(
        sequence["x2_robust_centroid_px"].to_numpy(float)
    )
    sequence["x4_centroid_isolated_jump_corrected_px"] = corrected
    sequence["x4_jump_corrected"] = changed
    output_path.parent.mkdir(parents=True, exist_ok=True)
    sequence.to_csv(output_path, index=False, encoding="utf-8-sig")
    metadata = {
        "status": "success",
        "scan_id": scan_id,
        "rows": int(len(sequence)),
        "primary_alpha": PRIMARY_ALPHA,
        "x_method": "X4_centroid_isolated_jump_corrected",
        "assessability_rule": {
            "assessable": ">=0.60",
            "uncertain": ">=0.40 and <0.60",
            "not_assessable": "<0.40",
        },
        "source_module": "10_xroi_refinement/develop_xroi_assessability.py",
    }
    output_path.with_suffix(output_path.suffix + ".metadata.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    return sequence


def _make_metrics(
    flow: np.ndarray,
    source_frames: np.ndarray,
    features: pd.DataFrame,
    tracking: pd.DataFrame,
    scan_id: str,
    diameter_um: float,
    velocity_mm_s: float,
    config: dict[str, Any],
    freeze: Any,
    qc_dir: Path,
) -> dict[str, Any]:
    metrics_cfg = dict(config["metrics"])
    metrics_cfg["central_fraction"] = float(freeze.PRIMARY_ROI_FRACTION)
    metrics_cfg["primary_guard_px"] = int(freeze.PRIMARY_GUARD)
    metrics_cfg["guards_px"] = list(freeze.GUARDS)
    metrics_cfg["primary_window_um"] = int(freeze.PRIMARY_WINDOW_UM[1])
    frame_rows: list[dict[str, Any]] = []
    tracking_by_frame = tracking.set_index("frame_index", drop=False)
    for frame_index in range(FRAMES):
        feature = features.iloc[frame_index].to_dict()
        track = tracking_by_frame.loc[frame_index]
        row: dict[str, Any] = dict(feature)
        for key, value in track.to_dict().items():
            if key in {"frame_index", "scan_id", "diameter_um", "alpha"}:
                continue
            row[f"tracking_{key}"] = value
        x_x4 = float(feature.get("x4_centroid_isolated_jump_corrected_px", np.nan))
        z_upper = float(track.get("z_upper_px", np.nan))
        metric = freeze._frame_metrics(
            flow[frame_index],
            x_x4=x_x4,
            z_upper=z_upper,
            diameter_um=float(diameter_um),
            metrics_cfg=metrics_cfg,
        )
        row.update(
            {
                "scan_id": scan_id,
                "dataset_id": scan_id,
                "frame_index_0based": int(frame_index),
                "frame_index_1based": int(frame_index + 1),
                "source_frame_index": int(source_frames[frame_index]),
                "diameter_um": float(diameter_um),
                "velocity_mm_s": float(velocity_mm_s),
                "alpha": float(freeze.PRIMARY_ALPHA),
                "roi_central_fraction": float(freeze.PRIMARY_ROI_FRACTION),
                "primary_guard_px": int(freeze.PRIMARY_GUARD),
                "raw_auc_window_start_um": int(freeze.PRIMARY_WINDOW_UM[0]),
                "raw_auc_window_stop_um": int(freeze.PRIMARY_WINDOW_UM[1]),
                "axial_um_per_px": float(freeze.AXIAL_UM_PER_PX),
                "lateral_um_per_px": float(freeze.LATERAL_UM_PER_PX),
            }
        )
        row.update(metric)
        row["primary_frame_included"] = bool(
            row["assessability_class"] == "assessable" and row["raw_auc_0_200_guard2_valid"]
        )
        row["sensitivity_frame_included"] = bool(
            row["assessability_class"] in ("assessable", "uncertain") and row["raw_auc_0_200_guard2_valid"]
        )
        row["primary_exclusion_reason"] = (
            "included"
            if row["primary_frame_included"]
            else ("not_assessable_or_uncertain" if row["assessability_class"] != "assessable" else row["invalid_reason"])
        )
        row["sensitivity_exclusion_reason"] = (
            "included"
            if row["sensitivity_frame_included"]
            else ("not_assessable" if row["assessability_class"] == "not_assessable" else row["invalid_reason"])
        )
        frame_rows.append(row)

    frame_table = pd.DataFrame(frame_rows).sort_values("frame_index_0based").reset_index(drop=True)
    qc_dir.mkdir(parents=True, exist_ok=True)
    frame_path = qc_dir / f"{scan_id}_frame_level_metrics.csv"
    frame_table.to_csv(frame_path, index=False, encoding="utf-8-sig")
    scan_row = freeze._scan_row(frame_table, scan_id, float(diameter_um), float(velocity_mm_s), 1)
    scan_table = pd.DataFrame([scan_row])
    scan_table.to_csv(qc_dir / f"{scan_id}_scan_level_primary.csv", index=False, encoding="utf-8-sig")
    assess_sens = freeze._assessability_sensitivity(frame_table)
    assess_sens.to_csv(qc_dir / f"{scan_id}_assessability_sensitivity.csv", index=False, encoding="utf-8-sig")
    guard_sens = freeze._guard_sensitivity(frame_table)
    guard_sens.to_csv(qc_dir / f"{scan_id}_guard_sensitivity.csv", index=False, encoding="utf-8-sig")
    spatial = pd.DataFrame(freeze._spatial_rows(frame_table, scan_id, float(diameter_um), float(velocity_mm_s)))
    spatial.to_csv(qc_dir / f"{scan_id}_spatial_front_middle_rear.csv", index=False, encoding="utf-8-sig")
    p95_qc = freeze._p95_qc(frame_table)
    p95_qc.to_csv(qc_dir / f"{scan_id}_p95_denominator_qc.csv", index=False, encoding="utf-8-sig")

    shape_pass = tuple(flow.shape) == (FRAMES, N_Z, N_X)
    finite_pass = bool(np.isfinite(flow).all())
    tracking_pass = len(tracking) == FRAMES and set(tracking["frame_index"].astype(int)) == set(range(FRAMES))
    valid_rate = float(scan_row["primary_valid_fraction_of_assessable"])
    rate_pass = bool(np.isfinite(valid_rate) and valid_rate >= float(getattr(freeze, "VALID_FRACTION_FLOOR", 0.5)))
    single_status = "PASS_SINGLE_SCAN_QC" if all((shape_pass, finite_pass, tracking_pass, rate_pass)) else "CHECK_SINGLE_SCAN_QC"
    qc_summary = dict(scan_row)
    qc_summary.update(
        {
            "flow_shape": "x".join(str(int(x)) for x in flow.shape),
            "flow_finite": finite_pass,
            "tracking_alpha015_complete": tracking_pass,
            "primary_valid_fraction_threshold": float(getattr(freeze, "VALID_FRACTION_FLOOR", 0.5)),
            "status": single_status,
            "global_formal15_gate": "NOT_ASSESSED_SINGLE_SCAN",
        }
    )
    pd.DataFrame([qc_summary]).to_csv(qc_dir / f"{scan_id}_single_scan_qc_summary.csv", index=False, encoding="utf-8-sig")
    gate_rows = [
        {"criterion": "Flow形状", "observed": str(tuple(flow.shape)), "threshold": "(500,351,500)", "status": "PASS" if shape_pass else "FAIL"},
        {"criterion": "Flow有限值", "observed": str(finite_pass), "threshold": "全部有限", "status": "PASS" if finite_pass else "FAIL"},
        {"criterion": "alpha=0.15跟踪帧数", "observed": int(len(tracking)), "threshold": "500帧且帧号0..499", "status": "PASS" if tracking_pass else "FAIL"},
        {"criterion": "主指标有效比例", "observed": valid_rate, "threshold": f">={getattr(freeze, 'VALID_FRACTION_FLOOR', 0.5):.2f}", "status": "PASS" if rate_pass else "FAIL"},
        {"criterion": "Formal15全局Gate A", "observed": "单卷未评估", "threshold": "需要完整冻结批次", "status": "NOT_ASSESSED"},
    ]
    pd.DataFrame(gate_rows).to_csv(qc_dir / f"{scan_id}_single_scan_gate_status.csv", index=False, encoding="utf-8-sig")
    return {
        "frame_table": frame_table,
        "scan_row": scan_row,
        "single_status": single_status,
        "paths": {
            "frame_metrics": str(frame_path),
            "scan_primary": str(qc_dir / f"{scan_id}_scan_level_primary.csv"),
            "qc_summary": str(qc_dir / f"{scan_id}_single_scan_qc_summary.csv"),
        },
    }


def run(args: argparse.Namespace) -> dict[str, Any]:
    raw = args.input_oct.expanduser().resolve()
    if not raw.is_file():
        raise FileNotFoundError(f"找不到原始 .oct 文件：{raw}")
    config_path = args.config.expanduser().resolve()
    config = _read_config(config_path)
    output = args.output_dir.expanduser().resolve()
    _prepare_output(output)
    scan_id = args.scan_id or raw.stem
    if not scan_id.replace("_", "").replace("-", "").isalnum():
        raise ValueError("scan-id 只允许字母、数字、下划线和短横线")
    velocity = float(args.velocity_mm_s) if args.velocity_mm_s is not None else float("nan")
    matlab_exe = _resolve_matlab(args.matlab)
    omag_dir = output / "01_omag"
    flow_path = omag_dir / "p_bld_ed.mat"
    export_log = omag_dir / "export.log"
    _run_matlab(raw, flow_path, export_log, scan_id, matlab_exe)
    if not flow_path.is_file():
        raise FileNotFoundError(f"MATLAB 返回成功但未找到输出：{flow_path}")

    flow, source_frames, flow_meta = load_flow(flow_path)
    if tuple(flow.shape) != (FRAMES, N_Z, N_X):
        raise ValueError(f"导出 Flow 形状错误：{flow.shape}，期望 {(FRAMES, N_Z, N_X)}")
    if not np.isfinite(flow).all():
        raise ValueError("导出 Flow 含非有限值")
    tracking_path = _run_tracking(flow_path, output / "02_tracking", scan_id, float(args.diameter_um), config_path)
    tracking = pd.read_csv(tracking_path)
    dev = _load_module("delivery_xroi_runtime", PYTHON_ROOT / "10_xroi_refinement" / "develop_xroi_assessability.py")
    freeze = _load_module("delivery_freeze_runtime", PYTHON_ROOT / "13_formal_9scan" / "formal_9scan_freeze.py")
    features_path = output / "03_features" / f"{scan_id}_features.csv"
    features = _make_features(
        flow,
        source_frames,
        tracking_path,
        scan_id,
        float(args.diameter_um),
        velocity,
        dev,
        features_path,
    )
    metric_result = _make_metrics(
        flow,
        source_frames,
        features,
        tracking,
        scan_id,
        float(args.diameter_um),
        velocity,
        config,
        freeze,
        output / "04_qc",
    )
    summary = {
        "status": "success",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "package": "Bscan_tail_auto_quantification_for_junior_20260903",
        "scan_id": scan_id,
        "input_oct": str(raw),
        "input_oct_size_bytes": int(raw.stat().st_size),
        "diameter_um": float(args.diameter_um),
        "velocity_mm_s": velocity,
        "matlab_executable": str(matlab_exe),
        "flow": {
            "path": str(flow_path),
            "shape_python_frame_z_x": list(flow.shape),
            "dtype": str(flow.dtype),
            "finite": bool(np.isfinite(flow).all()),
            "source_frame_mapping": [int(source_frames[0]), int(source_frames[-1])],
            "metadata_keys": sorted(flow_meta.keys()),
        },
        "tracking": {"path": str(tracking_path), "rows": int(len(tracking)), "primary_alpha": PRIMARY_ALPHA},
        "features": {"path": str(features_path), "rows": int(len(features)), "x_method": "X4_centroid_isolated_jump_corrected"},
        "qc": metric_result["paths"],
        "single_scan_qc_status": metric_result["single_status"],
        "formal_global_gate": "NOT_ASSESSED_SINGLE_SCAN",
        "frozen_parameters": {
            "central_roi_fraction": float(freeze.PRIMARY_ROI_FRACTION),
            "primary_guard_px": int(freeze.PRIMARY_GUARD),
            "guard_sensitivity_px": list(freeze.GUARDS),
            "raw_auc_window_um": list(freeze.PRIMARY_WINDOW_UM),
            "assessability": "primary=assessable only; sensitivity=assessable+uncertain; not_assessable excluded",
            "aggregation": "pooled median of valid assessable frames; no block medians",
        },
        "delivery_code_sha256": {
            "python/run_single_scan.py": _sha256(PYTHON_ROOT / "run_single_scan.py"),
            "python/13_formal_9scan/formal_9scan_freeze.py": _sha256(PYTHON_ROOT / "13_formal_9scan" / "formal_9scan_freeze.py"),
            "python/10_xroi_refinement/develop_xroi_assessability.py": _sha256(PYTHON_ROOT / "10_xroi_refinement" / "develop_xroi_assessability.py"),
            "matlab/export_p_bld_ed_volume.m": _sha256(MATLAB_ROOT / "export_p_bld_ed_volume.m"),
        },
    }
    (output / "run_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=True), encoding="utf-8")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-oct", required=True, type=Path, help="原始 OCT .oct 文件；只读")
    parser.add_argument("--diameter-um", required=True, type=float, help="血管仿真内径，单位 um")
    parser.add_argument("--velocity-mm-s", type=float, help="速度标签，单位 mm/s；不填则输出 NaN")
    parser.add_argument("--scan-id", help="输出标识；默认使用 .oct 文件名")
    parser.add_argument("--output-dir", required=True, type=Path, help="新建的输出目录，必须为空或不存在")
    parser.add_argument("--matlab", help="matlab.exe 路径；也可设置 MATLAB_EXE")
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG, help="交付包配置 JSON")
    args = parser.parse_args()
    summary = run(args)
    print(json.dumps({"status": summary["status"], "output_dir": str(args.output_dir.resolve()), "scan_id": summary["scan_id"]}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
