#!/usr/bin/env python3
"""Task 2: X4-centred fixed-128-um lateral-width sensitivity audit.

The frozen run001 remains the primary analysis. This script downloads the released
2-D numerical arrays, verifies package and per-frame SHA-256 hashes, replays the
frozen X4+X1 geometry as an input/implementation check, and then changes only the
lateral source/tail width to 128 um while retaining X4, z_top, vertical diameter,
500 um tail length, guard=0, background skip-3/take-5, linear SV, signed residuals,
and all other frozen quantification rules.

This is geometric QC. It is not a new scientific endpoint and no p-values are
computed.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

EXPECTED_VALID = 2422
EXPECTED_PACKAGES = 25
RELEASE_TAG = "formal-sv-d128-v21-run001"
RELEASE_BASE = f"https://github.com/bulbel-magnolia/OCTA-/releases/download/{RELEASE_TAG}"
FIXED_WIDTH_UM = 128.0


def first_present(columns: Iterable[str], candidates: Sequence[str], *, required: bool = True) -> str | None:
    cols = set(columns)
    for name in candidates:
        if name in cols:
            return name
    if required:
        raise KeyError(f"Required column absent; tried {candidates}")
    return None


def as_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False).astype(bool)
    if pd.api.types.is_numeric_dtype(series):
        return series.fillna(0).astype(float).ne(0)
    x = series.astype(str).str.strip().str.lower()
    true = {"true", "1", "yes", "y", "t"}
    false = {"false", "0", "no", "n", "f", "", "nan", "none", "na", "<na>"}
    bad = ~x.isin(true | false)
    if bad.any():
        raise ValueError(f"Unrecognized boolean values: {sorted(x[bad].unique())[:10]}")
    return x.isin(true)


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def safe_divide(a, b) -> float:
    try:
        aa, bb = float(a), float(b)
    except (TypeError, ValueError):
        return float("nan")
    if not np.isfinite(aa) or not np.isfinite(bb) or bb == 0.0:
        return float("nan")
    return aa / bb


def pct_delta(new, old) -> float:
    return 100.0 * safe_divide(float(new) - float(old), old)


def abs_rel_error(new: float, old: float) -> float:
    if not np.isfinite(new) or not np.isfinite(old):
        return float("nan")
    scale = max(abs(old), 1.0)
    return abs(new - old) / scale


def quantiles(values: Iterable[float]) -> dict[str, float | int]:
    a = np.asarray([float(x) for x in values if np.isfinite(float(x))], dtype=float)
    if a.size == 0:
        return {"n": 0, "mean": np.nan, "sd": np.nan, "min": np.nan, "q1": np.nan,
                "median": np.nan, "q3": np.nan, "iqr": np.nan, "p95": np.nan, "max": np.nan}
    q1, med, q3, p95 = np.quantile(a, [0.25, 0.5, 0.75, 0.95])
    return {"n": int(a.size), "mean": float(np.mean(a)),
            "sd": float(np.std(a, ddof=1)) if a.size > 1 else np.nan,
            "min": float(np.min(a)), "q1": float(q1), "median": float(med),
            "q3": float(q3), "iqr": float(q3-q1), "p95": float(p95), "max": float(np.max(a))}


def parse_package_name(name: str) -> tuple[str, int, int]:
    m = re.search(r"_(flow\d+)_(\d{3})_(\d{3})\.zip$", name)
    if not m:
        raise ValueError(f"Cannot parse package name: {name}")
    return m.group(1), int(m.group(2)), int(m.group(3))


def download_file(url: str, destination: Path) -> None:
    command = [
        "curl", "--location", "--fail", "--retry", "4", "--retry-all-errors",
        "--connect-timeout", "30", "--silent", "--show-error",
        "--output", str(destination), url,
    ]
    subprocess.run(command, check=True)


def prepare_inputs(root: Path, formal_dir: Path, task1_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    frame = pd.read_csv(formal_dir / "frame_results.csv", low_memory=False)
    task1 = pd.read_csv(task1_dir / "relative_tail_intensity_framewise.csv", low_memory=False)

    scan_col = first_present(frame.columns, ["scan_id", "scan"])
    frame_col = first_present(frame.columns, ["frame_index_0based", "frame_index", "bscan_index"])
    valid_col = first_present(frame.columns, ["valid"])
    x4_col = first_present(frame.columns, ["x4_centroid_isolated_jump_corrected_px"])
    xleft_col = first_present(frame.columns, ["x_left_edge_px"])
    xright_col = first_present(frame.columns, ["x_right_edge_px"])
    ztop_col = first_present(frame.columns, ["z_top_edge_px"])
    width_col = first_present(frame.columns, ["lateral_width_um"])
    zupper_col = first_present(frame.columns, ["z_upper_px"], required=False)

    geom = pd.DataFrame({
        "scan_id": frame[scan_col].astype(str),
        "frame_index": pd.to_numeric(frame[frame_col], errors="raise").astype(int),
        "formal_valid": as_bool(frame[valid_col]),
        "x4_px": pd.to_numeric(frame[x4_col], errors="coerce"),
        "baseline_x_left_edge_px": pd.to_numeric(frame[xleft_col], errors="coerce"),
        "baseline_x_right_edge_px": pd.to_numeric(frame[xright_col], errors="coerce"),
        "z_top_edge_px": pd.to_numeric(frame[ztop_col], errors="coerce"),
        "baseline_lateral_width_um": pd.to_numeric(frame[width_col], errors="coerce"),
    })
    if zupper_col is not None:
        geom["z_upper_px"] = pd.to_numeric(frame[zupper_col], errors="coerce")

    geom = geom.loc[geom.formal_valid].copy()
    if len(geom) != EXPECTED_VALID:
        raise AssertionError(f"Expected {EXPECTED_VALID} frozen-valid geometry rows, found {len(geom)}")
    if geom.duplicated(["scan_id", "frame_index"]).any():
        raise AssertionError("Duplicate frozen geometry keys")

    task1["scan_id"] = task1["scan_id"].astype(str)
    frame_task_col = first_present(task1.columns, ["frame_index", "frame_index_0based"])
    task1["frame_index"] = pd.to_numeric(task1[frame_task_col], errors="raise").astype(int)
    if len(task1) != EXPECTED_VALID:
        raise AssertionError(f"Task1 framewise row count mismatch: {len(task1)}")

    merged = task1.merge(geom.drop(columns=["formal_valid"]), on=["scan_id", "frame_index"], how="inner", validate="one_to_one")
    if len(merged) != EXPECTED_VALID:
        raise AssertionError(f"Task1/frozen geometry key mismatch: {len(merged)}")

    centre_from_edges = (merged["baseline_x_left_edge_px"] + merged["baseline_x_right_edge_px"]) / 2.0
    centre_error = np.abs(centre_from_edges.to_numpy(float) - merged["x4_px"].to_numpy(float))
    if not np.all(np.isfinite(centre_error)) or float(np.max(centre_error)) > 1e-10:
        raise AssertionError(f"Frozen source geometry is not exactly recentered on X4; max error {np.nanmax(centre_error)}")

    config = json.loads((formal_dir / "run_config.json").read_text(encoding="utf-8"))
    cal = config["calibration"]
    geo = config["geometry"]
    bg = config["background"]
    frozen_expected = {
        "diameter_um": 128.0, "dx_um": 12.7, "dz_um": 6.7,
        "tail_gap_um": 0.0, "tail_length_um": 500.0,
        "ellipse_supersample": 16, "skip_columns": 3, "strip_width_columns": 5,
    }
    observed = {
        "diameter_um": float(cal["diameter_um"]), "dx_um": float(cal["dx_um"]), "dz_um": float(cal["dz_um"]),
        "tail_gap_um": float(geo["tail_gap_um"]), "tail_length_um": float(geo["tail_length_um"]),
        "ellipse_supersample": int(geo["ellipse_supersample"]),
        "skip_columns": int(bg["skip_columns"]), "strip_width_columns": int(bg["strip_width_columns"]),
    }
    if observed != frozen_expected:
        raise AssertionError(f"Frozen configuration differs from Task2 specification: {observed}")
    return merged, frame, observed


def summarize_changes(framewise: pd.DataFrame) -> pd.DataFrame:
    rows = []
    scopes = [("overall", "all", framewise)]
    for scan in sorted(framewise.scan_id.unique(), key=lambda x: float(framewise.loc[framewise.scan_id.eq(x), "flow_mm_s"].iloc[0])):
        scopes.append(("scan", scan, framewise.loc[framewise.scan_id.eq(scan)]))
    metrics = [
        ("q_vessel", "q_vessel_base", "q_vessel_fixed128"),
        ("q_tail", "q_tail_base", "q_tail_fixed128"),
        ("relative_tail_burden", "ratio_tail_to_vessel_base", "ratio_tail_to_vessel_fixed128"),
        ("relative_tail_intensity", "ri_tail_base", "ri_tail_fixed128"),
    ]
    for scope_type, scope_value, group in scopes:
        fixed_valid = group.loc[group.fixed128_valid]
        for metric, base_col, fixed_col in metrics:
            base_stats = quantiles(fixed_valid[base_col])
            fixed_stats = quantiles(fixed_valid[fixed_col])
            delta = fixed_valid[fixed_col] - fixed_valid[base_col]
            delta_pct = 100.0 * delta / fixed_valid[base_col]
            d = quantiles(delta)
            ad = quantiles(np.abs(delta))
            dp = quantiles(delta_pct)
            adp = quantiles(np.abs(delta_pct))
            rows.append({
                "scope_type": scope_type,
                "scope_value": scope_value,
                "metric": metric,
                "baseline_valid_n": int(len(group)),
                "fixed128_valid_n": int(group.fixed128_valid.sum()),
                "fixed128_invalid_n": int((~group.fixed128_valid).sum()),
                "baseline_median": base_stats["median"],
                "baseline_q1": base_stats["q1"],
                "baseline_q3": base_stats["q3"],
                "fixed128_median": fixed_stats["median"],
                "fixed128_q1": fixed_stats["q1"],
                "fixed128_q3": fixed_stats["q3"],
                "delta_median": d["median"],
                "abs_delta_median": ad["median"],
                "abs_delta_p95": ad["p95"],
                "delta_pct_median": dp["median"],
                "abs_delta_pct_median": adp["median"],
                "abs_delta_pct_p95": adp["p95"],
            })
    return pd.DataFrame(rows)


def width_summary(framewise: pd.DataFrame) -> pd.DataFrame:
    records = []
    for scan, group in framewise.groupby("scan_id", sort=False):
        stats = quantiles(group["baseline_lateral_width_um"])
        delta = group["baseline_lateral_width_um"] - FIXED_WIDTH_UM
        ds = quantiles(delta)
        records.append({
            "scan_id": scan,
            "flow_mm_s": float(group["flow_mm_s"].iloc[0]),
            "n_valid_frames": int(len(group)),
            "apparent_width_median_um": stats["median"],
            "apparent_width_q1_um": stats["q1"],
            "apparent_width_q3_um": stats["q3"],
            "apparent_width_iqr_um": stats["iqr"],
            "apparent_width_min_um": stats["min"],
            "apparent_width_max_um": stats["max"],
            "fixed_width_um": FIXED_WIDTH_UM,
            "apparent_minus_fixed_median_um": ds["median"],
        })
    return pd.DataFrame(records).sort_values("flow_mm_s").reset_index(drop=True)


def format_pct(x) -> str:
    return "NA" if not np.isfinite(float(x)) else f"{float(x):.3f}%"


def format_num(x) -> str:
    return "NA" if not np.isfinite(float(x)) else f"{float(x):.6g}"


def write_readme(out: Path, framewise: pd.DataFrame, summary: pd.DataFrame, widths: pd.DataFrame,
                 validation: dict, frozen_sha: str, task1_sha: str, workflow_sha: str) -> None:
    overall = summary.loc[summary.scope_type.eq("overall")].set_index("metric")
    lines = [
        "# Fixed 128 µm lateral-width geometry sensitivity — Task 2",
        "",
        "This directory is a **geometry QC / sensitivity analysis**. Frozen run001 remains the primary analysis and is not modified.",
        "",
        "## Geometry contrast",
        "",
        "- Baseline: X4 centre + X1 apparent lateral width + frozen v2.1 z geometry.",
        "- Variant: X4 centre + **fixed 128 µm lateral width** + the same frozen z geometry.",
        "- Source vertical diameter remains 128 µm.",
        "- Tail width follows the corresponding source width; tail length remains 500 µm; guard remains 0.",
        "- Background rule remains skip 3 / take 5 per side and is recomputed from the **variant vessel edges**.",
        "- `sv_raw`, signed background subtraction, 16×16 ellipse fractional area, calibration, and all other quantification rules are unchanged.",
        "",
        "No parameter is chosen according to flow monotonicity or statistical significance. B-scans are spatial positions within a scan volume, not independent biological replicates; no p-values are reported.",
        "",
        "## Input and replay validation",
        "",
        f"- Frozen-valid frames requested: **{validation['baseline_valid_frames']}**.",
        f"- Fixed-128 valid frames: **{validation['fixed128_valid_frames']}**.",
        f"- Release packages verified: **{validation['release_packages_verified']}/{EXPECTED_PACKAGES}**.",
        f"- Valid-frame NPZ SHA-256 checks passed: **{validation['valid_frame_npz_sha256_verified']}**.",
        f"- Baseline geometry replay max relative error across audited scalar quantities: **{validation['baseline_replay_max_relative_error']:.3e}**.",
        f"- X4 vs baseline geometry-centre max difference: **{validation['x4_vs_baseline_geometry_center_max_abs_px']:.3e} px**.",
        "",
        "## Overall sensitivity",
        "",
        "| metric | baseline median | fixed-128 median | median |Δ| (%) | P95 |Δ| (%) |",
        "|---|---:|---:|---:|---:|",
    ]
    labels = {
        "q_vessel": "Q_V",
        "q_tail": "Q_T",
        "relative_tail_burden": "R_B (relative tail burden)",
        "relative_tail_intensity": "RI_tail (Relative Tail Intensity)",
    }
    for metric in ["q_vessel", "q_tail", "relative_tail_burden", "relative_tail_intensity"]:
        row = overall.loc[metric]
        lines.append(f"| {labels[metric]} | {format_num(row.baseline_median)} | {format_num(row.fixed128_median)} | {format_pct(row.abs_delta_pct_median)} | {format_pct(row.abs_delta_pct_p95)} |")
    lines.extend([
        "",
        "## Apparent-width context",
        "",
        "| scan | flow (mm/s) | valid frames | X1 apparent width median [Q1,Q3] µm | fixed width (µm) |",
        "|---|---:|---:|---:|---:|",
    ])
    for row in widths.itertuples(index=False):
        lines.append(f"| {row.scan_id} | {row.flow_mm_s:g} | {row.n_valid_frames} | {row.apparent_width_median_um:.3f} [{row.apparent_width_q1_um:.3f}, {row.apparent_width_q3_um:.3f}] | {row.fixed_width_um:.1f} |")
    lines.extend([
        "",
        "## Files",
        "",
        "- `fixed128_framewise.csv`: paired baseline/fixed-128 values for every frozen-valid frame.",
        "- `fixed128_sensitivity_summary.csv`: overall and per-scan descriptive sensitivity statistics.",
        "- `apparent_width_summary.csv`: frozen X1 apparent-width distribution by scan.",
        "- `release_package_audit.csv`: Release package and per-package array-verification accounting.",
        "- `input_sha256.csv`: repository input hashes used by the audit.",
        "- `validation.json`: inclusion, replay, array-hash, and variant-validity checks.",
        "- `provenance.json`: frozen/input/workflow provenance and fixed geometry definition.",
        "- `audit_fixed128_width.py`: reproducible script.",
        "",
        "## Interpretation boundary",
        "",
        "This comparison quantifies how much the frozen endpoints change when X1-dependent apparent width is replaced by the known 128 µm phantom diameter. It is not used to choose the version with a more favorable flow pattern. Scientific flow interpretation follows only after this QC is complete.",
        "",
        "## Provenance",
        "",
        f"- Frozen formal handoff SHA: `{frozen_sha}`.",
        f"- Task 1 result SHA: `{task1_sha}`.",
        f"- Workflow-trigger SHA: `{workflow_sha}`.",
        f"- 2-D array Release: `{RELEASE_TAG}`.",
    ])
    (out / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--formal-dir", default="results/formal_sv_d128_v21_full2500_run001")
    parser.add_argument("--task1-dir", default="analysis/formal_sv_d128_v21_run001/relative_tail_intensity_full2422")
    parser.add_argument("--output-dir", default="analysis/formal_sv_d128_v21_run001/fixed128_width_sensitivity_full2422")
    parser.add_argument("--frozen-source-sha", required=True)
    parser.add_argument("--task1-source-sha", required=True)
    parser.add_argument("--workflow-source-sha", required=True)
    args = parser.parse_args()

    root = args.root.resolve()
    formal_dir = (root / args.formal_dir).resolve()
    task1_dir = (root / args.task1_dir).resolve()
    out = (root / args.output_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(root / "src"))
    from svrecttail.geometry import VesselGeometry
    from svrecttail.quantification import quantify_frame

    base, frozen_frames, cfg = prepare_inputs(root, formal_dir, task1_dir)
    if len(base) != EXPECTED_VALID:
        raise AssertionError("Unexpected baseline-valid count")

    packages = pd.read_csv(formal_dir / "download_packages.csv")
    arrays = pd.read_csv(formal_dir / "arrays_sha256.csv")
    if len(packages) != EXPECTED_PACKAGES:
        raise AssertionError(f"Expected {EXPECTED_PACKAGES} packages; found {len(packages)}")
    arrays["scan_id"] = arrays["scan_id"].astype(str)
    arrays["frame_index_0based"] = pd.to_numeric(arrays["frame_index_0based"], errors="raise").astype(int)
    array_hash = {(r.scan_id, int(r.frame_index_0based)): str(r.sha256) for r in arrays.itertuples(index=False)}

    base_index = base.set_index(["scan_id", "frame_index"], drop=False)
    records: list[dict] = []
    package_records: list[dict] = []
    replay_max = 0.0
    replay_metric_max = {k: 0.0 for k in ["q_vessel", "q_tail", "ratio_tail_to_vessel", "source_area_um2", "tail_area_um2", "source_mean"]}
    verified_npz = 0
    fixed_invalid: list[dict] = []

    with tempfile.TemporaryDirectory(prefix="sv-fixed128-") as temp_name:
        temp = Path(temp_name)
        for pkg in packages.itertuples(index=False):
            filename = str(pkg.file)
            scan, lo, hi = parse_package_name(filename)
            expected_zip_sha = str(pkg.sha256)
            package_path = temp / filename
            url = f"{RELEASE_BASE}/{filename}"
            print(f"DOWNLOAD {filename}", flush=True)
            download_file(url, package_path)
            observed_zip_sha = sha256_file(package_path)
            observed_bytes = package_path.stat().st_size
            if observed_zip_sha != expected_zip_sha or observed_bytes != int(pkg.bytes):
                raise AssertionError(f"Release package verification failed: {filename}")

            wanted = base.loc[(base.scan_id.eq(scan)) & base.frame_index.between(lo, hi)].copy()
            verified_in_package = 0
            with zipfile.ZipFile(package_path, "r") as zf:
                names = set(zf.namelist())
                for row in wanted.itertuples(index=False):
                    frame = int(row.frame_index)
                    member = f"arrays/{scan}/frame_{frame:03d}.npz"
                    if member not in names:
                        raise FileNotFoundError(f"Missing {member} in {filename}")
                    npz_bytes = zf.read(member)
                    expected_npz_sha = array_hash[(scan, frame)]
                    observed_npz_sha = sha256_bytes(npz_bytes)
                    if observed_npz_sha != expected_npz_sha:
                        raise AssertionError(f"NPZ SHA mismatch for {scan}/{frame}")
                    verified_npz += 1
                    verified_in_package += 1
                    with np.load(io.BytesIO(npz_bytes), allow_pickle=False) as loaded:
                        if "sv_raw" not in loaded.files:
                            raise KeyError(f"sv_raw absent in {member}")
                        sv_raw = np.asarray(loaded["sv_raw"], dtype=np.float64)

                    if sv_raw.ndim != 2:
                        raise AssertionError(f"sv_raw not 2-D for {scan}/{frame}")
                    x4 = float(row.x4_px)
                    ztop = float(row.z_top_edge_px)
                    baseline_geom = VesselGeometry(
                        x_left_edge_px=float(row.baseline_x_left_edge_px),
                        x_right_edge_px=float(row.baseline_x_right_edge_px),
                        z_top_edge_px=ztop,
                        diameter_um=cfg["diameter_um"], dx_um=cfg["dx_um"], dz_um=cfg["dz_um"],
                    )
                    fixed_width_px = FIXED_WIDTH_UM / cfg["dx_um"]
                    fixed_geom = VesselGeometry(
                        x_left_edge_px=x4 - fixed_width_px / 2.0,
                        x_right_edge_px=x4 + fixed_width_px / 2.0,
                        z_top_edge_px=ztop,
                        diameter_um=cfg["diameter_um"], dx_um=cfg["dx_um"], dz_um=cfg["dz_um"],
                    )

                    qkwargs = dict(
                        tail_gap_um=cfg["tail_gap_um"], tail_length_um=cfg["tail_length_um"],
                        background_skip_columns=cfg["skip_columns"],
                        background_strip_width_columns=cfg["strip_width_columns"],
                        background_excluded_side=None,
                        ellipse_supersample=cfg["ellipse_supersample"], source_qc_valid=True,
                    )
                    replay = quantify_frame(sv_raw, baseline_geom, **qkwargs)
                    fixed = quantify_frame(sv_raw, fixed_geom, **qkwargs)

                    frozen_values = {
                        "q_vessel": float(row.q_vessel),
                        "q_tail": float(row.q_tail),
                        "ratio_tail_to_vessel": float(row.ratio_tail_to_vessel),
                        "source_area_um2": float(row.source_area_um2),
                        "tail_area_um2": float(row.tail_area_um2),
                        "source_mean": float(row.source_mean),
                    }
                    replay_values = {
                        "q_vessel": replay.q_vessel, "q_tail": replay.q_tail,
                        "ratio_tail_to_vessel": replay.ratio_tail_to_vessel,
                        "source_area_um2": replay.source_area_um2, "tail_area_um2": replay.tail_area_um2,
                        "source_mean": replay.source_mean,
                    }
                    for metric in frozen_values:
                        err = abs_rel_error(replay_values[metric], frozen_values[metric])
                        if not np.isfinite(err):
                            raise AssertionError(f"Nonfinite baseline replay comparison {metric} at {scan}/{frame}")
                        replay_metric_max[metric] = max(replay_metric_max[metric], err)
                        replay_max = max(replay_max, err)
                    if replay_max > 1e-8:
                        raise AssertionError(f"Baseline replay exceeded tolerance; max relative error {replay_max}")
                    if not replay.valid:
                        raise AssertionError(f"Frozen-valid frame failed baseline replay: {scan}/{frame}: {replay.invalid_reason}")

                    fixed_valid = bool(fixed.valid)
                    if not fixed_valid:
                        fixed_invalid.append({"scan_id": scan, "frame_index": frame, "reason": fixed.invalid_reason})

                    qv1 = fixed.q_vessel if fixed_valid else np.nan
                    qt1 = fixed.q_tail if fixed_valid else np.nan
                    rb1 = fixed.ratio_tail_to_vessel if fixed_valid else np.nan
                    av1 = fixed.source_area_um2 if fixed_valid else np.nan
                    at1 = fixed.tail_area_um2 if fixed_valid else np.nan
                    sv1 = fixed.source_mean if fixed_valid else np.nan
                    st1 = safe_divide(qt1, at1)
                    ri1 = safe_divide(st1, sv1)

                    rec = {
                        "scan_id": scan,
                        "frame_index": frame,
                        "frame_index_0based": frame,
                        "flow_mm_s": float(row.flow_mm_s),
                        "localization_source": row.localization_source,
                        "baseline_valid": True,
                        "fixed128_valid": fixed_valid,
                        "fixed128_invalid_reason": "" if fixed_valid else fixed.invalid_reason,
                        "x4_px": x4,
                        "z_top_edge_px": ztop,
                        "baseline_lateral_width_um": float(row.baseline_lateral_width_um),
                        "fixed128_lateral_width_um": FIXED_WIDTH_UM,
                        "q_vessel_base": float(row.q_vessel), "q_vessel_fixed128": qv1,
                        "q_tail_base": float(row.q_tail), "q_tail_fixed128": qt1,
                        "source_area_um2_base": float(row.source_area_um2), "source_area_um2_fixed128": av1,
                        "tail_area_um2_base": float(row.tail_area_um2), "tail_area_um2_fixed128": at1,
                        "source_mean_base": float(row.source_mean), "source_mean_fixed128": sv1,
                        "tail_mean_base": float(row.tail_mean), "tail_mean_fixed128": st1,
                        "ratio_tail_to_vessel_base": float(row.ratio_tail_to_vessel),
                        "ratio_tail_to_vessel_fixed128": rb1,
                        "ri_tail_base": float(row.ri_tail), "ri_tail_fixed128": ri1,
                        "baseline_background_left_columns": ";".join(map(str, replay.background.left_columns.tolist())),
                        "baseline_background_right_columns": ";".join(map(str, replay.background.right_columns.tolist())),
                        "fixed128_background_left_columns": ";".join(map(str, fixed.background.left_columns.tolist())),
                        "fixed128_background_right_columns": ";".join(map(str, fixed.background.right_columns.tolist())),
                    }
                    for metric, bcol, fcol in [
                        ("q_vessel", "q_vessel_base", "q_vessel_fixed128"),
                        ("q_tail", "q_tail_base", "q_tail_fixed128"),
                        ("ratio_tail_to_vessel", "ratio_tail_to_vessel_base", "ratio_tail_to_vessel_fixed128"),
                        ("ri_tail", "ri_tail_base", "ri_tail_fixed128"),
                    ]:
                        b, f = rec[bcol], rec[fcol]
                        rec[f"{metric}_delta"] = f - b if np.isfinite(f) else np.nan
                        rec[f"{metric}_delta_pct"] = pct_delta(f, b) if np.isfinite(f) else np.nan
                    records.append(rec)

            package_records.append({
                "package": filename,
                "scan_id": scan, "frame_lo": lo, "frame_hi": hi,
                "expected_bytes": int(pkg.bytes), "observed_bytes": int(observed_bytes),
                "expected_sha256": expected_zip_sha, "observed_sha256": observed_zip_sha,
                "package_verified": True,
                "valid_frames_in_package": int(len(wanted)),
                "valid_frame_npz_sha256_verified": int(verified_in_package),
            })
            package_path.unlink()
            print(f"PROCESSED {filename}: valid={len(wanted)}", flush=True)

    framewise = pd.DataFrame.from_records(records).sort_values(["flow_mm_s", "scan_id", "frame_index"]).reset_index(drop=True)
    if len(framewise) != EXPECTED_VALID:
        raise AssertionError(f"Fixed128 framewise count mismatch: {len(framewise)}")
    if verified_npz != EXPECTED_VALID:
        raise AssertionError(f"Expected {EXPECTED_VALID} verified valid NPZs; got {verified_npz}")
    if framewise.duplicated(["scan_id", "frame_index"]).any():
        raise AssertionError("Duplicate fixed128 output frame keys")

    summary = summarize_changes(framewise)
    widths = width_summary(framewise)
    x4_centre_error = np.abs((base.baseline_x_left_edge_px + base.baseline_x_right_edge_px)/2 - base.x4_px)

    validation = {
        "baseline_valid_frames": int(len(framewise)),
        "fixed128_valid_frames": int(framewise.fixed128_valid.sum()),
        "fixed128_invalid_frames": int((~framewise.fixed128_valid).sum()),
        "fixed128_invalid_examples": fixed_invalid[:20],
        "release_packages_expected": EXPECTED_PACKAGES,
        "release_packages_verified": int(sum(r["package_verified"] for r in package_records)),
        "valid_frame_npz_sha256_verified": int(verified_npz),
        "baseline_replay_max_relative_error": float(replay_max),
        "baseline_replay_max_relative_error_by_metric": {k: float(v) for k, v in replay_metric_max.items()},
        "x4_vs_baseline_geometry_center_max_abs_px": float(np.max(x4_centre_error)),
        "fixed_width_um": FIXED_WIDTH_UM,
        "background_rule": "recomputed from fixed geometry edges; skip 3, take 5 each side",
        "negative_residual_clipping": False,
        "parameter_selection_by_flow_pattern": False,
        "p_values_computed": False,
        "biological_replicate_claim": False,
    }

    framewise.to_csv(out / "fixed128_framewise.csv", index=False, float_format="%.17g")
    summary.to_csv(out / "fixed128_sensitivity_summary.csv", index=False, float_format="%.17g")
    widths.to_csv(out / "apparent_width_summary.csv", index=False, float_format="%.17g")
    pd.DataFrame(package_records).to_csv(out / "release_package_audit.csv", index=False)

    input_paths = [
        formal_dir / "frame_results.csv", formal_dir / "download_packages.csv",
        formal_dir / "arrays_sha256.csv", formal_dir / "run_config.json",
        task1_dir / "relative_tail_intensity_framewise.csv",
        root / "src/svrecttail/geometry.py", root / "src/svrecttail/background.py",
        root / "src/svrecttail/quantification.py",
    ]
    input_hashes = []
    for path in input_paths:
        input_hashes.append({
            "input_path": path.relative_to(root).as_posix(),
            "sha256": sha256_file(path), "size_bytes": int(path.stat().st_size),
        })
    pd.DataFrame(input_hashes).to_csv(out / "input_sha256.csv", index=False)

    provenance = {
        "analysis": "Task 2: fixed-128-um lateral geometry sensitivity",
        "status": "geometry_QC_not_scientific_endpoint",
        "frozen_input_directory": args.formal_dir,
        "task1_input_directory": args.task1_dir,
        "output_directory": args.output_dir,
        "formal_release_tag": RELEASE_TAG,
        "frozen_formal_source_sha": args.frozen_source_sha,
        "task1_result_sha": args.task1_source_sha,
        "workflow_source_sha": args.workflow_source_sha,
        "baseline_geometry": "X4-centered + X1 apparent lateral width + frozen v2.1 z_top",
        "variant_geometry": "X4-centered + fixed 128 um lateral width + frozen v2.1 z_top",
        "unchanged": {
            "vertical_source_diameter_um": 128.0,
            "tail_length_um": 500.0,
            "tail_gap_um": 0.0,
            "background_skip_alines": 3,
            "background_take_alines_each_side": 5,
            "ellipse_supersample": 16,
            "sv_signal": "linear sv_raw = Var_t(|E|), N denominator",
            "negative_background_residuals": "retained",
        },
        "background_edge_semantics": "same rule applied to each geometry; fixed128 background columns are selected relative to fixed128 vessel edges",
        "formal_selection": "only baseline valid == True frames enter paired sensitivity",
        "invalid_variant_handling": "variant-invalid frames are retained with validity/reason and excluded from numeric delta summaries; no zero fill",
        "statistical_note": "descriptive spatial sensitivity only; B-scans are not independent biological replicates",
    }
    (out / "provenance.json").write_text(json.dumps(provenance, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    (out / "validation.json").write_text(json.dumps(validation, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    write_readme(out, framewise, summary, widths, validation,
                 args.frozen_source_sha, args.task1_source_sha, args.workflow_source_sha)

    overall = summary.loc[summary.scope_type.eq("overall"), ["metric", "abs_delta_pct_median", "abs_delta_pct_p95"]]
    print(json.dumps({
        "baseline_valid_frames": validation["baseline_valid_frames"],
        "fixed128_valid_frames": validation["fixed128_valid_frames"],
        "packages_verified": validation["release_packages_verified"],
        "npz_verified": validation["valid_frame_npz_sha256_verified"],
        "baseline_replay_max_relative_error": validation["baseline_replay_max_relative_error"],
        "overall_abs_change_pct": overall.to_dict("records"),
        "output_dir": args.output_dir,
    }, indent=2))


if __name__ == "__main__":
    main()
