#!/usr/bin/env python3
"""Descriptive five-flow analysis for the frozen 128-um SV run001.

Scientific endpoints follow the frozen order:
    Flow -> Q_V, Sbar_V -> Q_T -> RI_tail -> RI(r)

The primary analysis always uses frozen run001 (X4 centre + X1 apparent width).
The fixed-128 result is carried only as geometry-sensitivity context. B-scans are
spatial positions within one scan volume per flow, not independent biological
replicates. No p-values, inferential tests, curve fitting, threshold detection,
zero filling, or depth interpolation are performed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

EXPECTED_VALID = 2422
EXPECTED_SCANS = ["flow01", "flow03", "flow05", "flow07", "flow10"]
DEPTH_TARGETS_UM = np.arange(0.0, 500.0 + 0.1, 10.0)
DZ_UM = 6.7
MAX_NEAREST_ERROR_UM = DZ_UM / 2.0 + 1e-9
PRIMARY_METRICS = ["q_vessel", "source_mean", "q_tail", "ri_tail"]
SUPPORT_METRICS = ["tail_mean", "ratio_tail_to_vessel"]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def finite_array(values: Iterable[float]) -> np.ndarray:
    a = pd.to_numeric(pd.Series(list(values)), errors="coerce").to_numpy(float)
    return a[np.isfinite(a)]


def describe(values: Iterable[float]) -> dict[str, float | int]:
    a = finite_array(values)
    if a.size == 0:
        return {k: np.nan for k in ["mean", "sd", "min", "q1", "median", "q3", "iqr", "max"]} | {"n": 0}
    q1, med, q3 = np.quantile(a, [0.25, 0.5, 0.75])
    return {
        "n": int(a.size),
        "mean": float(np.mean(a)),
        "sd": float(np.std(a, ddof=1)) if a.size > 1 else np.nan,
        "min": float(np.min(a)),
        "q1": float(q1),
        "median": float(med),
        "q3": float(q3),
        "iqr": float(q3 - q1),
        "max": float(np.max(a)),
    }


def scan_metric_summary(framewise: pd.DataFrame) -> pd.DataFrame:
    records: list[dict] = []
    metrics = PRIMARY_METRICS + SUPPORT_METRICS
    for (scan, flow), group in framewise.groupby(["scan_id", "flow_mm_s"], sort=True):
        rec: dict[str, object] = {
            "scan_id": scan,
            "flow_mm_s": float(flow),
            "n_valid_frames": int(len(group)),
            "frame_min": int(group.frame_index.min()),
            "frame_max": int(group.frame_index.max()),
        }
        for metric in metrics:
            stats = describe(group[metric])
            for name, value in stats.items():
                rec[f"{metric}_{name}"] = value
            med = stats["median"]
            rec[f"{metric}_relative_iqr"] = (
                float(stats["iqr"] / abs(med)) if np.isfinite(med) and med != 0 else np.nan
            )
        records.append(rec)
    return pd.DataFrame(records).sort_values("flow_mm_s").reset_index(drop=True)


def metric_pattern(summary: pd.DataFrame) -> pd.DataFrame:
    records: list[dict] = []
    metric_labels = {
        "q_vessel": "Q_V",
        "source_mean": "Sbar_V",
        "q_tail": "Q_T",
        "ri_tail": "RI_tail",
        "tail_mean": "Sbar_T",
        "ratio_tail_to_vessel": "R_B_relative_tail_burden",
    }
    ref = summary.loc[summary.flow_mm_s.eq(1.0)].iloc[0]
    for metric in PRIMARY_METRICS + SUPPORT_METRICS:
        ref_median = float(ref[f"{metric}_median"])
        for row in summary.itertuples(index=False):
            median = float(getattr(row, f"{metric}_median"))
            records.append({
                "metric": metric,
                "metric_label": metric_labels[metric],
                "endpoint_role": "primary" if metric in PRIMARY_METRICS else "supporting",
                "scan_id": row.scan_id,
                "flow_mm_s": float(row.flow_mm_s),
                "n_valid_frames": int(row.n_valid_frames),
                "median": median,
                "q1": float(getattr(row, f"{metric}_q1")),
                "q3": float(getattr(row, f"{metric}_q3")),
                "iqr": float(getattr(row, f"{metric}_iqr")),
                "relative_iqr": float(getattr(row, f"{metric}_relative_iqr")),
                "median_ratio_vs_flow01": median / ref_median,
                "median_change_pct_vs_flow01": 100.0 * (median - ref_median) / ref_median,
            })
    return pd.DataFrame(records)


def geometry_qc_context(fixed: pd.DataFrame) -> pd.DataFrame:
    fixed = fixed.copy()
    fixed["scan_id"] = fixed["scan_id"].astype(str)
    fixed["frame_index"] = pd.to_numeric(fixed["frame_index"], errors="raise").astype(int)
    fixed_valid = fixed.loc[fixed["fixed128_valid"].astype(bool)].copy()
    if len(fixed_valid) != EXPECTED_VALID:
        raise AssertionError(f"Expected {EXPECTED_VALID} paired fixed128-valid rows; found {len(fixed_valid)}")
    specs = [
        ("Q_V", "q_vessel_base", "q_vessel_fixed128"),
        ("Q_T", "q_tail_base", "q_tail_fixed128"),
        ("R_B_relative_tail_burden", "ratio_tail_to_vessel_base", "ratio_tail_to_vessel_fixed128"),
        ("RI_tail", "ri_tail_base", "ri_tail_fixed128"),
        ("Sbar_V", "source_mean_base", "source_mean_fixed128"),
        ("Sbar_T", "tail_mean_base", "tail_mean_fixed128"),
    ]
    records = []
    for (scan, flow), group in fixed_valid.groupby(["scan_id", "flow_mm_s"], sort=True):
        for label, base_col, var_col in specs:
            base_stats = describe(group[base_col])
            var_stats = describe(group[var_col])
            records.append({
                "scan_id": scan,
                "flow_mm_s": float(flow),
                "metric_label": label,
                "n_paired": int(len(group)),
                "primary_run001_median": base_stats["median"],
                "primary_run001_q1": base_stats["q1"],
                "primary_run001_q3": base_stats["q3"],
                "fixed128_qc_median": var_stats["median"],
                "fixed128_qc_q1": var_stats["q1"],
                "fixed128_qc_q3": var_stats["q3"],
                "fixed128_vs_primary_median_change_pct": 100.0 * (var_stats["median"] - base_stats["median"]) / base_stats["median"],
            })
    return pd.DataFrame(records).sort_values(["flow_mm_s", "metric_label"]).reset_index(drop=True)


def nearest_depth_records(task1_dir: Path, framewise: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    profile_paths = sorted((task1_dir / "profiles").glob("*_relative_intensity.csv.gz"))
    if len(profile_paths) != 25:
        raise AssertionError(f"Expected 25 Task1 profile chunks, found {len(profile_paths)}")

    valid_keys = set(zip(framewise.scan_id.astype(str), framewise.frame_index.astype(int)))
    out_records: list[dict] = []
    source_rows = 0
    tail_rows = 0

    for path in profile_paths:
        profile = pd.read_csv(path, compression="gzip", low_memory=False)
        source_rows += len(profile)
        required = {"scan_id", "frame_index", "flow_mm_s", "localization_source", "source_mean", "r_um", "T", "ri_r", "tail_z_fraction"}
        missing = required.difference(profile.columns)
        if missing:
            raise KeyError(f"{path.name} missing Task1 profile columns {sorted(missing)}")
        profile["scan_id"] = profile["scan_id"].astype(str)
        profile["frame_index"] = pd.to_numeric(profile["frame_index"], errors="raise").astype(int)
        profile["r_um"] = pd.to_numeric(profile["r_um"], errors="coerce")
        profile["ri_r"] = pd.to_numeric(profile["ri_r"], errors="coerce")
        profile["T"] = pd.to_numeric(profile["T"], errors="coerce")
        profile["tail_z_fraction"] = pd.to_numeric(profile["tail_z_fraction"], errors="coerce")

        profile = profile.loc[profile.tail_z_fraction.gt(0)].copy()
        tail_rows += len(profile)
        for (scan, frame), group in profile.groupby(["scan_id", "frame_index"], sort=False):
            key = (str(scan), int(frame))
            if key not in valid_keys:
                raise AssertionError(f"Unexpected non-valid frame in Task1 valid profile output: {key}")
            g = group.sort_values("r_um")
            r = g.r_um.to_numpy(float)
            if not np.isfinite(r).all() or r.size == 0:
                raise AssertionError(f"Invalid r grid at {key}")
            for target in DEPTH_TARGETS_UM:
                idx = int(np.argmin(np.abs(r - target)))
                selected = g.iloc[idx]
                error = float(selected.r_um - target)
                if abs(error) > MAX_NEAREST_ERROR_UM:
                    raise AssertionError(f"Nearest depth error {error} um exceeds half axial pixel for {key} at {target}")
                out_records.append({
                    "scan_id": str(scan),
                    "frame_index": int(frame),
                    "flow_mm_s": float(selected.flow_mm_s),
                    "localization_source": selected.localization_source,
                    "target_r_um": float(target),
                    "actual_r_um": float(selected.r_um),
                    "depth_error_um": error,
                    "abs_depth_error_um": abs(error),
                    "tail_z_fraction": float(selected.tail_z_fraction),
                    "source_mean": float(selected.source_mean),
                    "T": float(selected["T"]),
                    "ri_r": float(selected.ri_r),
                })

    anchors = pd.DataFrame.from_records(out_records)
    expected_rows = EXPECTED_VALID * len(DEPTH_TARGETS_UM)
    if len(anchors) != expected_rows:
        counts = anchors.groupby(["scan_id", "frame_index"]).size()
        raise AssertionError(f"Expected {expected_rows} anchor records; got {len(anchors)}. Per-frame count range {counts.min()}..{counts.max()}")
    if anchors.duplicated(["scan_id", "frame_index", "target_r_um"]).any():
        raise AssertionError("Duplicate depth anchor keys")
    if anchors.ri_r.isna().any() or (~np.isfinite(anchors.ri_r.to_numpy(float))).any():
        raise AssertionError("Nonfinite RI(r) selected at a depth anchor")

    key_count = anchors[["scan_id", "frame_index"]].drop_duplicates().shape[0]
    if key_count != EXPECTED_VALID:
        raise AssertionError(f"Depth anchors cover {key_count}, not {EXPECTED_VALID}, valid frames")

    validation = {
        "profile_chunks": len(profile_paths),
        "task1_profile_rows_read": int(source_rows),
        "task1_tail_support_rows": int(tail_rows),
        "depth_targets_count": int(len(DEPTH_TARGETS_UM)),
        "depth_target_min_um": float(DEPTH_TARGETS_UM.min()),
        "depth_target_max_um": float(DEPTH_TARGETS_UM.max()),
        "depth_target_step_um": 10.0,
        "depth_anchor_records": int(len(anchors)),
        "depth_anchor_unique_frames": int(key_count),
        "depth_anchor_max_abs_error_um": float(anchors.abs_depth_error_um.max()),
        "depth_anchor_median_abs_error_um": float(anchors.abs_depth_error_um.median()),
        "depth_interpolation_performed": False,
        "depth_curve_fitting_performed": False,
        "negative_ri_anchor_records": int((anchors.ri_r < 0).sum()),
    }
    return anchors.sort_values(["flow_mm_s", "scan_id", "frame_index", "target_r_um"]).reset_index(drop=True), validation


def depth_summary(anchors: pd.DataFrame) -> pd.DataFrame:
    records = []
    for (scan, flow, target), group in anchors.groupby(["scan_id", "flow_mm_s", "target_r_um"], sort=True):
        stats = describe(group.ri_r)
        tstats = describe(group.T)
        records.append({
            "scan_id": scan,
            "flow_mm_s": float(flow),
            "target_r_um": float(target),
            "n_frames": int(len(group)),
            "ri_r_mean": stats["mean"],
            "ri_r_q1": stats["q1"],
            "ri_r_median": stats["median"],
            "ri_r_q3": stats["q3"],
            "ri_r_iqr": stats["iqr"],
            "ri_r_min": stats["min"],
            "ri_r_max": stats["max"],
            "ri_r_negative_n": int((group.ri_r < 0).sum()),
            "ri_r_negative_fraction": float((group.ri_r < 0).mean()),
            "T_median": tstats["median"],
            "actual_r_median_um": float(group.actual_r_um.median()),
            "abs_depth_error_median_um": float(group.abs_depth_error_um.median()),
            "abs_depth_error_max_um": float(group.abs_depth_error_um.max()),
        })
    return pd.DataFrame(records).sort_values(["flow_mm_s", "target_r_um"]).reset_index(drop=True)


def selected_depth_table(depth: pd.DataFrame) -> pd.DataFrame:
    selected = {0.0, 50.0, 100.0, 200.0, 300.0, 400.0, 500.0}
    return depth.loc[depth.target_r_um.isin(selected)].copy().reset_index(drop=True)


def fmt_interval(med, q1, q3, *, scale=1.0, digits=3) -> str:
    values = np.asarray([med, q1, q3], dtype=float) / scale
    if not np.isfinite(values).all():
        return "NA"
    return f"{values[0]:.{digits}f} [{values[1]:.{digits}f}, {values[2]:.{digits}f}]"


def write_readme(out: Path, summary: pd.DataFrame, pattern: pd.DataFrame, depth_selected: pd.DataFrame,
                 geom_qc: pd.DataFrame, validation: dict, source_sha: str, task2_sha: str, workflow_sha: str) -> None:
    lines = [
        "# 128 µm SV-OCTA: five-flow descriptive analysis",
        "",
        "This is the first scientific flow-pattern analysis after completion of the frozen localization/quantification work, Relative Tail Intensity derivation, and fixed-128 geometry sensitivity QC.",
        "",
        "## Statistical unit and scope",
        "",
        "Each flow condition currently corresponds to one complete scan volume. The included B-scans are spatial positions within that volume, not independent vessels or independent experimental replicates. Results below describe volume-level spatial distributions and effect patterns. No ANOVA, t-test, p-value, or population-level flow-effect claim is made.",
        "",
        "The primary endpoint values are the frozen run001 geometry (X4 centre + X1 apparent width). Fixed-128 values are shown only as geometry-sensitivity context and do not replace the primary analysis.",
        "",
        "## Primary scan-level results",
        "",
        "| flow (mm/s) | valid frames | Q_V median [Q1,Q3] ×10^12 | Sbar_V median [Q1,Q3] ×10^8 | Q_T median [Q1,Q3] ×10^12 | RI_tail median [Q1,Q3] | Sbar_T median [Q1,Q3] ×10^7 |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in summary.itertuples(index=False):
        lines.append(
            f"| {row.flow_mm_s:g} | {row.n_valid_frames} | "
            f"{fmt_interval(row.q_vessel_median,row.q_vessel_q1,row.q_vessel_q3,scale=1e12)} | "
            f"{fmt_interval(row.source_mean_median,row.source_mean_q1,row.source_mean_q3,scale=1e8)} | "
            f"{fmt_interval(row.q_tail_median,row.q_tail_q1,row.q_tail_q3,scale=1e12)} | "
            f"{fmt_interval(row.ri_tail_median,row.ri_tail_q1,row.ri_tail_q3,digits=4)} | "
            f"{fmt_interval(row.tail_mean_median,row.tail_mean_q1,row.tail_mean_q3,scale=1e7)} |"
        )

    patt = pattern.set_index(["metric", "flow_mm_s"])
    lines.extend([
        "",
        "## Descriptive flow pattern",
        "",
        "Using the 1 mm/s scan median only as a descriptive reference:",
        "",
        "| flow (mm/s) | Δ median Q_V | Δ median Sbar_V | Δ median Q_T | Δ median RI_tail |",
        "|---:|---:|---:|---:|---:|",
    ])
    for flow in [1.0, 3.0, 5.0, 7.0, 10.0]:
        vals = [float(patt.loc[(metric, flow), "median_change_pct_vs_flow01"]) for metric in PRIMARY_METRICS]
        lines.append(f"| {flow:g} | {vals[0]:+.2f}% | {vals[1]:+.2f}% | {vals[2]:+.2f}% | {vals[3]:+.2f}% |")

    # Fixed128 context for RI specifically.
    fixed_ri = geom_qc.loc[geom_qc.metric_label.eq("RI_tail")].sort_values("flow_mm_s")
    lines.extend([
        "",
        "### Geometry-QC context",
        "",
        "The fixed-128 sensitivity changed absolute magnitudes, especially Q_T, but did not remove the distinct flow07 relative-intensity pattern. RI_tail medians under fixed-128 geometry are:",
        "",
        "| flow (mm/s) | primary RI_tail median | fixed-128 RI_tail median | fixed128 vs primary |",
        "|---:|---:|---:|---:|",
    ])
    for row in fixed_ri.itertuples(index=False):
        lines.append(f"| {row.flow_mm_s:g} | {row.primary_run001_median:.4f} | {row.fixed128_qc_median:.4f} | {row.fixed128_vs_primary_median_change_pct:+.2f}% |")

    lines.extend([
        "",
        "## Relative intensity with depth RI(r)",
        "",
        "For cross-volume description, every valid frame is sampled at fixed physical depth targets from 0 to 500 µm in 10 µm steps. For each target, the nearest original axial profile sample is selected and retained only within half one axial pixel (3.35 µm). This is nearest-sample reporting: **no interpolation and no curve fitting**. Negative background-corrected RI(r) values are retained.",
        "",
        f"Observed maximum nearest-sample depth mismatch: **{validation['depth_anchor_max_abs_error_um']:.3f} µm**.",
        "",
        "Selected depths are shown below; `ri_depth_summary_10um.csv` contains all 51 depth targets.",
        "",
        "| flow (mm/s) | depth (µm) | RI(r) median [Q1,Q3] | negative-frame fraction |",
        "|---:|---:|---:|---:|",
    ])
    for row in depth_selected.itertuples(index=False):
        lines.append(f"| {row.flow_mm_s:g} | {row.target_r_um:.0f} | {row.ri_r_median:.4f} [{row.ri_r_q1:.4f}, {row.ri_r_q3:.4f}] | {row.ri_r_negative_fraction:.3f} |")

    lines.extend([
        "",
        "## Interpretation",
        "",
        "Across the current five volumes, Q_V and Sbar_V are similar for 1, 3, 5 and 10 mm/s and lower in the 7 mm/s volume. Q_T varies less across the five volumes. Sbar_T is also comparatively stable, so the elevated RI_tail in flow07 is driven primarily by the lower vessel-reference intensity Sbar_V rather than by a proportional rise in absolute tail signal. The fixed-128 geometry QC preserves this RI_tail ordering, so the pattern is not explained solely by the narrower X1 apparent width in flow07.",
        "",
        "These statements describe the present scan volumes. Independent acquisitions are required before assigning the observed volume-to-volume pattern to a general flow-speed effect.",
        "",
        "RI(r) is reported as a signed relative signal at fixed physical depths. No detection depth is inferred because a matched blank is unavailable; detection depth remains `not_evaluated / NA`.",
        "",
        "## Outputs",
        "",
        "- `scan_metrics_summary.csv`: spatial descriptive distributions for Q_V, Sbar_V, Q_T, RI_tail and supporting Sbar_T/R_B.",
        "- `flow_metric_pattern.csv`: scan medians and descriptive changes relative to the 1 mm/s scan.",
        "- `geometry_qc_context.csv`: primary run001 vs fixed-128 scan medians; QC only.",
        "- `ri_depth_anchor_framewise.csv.gz`: nearest original RI(r) sample for every valid frame × 51 physical depth targets.",
        "- `ri_depth_summary_10um.csv`: scan-wise RI(r) median/IQR and signed-value counts at all depth targets.",
        "- `ri_depth_selected.csv`: 0/50/100/200/300/400/500 µm subset for compact review.",
        "- `input_sha256.csv`, `validation.json`, `provenance.json`, and this reproducible script.",
        "",
        "## Provenance",
        "",
        f"- Task 1 / primary derived metric source SHA: `{source_sha}`.",
        f"- Task 2 geometry-QC result SHA: `{task2_sha}`.",
        f"- Workflow-trigger SHA: `{workflow_sha}`.",
    ])
    (out / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--task1-dir", default="analysis/formal_sv_d128_v21_run001/relative_tail_intensity_full2422")
    parser.add_argument("--task2-dir", default="analysis/formal_sv_d128_v21_run001/fixed128_width_sensitivity_full2422")
    parser.add_argument("--output-dir", default="analysis/formal_sv_d128_v21_run001/flow_descriptive_d128_run001")
    parser.add_argument("--task1-source-sha", required=True)
    parser.add_argument("--task2-source-sha", required=True)
    parser.add_argument("--workflow-source-sha", required=True)
    args = parser.parse_args()

    root = args.root.resolve()
    task1_dir = (root / args.task1_dir).resolve()
    task2_dir = (root / args.task2_dir).resolve()
    out = (root / args.output_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)

    framewise = pd.read_csv(task1_dir / "relative_tail_intensity_framewise.csv", low_memory=False)
    framewise["scan_id"] = framewise["scan_id"].astype(str)
    framewise["frame_index"] = pd.to_numeric(framewise["frame_index"], errors="raise").astype(int)
    if len(framewise) != EXPECTED_VALID or not framewise.valid.astype(bool).all():
        raise AssertionError("Primary Task1 input is not exactly the 2422 frozen-valid frames")
    scans = framewise.sort_values("flow_mm_s").scan_id.drop_duplicates().tolist()
    if scans != EXPECTED_SCANS:
        raise AssertionError(f"Unexpected scan set/order: {scans}")

    fixed = pd.read_csv(task2_dir / "fixed128_framewise.csv", low_memory=False)
    summary = scan_metric_summary(framewise)
    pattern = metric_pattern(summary)
    geom_qc = geometry_qc_context(fixed)
    anchors, depth_validation = nearest_depth_records(task1_dir, framewise)
    depth = depth_summary(anchors)
    selected = selected_depth_table(depth)

    # Cross-file consistency checks.
    if not np.array_equal(summary.n_valid_frames.to_numpy(int), np.array([486, 468, 491, 493, 484])):
        raise AssertionError(f"Unexpected per-scan valid counts: {summary[['scan_id','n_valid_frames']].to_dict('records')}")
    expected_depth_rows = len(EXPECTED_SCANS) * len(DEPTH_TARGETS_UM)
    if len(depth) != expected_depth_rows:
        raise AssertionError(f"Expected {expected_depth_rows} scan-depth summary rows, found {len(depth)}")

    validation = {
        "primary_valid_frames": int(len(framewise)),
        "scan_valid_counts": {r.scan_id: int(r.n_valid_frames) for r in summary.itertuples(index=False)},
        "scans": EXPECTED_SCANS,
        "flow_mm_s": [1, 3, 5, 7, 10],
        "primary_geometry": "frozen run001 X4 center + X1 apparent width",
        "fixed128_role": "geometry QC only",
        "p_values_computed": False,
        "inferential_test_performed": False,
        "bscan_independent_replicate_claim": False,
        "zero_fill_or_invalid_frame_interpolation": False,
        "detection_depth": "not_evaluated / NA",
    }
    validation.update(depth_validation)

    summary.to_csv(out / "scan_metrics_summary.csv", index=False, float_format="%.17g")
    pattern.to_csv(out / "flow_metric_pattern.csv", index=False, float_format="%.17g")
    geom_qc.to_csv(out / "geometry_qc_context.csv", index=False, float_format="%.17g")
    anchors.to_csv(
        out / "ri_depth_anchor_framewise.csv.gz", index=False, float_format="%.17g",
        compression={"method": "gzip", "compresslevel": 9, "mtime": 0},
    )
    depth.to_csv(out / "ri_depth_summary_10um.csv", index=False, float_format="%.17g")
    selected.to_csv(out / "ri_depth_selected.csv", index=False, float_format="%.17g")

    input_paths = [
        task1_dir / "relative_tail_intensity_framewise.csv",
        task1_dir / "validation.json",
        task1_dir / "provenance.json",
        task2_dir / "fixed128_framewise.csv",
        task2_dir / "validation.json",
        task2_dir / "provenance.json",
    ] + sorted((task1_dir / "profiles").glob("*_relative_intensity.csv.gz"))
    pd.DataFrame([
        {"input_path": p.relative_to(root).as_posix(), "sha256": sha256_file(p), "size_bytes": int(p.stat().st_size)}
        for p in input_paths
    ]).to_csv(out / "input_sha256.csv", index=False)

    provenance = {
        "analysis": "128-um SV five-flow descriptive analysis",
        "scientific_endpoint_order": ["Q_V", "Sbar_V", "Q_T", "RI_tail", "RI(r)"],
        "primary_source": args.task1_dir,
        "primary_source_sha": args.task1_source_sha,
        "geometry_qc_source": args.task2_dir,
        "geometry_qc_source_sha": args.task2_source_sha,
        "workflow_source_sha": args.workflow_source_sha,
        "primary_geometry": "X4-centered + X1 apparent lateral width + frozen v2.1 z_upper",
        "fixed128_geometry_role": "sensitivity/QC only; not substituted for primary endpoints",
        "ri_depth_reporting": {
            "targets_um": DEPTH_TARGETS_UM.tolist(),
            "selection": "nearest original axial profile sample within half one dz pixel",
            "dz_um": DZ_UM,
            "maximum_allowed_error_um": MAX_NEAREST_ERROR_UM,
            "interpolation": False,
            "curve_fitting": False,
        },
        "statistics": "descriptive spatial distributions only; one scan volume per flow condition",
        "p_values": False,
        "detection_depth": "not_evaluated / NA without matched blank",
    }
    (out / "provenance.json").write_text(json.dumps(provenance, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    (out / "validation.json").write_text(json.dumps(validation, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    write_readme(out, summary, pattern, selected, geom_qc, validation,
                 args.task1_source_sha, args.task2_source_sha, args.workflow_source_sha)

    primary = pattern.loc[pattern.endpoint_role.eq("primary"),
                          ["metric_label", "flow_mm_s", "median", "median_change_pct_vs_flow01"]]
    print(json.dumps({
        "valid_frames": int(len(framewise)),
        "depth_anchor_rows": int(len(anchors)),
        "depth_anchor_max_abs_error_um": validation["depth_anchor_max_abs_error_um"],
        "primary_metric_pattern": primary.to_dict("records"),
        "output_dir": args.output_dir,
    }, indent=2))


if __name__ == "__main__":
    main()
