#!/usr/bin/env python3
"""Spatial autocorrelation analysis for the formal d128 SV run001 volumes.

This is a volume-internal spatial dependence analysis. It uses only formally valid
frames from Task 1, preserves original frame indices, and never treats B-scans as
independent biological replicates. Missing/invalid frames remain missing and are
not zero-filled or interpolated.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

EXPECTED_COUNTS = {
    "flow01": 486,
    "flow03": 468,
    "flow05": 491,
    "flow07": 493,
    "flow10": 484,
}
METRICS = {
    "q_vessel": "Q_V",
    "source_mean": "Sbar_V",
    "q_tail": "Q_T",
    "ri_tail": "RI_tail",
}
SELECTED_LAGS = [0, 1, 2, 3, 5, 10, 20, 30, 50, 75, 100, 150, 200]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def json_write(path: Path, data) -> None:
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def threshold_crossing(lags: np.ndarray, acf: np.ndarray, threshold: float):
    """First downward crossing, linearly interpolated between integer lags."""
    for i in range(1, len(lags)):
        a0, a1 = acf[i - 1], acf[i]
        if not np.isfinite(a0) or not np.isfinite(a1):
            continue
        if a0 > threshold and a1 <= threshold:
            if a0 == a1:
                return float(lags[i])
            fraction = (a0 - threshold) / (a0 - a1)
            return float(lags[i - 1] + fraction * (lags[i] - lags[i - 1]))
    return None


def first_nonpositive(lags: np.ndarray, acf: np.ndarray):
    for lag, value in zip(lags[1:], acf[1:]):
        if np.isfinite(value) and value <= 0:
            return int(lag)
    return None


def acf_pairwise_by_lag(series: np.ndarray, max_lag: int):
    """Global-mean/global-variance ACF with exact frame-index pair separation.

    ACF(h) = mean[(x_i-mu)(x_{i+h}-mu)] / mean[(x_i-mu)^2]
    over pairs where both positions are finite. Invalid positions stay NaN.
    """
    finite = np.isfinite(series)
    values = series[finite]
    if values.size < 3:
        raise ValueError("too few valid values")
    mu = float(np.mean(values))
    variance = float(np.mean((values - mu) ** 2))
    if not np.isfinite(variance) or variance <= 0:
        raise ValueError("non-positive variance")

    rows = []
    acfs = []
    lags = np.arange(max_lag + 1, dtype=int)
    for lag in lags:
        if lag == 0:
            left = series
            right = series
        else:
            left = series[:-lag]
            right = series[lag:]
        mask = np.isfinite(left) & np.isfinite(right)
        n_pairs = int(mask.sum())
        if n_pairs == 0:
            value = np.nan
        else:
            value = float(np.mean((left[mask] - mu) * (right[mask] - mu)) / variance)
        acfs.append(value)
        rows.append((lag, n_pairs, value))
    return np.asarray(lags), np.asarray(acfs, dtype=float), rows, mu, variance


def initial_positive_tau(acf: np.ndarray):
    """1 + 2*sum of the initial strictly-positive ACF sequence from lag 1."""
    positive = []
    stop_lag = None
    for lag in range(1, len(acf)):
        value = acf[lag]
        if not np.isfinite(value) or value <= 0:
            stop_lag = lag
            break
        positive.append(float(value))
    truncated = stop_lag is None
    tau = float(1.0 + 2.0 * sum(positive))
    return tau, stop_lag, truncated


def fmt(x, digits=2):
    if x is None or not np.isfinite(x):
        return "NA"
    return f"{x:.{digits}f}"


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=Path("."))
    p.add_argument("--task1-dir", type=Path, required=True)
    p.add_argument("--formal-dir", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--max-lag", type=int, default=200)
    p.add_argument("--source-sha", required=True)
    p.add_argument("--workflow-source-sha", required=True)
    args = p.parse_args(argv)

    root = args.root.resolve()
    task1 = (root / args.task1_dir).resolve()
    formal = (root / args.formal_dir).resolve()
    out = (root / args.output_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)

    frame_path = task1 / "relative_tail_intensity_framewise.csv"
    manifest_path = formal / "manifest.csv"
    frames = pd.read_csv(frame_path)
    manifest = pd.read_csv(manifest_path)

    if len(frames) != 2422:
        raise RuntimeError(f"expected 2422 Task-1 valid rows, got {len(frames)}")
    if not frames["valid"].astype(bool).all():
        raise RuntimeError("Task-1 framewise input contains non-valid rows")
    observed = frames.groupby("scan_id").size().to_dict()
    if observed != EXPECTED_COUNTS:
        raise RuntimeError(f"scan counts differ from frozen run: {observed}")
    if frames.duplicated(["scan_id", "frame_index_0based"]).any():
        raise RuntimeError("duplicate scan/frame keys")
    if not frames["frame_index_0based"].between(0, 499).all():
        raise RuntimeError("frame index outside 0..499")

    # Physical slow-axis spacing is intentionally not inferred.
    slow_spacing_available = bool(pd.to_numeric(manifest["slow_axis_spacing_um"], errors="coerce").notna().any())
    slow_position_available = bool(pd.to_numeric(manifest["slow_axis_position_um"], errors="coerce").notna().any())
    if slow_spacing_available or slow_position_available:
        raise RuntimeError("slow-axis physical metadata unexpectedly became available; review conversion policy")

    acf_records = []
    summary_records = []
    selected_records = []

    for scan_id in sorted(EXPECTED_COUNTS):
        scan_rows = frames.loc[frames.scan_id == scan_id].copy()
        flow = float(scan_rows["flow_mm_s"].iloc[0])
        n_valid = len(scan_rows)
        for metric, label in METRICS.items():
            if metric not in scan_rows:
                raise RuntimeError(f"missing metric: {metric}")
            values = pd.to_numeric(scan_rows[metric], errors="coerce")
            if not np.isfinite(values).all():
                raise RuntimeError(f"non-finite {metric} in {scan_id}")

            series = np.full(500, np.nan, dtype=float)
            idx = scan_rows["frame_index_0based"].astype(int).to_numpy()
            series[idx] = values.to_numpy(float)
            lags, acf, raw_rows, mean_value, variance = acf_pairwise_by_lag(series, args.max_lag)

            for lag, n_pairs, value in raw_rows:
                acf_records.append({
                    "scan_id": scan_id,
                    "flow_mm_s": flow,
                    "metric": metric,
                    "metric_label": label,
                    "lag_frames": int(lag),
                    "pair_count": int(n_pairs),
                    "acf": value,
                })
                if int(lag) in SELECTED_LAGS:
                    selected_records.append({
                        "scan_id": scan_id,
                        "flow_mm_s": flow,
                        "metric": metric,
                        "metric_label": label,
                        "lag_frames": int(lag),
                        "pair_count": int(n_pairs),
                        "acf": value,
                    })

            lag_1e = threshold_crossing(lags, acf, math.exp(-1.0))
            lag_0p1 = threshold_crossing(lags, acf, 0.1)
            lag_zero = first_nonpositive(lags, acf)
            tau, tau_stop_lag, tau_truncated = initial_positive_tau(acf)
            spatial_ess = float(n_valid / tau) if tau > 0 else None

            summary_records.append({
                "scan_id": scan_id,
                "flow_mm_s": flow,
                "metric": metric,
                "metric_label": label,
                "n_valid_frames": n_valid,
                "series_mean": mean_value,
                "series_variance": variance,
                "acf_lag1": float(acf[1]),
                "acf_lag2": float(acf[2]),
                "acf_lag5": float(acf[5]),
                "acf_lag10": float(acf[10]),
                "acf_lag20": float(acf[20]),
                "acf_lag50": float(acf[50]),
                "acf_lag100": float(acf[100]),
                "lag_at_acf_1_over_e_frames": lag_1e,
                "lag_at_acf_0p1_frames": lag_0p1,
                "first_nonpositive_lag_frames": lag_zero,
                "initial_positive_tau_frames": tau,
                "tau_stop_lag_frames": tau_stop_lag,
                "tau_truncated_at_max_lag": tau_truncated,
                "descriptive_spatial_ess": spatial_ess,
                "max_lag_evaluated": args.max_lag,
            })

    acf_df = pd.DataFrame(acf_records)
    summary_df = pd.DataFrame(summary_records)
    selected_df = pd.DataFrame(selected_records)
    acf_df.to_csv(out / "acf_by_scan_metric.csv", index=False)
    summary_df.to_csv(out / "correlation_length_summary.csv", index=False)
    selected_df.to_csv(out / "acf_selected_lags.csv", index=False)

    # Compact cross-scan summary: median range across the five volumes for each metric.
    aggregate = []
    for metric, label in METRICS.items():
        sub = summary_df.loc[summary_df.metric == metric]
        for field in [
            "acf_lag1", "acf_lag5", "acf_lag10", "acf_lag20",
            "lag_at_acf_1_over_e_frames", "lag_at_acf_0p1_frames",
            "first_nonpositive_lag_frames", "initial_positive_tau_frames",
            "descriptive_spatial_ess",
        ]:
            vals = pd.to_numeric(sub[field], errors="coerce").dropna().to_numpy(float)
            aggregate.append({
                "metric": metric,
                "metric_label": label,
                "quantity": field,
                "n_scans_with_value": int(len(vals)),
                "median_across_scans": float(np.median(vals)) if len(vals) else np.nan,
                "min_across_scans": float(np.min(vals)) if len(vals) else np.nan,
                "max_across_scans": float(np.max(vals)) if len(vals) else np.nan,
            })
    aggregate_df = pd.DataFrame(aggregate)
    aggregate_df.to_csv(out / "cross_scan_autocorrelation_summary.csv", index=False)

    input_sha = pd.DataFrame([
        {"file": str(frame_path.relative_to(root)), "sha256": sha256(frame_path)},
        {"file": str(manifest_path.relative_to(root)), "sha256": sha256(manifest_path)},
    ])
    input_sha.to_csv(out / "input_sha256.csv", index=False)

    validation = {
        "valid_frames": int(len(frames)),
        "scan_valid_counts": {k: int(v) for k, v in observed.items()},
        "metrics": list(METRICS.keys()),
        "frame_axis_length_per_scan": 500,
        "max_lag_evaluated_frames": int(args.max_lag),
        "invalid_frames_interpolated": False,
        "invalid_frames_zero_filled": False,
        "acf_uses_exact_frame_index_separation": True,
        "slow_axis_spacing_um_available": slow_spacing_available,
        "slow_axis_position_um_available": slow_position_available,
        "physical_length_conversion_performed": False,
        "p_values_computed": False,
        "bscan_independent_biological_replicate_claim": False,
        "descriptive_spatial_ess_is_not_biological_n": True,
        "acf_lag0_max_abs_error_from_1": float(np.max(np.abs(acf_df.loc[acf_df.lag_frames == 0, "acf"] - 1.0))),
    }
    json_write(out / "validation.json", validation)

    provenance = {
        "analysis": "volume-internal spatial autocorrelation for formal SV d128 run001",
        "input_task1": str(frame_path.relative_to(root)),
        "input_formal_manifest": str(manifest_path.relative_to(root)),
        "source_branch_sha_before_analysis": args.source_sha,
        "workflow_source_sha": args.workflow_source_sha,
        "selection": "Task-1 formal valid frames only (2422 total)",
        "metrics": METRICS,
        "acf_definition": "mean[(x_i-mu)(x_i+h-mu)] / mean[(x_i-mu)^2], using only finite exact-index pairs at lag h",
        "correlation_length_definitions": {
            "lag_at_acf_1_over_e_frames": "first downward crossing of ACF <= exp(-1), linearly interpolated between integer lags",
            "lag_at_acf_0p1_frames": "first downward crossing of ACF <= 0.1, linearly interpolated",
            "first_nonpositive_lag_frames": "first integer lag with ACF <= 0",
            "initial_positive_tau_frames": "1 + 2*sum ACF(h) over the initial strictly-positive sequence",
            "descriptive_spatial_ess": "n_valid / initial_positive_tau; descriptive spatial-information count only, never biological replicate n",
        },
        "slow_axis_physical_spacing": "not available in frozen manifest; no um conversion",
        "statistical_unit_note": "Each flow is one volume; B-scans are adjacent spatial samples, not independent experiments.",
    }
    json_write(out / "provenance.json", provenance)

    # Generate a compact README from the computed values.
    lines = [
        "# Spatial autocorrelation — formal SV d128 run001",
        "",
        "This analysis quantifies **within-volume spatial dependence** of the four frozen scalar metrics. It does not test a flow effect and does not treat B-scans as independent experimental replicates.",
        "",
        "## Method",
        "",
        "- Input: all 2422 Task-1 formal valid frames.",
        "- Frame index 0–499 is retained as the spatial coordinate. Invalid positions remain missing; no interpolation or zero fill is used.",
        f"- ACF is evaluated from lag 0 to {args.max_lag} B-scan positions.",
        "- The frozen manifest contains no slow-axis spacing/position in µm, so correlation lengths are reported in **B-scan-position units only**.",
        "- `descriptive_spatial_ess` is a volume-internal information-count approximation. It is **not** an independent biological/experimental sample size.",
        "",
        "## Per-volume correlation scales",
        "",
        "| scan | metric | ACF lag1 | ACF lag5 | ACF lag10 | 1/e lag | 0.1 lag | first <=0 lag | tau+ | spatial ESS* |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in summary_df.itertuples(index=False):
        lines.append(
            f"| {row.scan_id} | {row.metric_label} | {row.acf_lag1:.3f} | {row.acf_lag5:.3f} | {row.acf_lag10:.3f} | "
            f"{fmt(row.lag_at_acf_1_over_e_frames)} | {fmt(row.lag_at_acf_0p1_frames)} | "
            f"{fmt(row.first_nonpositive_lag_frames,0)} | {row.initial_positive_tau_frames:.2f} | {row.descriptive_spatial_ess:.1f} |"
        )
    lines += [
        "",
        "`tau+` = initial-positive integrated autocorrelation span. `spatial ESS*` = n_valid / tau+ and is descriptive only.",
        "",
        "## Cross-volume overview",
        "",
        "| metric | median ACF lag1 (range) | median 1/e lag (range) | median first <=0 lag (range) | median tau+ (range) |",
        "|---|---:|---:|---:|---:|",
    ]
    for metric, label in METRICS.items():
        sub = summary_df.loc[summary_df.metric == metric]
        def triplet(field):
            vals = pd.to_numeric(sub[field], errors="coerce").dropna().to_numpy(float)
            if len(vals) == 0:
                return "NA"
            return f"{np.median(vals):.2f} ({np.min(vals):.2f}–{np.max(vals):.2f})"
        lines.append(
            f"| {label} | {triplet('acf_lag1')} | {triplet('lag_at_acf_1_over_e_frames')} | "
            f"{triplet('first_nonpositive_lag_frames')} | {triplet('initial_positive_tau_frames')} |"
        )
    lines += [
        "",
        "## Interpretation boundary",
        "",
        "These correlation scales describe how rapidly the measured metric changes along the slow-axis positions **inside one acquisition**. They can later inform spacing/blocking and subsampling design. They do not create independent replicates for flow inference; independent acquisitions remain the experimental unit for that purpose.",
        "",
        "## Outputs",
        "",
        "- `acf_by_scan_metric.csv`: complete lag-by-lag ACF table.",
        "- `acf_selected_lags.csv`: compact selected-lag table.",
        "- `correlation_length_summary.csv`: per-scan/per-metric correlation scales.",
        "- `cross_scan_autocorrelation_summary.csv`: descriptive five-volume overview.",
        "- `validation.json`, `provenance.json`, `input_sha256.csv`.",
        "- `analyze_spatial_autocorrelation.py`: reproducible script.",
        "",
    ]
    (out / "README.md").write_text("\n".join(lines), encoding="utf-8")

    print(json.dumps({
        "valid_frames": len(frames),
        "scans": EXPECTED_COUNTS,
        "metrics": list(METRICS),
        "max_lag": args.max_lag,
        "slow_axis_um_conversion": "not_available",
        "summary": summary_df[["scan_id", "metric_label", "acf_lag1", "lag_at_acf_1_over_e_frames", "first_nonpositive_lag_frames", "initial_positive_tau_frames", "descriptive_spatial_ess"]].to_dict(orient="records"),
    }, ensure_ascii=False, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
