#!/usr/bin/env python3
"""Systematic subsampling stability for formal SV d128 run001.

This analysis asks how closely a spatially uniform subset of B-scans reproduces
full-volume median endpoints. It does not treat B-scans as independent biological
replicates and does not alter the frozen run001 quantification.

For each scan, stride s selects frames satisfying frame_index % s == phase. All
phases 0..s-1 are audited, so the result does not depend on a favorable start
position. Invalid frames remain absent; no zero filling or interpolation occurs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

METRICS = {
    "q_vessel": "Q_V",
    "source_mean": "Sbar_V",
    "q_tail": "Q_T",
    "ri_tail": "RI_tail",
}
PRIMARY = ("q_tail", "ri_tail")
DEFAULT_STRIDES = (2, 5, 10, 20, 25, 50)
THRESHOLDS = (1.0, 2.0, 5.0, 10.0)
EXPECTED_COUNTS = {"flow01": 486, "flow03": 468, "flow05": 491, "flow07": 493, "flow10": 484}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def q(values, p):
    a = np.asarray(values, dtype=float)
    a = a[np.isfinite(a)]
    return float(np.quantile(a, p)) if a.size else np.nan


def pct_error(value: float, reference: float) -> float:
    if not np.isfinite(value) or not np.isfinite(reference) or reference == 0:
        return np.nan
    return 100.0 * (value - reference) / reference


def summarize_group(g: pd.DataFrame) -> dict:
    abs_pct = g["abs_error_pct"].to_numpy(float)
    signed_pct = g["error_pct"].to_numpy(float)
    abs_raw = g["abs_error"].to_numpy(float)
    out = {
        "n_phase_records": int(len(g)),
        "selected_n_min": int(g["n_selected"].min()),
        "selected_n_median": float(g["n_selected"].median()),
        "selected_n_max": int(g["n_selected"].max()),
        "error_pct_median": q(signed_pct, 0.50),
        "error_pct_min": float(np.nanmin(signed_pct)),
        "error_pct_max": float(np.nanmax(signed_pct)),
        "abs_error_pct_median": q(abs_pct, 0.50),
        "abs_error_pct_p90": q(abs_pct, 0.90),
        "abs_error_pct_p95": q(abs_pct, 0.95),
        "abs_error_pct_max": float(np.nanmax(abs_pct)),
        "abs_error_median": q(abs_raw, 0.50),
        "abs_error_max": float(np.nanmax(abs_raw)),
    }
    for t in THRESHOLDS:
        tag = str(int(t))
        n = int(np.sum(abs_pct <= t))
        out[f"within_{tag}pct_n"] = n
        out[f"within_{tag}pct_fraction"] = n / len(g)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path, default=Path("."))
    ap.add_argument("--task1-dir", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--source-sha", required=True)
    ap.add_argument("--workflow-source-sha", required=True)
    ap.add_argument("--strides", type=int, nargs="*", default=list(DEFAULT_STRIDES))
    args = ap.parse_args()

    root = args.root.resolve()
    task1 = (root / args.task1_dir).resolve() if not args.task1_dir.is_absolute() else args.task1_dir.resolve()
    out = (root / args.output_dir).resolve() if not args.output_dir.is_absolute() else args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)

    input_file = task1 / "relative_tail_intensity_framewise.csv"
    df = pd.read_csv(input_file)
    required = {"scan_id", "frame_index_0based", "valid", *METRICS.keys()}
    missing = required - set(df.columns)
    if missing:
        raise RuntimeError(f"missing required columns: {sorted(missing)}")

    valid = df["valid"].astype(str).str.lower().isin(["true", "1", "yes"])
    df = df.loc[valid].copy()
    if len(df) != 2422:
        raise RuntimeError(f"expected 2422 valid rows, got {len(df)}")
    counts = df.groupby("scan_id").size().to_dict()
    if counts != EXPECTED_COUNTS:
        raise RuntimeError(f"scan valid counts changed: {counts}")
    if df.duplicated(["scan_id", "frame_index_0based"]).any():
        raise RuntimeError("duplicate scan/frame keys")
    for metric in METRICS:
        if not np.isfinite(df[metric].to_numpy(float)).all():
            raise RuntimeError(f"nonfinite valid values in {metric}")

    # Full-volume reference medians.
    refs = []
    ref_map = {}
    for scan, g in df.groupby("scan_id", sort=True):
        flow = float(g["flow_mm_s"].iloc[0]) if "flow_mm_s" in g.columns else float(scan.replace("flow", ""))
        for metric, label in METRICS.items():
            med = float(g[metric].median())
            ref_map[(scan, metric)] = med
            refs.append({
                "scan_id": scan,
                "flow_mm_s": flow,
                "metric": metric,
                "metric_label": label,
                "n_full_valid": int(len(g)),
                "full_median": med,
                "full_q1": float(g[metric].quantile(0.25)),
                "full_q3": float(g[metric].quantile(0.75)),
            })
    refs_df = pd.DataFrame(refs)
    refs_df.to_csv(out / "full_volume_reference.csv", index=False)

    # Exhaustive systematic phase audit.
    phase_rows = []
    coverage_checks = []
    strides = tuple(sorted(set(int(x) for x in args.strides)))
    if any(s < 2 or s > 100 for s in strides):
        raise ValueError("strides must be between 2 and 100")

    for scan, g in df.groupby("scan_id", sort=True):
        frame_index = g["frame_index_0based"].astype(int)
        flow = float(g["flow_mm_s"].iloc[0]) if "flow_mm_s" in g.columns else float(scan.replace("flow", ""))
        for stride in strides:
            assigned = np.zeros(len(g), dtype=int)
            for phase in range(stride):
                mask = (frame_index.to_numpy() % stride) == phase
                assigned += mask.astype(int)
                sub = g.loc[mask]
                if len(sub) < 5:
                    raise RuntimeError(f"too few selected valid frames: {scan}, stride={stride}, phase={phase}, n={len(sub)}")
                first_frame = int(sub["frame_index_0based"].min())
                last_frame = int(sub["frame_index_0based"].max())
                for metric, label in METRICS.items():
                    ref = ref_map[(scan, metric)]
                    est = float(sub[metric].median())
                    delta = est - ref
                    err_pct = pct_error(est, ref)
                    phase_rows.append({
                        "scan_id": scan,
                        "flow_mm_s": flow,
                        "stride_frames": stride,
                        "phase": phase,
                        "n_selected": int(len(sub)),
                        "first_selected_frame": first_frame,
                        "last_selected_frame": last_frame,
                        "metric": metric,
                        "metric_label": label,
                        "full_median": ref,
                        "subsample_median": est,
                        "error": delta,
                        "abs_error": abs(delta),
                        "error_pct": err_pct,
                        "abs_error_pct": abs(err_pct),
                    })
            coverage_checks.append({
                "scan_id": scan,
                "stride_frames": stride,
                "valid_frames": int(len(g)),
                "all_valid_frames_assigned_once": bool(np.all(assigned == 1)),
                "assignment_min": int(assigned.min()),
                "assignment_max": int(assigned.max()),
            })

    phase_df = pd.DataFrame(phase_rows)
    phase_df.to_csv(out / "systematic_phase_results.csv", index=False)
    coverage_df = pd.DataFrame(coverage_checks)
    coverage_df.to_csv(out / "phase_coverage_validation.csv", index=False)
    if not coverage_df["all_valid_frames_assigned_once"].all():
        raise RuntimeError("systematic phase partition failed")

    # Per-scan stability summary.
    per_scan_records = []
    for (scan, stride, metric), g in phase_df.groupby(["scan_id", "stride_frames", "metric"], sort=True):
        rec = {
            "scan_id": scan,
            "flow_mm_s": float(g["flow_mm_s"].iloc[0]),
            "stride_frames": int(stride),
            "metric": metric,
            "metric_label": METRICS[metric],
            "full_median": float(g["full_median"].iloc[0]),
        }
        rec.update(summarize_group(g))
        per_scan_records.append(rec)
    per_scan = pd.DataFrame(per_scan_records)
    per_scan.to_csv(out / "subsampling_summary_by_scan.csv", index=False)

    # Aggregate across all 5 volumes and all phases at each stride.
    cross_records = []
    for (stride, metric), g in phase_df.groupby(["stride_frames", "metric"], sort=True):
        rec = {
            "stride_frames": int(stride),
            "approx_positions_per_500": 500.0 / int(stride),
            "metric": metric,
            "metric_label": METRICS[metric],
            "n_scans": int(g["scan_id"].nunique()),
        }
        rec.update(summarize_group(g))
        worst = g.loc[g["abs_error_pct"].idxmax()]
        rec.update({
            "worst_scan_id": worst["scan_id"],
            "worst_phase": int(worst["phase"]),
            "worst_n_selected": int(worst["n_selected"]),
        })
        cross_records.append(rec)
    cross = pd.DataFrame(cross_records)
    cross.to_csv(out / "cross_scan_subsampling_summary.csv", index=False)

    # Compact planning table for primary endpoints. Thresholds are descriptive
    # planning criteria only; no inferential meaning is attached to them.
    planning = []
    for stride in strides:
        row = {"stride_frames": stride, "approx_positions_per_500": 500.0 / stride}
        primary_rows = phase_df[(phase_df["stride_frames"] == stride) & phase_df["metric"].isin(PRIMARY)]
        for metric in PRIMARY:
            g = primary_rows[primary_rows["metric"] == metric]
            row[f"{metric}_abs_error_pct_p95"] = q(g["abs_error_pct"], 0.95)
            row[f"{metric}_abs_error_pct_max"] = float(g["abs_error_pct"].max())
        row["joint_primary_abs_error_pct_max"] = float(primary_rows["abs_error_pct"].max())
        row["joint_primary_abs_error_pct_p95"] = q(primary_rows["abs_error_pct"], 0.95)
        for t in THRESHOLDS:
            row[f"all_primary_scan_phases_within_{int(t)}pct"] = bool((primary_rows["abs_error_pct"] <= t).all())
            row[f"primary_scan_phase_fraction_within_{int(t)}pct"] = float((primary_rows["abs_error_pct"] <= t).mean())
        planning.append(row)
    planning_df = pd.DataFrame(planning).sort_values("stride_frames")
    planning_df.to_csv(out / "primary_endpoint_planning_matrix.csv", index=False)

    # Identify the sparsest tested grid meeting each all-phase primary criterion.
    criterion_rows = []
    for t in THRESHOLDS:
        col = f"all_primary_scan_phases_within_{int(t)}pct"
        passing = planning_df[planning_df[col]]
        if len(passing):
            chosen = passing.sort_values("stride_frames").iloc[-1]
            criterion_rows.append({
                "criterion_abs_error_pct": t,
                "sparsest_tested_stride_meeting_all_primary_scan_phases": int(chosen["stride_frames"]),
                "approx_positions_per_500": float(chosen["approx_positions_per_500"]),
                "status": "descriptive_planning_candidate_only",
            })
        else:
            criterion_rows.append({
                "criterion_abs_error_pct": t,
                "sparsest_tested_stride_meeting_all_primary_scan_phases": np.nan,
                "approx_positions_per_500": np.nan,
                "status": "none_of_tested_grids_meet_criterion",
            })
    pd.DataFrame(criterion_rows).to_csv(out / "planning_criteria_candidates.csv", index=False)

    input_hashes = pd.DataFrame([
        {"file": str(input_file.relative_to(root)), "sha256": sha256(input_file)},
        {"file": "analysis/formal_sv_d128_v21_run001/spatial_autocorrelation_full2422/correlation_length_summary.csv",
         "sha256": sha256(root / "analysis/formal_sv_d128_v21_run001/spatial_autocorrelation_full2422/correlation_length_summary.csv")},
    ])
    input_hashes.to_csv(out / "input_sha256.csv", index=False)

    validation = {
        "valid_frames": int(len(df)),
        "scan_valid_counts": counts,
        "strides_frames": list(strides),
        "all_phase_partitions_assign_each_valid_frame_exactly_once": bool(coverage_df["all_valid_frames_assigned_once"].all()),
        "invalid_frames_interpolated": False,
        "invalid_frames_zero_filled": False,
        "subsample_uses_original_frame_index": True,
        "volume_estimator": "median",
        "p_values_computed": False,
        "inferential_test_performed": False,
        "bscan_independent_biological_replicate_claim": False,
        "slow_axis_physical_spacing_available": False,
        "planning_thresholds_are_inferential": False,
    }
    (out / "validation.json").write_text(json.dumps(validation, indent=2) + "\n", encoding="utf-8")

    provenance = {
        "analysis": "Systematic subsampling stability, formal SV d128 run001",
        "source_branch_sha": args.source_sha,
        "workflow_source_sha": args.workflow_source_sha,
        "input": str(input_file.relative_to(root)),
        "selection": "valid == True",
        "scheme": "for each stride s and every phase p=0..s-1, keep valid frames whose original frame_index_0based mod s equals p",
        "full_volume_reference": "median across all formally valid frames in that scan",
        "metrics": METRICS,
        "primary_endpoints_for_planning_table": list(PRIMARY),
        "planning_error_thresholds_pct": list(THRESHOLDS),
        "restriction": "descriptive within-volume sampling stability only; not biological sample-size estimation",
    }
    (out / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")

    # README with compact primary-endpoint table.
    lines = [
        "# Systematic subsampling stability — formal SV d128 run001",
        "",
        "This analysis tests whether spatially uniform subsets reproduce the **full-volume median** of the frozen valid-frame endpoints. It does not change run001, interpolate invalid frames, or treat B-scans as independent biological replicates.",
        "",
        "## Sampling design",
        "",
        "For each stride `s`, every phase `0..s-1` is evaluated: a phase keeps formally valid frames whose original `frame_index_0based % s == phase`. Thus no favorable starting offset is selected. Tested strides correspond approximately to 250, 100, 50, 25, 20 and 10 planned positions per 500-frame volume.",
        "",
        "## Primary endpoint stability",
        "",
        "| stride | approx positions / 500 | Q_T P95 abs error | Q_T max abs error | RI_tail P95 abs error | RI_tail max abs error | joint max |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for r in planning_df.itertuples(index=False):
        lines.append(
            f"| {int(r.stride_frames)} | {r.approx_positions_per_500:.1f} | "
            f"{r.q_tail_abs_error_pct_p95:.2f}% | {r.q_tail_abs_error_pct_max:.2f}% | "
            f"{r.ri_tail_abs_error_pct_p95:.2f}% | {r.ri_tail_abs_error_pct_max:.2f}% | "
            f"{r.joint_primary_abs_error_pct_max:.2f}% |"
        )
    lines += [
        "",
        "The error is the percentage difference between a systematic subset median and that scan's full-valid-volume median. P95 and maximum summarize **all five scans and all phase offsets**.",
        "",
        "## Interpretation boundary",
        "",
        "These results quantify spatial subsampling efficiency for one already-acquired volume. They do not estimate the number of independent acquisitions required per flow condition. The latter requires between-volume variance from true independent repeats.",
        "",
        "`planning_criteria_candidates.csv` reports which tested grids satisfy 1%, 2%, 5% or 10% absolute-error criteria for **every** Q_T and RI_tail scan-phase combination. Those thresholds are planning aids, not significance or validation thresholds.",
        "",
        "## Outputs",
        "",
        "- `full_volume_reference.csv`: full-valid-volume medians and quartiles.",
        "- `systematic_phase_results.csv`: every scan × stride × phase × metric result.",
        "- `subsampling_summary_by_scan.csv`: phase robustness within each scan.",
        "- `cross_scan_subsampling_summary.csv`: pooled descriptive phase-error summaries across scans.",
        "- `primary_endpoint_planning_matrix.csv`: compact Q_T / RI_tail planning view.",
        "- `planning_criteria_candidates.csv`: sparsest tested grid meeting each descriptive all-phase error criterion.",
        "- `phase_coverage_validation.csv`, `validation.json`, `provenance.json`, `input_sha256.csv`.",
        "- `analyze_subsampling_stability.py`: reproducible script.",
    ]
    (out / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(json.dumps({
        "valid_frames": len(df),
        "strides": list(strides),
        "phase_metric_records": len(phase_df),
        "primary_planning": planning_df.to_dict(orient="records"),
    }, indent=2))


if __name__ == "__main__":
    main()
