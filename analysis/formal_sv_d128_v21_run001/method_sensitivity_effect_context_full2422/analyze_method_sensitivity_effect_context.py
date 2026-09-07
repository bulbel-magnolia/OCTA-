#!/usr/bin/env python3
"""Compare observed scan-median flow contrasts with audited method sensitivities.

This is descriptive context only. It does not combine perturbations into a single
error term, does not compute p-values, and does not treat B-scans as biological
replicates. Background +2 is a local method perturbation, systematic subsampling
is a sampling-design perturbation, and fixed-128 width is a geometry stress test.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

EXPECTED_VALID = 2422
METRICS = [
    ("q_vessel", "Q_V"),
    ("source_mean", "Sbar_V"),
    ("q_tail", "Q_T"),
    ("ri_tail", "RI_tail"),
]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def qstats(values) -> dict:
    a = np.asarray(pd.to_numeric(pd.Series(values), errors="coerce"), dtype=float)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return {"n": 0, "median": np.nan, "p95": np.nan, "max": np.nan}
    return {
        "n": int(a.size),
        "median": float(np.quantile(a, .50)),
        "p95": float(np.quantile(a, .95)),
        "max": float(np.max(a)),
    }


def pct_delta(new, old):
    new = np.asarray(new, dtype=float)
    old = np.asarray(old, dtype=float)
    out = np.full_like(new, np.nan, dtype=float)
    m = np.isfinite(new) & np.isfinite(old) & (old != 0)
    out[m] = 100.0 * (new[m] - old[m]) / old[m]
    return out


def sensitivity_record(method, severity, metric, label, delta_pct, note):
    s = qstats(np.abs(delta_pct))
    return {
        "method": method,
        "perturbation_role": severity,
        "metric": metric,
        "metric_label": label,
        "n": s["n"],
        "abs_delta_pct_median": s["median"],
        "abs_delta_pct_p95": s["p95"],
        "abs_delta_pct_max": s["max"],
        "note": note,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path("."))
    ap.add_argument("--task1-dir", default="analysis/formal_sv_d128_v21_run001/relative_tail_intensity_full2422")
    ap.add_argument("--background-dir", default="analysis/formal_sv_d128_v21_run001/background_skip_plus2_full2500")
    ap.add_argument("--geometry-dir", default="analysis/formal_sv_d128_v21_run001/fixed128_width_sensitivity_full2422")
    ap.add_argument("--subsampling-dir", default="analysis/formal_sv_d128_v21_run001/subsampling_stability_full2422")
    ap.add_argument("--output-dir", default="analysis/formal_sv_d128_v21_run001/method_sensitivity_effect_context_full2422")
    ap.add_argument("--source-sha", required=True)
    ap.add_argument("--workflow-source-sha", required=True)
    args = ap.parse_args()

    root = args.root.resolve()
    t1dir = (root / args.task1_dir).resolve()
    bgdir = (root / args.background_dir).resolve()
    geodir = (root / args.geometry_dir).resolve()
    subdir = (root / args.subsampling_dir).resolve()
    out = (root / args.output_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)

    t1_path = t1dir / "relative_tail_intensity_framewise.csv"
    bg_path = bgdir / "background_plus2_framewise.csv"
    geo_path = geodir / "fixed128_framewise.csv"
    sub_path = subdir / "cross_scan_subsampling_summary.csv"

    t1 = pd.read_csv(t1_path, low_memory=False)
    bg = pd.read_csv(bg_path, low_memory=False)
    geo = pd.read_csv(geo_path, low_memory=False)
    sub = pd.read_csv(sub_path, low_memory=False)

    if len(t1) != EXPECTED_VALID:
        raise AssertionError(f"Task1 expected {EXPECTED_VALID} valid rows, got {len(t1)}")
    if len(bg) != EXPECTED_VALID:
        raise AssertionError(f"Background audit expected {EXPECTED_VALID} baseline-valid rows, got {len(bg)}")
    if len(geo) != EXPECTED_VALID:
        raise AssertionError(f"Geometry audit expected {EXPECTED_VALID} rows, got {len(geo)}")

    for df in (t1, bg, geo):
        df["scan_id"] = df["scan_id"].astype(str)
    if "frame_index" not in t1.columns:
        t1["frame_index"] = pd.to_numeric(t1["frame_index_0based"], errors="raise").astype(int)
    else:
        t1["frame_index"] = pd.to_numeric(t1["frame_index"], errors="raise").astype(int)
    bg["frame_index"] = pd.to_numeric(bg["frame_index_0based"], errors="raise").astype(int)
    if "frame_index" not in geo.columns:
        geo["frame_index"] = pd.to_numeric(geo["frame_index_0based"], errors="raise").astype(int)
    else:
        geo["frame_index"] = pd.to_numeric(geo["frame_index"], errors="raise").astype(int)

    keys = ["scan_id", "frame_index"]
    if t1.duplicated(keys).any() or bg.duplicated(keys).any() or geo.duplicated(keys).any():
        raise AssertionError("Duplicate scan/frame key in an input table")

    base = t1[keys + ["flow_mm_s", "q_vessel", "source_area_um2", "q_tail", "tail_area_um2", "source_mean", "ri_tail"]].copy()
    b = base.merge(bg[keys + ["variant_valid", "q_vessel_plus2", "q_tail_plus2", "ratio_plus2"]], on=keys, validate="one_to_one")
    if len(b) != EXPECTED_VALID or not b["variant_valid"].astype(bool).all():
        raise AssertionError("Background +2 does not cover all Task1-valid frames")
    b["source_mean_plus2"] = b["q_vessel_plus2"] / b["source_area_um2"]
    b["tail_mean_plus2"] = b["q_tail_plus2"] / b["tail_area_um2"]
    b["ri_tail_plus2"] = b["tail_mean_plus2"] / b["source_mean_plus2"]
    ri_check = b["ratio_plus2"] * b["source_area_um2"] / b["tail_area_um2"]
    if not np.allclose(b["ri_tail_plus2"], ri_check, rtol=2e-12, atol=1e-12):
        raise AssertionError("Background +2 RI identity failed")

    sens = []
    for metric, label in METRICS:
        variant_col = {
            "q_vessel": "q_vessel_plus2",
            "source_mean": "source_mean_plus2",
            "q_tail": "q_tail_plus2",
            "ri_tail": "ri_tail_plus2",
        }[metric]
        sens.append(sensitivity_record(
            "background_plus2", "local_method_perturbation", metric, label,
            pct_delta(b[variant_col], b[metric]),
            "Both bilateral background strips moved outward by 2 A-lines; source/tail geometry unchanged.",
        ))

    g = base.merge(geo[keys + ["fixed128_valid", "q_vessel_fixed128", "source_mean_fixed128", "q_tail_fixed128", "ri_tail_fixed128"]], on=keys, validate="one_to_one")
    if len(g) != EXPECTED_VALID or not g["fixed128_valid"].astype(bool).all():
        raise AssertionError("Fixed-128 audit does not cover all Task1-valid frames")
    for metric, label in METRICS:
        variant_col = {
            "q_vessel": "q_vessel_fixed128",
            "source_mean": "source_mean_fixed128",
            "q_tail": "q_tail_fixed128",
            "ri_tail": "ri_tail_fixed128",
        }[metric]
        sens.append(sensitivity_record(
            "fixed128_width", "geometry_stress_test", metric, label,
            pct_delta(g[variant_col], g[metric]),
            "X4-centered lateral width replaced by known 128 um phantom diameter; background recomputed from variant edges.",
        ))

    # Subsampling summaries are already across all scan x phase records for each metric/stride.
    for stride, role in [(2, "sampling_design_250_positions"), (5, "sampling_design_100_positions")]:
        for metric, label in METRICS:
            row = sub.loc[(pd.to_numeric(sub["stride_frames"]) == stride) & (sub["metric"].eq(metric))]
            if len(row) != 1:
                raise AssertionError(f"Expected one subsampling summary row for stride={stride}, metric={metric}; got {len(row)}")
            r = row.iloc[0]
            sens.append({
                "method": f"systematic_stride_{stride}",
                "perturbation_role": role,
                "metric": metric,
                "metric_label": label,
                "n": int(r["n_phase_records"]),
                "abs_delta_pct_median": float(r["abs_error_pct_median"]),
                "abs_delta_pct_p95": float(r["abs_error_pct_p95"]),
                "abs_delta_pct_max": float(r["abs_error_pct_max"]),
                "note": f"Systematic phase-complete subsampling, approximately {int(round(500/stride))} planned positions per 500-frame volume.",
            })

    sens_df = pd.DataFrame(sens)
    sens_df.to_csv(out / "method_sensitivity_summary.csv", index=False)

    # Scan medians and observed contrasts vs flow01.
    scan_rows = []
    for scan, grp in base.groupby("scan_id", sort=False):
        flow = float(grp["flow_mm_s"].iloc[0])
        for metric, label in METRICS:
            scan_rows.append({
                "scan_id": scan,
                "flow_mm_s": flow,
                "metric": metric,
                "metric_label": label,
                "n_valid_frames": int(len(grp)),
                "scan_median": float(np.median(grp[metric].to_numpy(float))),
            })
    scan_df = pd.DataFrame(scan_rows).sort_values(["metric", "flow_mm_s"])
    scan_df.to_csv(out / "scan_median_reference.csv", index=False)

    ref = scan_df.loc[scan_df.flow_mm_s.eq(1.0), ["metric", "scan_median"]].rename(columns={"scan_median":"flow01_median"})
    obs = scan_df.merge(ref, on="metric", validate="many_to_one")
    obs = obs.loc[~obs.flow_mm_s.eq(1.0)].copy()
    obs["observed_delta_pct_vs_flow01"] = 100.0 * (obs["scan_median"] - obs["flow01_median"]) / obs["flow01_median"]
    obs["observed_abs_delta_pct"] = obs["observed_delta_pct_vs_flow01"].abs()

    piv = sens_df.pivot(index="metric", columns="method", values=["abs_delta_pct_median", "abs_delta_pct_p95", "abs_delta_pct_max"])
    piv.columns = [f"{m}_{stat}" for stat, m in piv.columns]
    piv = piv.reset_index()
    cmp = obs.merge(piv, on="metric", validate="many_to_one")

    local_p95_cols = ["background_plus2_abs_delta_pct_p95", "systematic_stride_5_abs_delta_pct_p95"]
    local_med_cols = ["background_plus2_abs_delta_pct_median", "systematic_stride_5_abs_delta_pct_median"]
    cmp["local_method_p95_reference_pct"] = cmp[local_p95_cols].max(axis=1)
    cmp["local_method_typical_reference_pct"] = cmp[local_med_cols].max(axis=1)
    cmp["geometry_stress_p95_reference_pct"] = cmp["fixed128_width_abs_delta_pct_p95"]
    cmp["observed_to_local_p95_ratio"] = cmp["observed_abs_delta_pct"] / cmp["local_method_p95_reference_pct"]
    cmp["observed_to_geometry_p95_ratio"] = cmp["observed_abs_delta_pct"] / cmp["geometry_stress_p95_reference_pct"]

    def classify(r):
        x = float(r.observed_abs_delta_pct)
        lp = float(r.local_method_p95_reference_pct)
        gp = float(r.geometry_stress_p95_reference_pct)
        if x <= lp:
            return "within_local_method_sensitivity"
        if x <= gp:
            return "above_local_but_within_geometry_stress"
        return "above_local_and_geometry_stress"

    cmp["descriptive_context_class"] = cmp.apply(classify, axis=1)
    cmp["interpretation_boundary"] = "descriptive_scan_level_only; requires independent-volume replication"
    cmp.to_csv(out / "observed_effect_vs_method_sensitivity.csv", index=False)

    # Compact flow07 focus because it was the largest current volume-level deviation.
    flow07 = cmp.loc[cmp.flow_mm_s.eq(7.0), [
        "metric", "metric_label", "observed_delta_pct_vs_flow01", "observed_abs_delta_pct",
        "background_plus2_abs_delta_pct_median", "background_plus2_abs_delta_pct_p95",
        "fixed128_width_abs_delta_pct_median", "fixed128_width_abs_delta_pct_p95",
        "systematic_stride_5_abs_delta_pct_median", "systematic_stride_5_abs_delta_pct_p95",
        "systematic_stride_2_abs_delta_pct_median", "systematic_stride_2_abs_delta_pct_p95",
        "observed_to_local_p95_ratio", "observed_to_geometry_p95_ratio", "descriptive_context_class",
    ]].copy()
    flow07.to_csv(out / "flow07_context_summary.csv", index=False)

    validation = {
        "task1_valid_frames": int(len(t1)),
        "background_plus2_joined_frames": int(len(b)),
        "fixed128_joined_frames": int(len(g)),
        "metrics": [m for m, _ in METRICS],
        "observed_reference_scan": "flow01 / 1 mm/s",
        "observed_comparison_scans": ["flow03", "flow05", "flow07", "flow10"],
        "method_sensitivity_sources_are_not_combined_into_single_error": True,
        "background_role": "local_method_perturbation",
        "fixed128_role": "geometry_stress_test_not_random_error",
        "subsampling_role": "sampling_design_approximation",
        "p_values_computed": False,
        "inferential_test_performed": False,
        "bscan_independent_replicate_claim": False,
        "independent_volume_replication_required_for_flow_effect_claim": True,
    }
    (out / "validation.json").write_text(json.dumps(validation, indent=2) + "\n", encoding="utf-8")

    inputs = []
    for p in [t1_path, bg_path, geo_path, sub_path]:
        inputs.append({"file": str(p.relative_to(root)).replace("\\", "/"), "bytes": p.stat().st_size, "sha256": sha256_file(p)})
    pd.DataFrame(inputs).to_csv(out / "input_sha256.csv", index=False)

    provenance = {
        "analysis": "Observed volume-level flow contrast vs audited method-sensitivity context",
        "source_sha": args.source_sha,
        "workflow_source_sha": args.workflow_source_sha,
        "selection": "Task1 valid frames only (2422)",
        "observed_effect": "difference in per-volume medians relative to flow01",
        "no_composite_uncertainty": True,
        "reason_no_composite_uncertainty": "background, geometry, and sampling perturbations are not independent estimates of one random error distribution",
        "statistical_unit_note": "Each flow currently has one independent scan volume; B-scans are spatial positions, not independent experiments.",
    }
    (out / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")

    def f(x): return "NA" if not np.isfinite(float(x)) else f"{float(x):.2f}%"
    lines = [
        "# Observed flow contrasts vs method-sensitivity context",
        "",
        "This audit puts the current **volume-level** flow contrasts on the same percent scale as three audited perturbations. It does **not** estimate a single measurement-error SD and does not establish a flow-speed effect.",
        "",
        "- `background +2`: local method perturbation.",
        "- `fixed-128 width`: geometry **stress test**, not a random uncertainty estimate.",
        "- `systematic subsampling`: sampling-design approximation; 100-position (stride 5) is used as the local planning reference, while 250-position results are also retained.",
        "",
        "## Method sensitivity across all valid frames / scan phases",
        "",
        "| metric | background +2 median / P95 | fixed-128 median / P95 | ~100 positions median / P95 | ~250 positions median / P95 |",
        "|---|---:|---:|---:|---:|",
    ]
    for metric, label in METRICS:
        x = sens_df.loc[sens_df.metric.eq(metric)].set_index("method")
        lines.append(
            f"| {label} | {f(x.loc['background_plus2','abs_delta_pct_median'])} / {f(x.loc['background_plus2','abs_delta_pct_p95'])} | "
            f"{f(x.loc['fixed128_width','abs_delta_pct_median'])} / {f(x.loc['fixed128_width','abs_delta_pct_p95'])} | "
            f"{f(x.loc['systematic_stride_5','abs_delta_pct_median'])} / {f(x.loc['systematic_stride_5','abs_delta_pct_p95'])} | "
            f"{f(x.loc['systematic_stride_2','abs_delta_pct_median'])} / {f(x.loc['systematic_stride_2','abs_delta_pct_p95'])} |"
        )
    lines += [
        "",
        "## Current volume-level contrasts relative to 1 mm/s",
        "",
        "| flow | metric | observed median change | local-method P95 reference | geometry-stress P95 | context |",
        "|---:|---|---:|---:|---:|---|",
    ]
    for r in cmp.sort_values(["flow_mm_s", "metric"]).itertuples(index=False):
        lines.append(f"| {r.flow_mm_s:g} | {r.metric_label} | {r.observed_delta_pct_vs_flow01:+.2f}% | {r.local_method_p95_reference_pct:.2f}% | {r.geometry_stress_p95_reference_pct:.2f}% | {r.descriptive_context_class} |")
    lines += [
        "",
        "## Interpretation",
        "",
        "The classification is deliberately descriptive. `within_local_method_sensitivity` means the observed absolute scan-median contrast does not exceed the larger P95 of background +2 and ~100-position subsampling. `above_local_but_within_geometry_stress` means it exceeds those local perturbations but remains within the P95 of the stronger fixed-128 geometry stress test. `above_local_and_geometry_stress` means it exceeds both reference scales.",
        "",
        "None of these labels is a p-value or a causal flow-speed claim. Each flow still has one independent scan volume, so independent-volume replication is required.",
        "",
        "## Outputs",
        "",
        "- `method_sensitivity_summary.csv`",
        "- `observed_effect_vs_method_sensitivity.csv`",
        "- `flow07_context_summary.csv`",
        "- `scan_median_reference.csv`",
        "- `input_sha256.csv`",
        "- `validation.json`",
        "- `provenance.json`",
        "- `analyze_method_sensitivity_effect_context.py`",
    ]
    (out / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(json.dumps({
        "valid_frames": len(t1),
        "method_sensitivity_rows": len(sens_df),
        "observed_contrast_rows": len(cmp),
        "flow07": flow07[["metric_label", "observed_delta_pct_vs_flow01", "descriptive_context_class"]].to_dict("records"),
        "output_dir": str(out.relative_to(root)),
    }, indent=2))


if __name__ == "__main__":
    main()
