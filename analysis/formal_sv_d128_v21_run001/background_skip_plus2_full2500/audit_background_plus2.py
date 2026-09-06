#!/usr/bin/env python3
"""Full-volume sensitivity audit for moving bilateral SV background +2 A-lines.

Frozen run001 is not modified. Baseline: skip=3, width=5 each side. Variant:
skip=5, width=5. Only background location changes. Results are descriptive
because B-scans within one volume are spatial samples, not independent experiments.
"""
from __future__ import annotations

import csv
import math
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "results" / "formal_sv_d128_v21_full2500_run001"
OUT = Path(__file__).resolve().parent
BASE = RUN / "frame_results.csv"
SENS = RUN / "sensitivity_results.csv"
VARIANT = "background_skip_plus_2px"
SHIFT_ALINES = 2
SHIFT_UM = 25.4


def as_bool(v):
    return str(v).strip().lower() in {"true", "1", "yes"}


def num(v):
    try:
        x = float(v)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def get_frame_index(row):
    for key in ("frame_index_0based", "bscan_index"):
        if row.get(key, "") != "":
            x = num(row[key])
            if x is not None:
                return int(round(x))
    m = re.search(r"Bscan(\d+)", row.get("frame_id", ""), re.I)
    if m:
        return int(m.group(1))
    raise ValueError(f"cannot recover frame index from row: {row}")


def read_csv(path):
    with path.open(newline="", encoding="utf-8-sig") as h:
        return list(csv.DictReader(h))


def delta(new, old):
    return None if new is None or old is None else new - old


def pct(new, old):
    return None if new is None or old is None or old == 0 else 100.0 * (new - old) / old


def quantile(values, p):
    a = sorted(x for x in values if x is not None and math.isfinite(x))
    if not a:
        return None
    if len(a) == 1:
        return a[0]
    h = (len(a) - 1) * p
    lo, hi = int(math.floor(h)), int(math.ceil(h))
    if lo == hi:
        return a[lo]
    return a[lo] * (hi - h) + a[hi] * (h - lo)


def metric_stats(values, unit):
    a = [x for x in values if x is not None and math.isfinite(x)]
    aa = [abs(x) for x in a]
    out = {
        "n_metric": len(a),
        "signed_median": quantile(a, .50),
        "signed_q25": quantile(a, .25),
        "signed_q75": quantile(a, .75),
        "abs_median": quantile(aa, .50),
        "abs_p90": quantile(aa, .90),
        "abs_p95": quantile(aa, .95),
        "abs_max": max(aa) if aa else None,
        "positive_n": sum(x > 0 for x in a),
        "negative_n": sum(x < 0 for x in a),
        "zero_n": sum(x == 0 for x in a),
        "abs_gt_5_n": None, "abs_gt_5_pct": None,
        "abs_gt_10_n": None, "abs_gt_10_pct": None,
        "abs_gt_20_n": None, "abs_gt_20_pct": None,
        "abs_gt_0p02_n": None, "abs_gt_0p02_pct": None,
        "abs_gt_0p05_n": None, "abs_gt_0p05_pct": None,
        "abs_gt_0p10_n": None, "abs_gt_0p10_pct": None,
    }
    thresholds = ((5, "5"), (10, "10"), (20, "20")) if unit == "percent" else ((.02, "0p02"), (.05, "0p05"), (.10, "0p10"))
    for t, tag in thresholds:
        n = sum(x > t for x in aa)
        out[f"abs_gt_{tag}_n"] = n
        out[f"abs_gt_{tag}_pct"] = 100.0 * n / len(aa) if aa else None
    return out


def fmt(x, digits=3):
    return "NA" if x is None else f"{x:.{digits}f}"


base_rows = read_csv(BASE)
all_sens = read_csv(SENS)
sens_rows = [r for r in all_sens if r.get("variant") == VARIANT]
if len(base_rows) != 2500:
    raise RuntimeError(f"expected 2500 baseline rows, got {len(base_rows)}")

sens_by_key = {}
for row in sens_rows:
    key = (row["scan_id"], get_frame_index(row))
    if key in sens_by_key:
        raise RuntimeError(f"duplicate +2 row: {key}")
    sens_by_key[key] = row

framewise = []
for br in base_rows:
    if not as_bool(br.get("valid", "")):
        continue
    key = (br["scan_id"], get_frame_index(br))
    sr = sens_by_key.get(key)
    if sr is None:
        raise RuntimeError(f"missing +2 row for baseline-valid frame {key}")
    variant_valid = as_bool(sr.get("valid", ""))
    qv0, qt0, r0 = num(br.get("q_vessel")), num(br.get("q_tail")), num(br.get("ratio_tail_to_vessel"))
    qv1, qt1, r1 = num(sr.get("q_vessel")), num(sr.get("q_tail")), num(sr.get("ratio_tail_to_vessel"))
    if not variant_valid:
        qv1 = qt1 = r1 = None
    framewise.append({
        "scan_id": br["scan_id"],
        "frame_index_0based": get_frame_index(br),
        "flow_speed_mm_s": br.get("flow_speed_mm_s", ""),
        "localization_source": br.get("localization_source", ""),
        "baseline_valid": True,
        "variant_valid": variant_valid,
        "variant_invalid_reason": "" if variant_valid else sr.get("invalid_reason", ""),
        "q_vessel_base": qv0, "q_vessel_plus2": qv1,
        "q_vessel_delta": delta(qv1, qv0), "q_vessel_delta_pct": pct(qv1, qv0),
        "q_tail_base": qt0, "q_tail_plus2": qt1,
        "q_tail_delta": delta(qt1, qt0), "q_tail_delta_pct": pct(qt1, qt0),
        "ratio_base": r0, "ratio_plus2": r1,
        "ratio_delta": delta(r1, r0), "ratio_delta_pct": pct(r1, r0),
    })

expected_valid = sum(as_bool(r.get("valid", "")) for r in base_rows)
if len(framewise) != expected_valid:
    raise RuntimeError(f"baseline-valid mismatch: {len(framewise)} vs {expected_valid}")

frame_fields = list(framewise[0])
with (OUT / "background_plus2_framewise.csv").open("w", newline="", encoding="utf-8") as h:
    w = csv.DictWriter(h, fieldnames=frame_fields); w.writeheader(); w.writerows(framewise)

scopes = [("overall", "all", framewise)]
for scan in sorted({r["scan_id"] for r in framewise}):
    scopes.append(("flow", scan, [r for r in framewise if r["scan_id"] == scan]))
for src in sorted({r["localization_source"] for r in framewise}):
    scopes.append(("localization_source", src, [r for r in framewise if r["localization_source"] == src]))

metric_specs = [
    ("q_vessel_delta_pct", "percent"),
    ("q_tail_delta_pct", "percent"),
    ("ratio_delta", "absolute_ratio"),
    ("ratio_delta_pct", "percent"),
]
summary = []
for scope_type, scope_value, rows in scopes:
    n_base = len(rows)
    n_var = sum(r["variant_valid"] for r in rows)
    for metric, unit in metric_specs:
        rec = {
            "scope_type": scope_type, "scope_value": scope_value,
            "baseline_valid_n": n_base, "variant_valid_n": n_var,
            "variant_invalid_n": n_base - n_var,
            "variant_valid_pct": 100.0 * n_var / n_base if n_base else None,
            "metric": metric, "unit": unit,
        }
        rec.update(metric_stats([r[metric] for r in rows if r["variant_valid"]], unit))
        summary.append(rec)

summary_fields = list(summary[0])
with (OUT / "background_plus2_summary.csv").open("w", newline="", encoding="utf-8") as h:
    w = csv.DictWriter(h, fieldnames=summary_fields); w.writeheader(); w.writerows(summary)

extremes = sorted(
    [r for r in framewise if r["variant_valid"] and r["q_tail_delta_pct"] is not None],
    key=lambda r: abs(r["q_tail_delta_pct"]), reverse=True,
)[:50]
with (OUT / "background_plus2_extremes.csv").open("w", newline="", encoding="utf-8") as h:
    w = csv.DictWriter(h, fieldnames=frame_fields); w.writeheader(); w.writerows(extremes)

with (OUT / "audit_metadata.csv").open("w", newline="", encoding="utf-8") as h:
    w = csv.writer(h); w.writerow(["key", "value"])
    for k, v in [
        ("baseline_background_skip_alines", 3), ("variant_background_skip_alines", 5),
        ("strip_width_alines", 5), ("outward_shift_alines", SHIFT_ALINES),
        ("outward_shift_um", SHIFT_UM), ("baseline_valid_frames", len(framewise)),
        ("variant_valid_frames", sum(r["variant_valid"] for r in framewise)),
    ]:
        w.writerow([k, v])

by = {(r["scope_type"], r["scope_value"], r["metric"]): r for r in summary}
ov_qt = by[("overall", "all", "q_tail_delta_pct")]
ov_qv = by[("overall", "all", "q_vessel_delta_pct")]
ov_rd = by[("overall", "all", "ratio_delta")]
flow_qt = [by[("flow", s, "q_tail_delta_pct")] for s in sorted({r["scan_id"] for r in framewise})]
worst_flow = max(flow_qt, key=lambda r: r["abs_median"] if r["abs_median"] is not None else -1)
variant_n = sum(r["variant_valid"] for r in framewise)

readme = f"""# Background +2 A-line full-volume sensitivity audit

## What was changed

This is a **QC sensitivity audit**, not a replacement of frozen run001. Baseline background skips 3 A-lines from each vessel edge and then takes 5 A-lines per side. The variant skips 5 A-lines and still takes 5 per side, so both background strips move **{SHIFT_ALINES} A-lines = {SHIFT_UM:.1f} μm outward**. Source/tail geometry, SV definition, 500 μm tail window, signed residual handling and inclusion rules are unchanged.

Spatial B-scans in one volume are not independent experiments; this audit is descriptive and reports no p-values.

## Headline results

- Baseline-valid frames audited: **{len(framewise)}**.
- +2 variant valid: **{variant_n}/{len(framewise)} ({100*variant_n/len(framewise):.2f}%)**.
- |ΔQ_tail|: median **{fmt(ov_qt['abs_median'])}%**, P95 **{fmt(ov_qt['abs_p95'])}%**, max **{fmt(ov_qt['abs_max'])}%**.
- |ΔQ_vessel|: median **{fmt(ov_qv['abs_median'])}%**, P95 **{fmt(ov_qv['abs_p95'])}%**.
- |ΔR|: median **{fmt(ov_rd['abs_median'],4)}**, P95 **{fmt(ov_rd['abs_p95'],4)}**.
- Largest per-flow median |ΔQ_tail|: **{worst_flow['scope_value']} = {fmt(worst_flow['abs_median'])}%**.
- Frames with |ΔQ_tail| >5% / >10% / >20%: **{ov_qt['abs_gt_5_n']} / {ov_qt['abs_gt_10_n']} / {ov_qt['abs_gt_20_n']}**.

Large relative ΔQ_tail can occur when signed baseline Q_tail is close to zero. Therefore the framewise file retains absolute values/deltas, and the extremes file is diagnostic rather than an exclusion list.

## Files

- `background_plus2_framewise.csv`: all baseline-valid frames.
- `background_plus2_summary.csv`: overall, per-flow and localization-source robust summaries.
- `background_plus2_extremes.csv`: top 50 by |relative ΔQ_tail|.
- `audit_metadata.csv`: perturbation definition and row counts.
- `audit_background_plus2.py`: reproducible script.

## Interpretation rule

This audit tests whether a reasonable outward shift of local background materially changes Q_vessel, Q_tail or R. It must not be used to select whichever background gives a more monotonic or statistically favorable flow trend. The frozen directory `results/formal_sv_d128_v21_full2500_run001/` is untouched.
"""
(OUT / "README.md").write_text(readme, encoding="utf-8")
print(readme)
