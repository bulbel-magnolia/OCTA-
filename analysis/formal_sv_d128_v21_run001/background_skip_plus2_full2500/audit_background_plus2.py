#!/usr/bin/env python3
"""Audit sensitivity to shifting bilateral SV background strips 2 A-lines outward.

Frozen run001 remains unchanged. Baseline is skip=3, width=5; the existing
sensitivity variant background_skip_plus_2px is skip=5, width=5. Only the
background location changes. Spatial B-scans are descriptive samples, so this
script reports robust descriptive sensitivity and no inferential p-values.
"""
from __future__ import annotations

import csv
import math
import re
from pathlib import Path
from statistics import median

ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "results" / "formal_sv_d128_v21_full2500_run001"
OUT = Path(__file__).resolve().parent
BASE = RUN / "frame_results.csv"
SENS = RUN / "sensitivity_results.csv"
VARIANT = "background_skip_plus_2px"
DX_UM = 12.7
SHIFT_ALINES = 2
SHIFT_UM = DX_UM * SHIFT_ALINES


def b(v: str) -> bool:
    return str(v).strip().lower() in {"true", "1", "yes"}


def f(v: str):
    try:
        x = float(v)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def frame_index(row: dict[str, str]):
    for key in ("frame_index_0based", "bscan_index"):
        if row.get(key, "") != "":
            x = f(row[key])
            if x is not None:
                return int(round(x))
    m = re.search(r"Bscan(\d+)", row.get("frame_id", ""), re.I)
    if m:
        return int(m.group(1))
    raise ValueError(f"Cannot recover frame index: {row}")


def pct(new, old):
    if new is None or old is None or old == 0:
        return None
    return 100.0 * (new - old) / old


def delta(new, old):
    if new is None or old is None:
        return None
    return new - old


def q(values, p):
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


def fmt(x, nd=3):
    return "NA" if x is None else f"{x:.{nd}f}"


def metric_stats(values, thresholds=(5.0, 10.0, 20.0)):
    a = [x for x in values if x is not None and math.isfinite(x)]
    aa = [abs(x) for x in a]
    out = {
        "n_metric": len(a),
        "signed_median": q(a, .50),
        "signed_q25": q(a, .25),
        "signed_q75": q(a, .75),
        "abs_median": q(aa, .50),
        "abs_p90": q(aa, .90),
        "abs_p95": q(aa, .95),
        "abs_max": max(aa) if aa else None,
        "positive_n": sum(x > 0 for x in a),
        "negative_n": sum(x < 0 for x in a),
        "zero_n": sum(x == 0 for x in a),
    }
    for t in thresholds:
        out[f"abs_gt_{int(t)}_n"] = sum(x > t for x in aa)
        out[f"abs_gt_{int(t)}_pct"] = 100.0 * out[f"abs_gt_{int(t)}_n"] / len(aa) if aa else None
    return out


def read_csv(path):
    with path.open(newline="", encoding="utf-8-sig") as h:
        return list(csv.DictReader(h))


base_rows = read_csv(BASE)
sens_rows = [r for r in read_csv(SENS) if r.get("variant") == VARIANT]
if len(base_rows) != 2500:
    raise RuntimeError(f"Expected 2500 baseline rows, found {len(base_rows)}")

sens_by_key = {}
for r in sens_rows:
    k = (r["scan_id"], frame_index(r))
    if k in sens_by_key:
        raise RuntimeError(f"Duplicate sensitivity row: {k}")
    sens_by_key[k] = r

framewise = []
for br in base_rows:
    k = (br["scan_id"], frame_index(br))
    sr = sens_by_key.get(k)
    baseline_valid = b(br.get("valid", ""))
    variant_valid = bool(sr) and b(sr.get("valid", ""))
    if not baseline_valid:
        continue
    if sr is None:
        raise RuntimeError(f"Missing +2 sensitivity row for baseline-valid frame {k}")

    qv0, qt0, r0 = f(br.get("q_vessel")), f(br.get("q_tail")), f(br.get("ratio_tail_to_vessel"))
    qv1, qt1, r1 = f(sr.get("q_vessel")), f(sr.get("q_tail")), f(sr.get("ratio_tail_to_vessel"))
    if not variant_valid:
        qv1 = qt1 = r1 = None
    rec = {
        "scan_id": br["scan_id"],
        "frame_index_0based": frame_index(br),
        "flow_speed_mm_s": br.get("flow_speed_mm_s", ""),
        "localization_source": br.get("localization_source", ""),
        "baseline_valid": baseline_valid,
        "variant_valid": variant_valid,
        "variant_invalid_reason": "" if variant_valid else sr.get("invalid_reason", ""),
        "q_vessel_base": qv0,
        "q_vessel_plus2": qv1,
        "q_vessel_delta": delta(qv1, qv0),
        "q_vessel_delta_pct": pct(qv1, qv0),
        "q_tail_base": qt0,
        "q_tail_plus2": qt1,
        "q_tail_delta": delta(qt1, qt0),
        "q_tail_delta_pct": pct(qt1, qt0),
        "ratio_base": r0,
        "ratio_plus2": r1,
        "ratio_delta": delta(r1, r0),
        "ratio_delta_pct": pct(r1, r0),
    }
    framewise.append(rec)

if len(framewise) != sum(b(r.get("valid", "")) for r in base_rows):
    raise RuntimeError("Framewise audit row count does not match baseline-valid count")

fields = list(framewise[0].keys())
with (OUT / "background_plus2_framewise.csv").open("w", newline="", encoding="utf-8") as h:
    w = csv.DictWriter(h, fieldnames=fields)
    w.writeheader(); w.writerows(framewise)

scopes = [("overall", "all", framewise)]
for scan in sorted({r["scan_id"] for r in framewise}):
    scopes.append(("flow", scan, [r for r in framewise if r["scan_id"] == scan]))
for src in sorted({r["localization_source"] for r in framewise}):
    scopes.append(("localization_source", src, [r for r in framewise if r["localization_source"] == src]))

metric_cols = [
    ("q_vessel_delta_pct", "percent"),
    ("q_tail_delta_pct", "percent"),
    ("ratio_delta", "absolute_ratio"),
    ("ratio_delta_pct", "percent"),
]
summary = []
for scope_type, scope_value, rows in scopes:
    n_base = len(rows)
    n_var = sum(r["variant_valid"] for r in rows)
    for col, unit in metric_cols:
        st = metric_stats([r[col] for r in rows if r["variant_valid"]], thresholds=(5, 10, 20) if unit == "percent" else (.02, .05, .10))
        rec = {
            "scope_type": scope_type, "scope_value": scope_value,
            "baseline_valid_n": n_base, "variant_valid_n": n_var,
            "variant_invalid_n": n_base - n_var,
            "variant_valid_pct": 100.0 * n_var / n_base if n_base else None,
            "metric": col, "unit": unit,
        }
        rec.update(st)
        summary.append(rec)

sfields = list(summary[0].keys())
with (OUT / "background_plus2_summary.csv").open("w", newline="", encoding="utf-8") as h:
    w = csv.DictWriter(h, fieldnames=sfields)
    w.writeheader(); w.writerows(summary)

extreme_pool = [r for r in framewise if r["variant_valid"] and r["q_tail_delta_pct"] is not None]
extremes = sorted(extreme_pool, key=lambda r: abs(r["q_tail_delta_pct"]), reverse=True)[:50]
with (OUT / "background_plus2_extremes.csv").open("w", newline="", encoding="utf-8") as h:
    w = csv.DictWriter(h, fieldnames=fields)
    w.writeheader(); w.writerows(extremes)

# Compact machine-readable metadata.
with (OUT / "audit_metadata.csv").open("w", newline="", encoding="utf-8") as h:
    w = csv.writer(h)
    w.writerow(["key", "value"])
    w.writerow(["baseline_background_skip_alines", 3])
    w.writerow(["variant_background_skip_alines", 5])
    w.writerow(["strip_width_alines", 5])
    w.writerow(["outward_shift_alines", SHIFT_ALINES])
    w.writerow(["outward_shift_um", SHIFT_UM])
    w.writerow(["baseline_valid_frames", len(framewise)])
    w.writerow(["variant_valid_frames", sum(r["variant_valid"] for r in framewise)])

by_metric = {(r["scope_type"], r["scope_value"], r["metric"]): r for r in summary}
ov_qt = by_metric[("overall", "all", "q_tail_delta_pct")]
ov_qv = by_metric[("overall", "all", "q_vessel_delta_pct")]
ov_rd = by_metric[("overall", "all", "ratio_delta")]
flow_qt = [by_metric[("flow", s, "q_tail_delta_pct")] for s in sorted({r["scan_id"] for r in framewise})]
worst_flow = max(flow_qt, key=lambda r: r["abs_median"] if r["abs_median"] is not None else -1)

readme = f"""# Background +2 A-line full-volume sensitivity audit\n\n## Purpose\n\nThis is a **QC sensitivity audit**, not a replacement of the frozen run001 background. Baseline uses bilateral background strips after skipping 3 A-lines from each vessel edge; the variant skips 5 A-lines while keeping the strip width at 5 A-lines. Thus the sampled background moves {SHIFT_ALINES} A-lines = {SHIFT_UM:.1f} μm farther from the vessel. All source/tail geometry, SV definition, 500 μm tail window, signed residual handling and inclusion rules remain unchanged.\n\nSpatial B-scans within one scan volume are not independent experiments. This audit is descriptive only; no p-values are produced.\n\n## Headline results\n\n- Baseline-valid frames audited: **{len(framewise)}**.\n- +2 background variant valid: **{sum(r['variant_valid'] for r in framewise)} / {len(framewise)} ({100*sum(r['variant_valid'] for r in framewise)/len(framewise):.2f}%)**.\n- |ΔQ_tail| median: **{fmt(ov_qt['abs_median'])}%**; P95: **{fmt(ov_qt['abs_p95'])}%**; maximum: **{fmt(ov_qt['abs_max'])}%**.\n- |ΔQ_vessel| median: **{fmt(ov_qv['abs_median'])}%**; P95: **{fmt(ov_qv['abs_p95'])}%**.\n- |ΔR| median: **{fmt(ov_rd['abs_median'],4)}**; P95: **{fmt(ov_rd['abs_p95'],4)}**.\n- Flow with largest median |ΔQ_tail|: **{worst_flow['scope_value']}**, {fmt(worst_flow['abs_median'])}%.\n- Frames with |ΔQ_tail| > 5% / 10% / 20%: **{ov_qt['abs_gt_5_n']} / {ov_qt['abs_gt_10_n']} / {ov_qt['abs_gt_20_n']}**.\n\nA large relative ΔQ_tail can occur when baseline signed Q_tail is close to zero. For that reason `background_plus2_framewise.csv` retains the absolute Q_tail values and deltas, and `background_plus2_extremes.csv` is a diagnostic list rather than an exclusion list.\n\n## Files\n\n- `background_plus2_framewise.csv`: one row per baseline-valid frame.\n- `background_plus2_summary.csv`: overall, per-flow, and localization-source robust summaries.\n- `background_plus2_extremes.csv`: top 50 frames by absolute relative Q_tail change.\n- `audit_metadata.csv`: frozen perturbation definition and row counts.\n- `audit_background_plus2.py`: reproducible audit script.\n\n## Interpretation rule fixed before reading results\n\nThis audit asks whether moving a reasonable local background farther from the vessel materially changes Q_vessel, Q_tail or R. It must **not** be used to choose whichever background gives a more monotonic or statistically favorable flow pattern. The frozen run001 files under `results/formal_sv_d128_v21_full2500_run001/` are not modified.\n"""
(OUT / "README.md").write_text(readme, encoding="utf-8")
print(readme)
