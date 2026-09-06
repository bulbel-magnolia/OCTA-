from __future__ import annotations

import csv
import hashlib
import math
import os
from pathlib import Path
from statistics import mean

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
DX_UM = 12.7
DIAMETER_UM = 128.0
SCANS = ("flow01", "flow03", "flow05", "flow07", "flow10")
BASE = ROOT / "results" / "formal_sv_d128_v21_full2500_run001" / "tracking"
INPUTS = {scan: BASE / scan / f"{scan}_mentor_tracking.csv" for scan in SCANS}


def f(value: object) -> float:
    try:
        x = float(value)
    except (TypeError, ValueError):
        return math.nan
    return x if math.isfinite(x) else math.nan


def b(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1"}


def q(values: list[float], p: float) -> float:
    vals = sorted(x for x in values if math.isfinite(x))
    if not vals:
        return math.nan
    if len(vals) == 1:
        return vals[0]
    pos = (len(vals) - 1) * p
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return vals[lo]
    w = pos - lo
    return vals[lo] * (1.0 - w) + vals[hi] * w


def sd(values: list[float]) -> float:
    vals = [x for x in values if math.isfinite(x)]
    if len(vals) < 2:
        return math.nan
    m = mean(vals)
    return math.sqrt(sum((x - m) ** 2 for x in vals) / (len(vals) - 1))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def flow_value(scan_id: str) -> int:
    return int(scan_id.replace("flow", ""))


def fmt(x: float, digits: int = 4) -> str:
    return "NA" if not math.isfinite(x) else f"{x:.{digits}f}"


def summary(rows: list[dict[str, object]], scope: str, scan_id: str = "ALL") -> dict[str, object]:
    paired = [r for r in rows if bool(r["paired"])]
    deltas = [float(r["delta_px"]) for r in paired]
    absd = [float(r["abs_delta_px"]) for r in paired]
    n = len(paired)
    return {
        "scope": scope,
        "scan_id": scan_id,
        "n_frames": len(rows),
        "n_paired": n,
        "paired_fraction": n / len(rows) if rows else math.nan,
        "mean_delta_px": mean(deltas) if deltas else math.nan,
        "median_delta_px": q(deltas, 0.50),
        "sd_delta_px": sd(deltas),
        "q25_delta_px": q(deltas, 0.25),
        "q75_delta_px": q(deltas, 0.75),
        "mean_abs_delta_px": mean(absd) if absd else math.nan,
        "median_abs_delta_px": q(absd, 0.50),
        "p90_abs_delta_px": q(absd, 0.90),
        "p95_abs_delta_px": q(absd, 0.95),
        "p99_abs_delta_px": q(absd, 0.99),
        "max_abs_delta_px": max(absd) if absd else math.nan,
        "median_abs_delta_um": q(absd, 0.50) * DX_UM if absd else math.nan,
        "p95_abs_delta_um": q(absd, 0.95) * DX_UM if absd else math.nan,
        "max_abs_delta_um": max(absd) * DX_UM if absd else math.nan,
        "fraction_abs_le_0p25px": sum(x <= 0.25 for x in absd) / n if n else math.nan,
        "fraction_abs_le_0p5px": sum(x <= 0.50 for x in absd) / n if n else math.nan,
        "fraction_abs_le_1px": sum(x <= 1.0 for x in absd) / n if n else math.nan,
        "fraction_abs_gt_1px": sum(x > 1.0 for x in absd) / n if n else math.nan,
        "fraction_abs_gt_2px": sum(x > 2.0 for x in absd) / n if n else math.nan,
        "n_x1_fallback": sum(bool(r["x1_fallback"]) for r in rows),
        "n_x4_jump_corrected": sum(bool(r["x4_jump_corrected"]) for r in rows),
    }


def smoothness(rows: list[dict[str, object]], scan_id: str) -> list[dict[str, object]]:
    ordered = sorted(rows, key=lambda r: int(r["frame_index"]))
    out: list[dict[str, object]] = []
    for key, label in (("x1_midpoint_px", "X1_midpoint"), ("x4_px", "X4")):
        jumps: list[float] = []
        for left, right in zip(ordered, ordered[1:]):
            if int(right["frame_index"]) != int(left["frame_index"]) + 1:
                continue
            a = float(left[key])
            c = float(right[key])
            if math.isfinite(a) and math.isfinite(c):
                jumps.append(abs(c - a))
        out.append({
            "scan_id": scan_id,
            "trajectory": label,
            "n_adjacent_pairs": len(jumps),
            "mean_abs_step_px": mean(jumps) if jumps else math.nan,
            "median_abs_step_px": q(jumps, 0.50),
            "p95_abs_step_px": q(jumps, 0.95),
            "p99_abs_step_px": q(jumps, 0.99),
            "max_abs_step_px": max(jumps) if jumps else math.nan,
            "fraction_step_gt_1px": sum(x > 1.0 for x in jumps) / len(jumps) if jumps else math.nan,
            "fraction_step_gt_2px": sum(x > 2.0 for x in jumps) / len(jumps) if jumps else math.nan,
        })
    return out


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows and not fieldnames:
        raise RuntimeError(f"cannot infer columns for empty CSV: {path}")
    names = fieldnames or list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=names, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    frame_rows: list[dict[str, object]] = []
    metadata: list[dict[str, object]] = []

    for scan in SCANS:
        path = INPUTS[scan]
        if not path.is_file():
            raise FileNotFoundError(path)
        with path.open("r", encoding="utf-8", newline="") as fh:
            reader = csv.DictReader(fh)
            required = {
                "scan_id", "frame_index", "x1_local_geometry_px",
                "x4_centroid_isolated_jump_corrected_px", "x1_fallback",
                "x4_jump_corrected", "vessel_presence_prediction",
                "valid_local_body", "local_body_run_width_px",
                "tracking_class", "tracking_confidence", "assessability_score",
            }
            missing = sorted(required - set(reader.fieldnames or []))
            if missing:
                raise RuntimeError(f"{path} missing columns: {missing}")
            count = 0
            for raw in reader:
                count += 1
                x1 = f(raw.get("x1_local_geometry_px"))
                x4 = f(raw.get("x4_centroid_isolated_jump_corrected_px"))
                paired = math.isfinite(x1) and math.isfinite(x4)
                delta = x4 - x1 if paired else math.nan
                width = f(raw.get("local_body_run_width_px"))
                frame_rows.append({
                    "scan_id": scan,
                    "flow_mm_s": flow_value(scan),
                    "frame_index": int(float(raw["frame_index"])),
                    "x1_midpoint_px": x1,
                    "x4_px": x4,
                    "paired": paired,
                    "delta_px": delta,
                    "abs_delta_px": abs(delta) if paired else math.nan,
                    "delta_um": delta * DX_UM if paired else math.nan,
                    "abs_delta_um": abs(delta) * DX_UM if paired else math.nan,
                    "abs_delta_fraction_diameter": abs(delta) * DX_UM / DIAMETER_UM if paired else math.nan,
                    "local_body_run_width_px": width,
                    "abs_delta_fraction_local_width": abs(delta) / width if paired and math.isfinite(width) and width > 0 else math.nan,
                    "vessel_presence_prediction": raw.get("vessel_presence_prediction", ""),
                    "assessable": raw.get("vessel_presence_prediction", "") == "assessable",
                    "valid_local_body": b(raw.get("valid_local_body")),
                    "x1_fallback": b(raw.get("x1_fallback")),
                    "x4_jump_corrected": b(raw.get("x4_jump_corrected")),
                    "tracking_class": raw.get("tracking_class", ""),
                    "tracking_confidence": f(raw.get("tracking_confidence")),
                    "assessability_score": f(raw.get("assessability_score")),
                    "local_body_peak_cnr": f(raw.get("local_body_peak_cnr")),
                    "continuity_score": f(raw.get("continuity_score")),
                    "width_consistency_score": f(raw.get("width_consistency_score")),
                    "axial_completeness_score": f(raw.get("axial_completeness_score")),
                    "neighbor_support_score": f(raw.get("neighbor_support_score")),
                })
        if count != 500:
            raise RuntimeError(f"{scan}: expected 500 rows, got {count}")
        metadata.append({
            "scan_id": scan,
            "input_path": str(path.relative_to(ROOT)),
            "input_sha256": sha256(path),
            "row_count": count,
        })

    if len(frame_rows) != 2500:
        raise RuntimeError(f"expected 2500 total rows, got {len(frame_rows)}")

    all_summaries: list[dict[str, object]] = []
    all_summaries.append(summary(frame_rows, "all_frames", "ALL"))
    assessable_all = [r for r in frame_rows if bool(r["assessable"])]
    all_summaries.append(summary(assessable_all, "assessable_only", "ALL"))
    for scan in SCANS:
        subset = [r for r in frame_rows if r["scan_id"] == scan]
        all_summaries.append(summary(subset, "all_frames", scan))
        all_summaries.append(summary([r for r in subset if bool(r["assessable"])], "assessable_only", scan))

    smooth_rows: list[dict[str, object]] = []
    for scan in SCANS:
        smooth_rows.extend(smoothness([r for r in frame_rows if r["scan_id"] == scan], scan))

    extremes = sorted(
        [r for r in frame_rows if bool(r["paired"])],
        key=lambda r: (float(r["abs_delta_px"]), str(r["scan_id"]), int(r["frame_index"])),
        reverse=True,
    )[:50]

    write_csv(OUT / "x1_x4_framewise.csv", frame_rows)
    write_csv(OUT / "x1_x4_summary.csv", all_summaries)
    write_csv(OUT / "x1_x4_smoothness.csv", smooth_rows)
    write_csv(OUT / "x1_x4_extremes_top50.csv", extremes)
    write_csv(OUT / "audit_metadata.csv", metadata)

    overall = next(r for r in all_summaries if r["scope"] == "all_frames" and r["scan_id"] == "ALL")
    assess = next(r for r in all_summaries if r["scope"] == "assessable_only" and r["scan_id"] == "ALL")
    worst = extremes[0] if extremes else None
    smooth_lookup = {(r["scan_id"], r["trajectory"]): r for r in smooth_rows}

    lines = [
        "# X1 midpoint vs X4 full-volume stability audit",
        "",
        "This is a localization-only descriptive audit. It does **not** modify `results/formal_sv_d128_v21_full2500_run001`, the frozen formal geometry, or any Q_vessel/Q_tail/R endpoint.",
        "",
        "## Definitions",
        "",
        "- X1 midpoint: `x1_local_geometry_px`, the midpoint of the selected contiguous local-body run.",
        "- X4: `x4_centroid_isolated_jump_corrected_px`, i.e. robust X2 centroid after isolated 3-frame jump correction.",
        "- Signed difference: `delta_px = X4 - X1`; positive values mean X4 is to the larger-x side of X1.",
        f"- Lateral scale: {DX_UM} um/A-line; nominal vessel diameter: {DIAMETER_UM} um.",
        "- No B-scan is treated as an independent biological replicate; results are descriptive across spatial frames.",
        "",
        "## Whole 2500-frame result",
        "",
        f"- Paired X1/X4 coordinates: {overall['n_paired']}/{overall['n_frames']} ({100*float(overall['paired_fraction']):.2f}%).",
        f"- Signed X4-X1: mean {fmt(float(overall['mean_delta_px']))} px; median {fmt(float(overall['median_delta_px']))} px.",
        f"- Absolute difference: median {fmt(float(overall['median_abs_delta_px']))} px ({fmt(float(overall['median_abs_delta_um']),2)} um); p95 {fmt(float(overall['p95_abs_delta_px']))} px ({fmt(float(overall['p95_abs_delta_um']),2)} um); maximum {fmt(float(overall['max_abs_delta_px']))} px ({fmt(float(overall['max_abs_delta_um']),2)} um).",
        f"- Within 0.5 px: {100*float(overall['fraction_abs_le_0p5px']):.2f}%; within 1.0 px: {100*float(overall['fraction_abs_le_1px']):.2f}%; >2 px: {100*float(overall['fraction_abs_gt_2px']):.2f}%.",
        f"- X1 fallback frames: {overall['n_x1_fallback']}; X4 isolated-jump corrections: {overall['n_x4_jump_corrected']}.",
        "",
        "## Assessable subset",
        "",
        f"- Frames: {assess['n_frames']}.",
        f"- Signed median X4-X1: {fmt(float(assess['median_delta_px']))} px.",
        f"- Absolute difference: median {fmt(float(assess['median_abs_delta_px']))} px; p95 {fmt(float(assess['p95_abs_delta_px']))} px; maximum {fmt(float(assess['max_abs_delta_px']))} px.",
        f"- Within 0.5 px: {100*float(assess['fraction_abs_le_0p5px']):.2f}%; within 1.0 px: {100*float(assess['fraction_abs_le_1px']):.2f}%.",
        "",
        "## Adjacent-frame trajectory stability",
        "",
        "| Scan | X1 median step (px) | X1 p95 step | X4 median step (px) | X4 p95 step |",
        "|---|---:|---:|---:|---:|",
    ]
    for scan in SCANS:
        x1s = smooth_lookup[(scan, "X1_midpoint")]
        x4s = smooth_lookup[(scan, "X4")]
        lines.append(
            f"| {scan} | {fmt(float(x1s['median_abs_step_px']))} | {fmt(float(x1s['p95_abs_step_px']))} | {fmt(float(x4s['median_abs_step_px']))} | {fmt(float(x4s['p95_abs_step_px']))} |"
        )
    lines.extend([
        "",
        "## Largest discrepancy",
        "",
    ])
    if worst:
        lines.append(
            f"Largest |X4-X1| occurs at `{worst['scan_id']}` frame {int(worst['frame_index'])}: X1={fmt(float(worst['x1_midpoint_px']))} px, X4={fmt(float(worst['x4_px']))} px, delta={fmt(float(worst['delta_px']))} px ({fmt(float(worst['delta_um']),2)} um), assessability=`{worst['vessel_presence_prediction']}`, X1_fallback={worst['x1_fallback']}, X4_jump_corrected={worst['x4_jump_corrected']}."
        )
    lines.extend([
        "",
        "See `x1_x4_summary.csv` for flow-stratified values, `x1_x4_extremes_top50.csv` for the largest discrepancies, and `x1_x4_framewise.csv` for all 2500 frames.",
        "",
        "## Provenance",
        "",
        f"Workflow source SHA before generated commit: `{os.environ.get('GITHUB_SHA', 'local')}`.",
        "The five input tracking CSV SHA-256 values are recorded in `audit_metadata.csv`.",
        "",
    ])
    (OUT / "README.md").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
