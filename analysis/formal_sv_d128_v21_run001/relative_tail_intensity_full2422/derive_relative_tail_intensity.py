#!/usr/bin/env python3
"""Derive Relative Tail Intensity from frozen formal SV run001 outputs.

This is a post-processing analysis only. It never reconstructs OCT data and never
writes into the frozen formal results directory.

Formal definitions
------------------
A_T       = frozen fractional tail area (`tail_area_um2`)
Sbar_T    = Q_T / A_T
RI_tail   = Sbar_T / Sbar_V
RI(r)     = T(r) / Sbar_V

Only rows with frozen `valid == True` are included. Signed background-corrected
values are preserved; no clipping, normalization, gain, thresholding, or
interpolation is performed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

EXPECTED_VALID_FRAMES = 2422
EXPECTED_PROFILE_CHUNKS = 25
FORMAL_RELEASE_TAG = "formal-sv-d128-v21-run001"


def first_present(columns: Iterable[str], candidates: Sequence[str], *, required: bool = True) -> str | None:
    cols = set(columns)
    for name in candidates:
        if name in cols:
            return name
    if required:
        raise KeyError(f"None of the required columns are present: {candidates}")
    return None


def parse_bool_series(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False).astype(bool)
    if pd.api.types.is_numeric_dtype(series):
        return series.fillna(0).astype(float).ne(0)
    normalized = series.astype(str).str.strip().str.lower()
    true_values = {"true", "1", "yes", "y", "t"}
    false_values = {"false", "0", "no", "n", "f", "", "nan", "none", "na", "<na>"}
    unknown = ~normalized.isin(true_values | false_values)
    if unknown.any():
        vals = sorted(normalized[unknown].dropna().unique().tolist())[:10]
        raise ValueError(f"Unrecognized boolean values: {vals}")
    return normalized.isin(true_values)


def safe_divide(numerator: pd.Series | np.ndarray, denominator: pd.Series | np.ndarray) -> np.ndarray:
    num = np.asarray(numerator, dtype=float)
    den = np.asarray(denominator, dtype=float)
    out = np.full(num.shape, np.nan, dtype=float)
    ok = np.isfinite(num) & np.isfinite(den) & (den != 0.0)
    out[ok] = num[ok] / den[ok]
    return out


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def flow_from_scan_id(scan_id: str) -> float:
    match = re.search(r"(\d+(?:\.\d+)?)$", str(scan_id))
    if not match:
        return float("nan")
    return float(match.group(1))


def canonicalize_localization(localization: pd.DataFrame) -> pd.DataFrame:
    scan_col = first_present(localization.columns, ["scan_id", "scan"])
    frame_col = first_present(localization.columns, ["frame_index_0based", "frame_index"])
    loc_col = first_present(localization.columns, ["localization_source"], required=False)
    if loc_col is None:
        return pd.DataFrame(columns=["scan_id", "frame_index", "localization_source"])
    out = pd.DataFrame(
        {
            "scan_id": localization[scan_col].astype(str),
            "frame_index": pd.to_numeric(localization[frame_col], errors="raise").astype(int),
            "localization_source": localization[loc_col].astype(str),
        }
    )
    return out.drop_duplicates(["scan_id", "frame_index"], keep="last")


def prepare_framewise(formal_dir: Path) -> tuple[pd.DataFrame, dict]:
    frame_path = formal_dir / "frame_results.csv"
    raw = pd.read_csv(frame_path, low_memory=False)

    scan_col = first_present(raw.columns, ["scan_id", "scan"])
    frame_col = first_present(raw.columns, ["frame_index_0based", "frame_index"])
    valid_col = first_present(raw.columns, ["valid"])
    loc_col = first_present(raw.columns, ["localization_source"], required=False)
    flow_col = first_present(raw.columns, ["flow_mm_s", "flow_mm_per_s"], required=False)

    required_numeric = {
        "q_vessel": first_present(raw.columns, ["q_vessel"]),
        "q_tail": first_present(raw.columns, ["q_tail"]),
        "source_area_um2": first_present(raw.columns, ["source_area_um2"]),
        "tail_area_um2": first_present(raw.columns, ["tail_area_um2"]),
        "source_mean": first_present(raw.columns, ["source_mean"]),
        "ratio_tail_to_vessel": first_present(raw.columns, ["ratio_tail_to_vessel"]),
    }

    base = pd.DataFrame(
        {
            "scan_id": raw[scan_col].astype(str),
            "frame_index": pd.to_numeric(raw[frame_col], errors="raise").astype(int),
            "valid": parse_bool_series(raw[valid_col]),
        }
    )
    for canonical, source in required_numeric.items():
        base[canonical] = pd.to_numeric(raw[source], errors="coerce")

    if loc_col is not None:
        base["localization_source"] = raw[loc_col].astype(str)
    else:
        loc_path = formal_dir / "localization.csv"
        if not loc_path.exists():
            raise FileNotFoundError("localization_source absent from frame_results.csv and localization.csv missing")
        loc = canonicalize_localization(pd.read_csv(loc_path, low_memory=False))
        base = base.merge(loc, on=["scan_id", "frame_index"], how="left", validate="one_to_one")

    if flow_col is not None:
        base["flow_mm_s"] = pd.to_numeric(raw[flow_col], errors="coerce")
    else:
        base["flow_mm_s"] = base["scan_id"].map(flow_from_scan_id)

    detection_col = first_present(raw.columns, ["detection_status"], required=False)
    if detection_col is not None:
        base["detection_status"] = raw[detection_col].astype(str)

    if base.duplicated(["scan_id", "frame_index"]).any():
        dup = base.loc[base.duplicated(["scan_id", "frame_index"], keep=False), ["scan_id", "frame_index"]]
        raise ValueError(f"Duplicate frame keys in frozen frame_results.csv: {dup.head().to_dict('records')}")

    valid = base.loc[base["valid"]].copy()
    if len(valid) != EXPECTED_VALID_FRAMES:
        raise AssertionError(f"Expected {EXPECTED_VALID_FRAMES} valid frames; found {len(valid)}")

    if valid["localization_source"].isna().any() or valid["localization_source"].eq("nan").any():
        raise AssertionError("At least one valid frame lacks localization_source")

    if (valid["source_area_um2"] <= 0).any() or (valid["tail_area_um2"] <= 0).any():
        raise AssertionError("All valid frames must have positive frozen source and tail fractional areas")

    valid["tail_mean"] = safe_divide(valid["q_tail"], valid["tail_area_um2"])
    valid["ri_tail"] = safe_divide(valid["tail_mean"], valid["source_mean"])
    valid["frame_index_0based"] = valid["frame_index"]

    source_mean_recomputed = safe_divide(valid["q_vessel"], valid["source_area_um2"])
    burden_recomputed = safe_divide(valid["q_tail"], valid["q_vessel"])
    ri_from_integrals = safe_divide(
        safe_divide(valid["q_tail"], valid["tail_area_um2"]),
        safe_divide(valid["q_vessel"], valid["source_area_um2"]),
    )

    source_ok = np.allclose(valid["source_mean"].to_numpy(float), source_mean_recomputed, rtol=1e-10, atol=1e-12, equal_nan=True)
    burden_ok = np.allclose(valid["ratio_tail_to_vessel"].to_numpy(float), burden_recomputed, rtol=1e-10, atol=1e-12, equal_nan=True)
    ri_ok = np.allclose(valid["ri_tail"].to_numpy(float), ri_from_integrals, rtol=1e-12, atol=1e-14, equal_nan=True)
    if not source_ok:
        raise AssertionError("Frozen source_mean is inconsistent with q_vessel/source_area_um2")
    if not burden_ok:
        raise AssertionError("Frozen ratio_tail_to_vessel is inconsistent with q_tail/q_vessel")
    if not ri_ok:
        raise AssertionError("RI_tail implementation is inconsistent with its integral definition")

    def max_abs_delta(a: np.ndarray, b: np.ndarray) -> float | None:
        mask = np.isfinite(a) & np.isfinite(b)
        if not mask.any():
            return None
        return float(np.max(np.abs(a[mask] - b[mask])))

    validation = {
        "input_frame_count": int(len(base)),
        "valid_frame_count": int(len(valid)),
        "excluded_frame_count": int(len(base) - len(valid)),
        "source_area_nonpositive_valid": int((valid["source_area_um2"] <= 0).sum()),
        "tail_area_nonpositive_valid": int((valid["tail_area_um2"] <= 0).sum()),
        "source_mean_zero_valid": int((valid["source_mean"] == 0).sum()),
        "source_mean_nonfinite_valid": int((~np.isfinite(valid["source_mean"].to_numpy(float))).sum()),
        "ri_tail_nonfinite_valid": int((~np.isfinite(valid["ri_tail"].to_numpy(float))).sum()),
        "source_mean_identity_allclose": bool(source_ok),
        "relative_tail_burden_identity_allclose": bool(burden_ok),
        "ri_tail_identity_allclose": bool(ri_ok),
        "max_abs_source_mean_identity_error": max_abs_delta(valid["source_mean"].to_numpy(float), source_mean_recomputed),
        "max_abs_relative_tail_burden_identity_error": max_abs_delta(valid["ratio_tail_to_vessel"].to_numpy(float), burden_recomputed),
        "max_abs_ri_tail_identity_error": max_abs_delta(valid["ri_tail"].to_numpy(float), ri_from_integrals),
    }

    preferred = [
        "scan_id",
        "frame_index",
        "frame_index_0based",
        "flow_mm_s",
        "localization_source",
        "valid",
        "q_vessel",
        "q_tail",
        "source_area_um2",
        "tail_area_um2",
        "source_mean",
        "tail_mean",
        "ratio_tail_to_vessel",
        "ri_tail",
    ]
    if "detection_status" in valid.columns:
        preferred.append("detection_status")
    return valid[preferred].sort_values(["flow_mm_s", "scan_id", "frame_index"]).reset_index(drop=True), validation


def describe_series(series: pd.Series) -> dict[str, float | int]:
    values = pd.to_numeric(series, errors="coerce").to_numpy(float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return {"n": 0, "mean": np.nan, "sd": np.nan, "min": np.nan, "q1": np.nan, "median": np.nan, "q3": np.nan, "iqr": np.nan, "max": np.nan}
    q1, median, q3 = np.quantile(values, [0.25, 0.5, 0.75])
    return {
        "n": int(values.size),
        "mean": float(np.mean(values)),
        "sd": float(np.std(values, ddof=1)) if values.size > 1 else np.nan,
        "min": float(np.min(values)),
        "q1": float(q1),
        "median": float(median),
        "q3": float(q3),
        "iqr": float(q3 - q1),
        "max": float(np.max(values)),
    }


def make_scan_summary(framewise: pd.DataFrame) -> pd.DataFrame:
    metrics = [
        "q_vessel",
        "source_area_um2",
        "source_mean",
        "q_tail",
        "tail_area_um2",
        "tail_mean",
        "ratio_tail_to_vessel",
        "ri_tail",
    ]
    records: list[dict] = []
    for (scan_id, flow), group in framewise.groupby(["scan_id", "flow_mm_s"], sort=True, dropna=False):
        rec: dict[str, object] = {
            "scan_id": scan_id,
            "flow_mm_s": flow,
            "n_valid_frames": int(len(group)),
            "frame_min": int(group["frame_index"].min()),
            "frame_max": int(group["frame_index"].max()),
        }
        for metric in metrics:
            stats = describe_series(group[metric])
            for stat_name, value in stats.items():
                rec[f"{metric}_{stat_name}"] = value
        records.append(rec)
    return pd.DataFrame.from_records(records).sort_values(["flow_mm_s", "scan_id"]).reset_index(drop=True)


def canonical_profile_columns(profile: pd.DataFrame) -> dict[str, str | None]:
    return {
        "scan_id": first_present(profile.columns, ["scan_id", "scan"], required=False),
        "frame_index": first_present(profile.columns, ["frame_index_0based", "frame_index"]),
        "z_index_0based": first_present(profile.columns, ["z_index_0based", "row"], required=False),
        "z_um": first_present(profile.columns, ["z_um"], required=False),
        "r_px": first_present(profile.columns, ["r_px"], required=False),
        "r_um": first_present(profile.columns, ["r_um"]),
        "V": first_present(profile.columns, ["V", "vessel_mean", "profile_raw"], required=False),
        "B_left": first_present(profile.columns, ["B_left", "background_left"], required=False),
        "B_right": first_present(profile.columns, ["B_right", "background_right"], required=False),
        "B": first_present(profile.columns, ["B", "background_mean"], required=False),
        "T": first_present(profile.columns, ["T", "profile_bg_corrected"]),
        "P": first_present(profile.columns, ["P", "profile_integral"], required=False),
        "tail_z_fraction": first_present(profile.columns, ["tail_z_fraction"]),
        "validity": first_present(profile.columns, ["validity"], required=False),
        "frame_valid": first_present(profile.columns, ["frame_valid"], required=False),
        "effective_width_um": first_present(profile.columns, ["effective_width_um"], required=False),
    }


def process_profiles(formal_dir: Path, output_dir: Path, framewise: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    source_dir = formal_dir / "profiles"
    input_paths = sorted(source_dir.glob("*.csv.gz"))
    if len(input_paths) != EXPECTED_PROFILE_CHUNKS:
        raise AssertionError(f"Expected {EXPECTED_PROFILE_CHUNKS} frozen profile chunks; found {len(input_paths)}")

    profile_out_dir = output_dir / "profiles"
    profile_out_dir.mkdir(parents=True, exist_ok=True)

    key_lookup = framewise[["scan_id", "frame_index", "flow_mm_s", "localization_source", "valid", "source_mean"]].copy()
    key_lookup["scan_id"] = key_lookup["scan_id"].astype(str)
    key_lookup["frame_index"] = key_lookup["frame_index"].astype(int)
    valid_keys = set(zip(key_lookup["scan_id"], key_lookup["frame_index"]))

    profile_keys_seen: set[tuple[str, int]] = set()
    index_records: list[dict] = []
    total_input_rows = 0
    total_output_rows = 0
    total_negative_t = 0
    total_negative_ri = 0
    total_nonfinite_ri = 0

    for source_path in input_paths:
        profile = pd.read_csv(source_path, compression="gzip", low_memory=False)
        cols = canonical_profile_columns(profile)
        inferred_scan_match = re.match(r"(flow\d+)_frames_", source_path.name)
        inferred_scan = inferred_scan_match.group(1) if inferred_scan_match else None

        if cols["scan_id"] is None:
            if inferred_scan is None:
                raise ValueError(f"Cannot infer scan_id for {source_path.name}")
            scan_values = pd.Series(inferred_scan, index=profile.index, dtype="object")
        else:
            scan_values = profile[cols["scan_id"]].astype(str)

        frame_values = pd.to_numeric(profile[cols["frame_index"]], errors="raise").astype(int)
        working = pd.DataFrame({"scan_id": scan_values.astype(str), "frame_index": frame_values})
        mask = np.fromiter(((s, int(f)) in valid_keys for s, f in zip(working["scan_id"], working["frame_index"])), dtype=bool, count=len(working))
        filtered = profile.loc[mask].copy()
        working = working.loc[mask].copy()

        total_input_rows += int(len(profile))
        if filtered.empty:
            # A chunk may contain only formally excluded frames in principle; retain an empty indexed output.
            merged = working.copy()
        else:
            merged = working.merge(key_lookup, on=["scan_id", "frame_index"], how="left", validate="many_to_one")

        out = merged[["scan_id", "frame_index", "flow_mm_s", "localization_source", "valid", "source_mean"]].copy()
        out["frame_index_0based"] = out["frame_index"]

        for canonical in ["z_index_0based", "z_um", "r_px", "r_um", "V", "B_left", "B_right", "B", "P", "tail_z_fraction", "validity", "frame_valid", "effective_width_um"]:
            source_col = cols[canonical]
            if source_col is not None:
                out[canonical] = filtered[source_col].to_numpy()

        t_values = pd.to_numeric(filtered[cols["T"]], errors="coerce").to_numpy(float)
        out["T"] = t_values
        out["ri_r"] = safe_divide(t_values, out["source_mean"])

        # Formal tail coverage is required and remains fractional; no r-based replacement is allowed.
        if "tail_z_fraction" not in out.columns:
            raise AssertionError(f"Frozen profile {source_path.name} lacks tail_z_fraction")

        for s, f in zip(out["scan_id"], out["frame_index"]):
            profile_keys_seen.add((str(s), int(f)))

        total_output_rows += int(len(out))
        total_negative_t += int(np.sum(np.isfinite(t_values) & (t_values < 0)))
        ri_values = pd.to_numeric(out["ri_r"], errors="coerce").to_numpy(float)
        total_negative_ri += int(np.sum(np.isfinite(ri_values) & (ri_values < 0)))
        total_nonfinite_ri += int(np.sum(~np.isfinite(ri_values)))

        output_name = source_path.name.replace(".csv.gz", "_relative_intensity.csv.gz")
        output_path = profile_out_dir / output_name
        out.to_csv(
            output_path,
            index=False,
            float_format="%.17g",
            compression={"method": "gzip", "compresslevel": 9, "mtime": 0},
        )
        index_records.append(
            {
                "input_profile": str(source_path.as_posix()),
                "output_profile": str(output_path.as_posix()),
                "input_rows": int(len(profile)),
                "output_valid_rows": int(len(out)),
                "unique_valid_frames": int(out[["scan_id", "frame_index"]].drop_duplicates().shape[0]),
                "negative_T_rows": int(np.sum(np.isfinite(t_values) & (t_values < 0))),
                "negative_RI_rows": int(np.sum(np.isfinite(ri_values) & (ri_values < 0))),
                "nonfinite_RI_rows": int(np.sum(~np.isfinite(ri_values))),
            }
        )

    missing_keys = sorted(valid_keys - profile_keys_seen)
    unexpected_keys = sorted(profile_keys_seen - valid_keys)
    if missing_keys:
        raise AssertionError(f"Valid frames missing from frozen profiles: {missing_keys[:20]}")
    if unexpected_keys:
        raise AssertionError(f"Unexpected invalid keys survived profile filtering: {unexpected_keys[:20]}")

    profile_validation = {
        "profile_chunk_count": int(len(input_paths)),
        "profile_input_rows_all_frames": int(total_input_rows),
        "profile_output_rows_valid_frames": int(total_output_rows),
        "profile_unique_valid_frames": int(len(profile_keys_seen)),
        "profile_missing_valid_frames": int(len(missing_keys)),
        "profile_unexpected_frames_after_filter": int(len(unexpected_keys)),
        "profile_negative_T_rows_preserved": int(total_negative_t),
        "profile_negative_RI_rows_preserved": int(total_negative_ri),
        "profile_nonfinite_RI_rows": int(total_nonfinite_ri),
        "profile_interpolation_performed": False,
        "profile_clipping_performed": False,
    }
    return pd.DataFrame.from_records(index_records), profile_validation


def collect_input_hashes(formal_dir: Path) -> pd.DataFrame:
    paths = [
        formal_dir / "frame_results.csv",
        formal_dir / "localization.csv",
        formal_dir / "run_config.json",
        formal_dir / "tracking_config.json",
        formal_dir / "DATA_DICTIONARY.md",
    ]
    paths.extend(sorted((formal_dir / "profiles").glob("*.csv.gz")))
    records = []
    for path in paths:
        if not path.exists():
            continue
        records.append(
            {
                "input_path": path.as_posix(),
                "sha256": sha256_file(path),
                "size_bytes": int(path.stat().st_size),
            }
        )
    return pd.DataFrame.from_records(records)


def fmt_interval(median: float, q1: float, q3: float) -> str:
    if not all(np.isfinite([median, q1, q3])):
        return "NA"
    return f"{median:.6g} [{q1:.6g}, {q3:.6g}]"


def write_readme(output_dir: Path, summary: pd.DataFrame, validation: dict, frozen_source_sha: str, workflow_source_sha: str) -> None:
    lines = [
        "# Relative Tail Intensity — formal SV d128 v2.1 run001",
        "",
        "This directory is a post-processing derivative of the frozen formal run. The frozen directory `results/formal_sv_d128_v21_full2500_run001` is not modified, and no `.oct` reconstruction is performed here.",
        "",
        "## Formal inclusion",
        "",
        f"- Included frames: **{validation['valid_frame_count']}** with frozen `valid == True`.",
        f"- Excluded frames remain excluded: **{validation['excluded_frame_count']}**; no zero filling or interpolation is used.",
        "- Linear background-corrected SV values are used as stored; negative values are retained.",
        "- Flow/OMAG localization images are not used in the numerical Q_V/Q_T/RI integration.",
        "- Detection depth remains outside this analysis (`not_evaluated / NA` without a matched blank).",
        "",
        "## Definitions",
        "",
        "- `q_vessel` = Q_V: background-corrected source-ellipse integral.",
        "- `source_area_um2` = A_V: frozen fractional source area.",
        "- `source_mean` = Sbar_V = Q_V / A_V.",
        "- `q_tail` = Q_T: background-corrected 0–500 µm rectangular-tail integral.",
        "- `tail_area_um2` = A_T: **frozen fractional tail area** from the formal geometry; no ideal-width approximation is substituted.",
        "- `tail_mean` = Sbar_T = Q_T / A_T.",
        "- `ratio_tail_to_vessel` = R_B = Q_T / Q_V: **relative tail burden**, not tail intensity.",
        "- `ri_tail` = RI_tail = (Q_T/A_T)/(Q_V/A_V): **Relative Tail Intensity**.",
        "- `T` = T(r) = V(r) − B(r), retaining signed background-corrected depth signal.",
        "- `ri_r` = RI(r) = T(r)/Sbar_V.",
        "",
        "The depth-profile outputs preserve every available formal profile row belonging to a valid frame and keep `tail_z_fraction`. No cross-frame depth interpolation or curve fitting is introduced in Task 1.",
        "",
        "## Scan-level descriptive summary",
        "",
        "Frames are spatial samples within one scan volume per flow condition. The table is descriptive; the B-scans are **not** treated as independent biological replicates and no p-values are computed.",
        "",
        "| scan | flow (mm/s) | valid frames | Q_V median [Q1,Q3] | Sbar_V median [Q1,Q3] | Q_T median [Q1,Q3] | Sbar_T median [Q1,Q3] | R_B median [Q1,Q3] | RI_tail median [Q1,Q3] |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for _, row in summary.iterrows():
        lines.append(
            "| {scan} | {flow:g} | {n:d} | {qv} | {sv} | {qt} | {st} | {rb} | {ri} |".format(
                scan=row["scan_id"],
                flow=float(row["flow_mm_s"]),
                n=int(row["n_valid_frames"]),
                qv=fmt_interval(row["q_vessel_median"], row["q_vessel_q1"], row["q_vessel_q3"]),
                sv=fmt_interval(row["source_mean_median"], row["source_mean_q1"], row["source_mean_q3"]),
                qt=fmt_interval(row["q_tail_median"], row["q_tail_q1"], row["q_tail_q3"]),
                st=fmt_interval(row["tail_mean_median"], row["tail_mean_q1"], row["tail_mean_q3"]),
                rb=fmt_interval(row["ratio_tail_to_vessel_median"], row["ratio_tail_to_vessel_q1"], row["ratio_tail_to_vessel_q3"]),
                ri=fmt_interval(row["ri_tail_median"], row["ri_tail_q1"], row["ri_tail_q3"]),
            )
        )
    lines.extend(
        [
            "",
            "## Outputs",
            "",
            "- `relative_tail_intensity_framewise.csv`: one row per formally valid frame.",
            "- `relative_tail_intensity_scan_summary.csv`: spatial descriptive statistics by scan/flow.",
            "- `profiles/*_relative_intensity.csv.gz`: valid-frame depth profiles retaining T(r), fractional tail coverage, and adding RI(r).",
            "- `profile_index.csv`: row/frame accounting for every profile chunk.",
            "- `input_sha256.csv`: SHA-256 hashes of frozen scalar/profile inputs used by this derivation.",
            "- `validation.json`: identity checks, inclusion counts, profile coverage, and signed-value checks.",
            "- `provenance.json`: frozen source, workflow source, definitions, and analysis restrictions.",
            "- `derive_relative_tail_intensity.py`: reproducible derivation script.",
            "",
            "## Provenance",
            "",
            f"- Frozen formal handoff SHA: `{frozen_source_sha}`.",
            f"- Workflow-trigger source SHA: `{workflow_source_sha}`.",
            f"- Formal 2D-array release tag: `{FORMAL_RELEASE_TAG}`.",
        ]
    )
    (output_dir / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--formal-dir", default="results/formal_sv_d128_v21_full2500_run001")
    parser.add_argument("--output-dir", default="analysis/formal_sv_d128_v21_run001/relative_tail_intensity_full2422")
    parser.add_argument("--frozen-source-sha", required=True)
    parser.add_argument("--workflow-source-sha", required=True)
    args = parser.parse_args()

    root = args.root.resolve()
    formal_dir = (root / args.formal_dir).resolve()
    output_dir = (root / args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    if formal_dir == output_dir or formal_dir in output_dir.parents:
        raise AssertionError("Output directory must be independent of the frozen formal input directory")

    framewise, validation = prepare_framewise(formal_dir)
    summary = make_scan_summary(framewise)
    profile_index, profile_validation = process_profiles(formal_dir, output_dir, framewise)
    validation.update(profile_validation)
    validation["scan_valid_counts"] = {
        str(row.scan_id): int(row.n_valid_frames) for row in summary[["scan_id", "n_valid_frames"]].itertuples(index=False)
    }
    validation["no_p_values"] = True
    validation["no_parameter_tuning_by_flow_trend"] = True
    validation["negative_background_corrected_values_retained"] = True

    if validation["profile_unique_valid_frames"] != EXPECTED_VALID_FRAMES:
        raise AssertionError("Depth profiles do not cover all 2422 valid frames")

    framewise.to_csv(output_dir / "relative_tail_intensity_framewise.csv", index=False, float_format="%.17g")
    summary.to_csv(output_dir / "relative_tail_intensity_scan_summary.csv", index=False, float_format="%.17g")
    profile_index.to_csv(output_dir / "profile_index.csv", index=False)
    collect_input_hashes(formal_dir).to_csv(output_dir / "input_sha256.csv", index=False)

    provenance = {
        "analysis": "Task 1: Relative Tail Intensity, formal SV d128 v2.1 run001",
        "formal_release_tag": FORMAL_RELEASE_TAG,
        "frozen_input_directory": args.formal_dir,
        "output_directory": args.output_dir,
        "frozen_source_sha": args.frozen_source_sha,
        "workflow_source_sha": args.workflow_source_sha,
        "selection": "valid == True",
        "expected_and_observed_valid_frames": EXPECTED_VALID_FRAMES,
        "definitions": {
            "A_T": "tail_area_um2 from frozen fractional tail weights",
            "Sbar_T": "q_tail / tail_area_um2",
            "RI_tail": "(q_tail/tail_area_um2) / source_mean = (Q_T/A_T)/(Q_V/A_V)",
            "T_r": "frozen T(r) = V(r) - B(r)",
            "RI_r": "T(r) / source_mean",
            "R_B": "q_tail/q_vessel; relative tail burden, not tail intensity",
        },
        "frozen_signal_rules": {
            "signal": "linear sv_raw = Var_t(|E|), N denominator",
            "normalization": "none",
            "per_frame_gain": "none",
            "threshold_clipping": "none",
            "positive_truncation": "none",
            "negative_background_corrected_values": "retained",
        },
        "geometry": "X4 center + X1 apparent lateral width + frozen v2.1 z_upper; source ellipse; 500 um rectangular tail; guard 0",
        "background": "skip 3 A-lines outside each vessel edge, take 5 per side, combine left/right",
        "reconstruction_performed": False,
        "depth_interpolation_performed": False,
        "curve_fitting_performed": False,
        "p_values_computed": False,
        "detection_depth": "not_evaluated / NA without matched blank; not used here",
        "statistical_unit_note": "B-scans are adjacent spatial positions within one scan volume per flow, not independent biological replicates.",
    }
    (output_dir / "provenance.json").write_text(json.dumps(provenance, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    (output_dir / "validation.json").write_text(json.dumps(validation, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    write_readme(output_dir, summary, validation, args.frozen_source_sha, args.workflow_source_sha)

    print(json.dumps({
        "valid_frames": int(len(framewise)),
        "scan_valid_counts": validation["scan_valid_counts"],
        "profile_chunks": validation["profile_chunk_count"],
        "profile_rows": validation["profile_output_rows_valid_frames"],
        "ri_tail_nonfinite": validation["ri_tail_nonfinite_valid"],
        "ri_r_nonfinite_rows": validation["profile_nonfinite_RI_rows"],
        "output_dir": args.output_dir,
    }, indent=2))


if __name__ == "__main__":
    main()
