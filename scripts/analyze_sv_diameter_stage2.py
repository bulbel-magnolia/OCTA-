#!/usr/bin/env python3
"""SV diameter Stage 2 scientific analysis.

This script consumes the frozen no-background SV outputs already committed to the
repository. It does not modify any frozen D128 or multi-diameter result.

D128 did not historically expose 100/200/300 um no-background scalar windows.
To make the 4x5 common grid complete, the script performs a read-only replay from
the released D128 2-D sv_raw arrays using the same frozen geometry and the same
raw-SV quantifier. The replayed 500 um source/tail/RI values are checked against
the existing D128 no-background primary table before supplemental values are used.
"""
from __future__ import annotations

import hashlib
import io
import json
import math
import os
import re
import tempfile
import urllib.request
import zipfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis/sv_diameter_stage2_scientific_v1"
NEW = ROOT / "analysis/formal_sv_diameter_v1"
D128_NOBG = ROOT / "analysis/formal_sv_d128_v21_run001/no_background_relative_tail_v1_full2422"
D128_FORMAL = ROOT / "results/formal_sv_d128_v21_full2500_run001"

COMMON_FLOWS = [1, 3, 5, 7, 10]
DIAMETERS = [128, 235, 285, 500]
WINDOWS = [100, 200, 300, 500]
DEPTH_TARGETS = list(range(0, 501, 10))
SELECTED_DEPTHS = [0, 50, 100, 200, 300, 400, 500]
CONTRASTS = [(128, 235), (235, 285), (285, 500), (128, 500)]
SUBSETS = [
    "all_valid",
    "direct_candidate_supported_only",
    "central80_x1",
    "central80_z_top",
    "central80_x1_z_top",
    "direct_plus_central_geometry",
]
RELEASE_TAG = "formal-sv-d128-v21-run001"
RELEASE_BASE = f"https://github.com/bulbel-magnolia/OCTA-/releases/download/{RELEASE_TAG}"
D128_MAP = {
    "flow01": "D128_F01_V01",
    "flow03": "D128_F03_V01",
    "flow05": "D128_F05_V01",
    "flow07": "D128_F07_V01",
    "flow10": "D128_F10_V01",
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def write_csv(name: str, frame: pd.DataFrame) -> None:
    frame.to_csv(OUT / name, index=False, float_format="%.17g", lineterminator="\n")


def write_json(name: str, obj) -> None:
    with (OUT / name).open("w", encoding="utf-8", newline="\n") as f:
        json.dump(obj, f, ensure_ascii=False, indent=2, allow_nan=False)
        f.write("\n")


def safe_pct(new: float, old: float) -> float:
    if not (np.isfinite(new) and np.isfinite(old)) or old == 0:
        return np.nan
    return 100.0 * (new - old) / old


def rel_error(a: float, b: float) -> float:
    return abs(float(a) - float(b)) / max(abs(float(b)), np.finfo(float).tiny)


def stats(values) -> dict:
    a = pd.to_numeric(pd.Series(values), errors="coerce").to_numpy(float)
    a = a[np.isfinite(a)]
    if not len(a):
        return {k: np.nan for k in ("q1", "median", "q3", "mean", "sd", "min", "max")} | {"n": 0}
    q1, med, q3 = np.quantile(a, [0.25, 0.5, 0.75])
    return {
        "n": int(len(a)), "q1": float(q1), "median": float(med), "q3": float(q3),
        "mean": float(a.mean()), "sd": float(a.std(ddof=1)) if len(a) > 1 else np.nan,
        "min": float(a.min()), "max": float(a.max()),
    }


def spearman_rho(x, y) -> float:
    a = pd.to_numeric(pd.Series(x), errors="coerce")
    b = pd.to_numeric(pd.Series(y), errors="coerce")
    m = a.notna() & b.notna()
    if int(m.sum()) < 3:
        return np.nan
    ar = a[m].rank(method="average").to_numpy(float)
    br = b[m].rank(method="average").to_numpy(float)
    if np.std(ar) == 0 or np.std(br) == 0:
        return np.nan
    return float(np.corrcoef(ar, br)[0, 1])


def bool_series(s: pd.Series) -> pd.Series:
    if s.dtype == bool:
        return s
    return s.astype(str).str.lower().isin(["true", "1", "yes"])


def download(url: str, path: Path) -> None:
    req = urllib.request.Request(url, headers={"User-Agent": "OCTA-SV-Stage2/1.0"})
    with urllib.request.urlopen(req, timeout=180) as response, path.open("wb") as out:
        while True:
            block = response.read(4 * 1024 * 1024)
            if not block:
                break
            out.write(block)


def parse_package_name(name: str):
    m = re.search(r"_(flow\d{2})_(\d{3})_(\d{3})\.zip$", name)
    if not m:
        raise ValueError(f"Unrecognized D128 package name: {name}")
    return m.group(1), int(m.group(2)), int(m.group(3))


def load_new_framewise() -> pd.DataFrame:
    fw = pd.read_csv(NEW / "framewise_primary.csv")
    if len(fw) != 8509:
        raise AssertionError(f"Expected 8509 new-diameter valid frames, got {len(fw)}")
    required = [
        "scan_id", "diameter_um", "flow_mm_s", "frame_index_0based", "source_mean_raw",
        "source_area_um2", "apparent_width_um", "x4_px", "z_top_edge_px",
        "tail_mean_raw_100um", "tail_mean_raw_200um", "tail_mean_raw_300um", "tail_mean_raw_500um",
        "ri_tail_100um", "ri_tail_200um", "ri_tail_300um", "ri_tail_500um", "ri_tail",
        "localization_source", "candidate_status",
    ]
    missing = [c for c in required if c not in fw.columns]
    if missing:
        raise AssertionError(f"Missing new-diameter framewise fields: {missing}")
    fw = fw.copy()
    fw["diameter_um"] = fw["diameter_um"].astype(float)
    fw["flow_mm_s"] = fw["flow_mm_s"].astype(float)
    fw["direct_candidate_supported"] = fw["localization_source"].astype(str).eq("direct_candidate")
    fw["grid_scope"] = np.where(fw["flow_mm_s"].isin(COMMON_FLOWS), "common_grid", "d500_extension")
    return fw


def replay_d128() -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Read-only D128 replay for complete no-background window metrics and RI(r)."""
    import sys
    sys.path.insert(0, str(ROOT / "src"))
    from svrecttail.geometry import VesselGeometry
    from svrecttail.raw_sv_quantification import quantify_raw_sv

    base = pd.read_csv(D128_NOBG / "framewise_primary.csv")
    if len(base) != 2422:
        raise AssertionError(f"Expected 2422 D128 valid frames, got {len(base)}")
    loc = pd.read_csv(D128_FORMAL / "localization.csv")
    loc["frame_index_0based"] = pd.to_numeric(loc["frame_index_0based"]).astype(int)
    base["frame_index_0based"] = pd.to_numeric(base["frame_index_0based"]).astype(int)
    merged = base.merge(loc, on=["scan_id", "frame_index_0based"], how="left", validate="one_to_one", suffixes=("", "_loc"))
    if merged[["x_left_edge_px", "x_right_edge_px", "z_top_edge_px"]].isna().any().any():
        raise AssertionError("D128 valid rows failed to map to frozen localization geometry")

    arrays = pd.read_csv(D128_FORMAL / "arrays_sha256.csv")
    arrays["frame_index_0based"] = pd.to_numeric(arrays["frame_index_0based"]).astype(int)
    array_sha = {(str(r.scan_id), int(r.frame_index_0based)): str(r.sha256) for r in arrays.itertuples(index=False)}
    packages = pd.read_csv(D128_FORMAL / "download_packages.csv")

    rows = []
    depth_rows = []
    package_audit = []
    max_source_err = 0.0
    max_tail_err = 0.0
    max_ri_err = 0.0
    verified_npz = 0

    merged_keyed = merged.set_index(["scan_id", "frame_index_0based"], drop=False)
    with tempfile.TemporaryDirectory(prefix="sv-stage2-d128-") as td:
        td = Path(td)
        for pkg in packages.itertuples(index=False):
            name = str(pkg.file)
            scan, lo, hi = parse_package_name(name)
            path = td / name
            download(f"{RELEASE_BASE}/{name}", path)
            if sha256_file(path) != str(pkg.sha256):
                raise AssertionError(f"D128 package SHA mismatch: {name}")
            if hasattr(pkg, "bytes") and int(path.stat().st_size) != int(pkg.bytes):
                raise AssertionError(f"D128 package size mismatch: {name}")
            wanted = merged[(merged.scan_id.eq(scan)) & merged.frame_index_0based.between(lo, hi)]
            npz_count = 0
            with zipfile.ZipFile(path) as zf:
                for r in wanted.itertuples(index=False):
                    fi = int(r.frame_index_0based)
                    member = f"arrays/{scan}/frame_{fi:03d}.npz"
                    data = zf.read(member)
                    if sha256_bytes(data) != array_sha[(scan, fi)]:
                        raise AssertionError(f"D128 NPZ SHA mismatch: {scan}/{fi}")
                    npz_count += 1
                    verified_npz += 1
                    with np.load(io.BytesIO(data), allow_pickle=False) as d:
                        sv = np.asarray(d["sv_raw"], dtype=np.float64)
                    geom = VesselGeometry(
                        x_left_edge_px=float(r.x_left_edge_px),
                        x_right_edge_px=float(r.x_right_edge_px),
                        z_top_edge_px=float(r.z_top_edge_px),
                        diameter_um=128.0,
                        dx_um=12.7,
                        dz_um=6.7,
                    )
                    metrics, anchors = quantify_raw_sv(sv, geom, geometry_qc_valid=True)
                    max_source_err = max(max_source_err, rel_error(metrics["source_mean_raw"], r.source_mean_raw))
                    max_tail_err = max(max_tail_err, rel_error(metrics["tail_mean_raw_500um"], r.tail_mean_raw))
                    max_ri_err = max(max_ri_err, rel_error(metrics["ri_tail"], r.ri_tail))
                    direct = bool(getattr(r, "z_candidate_accepted"))
                    if direct:
                        loc_source = "direct_candidate"
                    elif bool(getattr(r, "z_short_gap_filled")):
                        loc_source = "short_gap_fill"
                    else:
                        loc_source = "other_valid"
                    item = {
                        "scan_id": D128_MAP[scan], "legacy_scan_id": scan,
                        "diameter_um": 128.0, "flow_mm_s": float(r.flow_mm_s),
                        "frame_index_0based": fi,
                        "source_mean_raw": float(metrics["source_mean_raw"]),
                        "source_area_um2": float(metrics["source_area_um2"]),
                        "apparent_width_um": float((r.x_right_edge_px - r.x_left_edge_px) * 12.7),
                        "x4_px": float((r.x_left_edge_px + r.x_right_edge_px) / 2.0),
                        "z_top_edge_px": float(r.z_top_edge_px),
                        "localization_source": loc_source,
                        "candidate_status": str(r.z_continuity_status),
                        "direct_candidate_supported": direct,
                        "grid_scope": "common_grid",
                    }
                    for w in WINDOWS:
                        item[f"tail_mean_raw_{w}um"] = float(metrics[f"tail_mean_raw_{w}um"])
                        item[f"ri_tail_{w}um"] = float(metrics[f"ri_tail_{w}um"])
                    item["ri_tail"] = float(metrics["ri_tail"])
                    rows.append(item)
                    for a in anchors:
                        depth_rows.append({
                            "scan_id": D128_MAP[scan], "legacy_scan_id": scan,
                            "diameter_um": 128.0, "flow_mm_s": float(r.flow_mm_s),
                            "frame_index_0based": fi, **a,
                        })
            package_audit.append({"package": name, "sha256_verified": True, "valid_frames": int(len(wanted)), "npz_verified": int(npz_count)})
            path.unlink()

    if verified_npz != 2422 or len(rows) != 2422:
        raise AssertionError(f"D128 replay coverage mismatch: {verified_npz}/{len(rows)}")
    if max(max_source_err, max_tail_err, max_ri_err) > 1e-10:
        raise AssertionError(f"D128 500um replay exceeded tolerance: source={max_source_err}, tail={max_tail_err}, ri={max_ri_err}")

    audit = {
        "frames_replayed": 2422,
        "release_packages_verified": int(len(package_audit)),
        "npz_sha256_verified": int(verified_npz),
        "source_mean_raw_max_relative_error_vs_frozen": float(max_source_err),
        "tail_mean_raw_500um_max_relative_error_vs_frozen": float(max_tail_err),
        "ri_tail_max_relative_error_vs_frozen": float(max_ri_err),
        "tolerance": 1e-10,
        "package_audit": package_audit,
    }
    return pd.DataFrame(rows), pd.DataFrame(depth_rows), audit


def build_volume_metrics(fw: pd.DataFrame) -> pd.DataFrame:
    rows = []
    metrics = ["source_mean_raw"] + [f"tail_mean_raw_{w}um" for w in WINDOWS] + [f"ri_tail_{w}um" for w in WINDOWS]
    for scan, g in fw.groupby("scan_id", sort=True):
        row = {
            "scan_id": scan,
            "diameter_um": float(g.diameter_um.iloc[0]),
            "flow_mm_s": float(g.flow_mm_s.iloc[0]),
            "nominal_frames": 500,
            "geometry_valid_frames": int(len(g)),
            "geometry_valid_fraction": float(len(g) / 500.0),
            "grid_scope": str(g.grid_scope.iloc[0]),
        }
        for m in metrics:
            s = stats(g[m])
            row[f"{m}_median"] = s["median"]
            row[f"{m}_q1"] = s["q1"]
            row[f"{m}_q3"] = s["q3"]
        row["ri_tail_median"] = row["ri_tail_500um_median"]
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["diameter_um", "flow_mm_s"]).reset_index(drop=True)


def validate_new_volume_replay(vol: pd.DataFrame) -> float:
    ref = pd.read_csv(NEW / "volume_summary.csv")
    merged = vol[vol.diameter_um.ne(128)].merge(ref, on=["scan_id", "diameter_um", "flow_mm_s"], validate="one_to_one", suffixes=("_stage2", "_formal"))
    errs = []
    for a, b in [
        ("source_mean_raw_median_stage2", "source_mean_raw_median_formal"),
        ("tail_mean_raw_500um_median_stage2", "tail_mean_raw_500um_median_formal"),
        ("ri_tail_median_stage2", "ri_tail_median_formal"),
    ]:
        errs.extend(rel_error(x, y) for x, y in zip(merged[a], merged[b]))
    return float(max(errs))


def build_contrasts(common: pd.DataFrame) -> pd.DataFrame:
    metrics = ["source_mean_raw"] + [f"tail_mean_raw_{w}um" for w in WINDOWS] + [f"ri_tail_{w}um" for w in WINDOWS]
    rows = []
    for flow in COMMON_FLOWS:
        g = common[common.flow_mm_s.eq(flow)].set_index("diameter_um")
        for metric in metrics:
            col = f"{metric}_median"
            for d0, d1 in CONTRASTS:
                v0, v1 = float(g.loc[d0, col]), float(g.loc[d1, col])
                rows.append({
                    "flow_mm_s": flow, "metric": metric,
                    "diameter_from_um": d0, "diameter_to_um": d1,
                    "value_from": v0, "value_to": v1,
                    "absolute_difference": v1 - v0,
                    "percent_difference": safe_pct(v1, v0),
                })
    return pd.DataFrame(rows)


def build_flow_variation(common: pd.DataFrame) -> pd.DataFrame:
    metrics = ["source_mean_raw"] + [f"tail_mean_raw_{w}um" for w in WINDOWS] + [f"ri_tail_{w}um" for w in WINDOWS]
    rows = []
    for d in DIAMETERS:
        g = common[common.diameter_um.eq(d)].sort_values("flow_mm_s")
        for metric in metrics:
            col = f"{metric}_median"
            vals = g[col].to_numpy(float)
            med = float(np.median(vals))
            i_min, i_max = int(np.argmin(vals)), int(np.argmax(vals))
            rows.append({
                "diameter_um": d, "metric": metric, "n_common_flow_volumes": len(vals),
                "median_across_flows": med, "min_across_flows": float(vals[i_min]),
                "max_across_flows": float(vals[i_max]), "absolute_range": float(vals.max() - vals.min()),
                "range_percent_of_median": float(100 * (vals.max() - vals.min()) / med) if med else np.nan,
                "flow_at_min_mm_s": float(g.iloc[i_min].flow_mm_s), "flow_at_max_mm_s": float(g.iloc[i_max].flow_mm_s),
                "flow1_value": float(g.loc[g.flow_mm_s.eq(1), col].iloc[0]),
                "flow10_value": float(g.loc[g.flow_mm_s.eq(10), col].iloc[0]),
                "flow1_to_flow10_percent_difference": safe_pct(
                    float(g.loc[g.flow_mm_s.eq(10), col].iloc[0]),
                    float(g.loc[g.flow_mm_s.eq(1), col].iloc[0]),
                ),
                "values_by_flow_json": json.dumps({str(int(r.flow_mm_s)): float(getattr(r, col)) for r in g.itertuples(index=False)}, separators=(",", ":")),
            })
    return pd.DataFrame(rows)


def aggregate_d128_depth(depth_fw: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (scan, target), g in depth_fw.groupby(["scan_id", "target_r_um"], sort=True):
        s = stats(g["RI_r"])
        rows.append({
            "scan_id": scan, "diameter_um": 128.0, "flow_mm_s": float(g.flow_mm_s.iloc[0]),
            "target_r_um": float(target), "n_frames": int(len(g)),
            "selected_r_um_median": float(g.selected_r_um.median()),
            "max_abs_sampling_error_um": float(g.abs_sampling_error_um.max()),
            "RI_r_q1": s["q1"], "RI_r_median": s["median"], "RI_r_q3": s["q3"],
            "raw_row_signal_median": float(g.S_tail_raw_r.median()),
        })
    return pd.DataFrame(rows)


def load_new_depth() -> pd.DataFrame:
    d = pd.read_csv(NEW / "depth_summary_10um.csv")
    needed = {"scan_id", "diameter_um", "flow_mm_s", "target_r_um", "RI_r_q1", "RI_r_median", "RI_r_q3", "raw_row_signal_median"}
    if not needed.issubset(d.columns):
        raise AssertionError(f"New depth summary missing fields: {sorted(needed - set(d.columns))}")
    return d


def build_depth_outputs(depth_all: pd.DataFrame):
    common = depth_all[depth_all.diameter_um.isin(DIAMETERS) & depth_all.flow_mm_s.isin(COMMON_FLOWS)].copy()
    common = common.sort_values(["diameter_um", "flow_mm_s", "target_r_um"])
    if len(common) != 20 * 51:
        raise AssertionError(f"Expected 1020 common-grid depth rows, got {len(common)}")
    selected = common[common.target_r_um.isin(SELECTED_DEPTHS)].copy()
    if len(selected) != 20 * 7:
        raise AssertionError("Selected-depth common-grid row count mismatch")

    contrast_rows = []
    for (flow, target), g in common.groupby(["flow_mm_s", "target_r_um"], sort=True):
        by = g.set_index("diameter_um")
        for d0, d1 in CONTRASTS:
            v0, v1 = float(by.loc[d0, "RI_r_median"]), float(by.loc[d1, "RI_r_median"])
            contrast_rows.append({
                "flow_mm_s": float(flow), "target_r_um": float(target),
                "diameter_from_um": d0, "diameter_to_um": d1,
                "ri_r_median_from": v0, "ri_r_median_to": v1,
                "absolute_difference": v1 - v0, "percent_difference": safe_pct(v1, v0),
            })
    contrasts = pd.DataFrame(contrast_rows)

    bands = [(0, 100), (100, 200), (200, 300), (300, 500)]
    band_rows = []
    for scan, g in common.groupby("scan_id", sort=True):
        for lo, hi in bands:
            if lo == 0:
                b = g[(g.target_r_um >= lo) & (g.target_r_um <= hi)]
            else:
                b = g[(g.target_r_um > lo) & (g.target_r_um <= hi)]
            band_rows.append({
                "scan_id": scan, "diameter_um": float(g.diameter_um.iloc[0]), "flow_mm_s": float(g.flow_mm_s.iloc[0]),
                "depth_band_um": f"{lo}-{hi}", "n_depth_anchors": int(len(b)),
                "median_of_volume_RI_r_medians": float(b.RI_r_median.median()),
                "mean_of_volume_RI_r_medians": float(b.RI_r_median.mean()),
                "near_edge_RI_r_median": float(b.sort_values("target_r_um").RI_r_median.iloc[0]),
                "deep_edge_RI_r_median": float(b.sort_values("target_r_um").RI_r_median.iloc[-1]),
            })
    return common, selected, contrasts, pd.DataFrame(band_rows)


def build_slow_axis(fw: pd.DataFrame) -> pd.DataFrame:
    f = fw.copy()
    f["slow_axis_segment"] = (f.frame_index_0based.astype(int) // 100).clip(0, 4)
    rows = []
    for (scan, seg), g in f.groupby(["scan_id", "slow_axis_segment"], sort=True):
        rows.append({
            "scan_id": scan, "diameter_um": float(g.diameter_um.iloc[0]), "flow_mm_s": float(g.flow_mm_s.iloc[0]),
            "segment_index": int(seg), "segment_frames_0based": f"{int(seg)*100}-{int(seg)*100+99}",
            "valid_frame_count": int(len(g)), "nominal_frame_count": 100,
            "valid_fraction": float(len(g) / 100.0),
            "source_mean_raw_median": float(g.source_mean_raw.median()),
            "tail_mean_raw_500um_median": float(g.tail_mean_raw_500um.median()),
            "ri_tail_median": float(g.ri_tail.median()),
            "grid_scope": str(g.grid_scope.iloc[0]),
        })
    return pd.DataFrame(rows).sort_values(["diameter_um", "flow_mm_s", "segment_index"])


def build_geometry_qc(fw: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for scan, g in fw.groupby("scan_id", sort=True):
        w, x, z, a = (stats(g[c]) for c in ["apparent_width_um", "x4_px", "z_top_edge_px", "source_area_um2"])
        direct = g.direct_candidate_supported.astype(bool)
        rows.append({
            "scan_id": scan, "diameter_um": float(g.diameter_um.iloc[0]), "flow_mm_s": float(g.flow_mm_s.iloc[0]),
            "nominal_frames": 500, "geometry_valid_frames": int(len(g)), "geometry_valid_fraction": float(len(g) / 500),
            "direct_candidate_supported_count": int(direct.sum()),
            "short_gap_or_other_valid_count": int((~direct).sum()),
            "direct_candidate_supported_fraction_of_valid": float(direct.mean()),
            "x1_apparent_width_um_q1": w["q1"], "x1_apparent_width_um_median": w["median"], "x1_apparent_width_um_q3": w["q3"],
            "x4_px_q1": x["q1"], "x4_px_median": x["median"], "x4_px_q3": x["q3"],
            "z_top_edge_px_q1": z["q1"], "z_top_edge_px_median": z["median"], "z_top_edge_px_q3": z["q3"],
            "source_area_um2_q1": a["q1"], "source_area_um2_median": a["median"], "source_area_um2_q3": a["q3"],
            "rho_x1_vs_source": spearman_rho(g.apparent_width_um, g.source_mean_raw),
            "rho_x1_vs_tail500": spearman_rho(g.apparent_width_um, g.tail_mean_raw_500um),
            "rho_x1_vs_ri": spearman_rho(g.apparent_width_um, g.ri_tail),
            "rho_z_top_vs_source": spearman_rho(g.z_top_edge_px, g.source_mean_raw),
            "rho_z_top_vs_tail500": spearman_rho(g.z_top_edge_px, g.tail_mean_raw_500um),
            "rho_z_top_vs_ri": spearman_rho(g.z_top_edge_px, g.ri_tail),
            "rho_source_area_vs_source": spearman_rho(g.source_area_um2, g.source_mean_raw),
            "rho_source_area_vs_tail500": spearman_rho(g.source_area_um2, g.tail_mean_raw_500um),
            "rho_source_area_vs_ri": spearman_rho(g.source_area_um2, g.ri_tail),
            "grid_scope": str(g.grid_scope.iloc[0]),
        })
    return pd.DataFrame(rows).sort_values(["diameter_um", "flow_mm_s"])


def subset_masks(g: pd.DataFrame) -> dict[str, pd.Series]:
    wlo, whi = g.apparent_width_um.quantile([0.1, 0.9])
    zlo, zhi = g.z_top_edge_px.quantile([0.1, 0.9])
    central_w = g.apparent_width_um.between(wlo, whi, inclusive="both")
    central_z = g.z_top_edge_px.between(zlo, zhi, inclusive="both")
    direct = g.direct_candidate_supported.astype(bool)
    return {
        "all_valid": pd.Series(True, index=g.index),
        "direct_candidate_supported_only": direct,
        "central80_x1": central_w,
        "central80_z_top": central_z,
        "central80_x1_z_top": central_w & central_z,
        "direct_plus_central_geometry": direct & central_w & central_z,
    }


def build_restricted(fw: pd.DataFrame):
    rows = []
    for scan, g in fw.groupby("scan_id", sort=True):
        masks = subset_masks(g)
        wlo, whi = g.apparent_width_um.quantile([0.1, 0.9])
        zlo, zhi = g.z_top_edge_px.quantile([0.1, 0.9])
        for name, mask in masks.items():
            s = g[mask]
            rows.append({
                "scan_id": scan, "diameter_um": float(g.diameter_um.iloc[0]), "flow_mm_s": float(g.flow_mm_s.iloc[0]),
                "subset": name, "n_frames": int(len(s)), "fraction_of_valid_frames": float(len(s) / len(g)),
                "x1_p10_um": float(wlo), "x1_p90_um": float(whi), "z_top_p10_px": float(zlo), "z_top_p90_px": float(zhi),
                "source_mean_raw_median": float(s.source_mean_raw.median()) if len(s) else np.nan,
                "tail_mean_raw_500um_median": float(s.tail_mean_raw_500um.median()) if len(s) else np.nan,
                "ri_tail_median": float(s.ri_tail.median()) if len(s) else np.nan,
                "grid_scope": str(g.grid_scope.iloc[0]),
            })
    full = pd.DataFrame(rows)
    common = full[full.diameter_um.isin(DIAMETERS) & full.flow_mm_s.isin(COMMON_FLOWS)].copy()
    diam_rows = []
    for (subset, d), g in common.groupby(["subset", "diameter_um"], sort=True):
        diam_rows.append({
            "subset": subset, "diameter_um": float(d), "n_common_flow_volumes": int(len(g)),
            "source_mean_raw_median_across_flows": float(g.source_mean_raw_median.median()),
            "tail_mean_raw_500um_median_across_flows": float(g.tail_mean_raw_500um_median.median()),
            "ri_tail_median_across_flows": float(g.ri_tail_median.median()),
        })
    diam = pd.DataFrame(diam_rows)
    contrast_rows = []
    for subset in SUBSETS:
        g = diam[diam.subset.eq(subset)].set_index("diameter_um")
        for metric in ["source_mean_raw", "tail_mean_raw_500um", "ri_tail"]:
            col = f"{metric}_median_across_flows"
            for d0, d1 in CONTRASTS:
                v0, v1 = float(g.loc[d0, col]), float(g.loc[d1, col])
                contrast_rows.append({
                    "subset": subset, "metric": metric, "diameter_from_um": d0, "diameter_to_um": d1,
                    "value_from": v0, "value_to": v1, "absolute_difference": v1 - v0,
                    "percent_difference": safe_pct(v1, v0),
                })
    return full, diam, pd.DataFrame(contrast_rows)


def build_d500_audit(restricted: pd.DataFrame, slow: pd.DataFrame, geometry: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, r in restricted[restricted.diameter_um.eq(500)].iterrows():
        flow, subset = float(r.flow_mm_s), str(r.subset)
        seg = slow[(slow.diameter_um.eq(500)) & (slow.flow_mm_s.eq(flow))]
        geom = geometry[(geometry.diameter_um.eq(500)) & (geometry.flow_mm_s.eq(flow))].iloc[0]
        item = {
            "flow_mm_s": flow, "subset": subset, "n_frames": int(r.n_frames),
            "source_mean_raw_median": float(r.source_mean_raw_median),
            "tail_mean_raw_500um_median": float(r.tail_mean_raw_500um_median),
            "ri_tail_median": float(r.ri_tail_median),
            "source_segment_min_median": float(seg.source_mean_raw_median.min()),
            "source_segment_max_median": float(seg.source_mean_raw_median.max()),
            "source_segment_range_pct_of_volume_subset_median": float(100 * (seg.source_mean_raw_median.max() - seg.source_mean_raw_median.min()) / r.source_mean_raw_median) if r.source_mean_raw_median else np.nan,
            "x1_apparent_width_um_median": float(geom.x1_apparent_width_um_median),
            "z_top_edge_px_median": float(geom.z_top_edge_px_median),
            "source_area_um2_median": float(geom.source_area_um2_median),
            "direct_candidate_supported_fraction_of_valid": float(geom.direct_candidate_supported_fraction_of_valid),
            "rho_x1_vs_source": float(geom.rho_x1_vs_source),
            "rho_z_top_vs_source": float(geom.rho_z_top_vs_source),
            "rho_source_area_vs_source": float(geom.rho_source_area_vs_source),
            "grid_scope": str(r.grid_scope),
        }
        if flow in COMMON_FLOWS:
            for comp in [235, 285]:
                c = restricted[(restricted.diameter_um.eq(comp)) & (restricted.flow_mm_s.eq(flow)) & (restricted.subset.eq(subset))]
                if len(c):
                    c = c.iloc[0]
                    item[f"source_vs_d{comp}_percent_difference"] = safe_pct(float(r.source_mean_raw_median), float(c.source_mean_raw_median))
                    item[f"tail500_vs_d{comp}_percent_difference"] = safe_pct(float(r.tail_mean_raw_500um_median), float(c.tail_mean_raw_500um_median))
                    item[f"ri_vs_d{comp}_percent_difference"] = safe_pct(float(r.ri_tail_median), float(c.ri_tail_median))
        rows.append(item)
    return pd.DataFrame(rows).sort_values(["flow_mm_s", "subset"])


def build_spatial_robustness(slow: pd.DataFrame) -> dict:
    common = slow[slow.diameter_um.isin(DIAMETERS) & slow.flow_mm_s.isin(COMMON_FLOWS)]
    ordered_128_235_285 = 0
    big_drop_then_small = 0
    d500_ri_rebound = 0
    d500_source_lower_d285 = 0
    total = 0
    details = []
    for (flow, seg), g in common.groupby(["flow_mm_s", "segment_index"], sort=True):
        by = g.set_index("diameter_um")
        if not all(d in by.index for d in DIAMETERS):
            continue
        r128, r235, r285, r500 = [float(by.loc[d, "ri_tail_median"]) for d in DIAMETERS]
        s285, s500 = float(by.loc[285, "source_mean_raw_median"]), float(by.loc[500, "source_mean_raw_median"])
        ordered = r128 > r235 > r285
        drop1 = abs(safe_pct(r235, r128))
        drop2 = abs(safe_pct(r285, r235))
        big_small = ordered and drop1 > drop2
        rebound = r500 > r285
        source_lower = s500 < s285
        ordered_128_235_285 += int(ordered)
        big_drop_then_small += int(big_small)
        d500_ri_rebound += int(rebound)
        d500_source_lower_d285 += int(source_lower)
        total += 1
        details.append({
            "flow_mm_s": float(flow), "segment_index": int(seg),
            "ri_128": r128, "ri_235": r235, "ri_285": r285, "ri_500": r500,
            "ordered_128_gt_235_gt_285": ordered,
            "abs_pct_drop_128_to_235": drop1, "abs_pct_drop_235_to_285": drop2,
            "larger_128_to_235_drop": big_small,
            "d500_ri_gt_d285": rebound,
            "d500_source_lt_d285": source_lower,
        })
    return {
        "total_common_flow_segments": total,
        "segments_with_ri_128_gt_235_gt_285": ordered_128_235_285,
        "segments_with_larger_128_to_235_than_235_to_285_drop": big_drop_then_small,
        "segments_with_d500_ri_gt_d285": d500_ri_rebound,
        "segments_with_d500_source_lt_d285": d500_source_lower_d285,
        "details": details,
    }


def build_restricted_robustness(contrasts: pd.DataFrame) -> dict:
    out = {}
    for subset in SUBSETS:
        g = contrasts[(contrasts.subset.eq(subset)) & contrasts.metric.eq("ri_tail")]
        c1 = float(g[(g.diameter_from_um.eq(128)) & (g.diameter_to_um.eq(235))].percent_difference.iloc[0])
        c2 = float(g[(g.diameter_from_um.eq(235)) & (g.diameter_to_um.eq(285))].percent_difference.iloc[0])
        c3 = float(g[(g.diameter_from_um.eq(285)) & (g.diameter_to_um.eq(500))].percent_difference.iloc[0])
        out[subset] = {
            "ri_128_to_235_percent": c1,
            "ri_235_to_285_percent": c2,
            "ri_285_to_500_percent": c3,
            "decline_128_to_235_to_285_preserved": bool(c1 < 0 and c2 < 0),
            "larger_first_drop_preserved": bool(abs(c1) > abs(c2)),
            "d500_rebound_preserved": bool(c3 > 0),
        }
    return out


def make_figures(common: pd.DataFrame, depth_common: pd.DataFrame) -> None:
    # 1. RI common grid: one trace per flow.
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for flow in COMMON_FLOWS:
        g = common[common.flow_mm_s.eq(flow)].sort_values("diameter_um")
        ax.plot(g.diameter_um, g.ri_tail_median, marker="o", label=f"{flow} mm/s")
    ax.set_xlabel("Physical vessel diameter (um)")
    ax.set_ylabel("Volume-median RI_tail")
    ax.set_xticks(DIAMETERS)
    ax.grid(True, alpha=0.3)
    ax.legend(title="Flow", frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "figure_common_grid_ri_tail.png", dpi=220)
    plt.close(fig)

    for metric, ylabel, filename in [
        ("source_mean_raw_median", "Volume-median vessel raw SV", "figure_source_by_diameter.png"),
        ("tail_mean_raw_500um_median", "Volume-median 500 um tail raw SV", "figure_tail500_by_diameter.png"),
        ("ri_tail_median", "Volume-median RI_tail", "figure_ri_tail_by_diameter.png"),
    ]:
        fig, ax = plt.subplots(figsize=(7.2, 4.8))
        for flow in COMMON_FLOWS:
            g = common[common.flow_mm_s.eq(flow)].sort_values("diameter_um")
            ax.plot(g.diameter_um, g[metric], marker="o", label=f"{flow} mm/s")
        ax.set_xlabel("Physical vessel diameter (um)")
        ax.set_ylabel(ylabel)
        ax.set_xticks(DIAMETERS)
        ax.grid(True, alpha=0.3)
        ax.legend(title="Flow", frameon=False)
        fig.tight_layout()
        fig.savefig(OUT / filename, dpi=220)
        plt.close(fig)

    # Across-common-flow descriptive median depth profile by diameter.
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for d in DIAMETERS:
        g = depth_common[depth_common.diameter_um.eq(d)]
        pivot = g.pivot(index="target_r_um", columns="flow_mm_s", values="RI_r_median").sort_index()
        center = pivot.median(axis=1)
        lo = pivot.min(axis=1)
        hi = pivot.max(axis=1)
        line = ax.plot(center.index, center.values, label=f"D{d}")[0]
        ax.fill_between(center.index.to_numpy(float), lo.to_numpy(float), hi.to_numpy(float), alpha=0.12, color=line.get_color())
    ax.set_xlabel("Distance below physical vessel bottom (um)")
    ax.set_ylabel("RI(r), median across common-flow volume medians")
    ax.grid(True, alpha=0.3)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "figure_depth_profiles_common_grid.png", dpi=220)
    plt.close(fig)


def build_readme(common, flowvar, depth_contrasts, spatial, restricted_robust, d500_audit, d128_audit, validation):
    ri_by_d = common.groupby("diameter_um").ri_tail_median.median().to_dict()
    source_by_d = common.groupby("diameter_um").source_mean_raw_median.median().to_dict()
    tail_by_d = common.groupby("diameter_um").tail_mean_raw_500um_median.median().to_dict()
    d1 = safe_pct(ri_by_d[235], ri_by_d[128])
    d2 = safe_pct(ri_by_d[285], ri_by_d[235])
    rebound = safe_pct(ri_by_d[500], ri_by_d[285])
    d500_source_vs_285 = safe_pct(source_by_d[500], source_by_d[285])
    d500_tail_vs_285 = safe_pct(tail_by_d[500], tail_by_d[285])

    # Depth: median fixed-flow contrast at each target.
    dep = depth_contrasts.groupby(["target_r_um", "diameter_from_um", "diameter_to_um"], as_index=False).percent_difference.median()
    c128235 = dep[(dep.diameter_from_um.eq(128)) & dep.diameter_to_um.eq(235)]
    c235285 = dep[(dep.diameter_from_um.eq(235)) & dep.diameter_to_um.eq(285)]
    c285500 = dep[(dep.diameter_from_um.eq(285)) & dep.diameter_to_um.eq(500)]
    near_128235 = float(c128235.loc[c128235.target_r_um.eq(0), "percent_difference"].iloc[0])
    deep_128235 = float(c128235.loc[c128235.target_r_um.eq(500), "percent_difference"].iloc[0])
    near_235285 = float(c235285.loc[c235285.target_r_um.eq(0), "percent_difference"].iloc[0])
    deep_235285 = float(c235285.loc[c235285.target_r_um.eq(500), "percent_difference"].iloc[0])
    near_285500 = float(c285500.loc[c285500.target_r_um.eq(0), "percent_difference"].iloc[0])
    deep_285500 = float(c285500.loc[c285500.target_r_um.eq(500), "percent_difference"].iloc[0])

    ri_flow = flowvar[flowvar.metric.eq("ri_tail_500um")].set_index("diameter_um")
    restricted_ok = sum(v["decline_128_to_235_to_285_preserved"] for v in restricted_robust.values())
    rebound_ok = sum(v["d500_rebound_preserved"] for v in restricted_robust.values())
    d500_common_all = d500_audit[(d500_audit.subset.eq("all_valid")) & d500_audit.flow_mm_s.isin(COMMON_FLOWS)]
    source_vs_d285_med = float(d500_common_all.source_vs_d285_percent_difference.median())
    tail_vs_d285_med = float(d500_common_all.tail500_vs_d285_percent_difference.median())
    ri_vs_d285_med = float(d500_common_all.ri_vs_d285_percent_difference.median())

    lines = [
        "# SV 人工血管拖尾直径分析 — Stage 2 scientific analysis v1",
        "",
        "本目录只分析已经冻结的 no-background raw-SV 结果，不修改 Stage 1 定量定义、连续性优先 v2.1 几何或 D128 冻结文件。独立实验单位仍然是 **scan volume**；B-scan 只表示同一 volume 内的 slow-axis 空间位置。",
        "",
        "## 1. 最主要的现象",
        "",
        f"在 1/3/5/7/10 mm/s 五个共同流速下，volume-median RI_tail 再取中位数依次为：D128 **{ri_by_d[128]:.6f}**、D235 **{ri_by_d[235]:.6f}**、D285 **{ri_by_d[285]:.6f}**、D500 **{ri_by_d[500]:.6f}**。",
        "",
        f"也就是说，128→235 um 时 RI_tail 下降 **{abs(d1):.2f}%**，235→285 um 再下降 **{abs(d2):.2f}%**。第一段下降明显更大。D500 相对 D285 又回升 **{rebound:.2f}%**。这个结果仍然是描述性扫描卷比较，不把它写成直径的因果效应。",
        "",
        "## 2. D500 为什么 RI 又高起来",
        "",
        f"D500 的 vessel raw SV 相比 D285 低 **{abs(d500_source_vs_285):.2f}%**，500 um absolute tail raw SV 也低 **{abs(d500_tail_vs_285):.2f}%**。因此 D500 并不是“绝对拖尾变强”。它的分母下降得更明显，使 tail/source 这个比值重新升高。",
        "",
        f"逐个共同流速直接比较时，D500 相对 D285 的 source 变化中位数为 **{source_vs_d285_med:+.2f}%**，tail 变化中位数为 **{tail_vs_d285_med:+.2f}%**，RI 变化中位数为 **{ri_vs_d285_med:+.2f}%**。这说明 D500 的 RI 回升主要要从低 vessel denominator 来理解，同时 numerator 本身没有升高。",
        "",
        "## 3. 这种直径排序是不是只出现在 volume 的某一小段",
        "",
        f"把每个 volume 的 500 个 slow-axis 位置固定分成五段后，共得到 **{spatial['total_common_flow_segments']}** 个“共同流速 × 空间段”比较。D128 > D235 > D285 的 RI 排序出现在 **{spatial['segments_with_ri_128_gt_235_gt_285']}/{spatial['total_common_flow_segments']}** 段；其中 128→235 的下降幅度大于 235→285 的段数为 **{spatial['segments_with_larger_128_to_235_than_235_to_285_drop']}/{spatial['total_common_flow_segments']}**。",
        "",
        f"D500 RI 高于 D285 出现在 **{spatial['segments_with_d500_ri_gt_d285']}/{spatial['total_common_flow_segments']}** 段；D500 vessel source 低于 D285 出现在 **{spatial['segments_with_d500_source_lt_d285']}/{spatial['total_common_flow_segments']}** 段。这里的‘段’仍然只是空间稳健性检查，不是新的独立重复。",
        "",
        "## 4. 限制几何和定位来源以后，主现象还在不在",
        "",
        f"预先定义的 6 个子集（全部有效帧、direct-only、X1 中央 80%、z_top 中央 80%、X1+z_top 中央 80%、direct+central geometry）中，128→235→285 的 RI 递减保留在 **{restricted_ok}/6** 个子集；D500 相对 D285 的 RI 回升保留在 **{rebound_ok}/6** 个子集。",
        "",
        "这一步的含义很直接：如果删掉几何极端位置、只保留直接候选支持帧后排序仍然存在，那么主要现象就不容易用‘某个直径组恰好收进了更多异常定位帧’来解释。具体每个子集的数值见 `restricted_subset_summary.csv` 和 `restricted_subset_contrasts.csv`。",
        "",
        "## 5. 直径差异随拖尾深度怎么变化",
        "",
        f"对每个固定流速先比较 RI(r)，再对五个共同流速的百分比差异取中位数：在 vessel bottom 附近（0 um），128→235 为 **{near_128235:+.2f}%**、235→285 为 **{near_235285:+.2f}%**；到 500 um 深度分别为 **{deep_128235:+.2f}%** 和 **{deep_235285:+.2f}%**。",
        "",
        f"D500 相对 D285 的 RI(r) 差异在 0 um 为 **{near_285500:+.2f}%**，在 500 um 为 **{deep_285500:+.2f}%**。深部仍然存在相对差异时，应理解为 raw-SV floor 与 vessel denominator 共同决定的 normalized floor；这里没有做背景扣除，所以不能把深部非零 RI(r) 当成拖尾终点或 detection depth。",
        "",
        "## 6. Flow 的描述性结构",
        "",
    ]
    for d in DIAMETERS:
        r = ri_flow.loc[d]
        lines.append(f"- D{d}: 五个共同流速的 RI_tail 范围为 **{r.min_across_flows:.6f}–{r.max_across_flows:.6f}**，相对于该直径跨流速中位数的范围约 **{r.range_percent_of_median:.1f}%**；最低点在 {r.flow_at_min_mm_s:g} mm/s，最高点在 {r.flow_at_max_mm_s:g} mm/s。")
    lines += [
        "",
        "不同直径的最高/最低流速位置并不统一，因此当前数据没有显示一个可以跨直径直接概括成‘随 flow 单调增加/降低’的共同结构。每个 diameter×flow 只有一个 independent volume，这里不进行 flow 的显著性检验或因果解释。",
        "",
        "## 7. D128 100/200/300 um 窗口如何补齐",
        "",
        "D128 历史 no-background 主表只保存了 500 um scalar，但冻结 release 保存了逐帧二维 `sv_raw`。本 Stage 2 仅在新目录中重新读取这些数组，用同一 X1/X4/z_top、真实 128 um 直径和同一 `quantify_raw_sv` 计算 100/200/300/500 um；不覆盖任何 D128 文件。",
        "",
        f"2422/2422 帧完成只读 replay。与冻结 D128 500 um 结果相比，source/tail/RI 最大相对误差分别为 **{d128_audit['source_mean_raw_max_relative_error_vs_frozen']:.3e} / {d128_audit['tail_mean_raw_500um_max_relative_error_vs_frozen']:.3e} / {d128_audit['ri_tail_max_relative_error_vs_frozen']:.3e}**，均低于 1e-10。",
        "",
        "## 8. 当前可以支持的 Stage 2 判断",
        "",
        "1. **128→235→285 的相对拖尾下降是一个跨共同流速、跨 slow-axis 空间段并对几何/QC 限制具有稳健性的描述性模式。**",
        "2. **主要下降集中在 128→235；235→285 是较小的进一步下降。**",
        "3. **D500 的 RI 回升不代表 absolute tail raw SV 更高；低 vessel raw SV denominator 是核心组成因素。**",
        "4. **直径差异在靠近 vessel 的区域更容易拉开，深部逐渐受到 raw-SV floor / vessel denominator 的影响。**",
        "5. **不同直径下的 flow 变化形态不完全一致，当前不建立统一 Flow×Diameter 模型。**",
        "6. **几何/QC 组成不能单独解释主要直径排序。D500 denominator 低值若在所有限制分析中持续存在，应把它视为需要进一步做 source-ellipse 内部信号分布审计的信号层特征，而不是通过重新设计 source ROI 去消除。**",
        "",
        "## 9. 输出文件",
        "",
        "- `common_grid_volume_metrics.csv`: 4×5 主网格的 source / 100–500 um tail / RI。",
        "- `diameter_contrasts_by_flow.csv`: 每个固定流速下的直径差值和百分比差。",
        "- `flow_variation_by_diameter.csv`: 每个固定直径在五个共同流速中的变化范围。",
        "- `source_tail_ratio_decomposition.csv`: source、absolute 500 um tail、RI 并排。",
        "- `depth_profile_summary.csv`, `depth_selected_common_grid.csv`, `depth_diameter_contrasts.csv`, `depth_band_summary.csv`: RI(r) 深度结果。",
        "- `slow_axis_segment_summary.csv`, `slow_axis_robustness_details.csv`: volume 内空间稳健性。",
        "- `geometry_qc_summary.csv`: X1/X4/z_top/source area 与 source/tail/RI 的描述性关系。",
        "- `restricted_subset_summary.csv`, `restricted_subset_diameter_summary.csv`, `restricted_subset_contrasts.csv`: 预定义限制子集。",
        "- `d500_denominator_audit.csv`: D500 denominator 的 flow/segment/geometry/QC 审计。",
        "- `d500_extra_flow_extension.csv`: D500 的 2/9/12 mm/s 次级扩展。",
        "- PNG figures、`validation.json`、`provenance.json`。",
        "",
        f"Validation status: **{validation['status']}**. No t-test, ANOVA, p-value, mixed model, regression surface, nonlinear fit or detection-depth estimate was computed.",
    ]
    (OUT / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    new_fw = load_new_framewise()
    d128_fw, d128_depth_fw, d128_audit = replay_d128()
    fw = pd.concat([d128_fw, new_fw], ignore_index=True, sort=False)
    if len(fw) != 2422 + 8509:
        raise AssertionError("Combined valid-frame count mismatch")

    volume = build_volume_metrics(fw)
    if len(volume) != 23:
        raise AssertionError(f"Expected 23 total volumes, got {len(volume)}")
    common = volume[volume.diameter_um.isin(DIAMETERS) & volume.flow_mm_s.isin(COMMON_FLOWS)].copy()
    if len(common) != 20 or common.duplicated(["diameter_um", "flow_mm_s"]).any():
        raise AssertionError("4x5 common grid is not complete and unique")
    new_volume_max_error = validate_new_volume_replay(volume)
    if new_volume_max_error > 1e-12:
        raise AssertionError(f"Stage2 new-diameter volume replay mismatch: {new_volume_max_error}")

    contrasts = build_contrasts(common)
    flowvar = build_flow_variation(common)
    source_tail_ratio = common[[
        "scan_id", "diameter_um", "flow_mm_s", "source_mean_raw_median",
        "tail_mean_raw_500um_median", "ri_tail_median", "geometry_valid_frames", "geometry_valid_fraction"
    ]].copy()
    for flow in COMMON_FLOWS:
        ref = source_tail_ratio[(source_tail_ratio.flow_mm_s.eq(flow)) & source_tail_ratio.diameter_um.eq(128)].iloc[0]
        idx = source_tail_ratio.flow_mm_s.eq(flow)
        source_tail_ratio.loc[idx, "source_vs_d128_same_flow_pct"] = source_tail_ratio.loc[idx, "source_mean_raw_median"].map(lambda v: safe_pct(v, ref.source_mean_raw_median))
        source_tail_ratio.loc[idx, "tail500_vs_d128_same_flow_pct"] = source_tail_ratio.loc[idx, "tail_mean_raw_500um_median"].map(lambda v: safe_pct(v, ref.tail_mean_raw_500um_median))
        source_tail_ratio.loc[idx, "ri_vs_d128_same_flow_pct"] = source_tail_ratio.loc[idx, "ri_tail_median"].map(lambda v: safe_pct(v, ref.ri_tail_median))

    new_depth = load_new_depth()
    d128_depth = aggregate_d128_depth(d128_depth_fw)
    depth_all = pd.concat([d128_depth, new_depth], ignore_index=True, sort=False)
    depth_common, depth_selected, depth_contrasts, depth_bands = build_depth_outputs(depth_all)

    slow = build_slow_axis(fw)
    geometry = build_geometry_qc(fw)
    restricted, restricted_diam, restricted_contrasts = build_restricted(fw)
    d500_audit = build_d500_audit(restricted, slow, geometry)
    spatial = build_spatial_robustness(slow)
    restricted_robust = build_restricted_robustness(restricted_contrasts)

    d500_extra = volume[(volume.diameter_um.eq(500)) & volume.flow_mm_s.isin([2, 9, 12])].copy()

    write_csv("common_grid_volume_metrics.csv", common)
    write_csv("all_volume_metrics_with_d500_extension.csv", volume)
    write_csv("diameter_contrasts_by_flow.csv", contrasts)
    write_csv("flow_variation_by_diameter.csv", flowvar)
    write_csv("source_tail_ratio_decomposition.csv", source_tail_ratio)
    write_csv("depth_profile_summary.csv", depth_common)
    write_csv("depth_selected_common_grid.csv", depth_selected)
    write_csv("depth_diameter_contrasts.csv", depth_contrasts)
    write_csv("depth_band_summary.csv", depth_bands)
    write_csv("slow_axis_segment_summary.csv", slow)
    write_csv("slow_axis_robustness_details.csv", pd.DataFrame(spatial["details"]))
    write_csv("geometry_qc_summary.csv", geometry)
    write_csv("restricted_subset_summary.csv", restricted)
    write_csv("restricted_subset_diameter_summary.csv", restricted_diam)
    write_csv("restricted_subset_contrasts.csv", restricted_contrasts)
    write_csv("d500_denominator_audit.csv", d500_audit)
    write_csv("d500_extra_flow_extension.csv", d500_extra)

    make_figures(common, depth_common)

    validation = {
        "status": "passed",
        "common_grid_volumes": 20,
        "total_volumes_with_d500_extension": 23,
        "new_diameter_valid_frames": 8509,
        "d128_valid_frames": 2422,
        "combined_valid_spatial_frames": int(len(fw)),
        "common_grid_depth_rows": int(len(depth_common)),
        "selected_depth_rows": int(len(depth_selected)),
        "slow_axis_segment_rows": int(len(slow)),
        "restricted_subset_rows": int(len(restricted)),
        "new_diameter_volume_reaggregation_max_relative_error": float(new_volume_max_error),
        "d128_supplemental_replay": d128_audit,
        "spatial_robustness": {k: v for k, v in spatial.items() if k != "details"},
        "restricted_subset_robustness": restricted_robust,
        "statistical_boundary": {
            "experimental_unit": "scan volume",
            "bscan": "slow-axis spatial sample; not independent replicate",
            "inferential_statistics_computed": False,
            "p_values_computed": False,
            "regression_or_curve_fit_computed": False,
            "detection_depth_computed": False,
        },
    }
    write_json("validation.json", validation)

    input_paths = [
        NEW / "framewise_primary.csv", NEW / "volume_summary.csv", NEW / "depth_summary_10um.csv",
        NEW / "depth_selected.csv", NEW / "validation.json", NEW / "provenance.json",
        D128_NOBG / "framewise_primary.csv", D128_NOBG / "scan_flow_summary.csv",
        D128_NOBG / "depth_full_10um.csv", D128_NOBG / "depth_selected.csv",
        D128_FORMAL / "localization.csv", D128_FORMAL / "arrays_sha256.csv", D128_FORMAL / "download_packages.csv",
    ]
    provenance = {
        "repository": "bulbel-magnolia/OCTA-",
        "branch": "analysis/sv-diameter-stage2",
        "workflow_source_sha": os.environ.get("GITHUB_SHA", "local"),
        "script": "scripts/analyze_sv_diameter_stage2.py",
        "script_sha256": sha256_file(Path(__file__)),
        "frozen_input_sha256": {str(p.relative_to(ROOT)): sha256_file(p) for p in input_paths},
        "signal_definition": {
            "sv_raw": "var(abs(E), 1, 3)", "background_subtraction": False,
            "normalization": False, "log": False, "gain": False,
            "clipping": False, "positive_truncation": False,
        },
        "geometry_definition": {
            "lateral_center": "X4", "lateral_span": "X1 apparent width",
            "axial_source_height": "physical vessel diameter", "tail_guard_um": 0,
            "tail_windows_um": WINDOWS, "primary_tail_window_um": 500,
            "dx_um": 12.7, "dz_um": 6.7,
        },
        "d128_supplemental_note": "Read-only replay of released 2-D sv_raw arrays in Stage2 output only; frozen D128 files unchanged. 500um values validated against existing no-background D128 primary output before supplemental 100/200/300um values were accepted.",
        "experimental_unit": "scan volume",
    }
    write_json("provenance.json", provenance)

    build_readme(common, flowvar, depth_contrasts, spatial, restricted_robust, d500_audit, d128_audit, validation)

    # A small machine-readable finding summary is useful for review without parsing prose.
    finding_summary = {
        "ri_tail_median_across_common_flows_by_diameter": {str(int(d)): float(v) for d, v in common.groupby("diameter_um").ri_tail_median.median().items()},
        "spatial_robustness": {k: v for k, v in spatial.items() if k != "details"},
        "restricted_subset_robustness": restricted_robust,
    }
    write_json("findings_summary.json", finding_summary)
    print(json.dumps({"status": "passed", "output_dir": str(OUT), "common_grid_volumes": 20}, ensure_ascii=False))


if __name__ == "__main__":
    main()
