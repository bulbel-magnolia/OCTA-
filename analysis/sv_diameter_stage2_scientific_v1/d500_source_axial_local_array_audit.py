#!/usr/bin/env python3
"""Secondary local-array audit for the low D500 frozen source denominator.

This script does NOT redefine the frozen source ROI. It subdivides the same
fractional source ellipse into three relative axial thirds (upper/middle/lower)
and reports raw-SV means and contributions. The three band weights are built
with the same 16x16 deterministic supersampling convention as the frozen
ellipse implementation, and must exactly reconstruct the full frozen ellipse
at the supersample level.

Required input is the original per-frame MAT export containing `sv_raw` for the
new-diameter scans. Those MAT files were deliberately not committed during the
intake/formal-quantification stage; point --mat-root to the retained local
interim directory when available.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from svrecttail.geometry import VesselGeometry, ellipse_weights, interval_overlap_weights

FORMAL = ROOT / "analysis/formal_sv_diameter_v1/framewise_primary.csv"
OUT_DEFAULT = ROOT / "analysis/sv_diameter_stage2_scientific_v1/d500_source_axial_audit"
DX_UM = 12.7
DZ_UM = 6.7
SUPERSAMPLE = 16
BANDS = (("upper", 0.0, 1/3), ("middle", 1/3, 2/3), ("lower", 2/3, 1.0))
TOL = 1e-10


def sha256(path: Path) -> str:
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def relerr(a: float, b: float) -> float:
    return abs(float(a) - float(b)) / max(abs(float(b)), np.finfo(float).tiny)


def find_mat(mat_root: Path, scan_id: str, frame: int) -> Path:
    name = f"frame_{frame:03d}.mat"
    candidates = [
        mat_root / scan_id / name,
        mat_root / f"{scan_id}_interim" / name,
        mat_root / name,
    ]
    hits = [p for p in candidates if p.is_file()]
    if not hits:
        hits = list(mat_root.glob(f"**/{scan_id}/**/{name}"))
    if len(hits) != 1:
        raise FileNotFoundError(f"Expected one MAT for {scan_id}/{name}; found {len(hits)}")
    return hits[0]


def geometry_from_row(r) -> VesselGeometry:
    return VesselGeometry(
        x_left_edge_px=float(r.x_left_edge_px),
        x_right_edge_px=float(r.x_right_edge_px),
        z_top_edge_px=float(r.z_top_edge_px),
        diameter_um=float(r.physical_diameter_um),
        dx_um=DX_UM,
        dz_um=DZ_UM,
    )


def band_weights(shape: tuple[int, int], g: VesselGeometry, lo: float, hi: float) -> np.ndarray:
    """Fractional ellipse weights restricted to one relative axial band.

    The same 16x16 subpixel samples used by ellipse_weights are partitioned by
    normalized physical axial coordinate u=(z-z_top)/diameter. Thus the three
    band arrays sum exactly to the full source ellipse at the supersample level.
    """
    nz, nx = shape
    w = np.zeros((nz, nx), dtype=np.float64)
    x_overlap = interval_overlap_weights(nx, g.x_left_edge_px, g.x_right_edge_px)
    z_overlap = interval_overlap_weights(nz, g.z_top_edge_px, g.z_bottom_edge_px)
    xi = np.flatnonzero(x_overlap > 0)
    zi = np.flatnonzero(z_overlap > 0)
    if not len(xi) or not len(zi):
        return w
    offsets = (np.arange(SUPERSAMPLE, dtype=float) + 0.5) / SUPERSAMPLE - 0.5
    xs = xi[:, None] + offsets[None, :]
    zs = zi[:, None] + offsets[None, :]
    xr = g.lateral_width_um / 2.0
    zr = g.diameter_um / 2.0
    xterm = ((xs - g.x_center_px) * g.dx_um / xr) ** 2
    zterm = ((zs - g.z_center_px) * g.dz_um / zr) ** 2
    inside = zterm[:, None, :, None] + xterm[None, :, None, :] <= 1.0
    u = ((zs - g.z_top_edge_px) * g.dz_um / g.diameter_um)
    if hi == 1.0:
        inband_z = (u >= lo) & (u <= hi)
    else:
        inband_z = (u >= lo) & (u < hi)
    local = (inside & inband_z[:, None, :, None]).mean(axis=(2, 3), dtype=np.float64)
    w[np.ix_(zi, xi)] = local
    return w


def weighted_stats(sv: np.ndarray, w: np.ndarray, pixel_area: float) -> tuple[float, float, float]:
    m = w > 0
    area = float(w.sum() * pixel_area)
    q = float(np.sum(sv[m] * w[m], dtype=np.float64) * pixel_area)
    return area, q, q / area


def qtls(x):
    a = np.asarray(x, dtype=float)
    q1, med, q3 = np.quantile(a, [0.25, 0.5, 0.75])
    return float(q1), float(med), float(q3)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mat-root", required=True, help="Retained local per-frame MAT export root")
    ap.add_argument("--output-dir", default=str(OUT_DEFAULT))
    ap.add_argument("--diameters", default="285,500", help="Comparator diameters; default D285+D500")
    args = ap.parse_args()

    mat_root = Path(args.mat_root).expanduser().resolve()
    out = Path(args.output_dir).expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)
    diameters = {float(x) for x in args.diameters.split(",")}

    fw = pd.read_csv(FORMAL)
    fw = fw[fw.diameter_um.isin(diameters)].copy()
    if fw.empty:
        raise RuntimeError("No requested diameters in frozen framewise table")

    # Recreate the same central-geometry restrictions used in Stage 2.
    fw["central80_x1"] = False
    fw["central80_z_top"] = False
    for scan, g in fw.groupby("scan_id"):
        xlo, xhi = np.quantile(g.apparent_width_um, [0.1, 0.9])
        zlo, zhi = np.quantile(g.z_top_edge_px, [0.1, 0.9])
        fw.loc[g.index, "central80_x1"] = g.apparent_width_um.between(xlo, xhi, inclusive="both")
        fw.loc[g.index, "central80_z_top"] = g.z_top_edge_px.between(zlo, zhi, inclusive="both")
    fw["direct_only"] = fw.localization_source.eq("direct_candidate")
    fw["central_geometry"] = fw.central80_x1 & fw.central80_z_top
    fw["direct_plus_central"] = fw.direct_only & fw.central_geometry

    pixel_area = DX_UM * DZ_UM
    rows = []
    reconstruction_errors = []
    sha_verified = 0

    for r in fw.itertuples(index=False):
        p = find_mat(mat_root, str(r.scan_id), int(r.frame_index_0based))
        digest = sha256(p)
        if digest != str(r.input_mat_sha256):
            raise RuntimeError(f"MAT SHA mismatch: {r.scan_id}/{p.name}")
        sha_verified += 1
        mat = loadmat(p, variable_names=["sv_raw"], simplify_cells=True)
        sv = np.asarray(mat["sv_raw"], dtype=np.float64)
        if sv.shape != (351, 500) or not np.isfinite(sv).all():
            raise RuntimeError(f"Invalid sv_raw: {r.scan_id}/{p.name}")

        geom = geometry_from_row(r)
        full = ellipse_weights(sv.shape, geom, supersample=SUPERSAMPLE)
        bw = {name: band_weights(sv.shape, geom, lo, hi) for name, lo, hi in BANDS}
        summed = bw["upper"] + bw["middle"] + bw["lower"]
        if not np.array_equal(summed, full):
            # Floating means of identical boolean partitions should be exact;
            # keep a numeric fallback with a strict machine-level gate.
            if float(np.max(np.abs(summed - full))) > 1e-15:
                raise RuntimeError("Axial band weights do not reconstruct frozen ellipse")

        full_area, full_q, full_mean = weighted_stats(sv, full, pixel_area)
        errs = {
            "source_area": relerr(full_area, r.source_area_um2),
            "source_q": relerr(full_q, r.source_q_raw),
            "source_mean": relerr(full_mean, r.source_mean_raw),
        }
        reconstruction_errors.append(max(errs.values()))
        if max(errs.values()) > TOL:
            raise RuntimeError(f"Frozen source replay mismatch {r.scan_id}/{p.name}: {errs}")

        rec = dict(
            scan_id=r.scan_id, diameter_um=float(r.diameter_um), flow_mm_s=float(r.flow_mm_s),
            frame_index_0based=int(r.frame_index_0based), slow_axis_segment=int(r.frame_index_0based)//100,
            source_mean_raw_frozen=float(r.source_mean_raw), source_q_raw_frozen=float(r.source_q_raw),
            source_area_um2_frozen=float(r.source_area_um2), localization_source=r.localization_source,
            central80_x1=bool(r.central80_x1), central80_z_top=bool(r.central80_z_top),
            direct_only=bool(r.direct_only), central_geometry=bool(r.central_geometry),
            direct_plus_central=bool(r.direct_plus_central), input_mat_sha256=digest,
        )
        qsum = 0.0
        for name, _, _ in BANDS:
            area, q, mean = weighted_stats(sv, bw[name], pixel_area)
            rec[f"{name}_area_um2"] = area
            rec[f"{name}_q_raw"] = q
            rec[f"{name}_mean_raw"] = mean
            rec[f"{name}_q_fraction"] = q / full_q
            qsum += q
        rec["lower_to_upper_mean_ratio"] = rec["lower_mean_raw"] / rec["upper_mean_raw"]
        rec["middle_to_upper_mean_ratio"] = rec["middle_mean_raw"] / rec["upper_mean_raw"]
        rec["band_q_sum_relative_error"] = relerr(qsum, full_q)
        rows.append(rec)

    fr = pd.DataFrame(rows).sort_values(["diameter_um", "flow_mm_s", "frame_index_0based"])
    fr.to_csv(out / "source_axial_band_framewise.csv", index=False, float_format="%.17g")

    summary_rows = []
    subset_defs = {
        "all_valid": np.ones(len(fr), dtype=bool),
        "direct_only": fr.direct_only.to_numpy(bool),
        "central_geometry": fr.central_geometry.to_numpy(bool),
        "direct_plus_central": fr.direct_plus_central.to_numpy(bool),
    }
    for subset, mask in subset_defs.items():
        sg = fr.loc[mask]
        for (scan, diameter, flow), g in sg.groupby(["scan_id", "diameter_um", "flow_mm_s"], sort=True):
            row = dict(subset=subset, scan_id=scan, diameter_um=diameter, flow_mm_s=flow, n_frames=len(g))
            for band in ("upper", "middle", "lower"):
                q1, med, q3 = qtls(g[f"{band}_mean_raw"])
                row.update({f"{band}_mean_raw_q1": q1, f"{band}_mean_raw_median": med, f"{band}_mean_raw_q3": q3,
                            f"{band}_q_fraction_median": float(g[f"{band}_q_fraction"].median())})
            row["source_mean_raw_median"] = float(g.source_mean_raw_frozen.median())
            row["lower_to_upper_mean_ratio_median"] = float(g.lower_to_upper_mean_ratio.median())
            summary_rows.append(row)
    sm = pd.DataFrame(summary_rows)
    sm.to_csv(out / "source_axial_band_volume_summary.csv", index=False, float_format="%.17g")

    # D500-vs-D285 same-flow descriptive band contrasts.
    contrasts = []
    common = [1.0, 3.0, 5.0, 7.0, 10.0]
    for subset in sm.subset.unique():
        s = sm[sm.subset.eq(subset)]
        for flow in common:
            a = s[(s.diameter_um.eq(285)) & (s.flow_mm_s.eq(flow))]
            b = s[(s.diameter_um.eq(500)) & (s.flow_mm_s.eq(flow))]
            if len(a) != 1 or len(b) != 1:
                continue
            ar, br = a.iloc[0], b.iloc[0]
            for band in ("upper", "middle", "lower"):
                va, vb = ar[f"{band}_mean_raw_median"], br[f"{band}_mean_raw_median"]
                contrasts.append(dict(subset=subset, flow_mm_s=flow, band=band,
                    d285_mean_raw=va, d500_mean_raw=vb,
                    absolute_difference=vb-va, percent_difference=100*(vb-va)/va))
    pd.DataFrame(contrasts).to_csv(out / "d500_vs_d285_axial_band_contrasts.csv", index=False, float_format="%.17g")

    # D500 spatial robustness by slow-axis fifth.
    seg_rows = []
    d500 = fr[fr.diameter_um.eq(500)]
    for (scan, flow, seg), g in d500.groupby(["scan_id", "flow_mm_s", "slow_axis_segment"], sort=True):
        row = dict(scan_id=scan, flow_mm_s=flow, slow_axis_segment=int(seg), n_frames=len(g),
                   source_mean_raw_median=float(g.source_mean_raw_frozen.median()))
        for band in ("upper", "middle", "lower"):
            row[f"{band}_mean_raw_median"] = float(g[f"{band}_mean_raw"].median())
        row["lower_to_upper_mean_ratio_median"] = float(g.lower_to_upper_mean_ratio.median())
        seg_rows.append(row)
    pd.DataFrame(seg_rows).to_csv(out / "d500_axial_band_slow_axis_summary.csv", index=False, float_format="%.17g")

    validation = {
        "status": "passed",
        "mat_root": str(mat_root),
        "frames_analyzed": int(len(fr)),
        "mat_sha256_verified": int(sha_verified),
        "diameters": sorted(diameters),
        "supersample": SUPERSAMPLE,
        "band_definition": "relative physical axial thirds of the unchanged frozen source ellipse",
        "max_frozen_source_replay_relative_error": float(max(reconstruction_errors, default=0.0)),
        "max_band_q_sum_relative_error": float(fr.band_q_sum_relative_error.max()),
        "tolerance": TOL,
        "source_roi_redefined": False,
        "background_subtraction": False,
        "normalization": False,
        "inferential_statistics": False,
    }
    (out / "validation.json").write_text(json.dumps(validation, indent=2) + "\n", encoding="utf-8")

    provenance = {
        "frozen_framewise_source": str(FORMAL.relative_to(ROOT)),
        "analysis_script": str(Path(__file__).relative_to(ROOT)),
        "signal": "sv_raw = var(abs(E),1,3)",
        "geometry": "frozen X4/X1/z_top + physical diameter; unchanged ellipse",
        "comparison": "D500 primary; D285 same-flow comparator when available",
    }
    (out / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")

    print(json.dumps(validation, indent=2))


if __name__ == "__main__":
    main()
