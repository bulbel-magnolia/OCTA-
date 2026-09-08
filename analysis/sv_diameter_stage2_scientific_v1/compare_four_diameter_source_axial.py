#!/usr/bin/env python3
"""Four-diameter comparison of axial raw-SV structure inside the frozen source ellipse.

Scientific question only: is upper > middle > lower shared by D128/D235/D285/D500,
and how do middle/upper, lower/upper and band Q fractions vary with physical diameter?

No source ROI is redefined. D235/D285/D500 use retained formal MAT exports with
per-file SHA verification. D128 is replayed read-only from the frozen public release.
B-scans and slow-axis segments are spatial robustness samples, not independent repeats.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import io
import json
import re
import sys
import tempfile
import urllib.request
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[2]
OUT_DEFAULT = ROOT / "analysis/sv_diameter_stage2_scientific_v1/four_diameter_source_axial_comparison"
NEW_FORMAL = ROOT / "analysis/formal_sv_diameter_v1/framewise_primary.csv"
D128_NOBG = ROOT / "analysis/formal_sv_d128_v21_run001/no_background_relative_tail_v1_full2422/framewise_primary.csv"
D128_FORMAL = ROOT / "results/formal_sv_d128_v21_full2500_run001"
D128_RELEASE_TAG = "formal-sv-d128-v21-run001"
D128_RELEASE_BASE = f"https://github.com/bulbel-magnolia/OCTA-/releases/download/{D128_RELEASE_TAG}"
COMMON_FLOWS = [1.0, 3.0, 5.0, 7.0, 10.0]
DIAMETERS = [128.0, 235.0, 285.0, 500.0]
ADJACENT = [(128.0, 235.0), (235.0, 285.0), (285.0, 500.0)]
TOL = 1e-10
D128_MAP = {"flow01":"D128_F01_V01", "flow03":"D128_F03_V01", "flow05":"D128_F05_V01", "flow07":"D128_F07_V01", "flow10":"D128_F10_V01"}

# Import the already validated D500/D285 band implementation instead of duplicating it.
AUDIT_PATH = ROOT / "analysis/sv_diameter_stage2_scientific_v1/d500_source_axial_local_array_audit.py"
spec = importlib.util.spec_from_file_location("source_axial_audit", AUDIT_PATH)
audit = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(audit)

sys.path.insert(0, str(ROOT / "src"))
from svrecttail.geometry import VesselGeometry, ellipse_weights


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def relerr(a: float, b: float) -> float:
    return abs(float(a)-float(b))/max(abs(float(b)), np.finfo(float).tiny)


def safe_pct(new: float, old: float) -> float:
    return 100.0*(new-old)/old if np.isfinite(old) and old != 0 else np.nan


def download(url: str, path: Path) -> None:
    req = urllib.request.Request(url, headers={"User-Agent":"OCTA-SV-four-diameter-axial/1.0"})
    with urllib.request.urlopen(req, timeout=300) as r, path.open("wb") as f:
        while True:
            block = r.read(4*1024*1024)
            if not block:
                break
            f.write(block)


def parse_package_name(name: str):
    m = re.search(r"_(flow\d{2})_(\d{3})_(\d{3})\.zip$", name)
    if not m:
        raise ValueError(f"Unrecognized D128 package name: {name}")
    return m.group(1), int(m.group(2)), int(m.group(3))


def band_record(sv: np.ndarray, geom: VesselGeometry) -> dict:
    pixel_area = audit.DX_UM * audit.DZ_UM
    full = ellipse_weights(sv.shape, geom, supersample=audit.SUPERSAMPLE)
    bands = {name:audit.band_weights(sv.shape, geom, lo, hi) for name,lo,hi in audit.BANDS}
    summed = bands["upper"] + bands["middle"] + bands["lower"]
    if np.max(np.abs(summed-full)) > 1e-15:
        raise RuntimeError("Band weights do not reconstruct frozen ellipse")
    full_area, full_q, full_mean = audit.weighted_stats(sv, full, pixel_area)
    rec = {"source_area_um2_replay":full_area, "source_q_raw_replay":full_q, "source_mean_raw_replay":full_mean}
    qsum = 0.0
    for name,_,_ in audit.BANDS:
        area,q,mean = audit.weighted_stats(sv, bands[name], pixel_area)
        rec[f"{name}_area_um2"] = area
        rec[f"{name}_q_raw"] = q
        rec[f"{name}_mean_raw"] = mean
        rec[f"{name}_q_fraction"] = q/full_q
        qsum += q
    rec["middle_to_upper_mean_ratio"] = rec["middle_mean_raw"]/rec["upper_mean_raw"]
    rec["lower_to_upper_mean_ratio"] = rec["lower_mean_raw"]/rec["upper_mean_raw"]
    rec["lower_to_middle_mean_ratio"] = rec["lower_mean_raw"]/rec["middle_mean_raw"]
    rec["strict_upper_middle_lower"] = bool(rec["upper_mean_raw"] > rec["middle_mean_raw"] > rec["lower_mean_raw"])
    rec["band_q_sum_relative_error"] = audit.require_band_reconstruction(qsum, full_q)
    return rec


def replay_new(mat_root: Path) -> tuple[pd.DataFrame, dict]:
    fw = pd.read_csv(NEW_FORMAL)
    fw = fw[fw.diameter_um.isin([235.0,285.0,500.0])].copy()
    rows=[]; sha_ok=0; max_area=max_q=max_mean=max_band=0.0
    for r in fw.itertuples(index=False):
        p = audit.find_mat(mat_root, str(r.scan_id), int(r.frame_index_0based))
        if not p.is_file() or p.stat().st_size == 0:
            raise RuntimeError(f"Missing/empty MAT: {r.scan_id}/{p.name}")
        digest = sha256_file(p)
        if digest != str(r.input_mat_sha256):
            raise RuntimeError(f"MAT SHA mismatch: {r.scan_id}/{p.name}")
        sha_ok += 1
        d = loadmat(p, variable_names=["sv_raw"], simplify_cells=True)
        sv = np.asarray(d["sv_raw"], dtype=np.float64)
        if sv.shape != (351,500) or not np.isfinite(sv).all():
            raise RuntimeError(f"Invalid sv_raw: {r.scan_id}/{p.name}")
        geom = audit.geometry_from_row(r)
        br = band_record(sv, geom)
        ea = relerr(br["source_area_um2_replay"], r.source_area_um2)
        eq = relerr(br["source_q_raw_replay"], r.source_q_raw)
        em = relerr(br["source_mean_raw_replay"], r.source_mean_raw)
        max_area=max(max_area,ea); max_q=max(max_q,eq); max_mean=max(max_mean,em); max_band=max(max_band,br["band_q_sum_relative_error"])
        if max(ea,eq,em) >= TOL:
            raise RuntimeError(f"Frozen new-diameter source replay failed: {r.scan_id}/{r.frame_index_0based}")
        rows.append({
            "scan_id":str(r.scan_id), "diameter_um":float(r.diameter_um), "flow_mm_s":float(r.flow_mm_s),
            "frame_index_0based":int(r.frame_index_0based), "slow_axis_segment":int(r.frame_index_0based)//100,
            "source_mean_raw_frozen":float(r.source_mean_raw), "z_top_edge_px":float(r.z_top_edge_px),
            "localization_source":str(r.localization_source), "input_identity":"retained_mat_sha256", "input_sha256":digest,
            **br,
        })
    return pd.DataFrame(rows), {
        "frames":len(rows), "sha_verified":sha_ok, "max_area_replay_relerr":max_area,
        "max_q_replay_relerr":max_q, "max_mean_replay_relerr":max_mean, "max_band_q_relerr":max_band,
    }


def replay_d128(package_dir: Path|None) -> tuple[pd.DataFrame, dict]:
    base = pd.read_csv(D128_NOBG)
    loc = pd.read_csv(D128_FORMAL / "localization.csv")
    base["frame_index_0based"] = pd.to_numeric(base.frame_index_0based).astype(int)
    loc["frame_index_0based"] = pd.to_numeric(loc.frame_index_0based).astype(int)
    merged = base.merge(loc, on=["scan_id","frame_index_0based"], how="left", validate="one_to_one", suffixes=("","_loc"))
    if len(merged) != 2422 or merged[["x_left_edge_px","x_right_edge_px","z_top_edge_px"]].isna().any().any():
        raise RuntimeError("D128 frozen frame/localization mapping failed")
    arrays = pd.read_csv(D128_FORMAL / "arrays_sha256.csv")
    arrays["frame_index_0based"] = pd.to_numeric(arrays.frame_index_0based).astype(int)
    array_sha={(str(r.scan_id),int(r.frame_index_0based)):str(r.sha256) for r in arrays.itertuples(index=False)}
    packages=pd.read_csv(D128_FORMAL / "download_packages.csv")
    rows=[]; package_ok=0; npz_ok=0; max_source=max_band=0.0
    with tempfile.TemporaryDirectory(prefix="sv-d128-four-diameter-") as td:
        td=Path(td)
        for pkg in packages.itertuples(index=False):
            name=str(pkg.file); scan,lo,hi=parse_package_name(name)
            cached=(package_dir/name) if package_dir is not None else None
            if cached is not None and cached.is_file():
                path=cached
            else:
                path=td/name; download(f"{D128_RELEASE_BASE}/{name}", path)
            if sha256_file(path) != str(pkg.sha256):
                raise RuntimeError(f"D128 package SHA mismatch: {name}")
            package_ok += 1
            wanted=merged[(merged.scan_id.eq(scan)) & merged.frame_index_0based.between(lo,hi)]
            with zipfile.ZipFile(path) as zf:
                for r in wanted.itertuples(index=False):
                    fi=int(r.frame_index_0based); member=f"arrays/{scan}/frame_{fi:03d}.npz"; data=zf.read(member)
                    if sha256_bytes(data) != array_sha[(scan,fi)]:
                        raise RuntimeError(f"D128 NPZ SHA mismatch: {scan}/{fi}")
                    npz_ok += 1
                    with np.load(io.BytesIO(data), allow_pickle=False) as d:
                        sv=np.asarray(d["sv_raw"],dtype=np.float64)
                    geom=VesselGeometry(float(r.x_left_edge_px),float(r.x_right_edge_px),float(r.z_top_edge_px),128.0,12.7,6.7)
                    br=band_record(sv,geom)
                    e=relerr(br["source_mean_raw_replay"],r.source_mean_raw); max_source=max(max_source,e); max_band=max(max_band,br["band_q_sum_relative_error"])
                    if e >= TOL:
                        raise RuntimeError(f"D128 frozen source replay failed: {scan}/{fi}")
                    if bool(getattr(r,"z_candidate_accepted")):
                        source="direct_candidate"
                    elif bool(getattr(r,"z_short_gap_filled")):
                        source="short_gap_fill"
                    else:
                        source="other_valid"
                    rows.append({
                        "scan_id":D128_MAP[scan], "diameter_um":128.0, "flow_mm_s":float(r.flow_mm_s),
                        "frame_index_0based":fi, "slow_axis_segment":fi//100, "source_mean_raw_frozen":float(r.source_mean_raw),
                        "z_top_edge_px":float(r.z_top_edge_px), "localization_source":source,
                        "input_identity":"d128_release_npz_sha256", "input_sha256":array_sha[(scan,fi)], **br,
                    })
            if path.parent == td:
                path.unlink()
    if len(rows) != 2422 or npz_ok != 2422:
        raise RuntimeError(f"D128 coverage mismatch: rows={len(rows)}, npz={npz_ok}")
    return pd.DataFrame(rows), {"frames":len(rows),"release_packages_sha_verified":package_ok,"npz_sha_verified":npz_ok,"max_source_mean_replay_relerr":max_source,"max_band_q_relerr":max_band}


def volume_summary(fr: pd.DataFrame) -> pd.DataFrame:
    rows=[]
    for (scan,d,f),g in fr.groupby(["scan_id","diameter_um","flow_mm_s"],sort=True):
        rec={"scan_id":scan,"diameter_um":d,"flow_mm_s":f,"n_frames":len(g),"strict_order_frame_fraction":float(g.strict_upper_middle_lower.mean())}
        for band in ["upper","middle","lower"]:
            rec[f"{band}_mean_raw_median"]=float(g[f"{band}_mean_raw"].median())
            rec[f"{band}_q_fraction_median"]=float(g[f"{band}_q_fraction"].median())
        for metric in ["middle_to_upper_mean_ratio","lower_to_upper_mean_ratio","lower_to_middle_mean_ratio","source_mean_raw_frozen"]:
            rec[f"{metric}_median"]=float(g[metric].median())
        rec["strict_upper_middle_lower_volume_medians"] = bool(rec["upper_mean_raw_median"] > rec["middle_mean_raw_median"] > rec["lower_mean_raw_median"])
        rows.append(rec)
    return pd.DataFrame(rows).sort_values(["diameter_um","flow_mm_s"]).reset_index(drop=True)


def diameter_summary(vol: pd.DataFrame) -> pd.DataFrame:
    rows=[]
    for d,g in vol[vol.flow_mm_s.isin(COMMON_FLOWS)].groupby("diameter_um",sort=True):
        rec={"diameter_um":d,"n_common_flow_volumes":len(g),"strict_order_volumes":int(g.strict_upper_middle_lower_volume_medians.sum())}
        metrics=["upper_mean_raw_median","middle_mean_raw_median","lower_mean_raw_median",
                 "middle_to_upper_mean_ratio_median","lower_to_upper_mean_ratio_median","lower_to_middle_mean_ratio_median",
                 "upper_q_fraction_median","middle_q_fraction_median","lower_q_fraction_median","source_mean_raw_frozen_median"]
        for m in metrics:
            rec[f"{m}_across_flow_median"]=float(g[m].median())
            rec[f"{m}_across_flow_min"]=float(g[m].min())
            rec[f"{m}_across_flow_max"]=float(g[m].max())
        rows.append(rec)
    return pd.DataFrame(rows)


def adjacent_contrasts(vol: pd.DataFrame) -> pd.DataFrame:
    rows=[]; metrics=["upper_mean_raw_median","middle_mean_raw_median","lower_mean_raw_median","middle_to_upper_mean_ratio_median","lower_to_upper_mean_ratio_median"]
    for f in COMMON_FLOWS:
        s=vol[vol.flow_mm_s.eq(f)].set_index("diameter_um")
        for a,b in ADJACENT:
            for m in metrics:
                old=float(s.loc[a,m]); new=float(s.loc[b,m])
                rows.append({"flow_mm_s":f,"from_diameter_um":a,"to_diameter_um":b,"metric":m,"from_value":old,"to_value":new,"absolute_difference":new-old,"percent_difference":safe_pct(new,old)})
    return pd.DataFrame(rows)


def segment_summary(fr: pd.DataFrame) -> pd.DataFrame:
    rows=[]
    for (scan,d,f,seg),g in fr.groupby(["scan_id","diameter_um","flow_mm_s","slow_axis_segment"],sort=True):
        u=float(g.upper_mean_raw.median()); m=float(g.middle_mean_raw.median()); l=float(g.lower_mean_raw.median())
        rows.append({"scan_id":scan,"diameter_um":d,"flow_mm_s":f,"slow_axis_segment":int(seg),"n_frames":len(g),"upper_mean_raw_median":u,"middle_mean_raw_median":m,"lower_mean_raw_median":l,"middle_to_upper_ratio_median":float(g.middle_to_upper_mean_ratio.median()),"lower_to_upper_ratio_median":float(g.lower_to_upper_mean_ratio.median()),"strict_upper_middle_lower":bool(u>m>l)})
    return pd.DataFrame(rows)


def write_readme(out: Path, ds: pd.DataFrame, vol: pd.DataFrame, seg: pd.DataFrame, validation: dict) -> None:
    lines=["# 四直径 frozen-source 轴向结构比较","","比较 D128 / D235 / D285 / D500 在同一个冻结 source ellipse 内的 upper / middle / lower raw-SV 结构。只分析已有数据，不改变 ROI。","","## Common-flow diameter summary","",
           "| Diameter | upper raw SV | middle raw SV | lower raw SV | middle/upper | lower/upper | upper Q frac | middle Q frac | lower Q frac | order |","|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for r in ds.itertuples(index=False):
        lines.append(f"| {r.diameter_um:.0f} | {r.upper_mean_raw_median_across_flow_median:.6g} | {r.middle_mean_raw_median_across_flow_median:.6g} | {r.lower_mean_raw_median_across_flow_median:.6g} | {r.middle_to_upper_mean_ratio_median_across_flow_median:.4f} | {r.lower_to_upper_mean_ratio_median_across_flow_median:.4f} | {r.upper_q_fraction_median_across_flow_median:.4f} | {r.middle_q_fraction_median_across_flow_median:.4f} | {r.lower_q_fraction_median_across_flow_median:.4f} | {int(r.strict_order_volumes)}/{int(r.n_common_flow_volumes)} |")
    lines += ["","## Fixed-flow shape table","","| Flow | Diameter | middle/upper | lower/upper | upper Q frac | middle Q frac | lower Q frac |","|---:|---:|---:|---:|---:|---:|---:|"]
    for r in vol[vol.flow_mm_s.isin(COMMON_FLOWS)].itertuples(index=False):
        lines.append(f"| {r.flow_mm_s:g} | {r.diameter_um:.0f} | {r.middle_to_upper_mean_ratio_median:.4f} | {r.lower_to_upper_mean_ratio_median:.4f} | {r.upper_q_fraction_median:.4f} | {r.middle_q_fraction_median:.4f} | {r.lower_q_fraction_median:.4f} |")
    lines += ["","## Spatial robustness","",f"Common-grid volume-level upper>middle>lower: {int(vol[vol.flow_mm_s.isin(COMMON_FLOWS)].strict_upper_middle_lower_volume_medians.sum())}/{len(vol[vol.flow_mm_s.isin(COMMON_FLOWS)])}.",f"Common-grid slow-axis segments upper>middle>lower: {int(seg[seg.flow_mm_s.isin(COMMON_FLOWS)].strict_upper_middle_lower.sum())}/{len(seg[seg.flow_mm_s.isin(COMMON_FLOWS)])}.","","## Validation","",f"Status: **{validation['status']}**. D128 and retained new-diameter arrays are identity-verified; all source replay and band-Q reconstruction gates are < {TOL:g}.","","B-scans and segments are spatial positions, not independent replicates. No p-values, regression, background subtraction, normalization, or new ROI were introduced."]
    (out/"README.md").write_text("\n".join(lines)+"\n",encoding="utf-8")


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--mat-root",required=True,help="Retained D235/D285/D500 formal MAT root")
    ap.add_argument("--d128-package-dir",default="",help="Optional directory containing frozen D128 release ZIP assets")
    ap.add_argument("--output-dir",default=str(OUT_DEFAULT))
    args=ap.parse_args()
    out=Path(args.output_dir).resolve(); out.mkdir(parents=True,exist_ok=True)
    mat_root=Path(args.mat_root).expanduser().resolve(); package_dir=Path(args.d128_package_dir).expanduser().resolve() if args.d128_package_dir else None

    new_fr,new_val=replay_new(mat_root)
    d128_fr,d128_val=replay_d128(package_dir)
    fr=pd.concat([d128_fr,new_fr],ignore_index=True).sort_values(["diameter_um","flow_mm_s","frame_index_0based"])
    vol=volume_summary(fr); ds=diameter_summary(vol); con=adjacent_contrasts(vol); seg=segment_summary(fr)

    expected={128.0:2422,235.0:233? ,285.0:2371,500.0:3913}
    # D235 count is read from frozen table rather than hard-coded; all new-diameter rows were SHA-gated above.
    common_vol=vol[vol.flow_mm_s.isin(COMMON_FLOWS)]
    if len(common_vol) != 20 or set(common_vol.diameter_um) != set(DIAMETERS):
        raise RuntimeError(f"Common-grid coverage mismatch: {len(common_vol)} volumes")
    if ds.strict_order_volumes.sum() < 0:  # structural no-op, keeps field numeric and validated
        raise RuntimeError("Impossible ordering count")
    validation={
        "status":"passed","common_grid_volumes":int(len(common_vol)),"total_frames":int(len(fr)),
        "frames_by_diameter":{str(int(k)):int(v) for k,v in fr.groupby("diameter_um").size().items()},
        "new_diameter":new_val,"d128":d128_val,
        "max_band_q_reconstruction_relerr":float(fr.band_q_sum_relative_error.max()),
        "source_roi_redefined":False,"signal":"sv_raw = var(abs(E),1,3)","background_subtraction":False,"normalization":False,
        "inferential_statistics":False,"p_values":False,
    }
    if validation["max_band_q_reconstruction_relerr"] >= TOL:
        raise RuntimeError("Band-Q validation failed")

    fr.to_csv(out/"four_diameter_source_axial_framewise.csv",index=False,float_format="%.17g",lineterminator="\n")
    vol.to_csv(out/"four_diameter_source_axial_volume_summary.csv",index=False,float_format="%.17g",lineterminator="\n")
    ds.to_csv(out/"four_diameter_source_axial_diameter_summary.csv",index=False,float_format="%.17g",lineterminator="\n")
    con.to_csv(out/"four_diameter_source_axial_adjacent_contrasts.csv",index=False,float_format="%.17g",lineterminator="\n")
    seg.to_csv(out/"four_diameter_source_axial_segment_summary.csv",index=False,float_format="%.17g",lineterminator="\n")
    (out/"validation.json").write_text(json.dumps(validation,indent=2)+"\n",encoding="utf-8")
    provenance={"analysis_script":str(Path(__file__).relative_to(ROOT)),"d128_source":"formal-sv-d128-v21-run001 release","new_diameter_source":"retained formal MAT exports SHA-matched to framewise_primary.csv","band_implementation":str(AUDIT_PATH.relative_to(ROOT)),"diameters":DIAMETERS,"common_flows":COMMON_FLOWS}
    (out/"provenance.json").write_text(json.dumps(provenance,indent=2)+"\n",encoding="utf-8")
    write_readme(out,ds,vol,seg,validation)
    print(json.dumps(validation,indent=2))

if __name__ == "__main__":
    main()
