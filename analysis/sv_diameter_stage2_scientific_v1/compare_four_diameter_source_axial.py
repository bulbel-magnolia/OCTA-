#!/usr/bin/env python3
"""Compare axial raw-SV structure of frozen source ellipses across four diameters."""
from __future__ import annotations

import argparse, hashlib, importlib.util, io, json, re, sys, tempfile, urllib.request, zipfile
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[2]
OUT_DEFAULT = ROOT / "analysis/sv_diameter_stage2_scientific_v1/four_diameter_source_axial_comparison"
NEW_FORMAL = ROOT / "analysis/formal_sv_diameter_v1/framewise_primary.csv"
D128_NOBG = ROOT / "analysis/formal_sv_d128_v21_run001/no_background_relative_tail_v1_full2422/framewise_primary.csv"
D128_FORMAL = ROOT / "results/formal_sv_d128_v21_full2500_run001"
RELEASE_TAG = "formal-sv-d128-v21-run001"
RELEASE_BASE = f"https://github.com/bulbel-magnolia/OCTA-/releases/download/{RELEASE_TAG}"
COMMON_FLOWS = [1.0, 3.0, 5.0, 7.0, 10.0]
DIAMETERS = [128.0, 235.0, 285.0, 500.0]
ADJACENT = [(128.0,235.0),(235.0,285.0),(285.0,500.0)]
TOL = 1e-10
D128_MAP = {"flow01":"D128_F01_V01","flow03":"D128_F03_V01","flow05":"D128_F05_V01","flow07":"D128_F07_V01","flow10":"D128_F10_V01"}

AUDIT_PATH = ROOT / "analysis/sv_diameter_stage2_scientific_v1/d500_source_axial_local_array_audit.py"
spec = importlib.util.spec_from_file_location("source_axial_audit", AUDIT_PATH)
audit = importlib.util.module_from_spec(spec); assert spec.loader is not None; spec.loader.exec_module(audit)
sys.path.insert(0, str(ROOT / "src"))
from svrecttail.geometry import VesselGeometry, ellipse_weights


def sha_file(p: Path) -> str:
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(4*1024*1024),b""): h.update(b)
    return h.hexdigest()

def sha_bytes(b: bytes) -> str: return hashlib.sha256(b).hexdigest()
def relerr(a,b): return abs(float(a)-float(b))/max(abs(float(b)),np.finfo(float).tiny)
def pct(new,old): return 100.0*(new-old)/old if old != 0 else np.nan

def download(url: str, path: Path):
    req=urllib.request.Request(url,headers={"User-Agent":"OCTA-four-diameter-axial/1.0"})
    with urllib.request.urlopen(req,timeout=300) as r, path.open("wb") as f:
        while True:
            b=r.read(4*1024*1024)
            if not b: break
            f.write(b)

def parse_pkg(name: str):
    m=re.search(r"_(flow\d{2})_(\d{3})_(\d{3})\.zip$",name)
    if not m: raise ValueError(name)
    return m.group(1),int(m.group(2)),int(m.group(3))


def band_record(sv: np.ndarray, geom: VesselGeometry) -> dict:
    pixel_area=audit.DX_UM*audit.DZ_UM
    full=ellipse_weights(sv.shape,geom,supersample=audit.SUPERSAMPLE)
    bands={n:audit.band_weights(sv.shape,geom,lo,hi) for n,lo,hi in audit.BANDS}
    if np.max(np.abs(bands["upper"]+bands["middle"]+bands["lower"]-full)) > 1e-15:
        raise RuntimeError("band weights do not reconstruct frozen ellipse")
    area,q,mean=audit.weighted_stats(sv,full,pixel_area)
    rec={"source_area_um2_replay":area,"source_q_raw_replay":q,"source_mean_raw_replay":mean}
    qsum=0.0
    for n,_,_ in audit.BANDS:
        a,qq,m=audit.weighted_stats(sv,bands[n],pixel_area)
        rec[f"{n}_area_um2"]=a; rec[f"{n}_q_raw"]=qq; rec[f"{n}_mean_raw"]=m; rec[f"{n}_q_fraction"]=qq/q; qsum+=qq
    rec["middle_to_upper_mean_ratio"]=rec["middle_mean_raw"]/rec["upper_mean_raw"]
    rec["lower_to_upper_mean_ratio"]=rec["lower_mean_raw"]/rec["upper_mean_raw"]
    rec["lower_to_middle_mean_ratio"]=rec["lower_mean_raw"]/rec["middle_mean_raw"]
    rec["strict_upper_middle_lower"]=bool(rec["upper_mean_raw"]>rec["middle_mean_raw"]>rec["lower_mean_raw"])
    rec["band_q_sum_relative_error"]=audit.require_band_reconstruction(qsum,q)
    return rec


def replay_new(mat_root: Path):
    fw=pd.read_csv(NEW_FORMAL); fw=fw[fw.diameter_um.isin([235.0,285.0,500.0])].copy()
    rows=[]; max_area=max_q=max_mean=max_band=0.0
    for r in fw.itertuples(index=False):
        p=audit.find_mat(mat_root,str(r.scan_id),int(r.frame_index_0based))
        if not p.is_file() or p.stat().st_size==0: raise RuntimeError(f"missing/empty MAT {r.scan_id}/{p.name}")
        digest=sha_file(p)
        if digest != str(r.input_mat_sha256): raise RuntimeError(f"MAT SHA mismatch {r.scan_id}/{p.name}")
        sv=np.asarray(loadmat(p,variable_names=["sv_raw"],simplify_cells=True)["sv_raw"],dtype=np.float64)
        if sv.shape!=(351,500) or not np.isfinite(sv).all(): raise RuntimeError(f"invalid sv_raw {r.scan_id}/{p.name}")
        br=band_record(sv,audit.geometry_from_row(r))
        ea=relerr(br["source_area_um2_replay"],r.source_area_um2); eq=relerr(br["source_q_raw_replay"],r.source_q_raw); em=relerr(br["source_mean_raw_replay"],r.source_mean_raw)
        max_area=max(max_area,ea); max_q=max(max_q,eq); max_mean=max(max_mean,em); max_band=max(max_band,br["band_q_sum_relative_error"])
        if max(ea,eq,em)>=TOL: raise RuntimeError(f"new-diameter replay failed {r.scan_id}/{r.frame_index_0based}")
        rows.append({"scan_id":str(r.scan_id),"diameter_um":float(r.diameter_um),"flow_mm_s":float(r.flow_mm_s),"frame_index_0based":int(r.frame_index_0based),"slow_axis_segment":int(r.frame_index_0based)//100,"source_mean_raw_frozen":float(r.source_mean_raw),"z_top_edge_px":float(r.z_top_edge_px),"localization_source":str(r.localization_source),"input_identity":"retained_mat_sha256","input_sha256":digest,**br})
    return pd.DataFrame(rows),{"frames":len(rows),"sha_verified":len(rows),"max_area_replay_relerr":max_area,"max_q_replay_relerr":max_q,"max_mean_replay_relerr":max_mean,"max_band_q_relerr":max_band}


def replay_d128(package_dir: Path|None):
    base=pd.read_csv(D128_NOBG); loc=pd.read_csv(D128_FORMAL/"localization.csv")
    base["frame_index_0based"]=pd.to_numeric(base.frame_index_0based).astype(int); loc["frame_index_0based"]=pd.to_numeric(loc.frame_index_0based).astype(int)
    merged=base.merge(loc,on=["scan_id","frame_index_0based"],how="left",validate="one_to_one",suffixes=("","_loc"))
    if len(merged)!=2422 or merged[["x_left_edge_px","x_right_edge_px","z_top_edge_px"]].isna().any().any(): raise RuntimeError("D128 mapping failed")
    arrays=pd.read_csv(D128_FORMAL/"arrays_sha256.csv"); arrays["frame_index_0based"]=pd.to_numeric(arrays.frame_index_0based).astype(int)
    arr_sha={(str(r.scan_id),int(r.frame_index_0based)):str(r.sha256) for r in arrays.itertuples(index=False)}
    packages=pd.read_csv(D128_FORMAL/"download_packages.csv")
    rows=[]; pkg_ok=npz_ok=0; max_source=max_band=0.0
    with tempfile.TemporaryDirectory(prefix="sv-d128-axial-") as td0:
        td=Path(td0)
        for pkg in packages.itertuples(index=False):
            name=str(pkg.file); scan,lo,hi=parse_pkg(name); cached=(package_dir/name) if package_dir else None
            path=cached if cached is not None and cached.is_file() else td/name
            if path.parent==td: download(f"{RELEASE_BASE}/{name}",path)
            if sha_file(path)!=str(pkg.sha256): raise RuntimeError(f"D128 package SHA mismatch {name}")
            pkg_ok+=1; wanted=merged[(merged.scan_id.eq(scan)) & merged.frame_index_0based.between(lo,hi)]
            with zipfile.ZipFile(path) as zf:
                for r in wanted.itertuples(index=False):
                    fi=int(r.frame_index_0based); member=f"arrays/{scan}/frame_{fi:03d}.npz"; data=zf.read(member)
                    if sha_bytes(data)!=arr_sha[(scan,fi)]: raise RuntimeError(f"D128 NPZ SHA mismatch {scan}/{fi}")
                    npz_ok+=1
                    with np.load(io.BytesIO(data),allow_pickle=False) as d: sv=np.asarray(d["sv_raw"],dtype=np.float64)
                    geom=VesselGeometry(float(r.x_left_edge_px),float(r.x_right_edge_px),float(r.z_top_edge_px),128.0,12.7,6.7); br=band_record(sv,geom)
                    e=relerr(br["source_mean_raw_replay"],r.source_mean_raw); max_source=max(max_source,e); max_band=max(max_band,br["band_q_sum_relative_error"])
                    if e>=TOL: raise RuntimeError(f"D128 source replay failed {scan}/{fi}")
                    source="direct_candidate" if bool(getattr(r,"z_candidate_accepted")) else ("short_gap_fill" if bool(getattr(r,"z_short_gap_filled")) else "other_valid")
                    rows.append({"scan_id":D128_MAP[scan],"diameter_um":128.0,"flow_mm_s":float(r.flow_mm_s),"frame_index_0based":fi,"slow_axis_segment":fi//100,"source_mean_raw_frozen":float(r.source_mean_raw),"z_top_edge_px":float(r.z_top_edge_px),"localization_source":source,"input_identity":"d128_release_npz_sha256","input_sha256":arr_sha[(scan,fi)],**br})
            if path.parent==td: path.unlink()
    if len(rows)!=2422 or npz_ok!=2422: raise RuntimeError(f"D128 coverage mismatch rows={len(rows)} npz={npz_ok}")
    return pd.DataFrame(rows),{"frames":len(rows),"release_packages_sha_verified":pkg_ok,"npz_sha_verified":npz_ok,"max_source_mean_replay_relerr":max_source,"max_band_q_relerr":max_band}


def summarize(fr: pd.DataFrame):
    vr=[]
    for (scan,d,f),g in fr.groupby(["scan_id","diameter_um","flow_mm_s"],sort=True):
        r={"scan_id":scan,"diameter_um":d,"flow_mm_s":f,"n_frames":len(g),"strict_order_frame_fraction":float(g.strict_upper_middle_lower.mean())}
        for b in ["upper","middle","lower"]:
            r[f"{b}_mean_raw_median"]=float(g[f"{b}_mean_raw"].median()); r[f"{b}_q_fraction_median"]=float(g[f"{b}_q_fraction"].median())
        for m in ["middle_to_upper_mean_ratio","lower_to_upper_mean_ratio","lower_to_middle_mean_ratio","source_mean_raw_frozen"]: r[f"{m}_median"]=float(g[m].median())
        r["strict_upper_middle_lower_volume_medians"]=bool(r["upper_mean_raw_median"]>r["middle_mean_raw_median"]>r["lower_mean_raw_median"]); vr.append(r)
    vol=pd.DataFrame(vr).sort_values(["diameter_um","flow_mm_s"]).reset_index(drop=True)
    dr=[]
    for d,g in vol[vol.flow_mm_s.isin(COMMON_FLOWS)].groupby("diameter_um",sort=True):
        r={"diameter_um":d,"n_common_flow_volumes":len(g),"strict_order_volumes":int(g.strict_upper_middle_lower_volume_medians.sum())}
        for m in ["upper_mean_raw_median","middle_mean_raw_median","lower_mean_raw_median","middle_to_upper_mean_ratio_median","lower_to_upper_mean_ratio_median","lower_to_middle_mean_ratio_median","upper_q_fraction_median","middle_q_fraction_median","lower_q_fraction_median","source_mean_raw_frozen_median"]:
            r[f"{m}_across_flow_median"]=float(g[m].median()); r[f"{m}_across_flow_min"]=float(g[m].min()); r[f"{m}_across_flow_max"]=float(g[m].max())
        dr.append(r)
    ds=pd.DataFrame(dr)
    cr=[]
    for f in COMMON_FLOWS:
        s=vol[vol.flow_mm_s.eq(f)].set_index("diameter_um")
        for a,b in ADJACENT:
            for m in ["upper_mean_raw_median","middle_mean_raw_median","lower_mean_raw_median","middle_to_upper_mean_ratio_median","lower_to_upper_mean_ratio_median"]:
                old=float(s.loc[a,m]); new=float(s.loc[b,m]); cr.append({"flow_mm_s":f,"from_diameter_um":a,"to_diameter_um":b,"metric":m,"from_value":old,"to_value":new,"absolute_difference":new-old,"percent_difference":pct(new,old)})
    sr=[]
    for (scan,d,f,seg),g in fr.groupby(["scan_id","diameter_um","flow_mm_s","slow_axis_segment"],sort=True):
        u=float(g.upper_mean_raw.median()); m=float(g.middle_mean_raw.median()); l=float(g.lower_mean_raw.median()); sr.append({"scan_id":scan,"diameter_um":d,"flow_mm_s":f,"slow_axis_segment":int(seg),"n_frames":len(g),"upper_mean_raw_median":u,"middle_mean_raw_median":m,"lower_mean_raw_median":l,"middle_to_upper_ratio_median":float(g.middle_to_upper_mean_ratio.median()),"lower_to_upper_ratio_median":float(g.lower_to_upper_mean_ratio.median()),"strict_upper_middle_lower":bool(u>m>l)})
    return vol,ds,pd.DataFrame(cr),pd.DataFrame(sr)


def write_readme(out: Path, ds: pd.DataFrame, vol: pd.DataFrame, seg: pd.DataFrame, validation: dict):
    L=["# 四直径 frozen-source 轴向结构比较","","只比较已有 D128 / D235 / D285 / D500 数据；source ellipse 与三等分规则保持冻结。","","| Diameter | upper | middle | lower | middle/upper | lower/upper | upper Q | middle Q | lower Q | order |","|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for r in ds.itertuples(index=False): L.append(f"| {r.diameter_um:.0f} | {r.upper_mean_raw_median_across_flow_median:.6g} | {r.middle_mean_raw_median_across_flow_median:.6g} | {r.lower_mean_raw_median_across_flow_median:.6g} | {r.middle_to_upper_mean_ratio_median_across_flow_median:.4f} | {r.lower_to_upper_mean_ratio_median_across_flow_median:.4f} | {r.upper_q_fraction_median_across_flow_median:.4f} | {r.middle_q_fraction_median_across_flow_median:.4f} | {r.lower_q_fraction_median_across_flow_median:.4f} | {int(r.strict_order_volumes)}/{int(r.n_common_flow_volumes)} |")
    cvol=vol[vol.flow_mm_s.isin(COMMON_FLOWS)]; cseg=seg[seg.flow_mm_s.isin(COMMON_FLOWS)]
    L += ["","## Spatial robustness",f"- common-grid volume order upper>middle>lower: {int(cvol.strict_upper_middle_lower_volume_medians.sum())}/{len(cvol)}",f"- common-grid slow-axis segment order: {int(cseg.strict_upper_middle_lower.sum())}/{len(cseg)}","","## Validation",f"- status: {validation['status']}",f"- band-Q max relative error: {validation['max_band_q_reconstruction_relerr']:.3e}","- no ROI redefinition, normalization, background subtraction, regression or p-values."]
    (out/"README.md").write_text("\n".join(L)+"\n",encoding="utf-8")


def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--mat-root",required=True); ap.add_argument("--d128-package-dir",default=""); ap.add_argument("--output-dir",default=str(OUT_DEFAULT)); args=ap.parse_args()
    out=Path(args.output_dir).resolve(); out.mkdir(parents=True,exist_ok=True); mat_root=Path(args.mat_root).expanduser().resolve(); package_dir=Path(args.d128_package_dir).expanduser().resolve() if args.d128_package_dir else None
    new_fr,new_val=replay_new(mat_root); d128_fr,d128_val=replay_d128(package_dir); fr=pd.concat([d128_fr,new_fr],ignore_index=True).sort_values(["diameter_um","flow_mm_s","frame_index_0based"])
    vol,ds,con,seg=summarize(fr); common_vol=vol[vol.flow_mm_s.isin(COMMON_FLOWS)]
    if len(common_vol)!=20 or set(common_vol.diameter_um)!=set(DIAMETERS): raise RuntimeError(f"common-grid coverage mismatch {len(common_vol)}")
    validation={"status":"passed","common_grid_volumes":20,"total_frames":int(len(fr)),"frames_by_diameter":{str(int(k)):int(v) for k,v in fr.groupby("diameter_um").size().items()},"new_diameter":new_val,"d128":d128_val,"max_band_q_reconstruction_relerr":float(fr.band_q_sum_relative_error.max()),"source_roi_redefined":False,"background_subtraction":False,"normalization":False,"inferential_statistics":False,"p_values":False}
    if validation["max_band_q_reconstruction_relerr"]>=TOL: raise RuntimeError("band-Q validation failed")
    fr.to_csv(out/"four_diameter_source_axial_framewise.csv",index=False,float_format="%.17g",lineterminator="\n"); vol.to_csv(out/"four_diameter_source_axial_volume_summary.csv",index=False,float_format="%.17g",lineterminator="\n"); ds.to_csv(out/"four_diameter_source_axial_diameter_summary.csv",index=False,float_format="%.17g",lineterminator="\n"); con.to_csv(out/"four_diameter_source_axial_adjacent_contrasts.csv",index=False,float_format="%.17g",lineterminator="\n"); seg.to_csv(out/"four_diameter_source_axial_segment_summary.csv",index=False,float_format="%.17g",lineterminator="\n")
    (out/"validation.json").write_text(json.dumps(validation,indent=2)+"\n",encoding="utf-8"); (out/"provenance.json").write_text(json.dumps({"analysis_script":str(Path(__file__).relative_to(ROOT)),"d128_source":RELEASE_TAG,"new_diameter_source":"retained formal MAT exports SHA-matched to framewise_primary.csv","band_implementation":str(AUDIT_PATH.relative_to(ROOT)),"diameters":DIAMETERS,"common_flows":COMMON_FLOWS},indent=2)+"\n",encoding="utf-8"); write_readme(out,ds,vol,seg,validation); print(json.dumps(validation,indent=2))

if __name__=="__main__": main()
