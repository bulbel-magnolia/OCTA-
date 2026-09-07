#!/usr/bin/env python3
"""Three descriptive follow-up audits for no-background SV Relative Tail Intensity.

1) source-tail-ratio coupling audit
2) exact 100/200/300/500 um tail-window sensitivity
3) within-volume 5x100-frame spatial-segment stability

Primary endpoint remains:
    RI_tail,L = mean(raw SV in tail ROI from vessel bottom to L) / mean(raw SV in vessel ROI)
No background is subtracted. B-scans are spatial positions, not independent experimental replicates.
No p-values or flow-effect inference are produced.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import io
import json
import math
import tempfile
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

EXPECTED_COUNTS = {"flow01":486,"flow03":468,"flow05":491,"flow07":493,"flow10":484}
EXPECTED_VALID = sum(EXPECTED_COUNTS.values())
EXPECTED_PACKAGES = 25
WINDOWS_UM = (100.0, 200.0, 300.0, 500.0)
SEGMENTS = ((0,99),(100,199),(200,299),(300,399),(400,499))
RELEASE_TAG = "formal-sv-d128-v21-run001"
RELEASE_BASE = f"https://github.com/bulbel-magnolia/OCTA-/releases/download/{RELEASE_TAG}"


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def safe_div(a, b):
    a, b = float(a), float(b)
    return a / b if np.isfinite(a) and np.isfinite(b) and b != 0 else np.nan


def weighted_integral(img: np.ndarray, weights: np.ndarray, pixel_area: float) -> float:
    m = weights > 0
    if not np.any(m) or not np.isfinite(img[m]).all():
        return np.nan
    return float(np.sum(img[m] * weights[m], dtype=np.float64) * pixel_area)


def pct_delta(new, old):
    return 100.0 * safe_div(float(new) - float(old), old)


def qstats(values):
    a = pd.to_numeric(pd.Series(values), errors="coerce").to_numpy(float)
    a = a[np.isfinite(a)]
    if not len(a):
        return dict(n=0,q1=np.nan,median=np.nan,q3=np.nan,iqr=np.nan,mean=np.nan,sd=np.nan,min=np.nan,max=np.nan)
    q1, med, q3 = np.quantile(a, [0.25,0.5,0.75])
    return dict(n=int(len(a)),q1=float(q1),median=float(med),q3=float(q3),iqr=float(q3-q1),mean=float(a.mean()),sd=float(a.std(ddof=1)) if len(a)>1 else np.nan,min=float(a.min()),max=float(a.max()))


def corr(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 3:
        return np.nan
    xx, yy = x[m], y[m]
    if np.std(xx) == 0 or np.std(yy) == 0:
        return np.nan
    return float(np.corrcoef(xx, yy)[0,1])


def spearman(x, y):
    xs = pd.Series(np.asarray(x,float)).rank(method="average").to_numpy(float)
    ys = pd.Series(np.asarray(y,float)).rank(method="average").to_numpy(float)
    return corr(xs, ys)


def coupling_audit(primary: pd.DataFrame, out: Path):
    summary = []
    deciles = []
    quintiles = []
    for scan, g0 in primary.groupby("scan_id", sort=False):
        g = g0.sort_values("frame_index_0based").copy()
        flow = float(g.flow_mm_s.iloc[0])
        source = g.source_mean_raw.to_numpy(float)
        tail = g.tail_mean_raw.to_numpy(float)
        ri = g.ri_tail.to_numpy(float)
        if not (np.all(source>0) and np.all(tail>0) and np.all(ri>0)):
            raise AssertionError(f"nonpositive raw metric in {scan}")
        qs, qt, qr = qstats(source), qstats(tail), qstats(ri)
        ls, lt, lr = np.log(source), np.log(tail), np.log(ri)
        identity_err = float(np.max(np.abs(lr - (lt-ls))))
        summary.append({
            "scan_id":scan,"flow_mm_s":flow,"n_valid":len(g),
            "source_mean_median":qs["median"],"source_iqr_over_median":safe_div(qs["iqr"],qs["median"]),
            "tail_mean_median":qt["median"],"tail_iqr_over_median":safe_div(qt["iqr"],qt["median"]),
            "ri_tail_median":qr["median"],"ri_iqr_over_median":safe_div(qr["iqr"],qr["median"]),
            "pearson_source_vs_ri":corr(source,ri),"spearman_source_vs_ri":spearman(source,ri),
            "pearson_tail_vs_ri":corr(tail,ri),"spearman_tail_vs_ri":spearman(tail,ri),
            "pearson_source_vs_tail":corr(source,tail),"spearman_source_vs_tail":spearman(source,tail),
            "var_log_source":float(np.var(ls,ddof=1)),"var_log_tail":float(np.var(lt,ddof=1)),
            "cov_log_tail_source":float(np.cov(lt,ls,ddof=1)[0,1]),"var_log_ri":float(np.var(lr,ddof=1)),
            "log_identity_max_abs_error":identity_err,
        })
        # Source-decile audit. qcut duplicates='drop' is deterministic for these continuous data.
        bins = pd.qcut(g.source_mean_raw, q=10, labels=False, duplicates="drop")
        for b, gb in g.groupby(bins, sort=True):
            s,t,r=qstats(gb.source_mean_raw),qstats(gb.tail_mean_raw),qstats(gb.ri_tail)
            deciles.append({"scan_id":scan,"flow_mm_s":flow,"source_decile":int(b)+1,"n":len(gb),
                            "source_median":s["median"],"tail_median":t["median"],"ri_median":r["median"],"ri_q1":r["q1"],"ri_q3":r["q3"]})
        q5 = pd.qcut(g.source_mean_raw, q=5, labels=False, duplicates="drop")
        for b, gb in g.groupby(q5, sort=True):
            s,t,r=qstats(gb.source_mean_raw),qstats(gb.tail_mean_raw),qstats(gb.ri_tail)
            quintiles.append({"scan_id":scan,"flow_mm_s":flow,"source_quintile":int(b)+1,"n":len(gb),
                              "source_median":s["median"],"tail_median":t["median"],"ri_median":r["median"]})
    sm = pd.DataFrame(summary).sort_values("flow_mm_s")
    de = pd.DataFrame(deciles).sort_values(["flow_mm_s","source_decile"])
    qu = pd.DataFrame(quintiles).sort_values(["flow_mm_s","source_quintile"])
    sm.to_csv(out/"coupling_summary.csv",index=False,float_format="%.17g")
    de.to_csv(out/"coupling_source_deciles.csv",index=False,float_format="%.17g")
    qu.to_csv(out/"coupling_source_quintiles.csv",index=False,float_format="%.17g")

    # Flow07 context relative to flow01 and the median of the other four volume medians.
    med = sm.set_index("scan_id")
    other = med.drop(index="flow07")
    rows=[]
    for metric in ("source_mean_median","tail_mean_median","ri_tail_median"):
        v7=float(med.loc["flow07",metric]); v1=float(med.loc["flow01",metric]); vo=float(other[metric].median())
        rows.append({"metric":metric,"flow07":v7,"flow01":v1,"other4_scan_median":vo,
                     "flow07_vs_flow01_pct":pct_delta(v7,v1),"flow07_vs_other4_median_pct":pct_delta(v7,vo)})
    ctx=pd.DataFrame(rows)
    ctx.to_csv(out/"flow07_coupling_context.csv",index=False,float_format="%.17g")
    return sm,de,ctx


def exact_windows(root: Path, primary: pd.DataFrame, out: Path):
    formal=root/"results/formal_sv_d128_v21_full2500_run001"
    task1=root/"analysis/formal_sv_d128_v21_run001/relative_tail_intensity_full2422"
    helper_path=root/"analysis/formal_sv_d128_v21_run001/fixed128_width_sensitivity_full2422/audit_fixed128_width.py"
    H=load_module("fixed128_helper_followup",helper_path)
    import sys
    sys.path.insert(0,str(root/"src"))
    from svrecttail.geometry import VesselGeometry, ellipse_weights, interval_overlap_weights

    base,_,cfg=H.prepare_inputs(root,formal,task1)
    p = primary.rename(columns={"frame_index_0based":"frame_index"})[["scan_id","frame_index","source_mean_raw","tail_mean_raw","ri_tail"]]
    merged=base.merge(p,on=["scan_id","frame_index"],how="inner",validate="one_to_one")
    if len(merged)!=EXPECTED_VALID:
        raise AssertionError("geometry/primary key mismatch")
    packages=pd.read_csv(formal/"download_packages.csv")
    arrays=pd.read_csv(formal/"arrays_sha256.csv")
    arrays["scan_id"]=arrays.scan_id.astype(str); arrays["frame_index_0based"]=pd.to_numeric(arrays.frame_index_0based).astype(int)
    array_hash={(r.scan_id,int(r.frame_index_0based)):str(r.sha256) for r in arrays.itertuples(index=False)}
    if len(packages)!=EXPECTED_PACKAGES:
        raise AssertionError("package count mismatch")

    records=[]; pkg_audit=[]; verified=0; source_replay_max=0.; ri500_replay_max=0.; pixel=cfg["dx_um"]*cfg["dz_um"]
    with tempfile.TemporaryDirectory(prefix="sv-three-followups-") as td:
        td=Path(td)
        for pkg in packages.itertuples(index=False):
            name=str(pkg.file); scan,lo,hi=H.parse_package_name(name); path=td/name
            H.download_file(f"{RELEASE_BASE}/{name}",path)
            zsha=H.sha256_file(path)
            if zsha!=str(pkg.sha256) or path.stat().st_size!=int(pkg.bytes):
                raise AssertionError(f"package verification failed {name}")
            wanted=merged[(merged.scan_id.eq(scan)) & merged.frame_index.between(lo,hi)].copy(); vin=0
            with zipfile.ZipFile(path,"r") as zf:
                names=set(zf.namelist())
                for row in wanted.itertuples(index=False):
                    fi=int(row.frame_index); member=f"arrays/{scan}/frame_{fi:03d}.npz"
                    if member not in names: raise FileNotFoundError(member)
                    data=zf.read(member)
                    if H.sha256_bytes(data)!=array_hash[(scan,fi)]: raise AssertionError(f"NPZ SHA mismatch {scan}/{fi}")
                    verified+=1; vin+=1
                    with np.load(io.BytesIO(data),allow_pickle=False) as d:
                        sv=np.asarray(d["sv_raw"],dtype=np.float64)
                    geom=VesselGeometry(x_left_edge_px=float(row.baseline_x_left_edge_px),x_right_edge_px=float(row.baseline_x_right_edge_px),
                                        z_top_edge_px=float(row.z_top_edge_px),diameter_um=cfg["diameter_um"],dx_um=cfg["dx_um"],dz_um=cfg["dz_um"])
                    sw=ellipse_weights(sv.shape,geom,supersample=cfg["ellipse_supersample"])
                    xw=interval_overlap_weights(sv.shape[1],geom.x_left_edge_px,geom.x_right_edge_px)
                    av=float(sw.sum()*pixel); qv=weighted_integral(sv,sw,pixel); source_mean=safe_div(qv,av)
                    source_replay_max=max(source_replay_max,abs(source_mean-float(row.source_mean_raw))/max(abs(float(row.source_mean_raw)),1.0))
                    for L in WINDOWS_UM:
                        z0=geom.z_bottom_edge_px
                        zw=interval_overlap_weights(sv.shape[0],z0,z0+L/cfg["dz_um"])
                        tw=np.multiply.outer(zw,xw); at=float(tw.sum()*pixel); qt=weighted_integral(sv,tw,pixel); tail_mean=safe_div(qt,at); ri=safe_div(tail_mean,source_mean)
                        records.append({"scan_id":scan,"frame_index_0based":fi,"flow_mm_s":float(row.flow_mm_s),"window_um":L,
                                        "source_mean_raw":source_mean,"tail_mean_raw_window":tail_mean,"ri_tail_window":ri,"tail_area_um2":at})
                        if L==500.0:
                            ri500_replay_max=max(ri500_replay_max,abs(ri-float(row.ri_tail))/max(abs(float(row.ri_tail)),1e-12))
            pkg_audit.append({"package":name,"scan_id":scan,"verified":True,"valid_frames":len(wanted),"npz_verified":vin}); path.unlink()
    wf=pd.DataFrame(records).sort_values(["flow_mm_s","frame_index_0based","window_um"])
    if len(wf)!=EXPECTED_VALID*len(WINDOWS_UM) or verified!=EXPECTED_VALID:
        raise AssertionError("window output coverage mismatch")
    if not np.isfinite(wf[["source_mean_raw","tail_mean_raw_window","ri_tail_window"]].to_numpy(float)).all():
        raise AssertionError("nonfinite exact-window output")
    wf.to_csv(out/"window_framewise.csv.gz",index=False,float_format="%.17g",compression="gzip")
    pd.DataFrame(pkg_audit).to_csv(out/"release_package_audit.csv",index=False)

    summary=[]
    for (scan,L),g in wf.groupby(["scan_id","window_um"],sort=False):
        rs=qstats(g.ri_tail_window); ts=qstats(g.tail_mean_raw_window); ss=qstats(g.source_mean_raw)
        summary.append({"scan_id":scan,"flow_mm_s":float(g.flow_mm_s.iloc[0]),"window_um":float(L),"n_valid":len(g),
                        "ri_q1":rs["q1"],"ri_median":rs["median"],"ri_q3":rs["q3"],
                        "tail_mean_median":ts["median"],"source_mean_median":ss["median"]})
    ws=pd.DataFrame(summary).sort_values(["flow_mm_s","window_um"])
    refs=ws[ws.flow_mm_s.eq(1)].set_index("window_um").ri_median.to_dict()
    ws["ri_delta_pct_vs_flow01"]=[pct_delta(r.ri_median,refs[r.window_um]) for r in ws.itertuples(index=False)]
    ws.to_csv(out/"window_scan_summary.csv",index=False,float_format="%.17g")
    # One row per window for the central question: where is flow07 elevated?
    rows=[]
    for L,g in ws.groupby("window_um"):
        v7=float(g[g.flow_mm_s.eq(7)].ri_median.iloc[0]); v1=float(g[g.flow_mm_s.eq(1)].ri_median.iloc[0])
        other=float(g[~g.flow_mm_s.eq(7)].ri_median.median())
        rows.append({"window_um":float(L),"flow07_ri_median":v7,"flow01_ri_median":v1,"other4_scan_median_ri":other,
                     "flow07_vs_flow01_pct":pct_delta(v7,v1),"flow07_vs_other4_median_pct":pct_delta(v7,other),
                     "flow07_rank_desc":int(g.ri_median.rank(method="min",ascending=False)[g.flow_mm_s.eq(7)].iloc[0])})
    wc=pd.DataFrame(rows).sort_values("window_um"); wc.to_csv(out/"window_flow07_context.csv",index=False,float_format="%.17g")
    validation={"release_packages_verified":int(len(pkg_audit)),"valid_frame_npz_sha256_verified":verified,
                "source_mean_replay_max_relative_error":float(source_replay_max),"ri500_replay_max_relative_error":float(ri500_replay_max)}
    return wf,ws,wc,validation


def segment_audit(primary: pd.DataFrame, out: Path):
    segrows=[]
    full={}
    for scan,g in primary.groupby("scan_id",sort=False):
        full[scan]={m:float(g[m].median()) for m in ["source_mean_raw","tail_mean_raw","ri_tail"]}
        flow=float(g.flow_mm_s.iloc[0])
        for sid,(lo,hi) in enumerate(SEGMENTS,1):
            sg=g[g.frame_index_0based.between(lo,hi)].copy()
            if len(sg)==0: raise AssertionError(f"empty segment {scan}/{sid}")
            ss,ts,rs=qstats(sg.source_mean_raw),qstats(sg.tail_mean_raw),qstats(sg.ri_tail)
            segrows.append({"scan_id":scan,"flow_mm_s":flow,"segment_id":sid,"frame_lo":lo,"frame_hi":hi,"n_valid":len(sg),
                            "source_median":ss["median"],"tail_median":ts["median"],"ri_median":rs["median"],"ri_q1":rs["q1"],"ri_q3":rs["q3"],
                            "source_delta_pct_vs_full":pct_delta(ss["median"],full[scan]["source_mean_raw"]),
                            "tail_delta_pct_vs_full":pct_delta(ts["median"],full[scan]["tail_mean_raw"]),
                            "ri_delta_pct_vs_full":pct_delta(rs["median"],full[scan]["ri_tail"])})
    seg=pd.DataFrame(segrows).sort_values(["flow_mm_s","segment_id"]); seg.to_csv(out/"segment_summary.csv",index=False,float_format="%.17g")
    stable=[]
    for scan,g in seg.groupby("scan_id",sort=False):
        stable.append({"scan_id":scan,"flow_mm_s":float(g.flow_mm_s.iloc[0]),"segments":len(g),"valid_min":int(g.n_valid.min()),"valid_max":int(g.n_valid.max()),
                       "ri_segment_min":float(g.ri_median.min()),"ri_segment_max":float(g.ri_median.max()),
                       "ri_segment_range_pct_of_full":100.0*(float(g.ri_median.max())-float(g.ri_median.min()))/full[scan]["ri_tail"],
                       "max_abs_ri_segment_delta_pct_vs_full":float(np.max(np.abs(g.ri_delta_pct_vs_full))),
                       "max_abs_source_segment_delta_pct_vs_full":float(np.max(np.abs(g.source_delta_pct_vs_full))),
                       "max_abs_tail_segment_delta_pct_vs_full":float(np.max(np.abs(g.tail_delta_pct_vs_full)))})
    st=pd.DataFrame(stable).sort_values("flow_mm_s"); st.to_csv(out/"segment_stability_by_scan.csv",index=False,float_format="%.17g")
    ctx=[]
    for sid,g in seg.groupby("segment_id"):
        v7=g[g.flow_mm_s.eq(7)].iloc[0]; others=g[~g.flow_mm_s.eq(7)]
        for metric in ("source_median","tail_median","ri_median"):
            other_med=float(others[metric].median()); val=float(v7[metric])
            ctx.append({"segment_id":int(sid),"frame_lo":int(v7.frame_lo),"frame_hi":int(v7.frame_hi),"metric":metric,
                        "flow07_value":val,"other4_scan_median":other_med,"flow07_vs_other4_median_pct":pct_delta(val,other_med),
                        "flow07_rank_desc":int(g[metric].rank(method="min",ascending=False)[g.flow_mm_s.eq(7)].iloc[0])})
    sc=pd.DataFrame(ctx).sort_values(["segment_id","metric"]); sc.to_csv(out/"segment_flow07_context.csv",index=False,float_format="%.17g")
    return seg,st,sc


def write_readme(out: Path, coupling: pd.DataFrame, ctx: pd.DataFrame, windows: pd.DataFrame, wctx: pd.DataFrame,
                 seg: pd.DataFrame, st: pd.DataFrame, sctx: pd.DataFrame, val: dict):
    c=coupling.set_index("scan_id"); f7ctx=ctx.set_index("metric")
    lines=["# No-background SV follow-up: three descriptive audits v1","",
           "Primary definition remains `RI_tail = mean(raw SV in tail ROI) / mean(raw SV in vessel ROI)`. No background subtraction is used.",
           "B-scans are spatial samples inside one scan volume. No p-values or independent-flow inference are reported.","",
           "## 1. Source-tail-ratio coupling","",
           "| scan | source vs RI Spearman | tail vs RI Spearman | source vs tail Spearman | source IQR/median | tail IQR/median | RI IQR/median |",
           "|---|---:|---:|---:|---:|---:|---:|"]
    for r in coupling.itertuples(index=False):
        lines.append(f"| {r.scan_id} | {r.spearman_source_vs_ri:.3f} | {r.spearman_tail_vs_ri:.3f} | {r.spearman_source_vs_tail:.3f} | {r.source_iqr_over_median:.3f} | {r.tail_iqr_over_median:.3f} | {r.ri_iqr_over_median:.3f} |")
    lines += ["", "Flow07 volume-level context:", ""]
    for m,label in [("source_mean_median","vessel raw mean"),("tail_mean_median","tail raw mean"),("ri_tail_median","RI_tail")]:
        r=f7ctx.loc[m]; lines.append(f"- {label}: flow07 vs flow01 **{r.flow07_vs_flow01_pct:+.2f}%**; vs median of other four scan medians **{r.flow07_vs_other4_median_pct:+.2f}%**.")
    lines += ["", "The negative source-vs-RI association is expected partly from the ratio definition itself. It is reported as a coupling audit, not as a biological association.","",
              "## 2. Exact tail-window sensitivity","", "| window (µm) | flow01 | flow03 | flow05 | flow07 | flow10 | flow07 vs flow01 |",
              "|---:|---:|---:|---:|---:|---:|---:|"]
    for L,g in windows.groupby("window_um"):
        d={int(r.flow_mm_s):r.ri_median for r in g.itertuples(index=False)}; cc=wctx[wctx.window_um.eq(L)].iloc[0]
        lines.append(f"| {L:.0f} | {d[1]:.4f} | {d[3]:.4f} | {d[5]:.4f} | **{d[7]:.4f}** | {d[10]:.4f} | **{cc.flow07_vs_flow01_pct:+.2f}%** |")
    lines += ["", "All four windows use exact native-pixel integration from the vessel bottom; no 10-µm anchor approximation is used.","",
              "## 3. Five-segment within-volume spatial stability","", "| scan | RI segment min | RI segment max | range / full median | max |segment-full| |",
              "|---|---:|---:|---:|---:|"]
    for r in st.itertuples(index=False):
        lines.append(f"| {r.scan_id} | {r.ri_segment_min:.4f} | {r.ri_segment_max:.4f} | {r.ri_segment_range_pct_of_full:.1f}% | {r.max_abs_ri_segment_delta_pct_vs_full:.1f}% |")
    f7seg=sctx[(sctx.metric=="ri_median")]
    lines += ["", "Flow07 RI relative to the median of the other four scans in the same 100-frame segment:",""]
    for r in f7seg.itertuples(index=False): lines.append(f"- frames {r.frame_lo}-{r.frame_hi}: **{r.flow07_vs_other4_median_pct:+.2f}%**, rank {r.flow07_rank_desc}/5.")
    lines += ["", "## Validation","", f"- valid frames: **{EXPECTED_VALID}**", f"- exact-window raw-array packages verified: **{val['release_packages_verified']}/{EXPECTED_PACKAGES}**",
              f"- valid-frame NPZ SHA verified: **{val['valid_frame_npz_sha256_verified']}**", f"- 500-µm RI replay max relative error: **{val['ri500_replay_max_relative_error']:.3e}**",
              f"- log identity `log(RI)=log(tail)-log(source)` max absolute numerical error: **{coupling.log_identity_max_abs_error.max():.3e}**",
              "- no background subtraction, interpolation, zero-fill, p-values, or B-scan pseudo-replication.","",
              "## Files","", "- `coupling_summary.csv`, `coupling_source_deciles.csv`, `coupling_source_quintiles.csv`, `flow07_coupling_context.csv`",
              "- `window_framewise.csv.gz`, `window_scan_summary.csv`, `window_flow07_context.csv`",
              "- `segment_summary.csv`, `segment_stability_by_scan.csv`, `segment_flow07_context.csv`",
              "- `validation.json`, `provenance.json`, `input_sha256.csv`, `release_package_audit.csv`"]
    (out/"README.md").write_text("\n".join(lines)+"\n",encoding="utf-8")


def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--root",default="."); ap.add_argument("--output-dir",default="analysis/formal_sv_d128_v21_run001/no_background_followup_three_audits_v1_full2422"); ap.add_argument("--source-sha",required=True); ap.add_argument("--workflow-source-sha",required=True); args=ap.parse_args()
    root=Path(args.root).resolve(); out=(root/args.output_dir).resolve(); out.mkdir(parents=True,exist_ok=True)
    primary_path=root/"analysis/formal_sv_d128_v21_run001/no_background_relative_tail_v1_full2422/framewise_primary.csv"
    primary=pd.read_csv(primary_path)
    if len(primary)!=EXPECTED_VALID or primary.groupby("scan_id").size().to_dict()!=EXPECTED_COUNTS: raise AssertionError("primary coverage mismatch")
    if primary.duplicated(["scan_id","frame_index_0based"]).any(): raise AssertionError("duplicate primary keys")
    for col in ["source_mean_raw","tail_mean_raw","ri_tail"]:
        if not np.isfinite(primary[col].to_numpy(float)).all(): raise AssertionError(f"nonfinite {col}")
    coupling,deciles,ctx=coupling_audit(primary,out)
    wf,ws,wctx,wval=exact_windows(root,primary,out)
    seg,st,sctx=segment_audit(primary,out)
    validation={"valid_frames":EXPECTED_VALID,"scan_valid_counts":EXPECTED_COUNTS,"window_frame_rows":len(wf),"window_lengths_um":list(WINDOWS_UM),
                "segments":[[a,b] for a,b in SEGMENTS],**wval,"p_values_computed":False,"background_subtracted":False,
                "invalid_frames_interpolated":False,"invalid_frames_zero_filled":False,"bscans_treated_as_independent_replicates":False}
    (out/"validation.json").write_text(json.dumps(validation,indent=2)+"\n",encoding="utf-8")
    provenance={"analysis":"three descriptive follow-ups for no-background RI_tail","source_sha":args.source_sha,"workflow_source_sha":args.workflow_source_sha,
                "primary_definition":"mean raw SV tail / mean raw SV vessel","audits":["source-tail-ratio coupling","exact 100/200/300/500 um windows","five 100-frame spatial segments"],
                "statistical_scope":"within-volume descriptive only; no flow p-values; B-scans not independent experimental units"}
    (out/"provenance.json").write_text(json.dumps(provenance,indent=2)+"\n",encoding="utf-8")
    inputs=[primary_path,root/"results/formal_sv_d128_v21_full2500_run001/download_packages.csv",root/"results/formal_sv_d128_v21_full2500_run001/arrays_sha256.csv",root/"results/formal_sv_d128_v21_full2500_run001/run_config.json",
            root/"analysis/formal_sv_d128_v21_run001/fixed128_width_sensitivity_full2422/audit_fixed128_width.py",root/"src/svrecttail/geometry.py"]
    pd.DataFrame([{"input_path":p.relative_to(root).as_posix(),"sha256":sha256(p),"size_bytes":p.stat().st_size} for p in inputs]).to_csv(out/"input_sha256.csv",index=False)
    write_readme(out,coupling,ctx,ws,wctx,seg,st,sctx,validation)
    print(json.dumps({"valid_frames":EXPECTED_VALID,"coupling":coupling[["scan_id","spearman_source_vs_ri","spearman_tail_vs_ri"]].to_dict("records"),
                      "window_flow07":wctx.to_dict("records"),"segment_flow07_ri":sctx[sctx.metric.eq("ri_median")].to_dict("records"),"validation":wval},indent=2))

if __name__=="__main__": main()
