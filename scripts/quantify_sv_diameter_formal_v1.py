"""Frozen diameter-stage identity validation and no-background raw-SV quantification.

Only volume-level descriptive aggregation is performed. Local paths locate inputs;
committed manifests and hashes define identity and experimental conditions.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import platform
import re
import subprocess
import sys
import time
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from svrecttail.geometry import VesselGeometry, ellipse_is_complete, interval_is_complete
from svrecttail.mentor_tracking import build_localization_from_tracking
from svrecttail.raw_sv_quantification import quantify_raw_sv, METRICS, relative_error

BASE = "d250615f84b371e0cc352abf5a22e4453de7b7b0"
OUT = ROOT / "analysis/formal_sv_diameter_v1"
LOCAL = ROOT / "outputs/formal_quant_v1"
INTAKE = ROOT / "results/diameter_intake_v1"
D128 = ROOT / "analysis/formal_sv_d128_v21_run001"
FORMAL128 = ROOT / "results/formal_sv_d128_v21_full2500_run001"
TRACK_CONFIG = ROOT / "config/tracking_config.continuity_first_v2_1.a020_n4.json"
KERNEL = ROOT / "src/svrecttail/raw_sv_quantification.py"
CONFIG = dict(
    schema_version="SV_Diameter_NoBackground_RI_v1",
    signal=dict(name="sv_raw", definition="var(abs(E), 1, 3)", variance_denominator="N",
                normalization=False, log=False, gain=False, clipping=False,
                positive_truncation=False, background_subtraction=False),
    geometry=dict(lateral_centre="X4", lateral_span="X1 apparent width",
                  axial_source="manifest physical diameter", dx_um=12.7, dz_um=6.7,
                  coordinates="zero-based pixel centres; pixel i spans [i-0.5,i+0.5]",
                  ellipse_supersample=16, tail_gap_um=0, windows_um=[100,200,300,500],
                  primary_window_um=500, bottom="z_top_edge_px + diameter_um / dz_um; unrounded"),
    depth=dict(targets_um=list(range(0,501,10)), selection="nearest native axial centre among pixels with positive 500um tail overlap",
               interpolation=False, tie_break="first ascending native row", max_error="dz_um/2 + 1e-9 um"),
    entry_rule="geometry_qc_valid == true", experimental_unit="scan volume",
    bscan="slow-axis spatial sample; not independent experiment",
    aggregation=dict(quantiles="numpy linear", sd_ddof=1),
    numerical_tolerances=dict(d128_relative=1e-10, direct_profile_relative=1e-12, area_relative=1e-12))


def sha(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def write_json(path, obj):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open("w", encoding="utf-8", newline="\n") as f:
        json.dump(obj, f, ensure_ascii=False, indent=2, allow_nan=False)
        f.write("\n")


def write_csv(path, frame):
    frame.to_csv(path, index=False, float_format="%.17g", lineterminator="\n",
                 compression={"method":"gzip", "mtime":0} if str(path).endswith(".gz") else None)


def protected_snapshot():
    paths = subprocess.check_output(["git", "ls-files", "analysis/formal_sv_d128_v21_run001",
        "results/formal_sv_d128_v21_full2500_run001", "results/diameter_intake_v1",
        "data/diameter_stage_v1"], cwd=ROOT, text=True).splitlines()
    return {p:sha(ROOT/p) for p in paths}


def prepare():
    if subprocess.check_output(["git","merge-base","HEAD",BASE],cwd=ROOT,text=True).strip() != BASE:
        raise ValueError("Wrong base")
    LOCAL.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    before = LOCAL/"protected_before.json"
    if not before.exists():
        write_json(before, protected_snapshot())
    if protected_snapshot() != read_json(before):
        raise ValueError("Protected inputs changed")
    if (OUT/"formal_config.json").exists() and read_json(OUT/"formal_config.json") != CONFIG:
        raise ValueError("Frozen formal config differs")
    write_json(OUT/"formal_config.json", CONFIG)


def validate_mat(path, audit, row, source_path):
    """Hash and parse the same bytes, before any numerical use, for every frame."""
    if not path.is_file() or not path.stat().st_size:
        raise ValueError(f"Missing/empty MAT: {row['scan_id']}/{path.name}")
    data = path.read_bytes()
    digest = hashlib.sha256(data).hexdigest()
    if digest != audit.sha256 or len(data) != int(audit.size_bytes):
        raise ValueError(f"Frozen MAT SHA/size mismatch: {row['scan_id']}/{path.name}")
    frame = int(audit.frame_index)
    if path.name != f"frame_{frame:03d}.mat":
        raise ValueError("Frame filename mismatch")
    mat = loadmat(io.BytesIO(data), variable_names=["sv_raw", "metadata"], simplify_cells=True)
    m, sv = mat["metadata"], mat["sv_raw"]
    if Path(m["source_file"]).resolve() != source_path.resolve():
        raise ValueError("MAT OCT identity mismatch")
    if m["bscan_index_matlab_1based"] != frame+1:
        raise ValueError("MAT B-scan identity mismatch")
    if (m["formal_signal_name"],m["formal_signal_definition"],m["variance_denominator"]) != (
        "sv_raw","var(abs(E), 1, 3)","N"):
        raise ValueError("MAT signal definition mismatch")
    rec = m["reconstruction"]
    if Path(rec["source_file"]).resolve() != source_path.resolve() or rec["bscan_index"] != frame+1:
        raise ValueError("MAT reconstruction identity mismatch")
    if rec["n_bscan_positions"] != int(row["bscan_count"]) or m["dimension_order"] != "depth x A-line":
        raise ValueError("MAT dimensions metadata mismatch")
    h = rec["header"]
    if (h["SPL"],h["nX"],h["frame_per_pos"]) != (int(row["spectral_samples"]),int(row["alines_per_bscan"]),int(row["repeat_count"])):
        raise ValueError("MAT acquisition metadata mismatch")
    if sv.shape != (351,500) or not np.isfinite(sv).all():
        raise ValueError("MAT shape/finite check failed")
    return sv, digest


def geometry_record(record, row):
    details = dict(geometry_qc_valid=False, invalid_reason="", physical_diameter_um=float(row["diameter_um"]))
    try:
        loc = build_localization_from_tracking(record, diameter_um=float(row["diameter_um"]),
                dx_um=float(row["dx_um"]), dz_um=float(row["dz_um"]))
        g = loc.geometry
        complete = ellipse_is_complete((351,500),g) and interval_is_complete(
            351,g.z_bottom_edge_px,g.z_bottom_edge_px+500/g.dz_um) and interval_is_complete(500,g.x_left_edge_px,g.x_right_edge_px)
        details.update(geometry_qc_valid=bool(loc.source_qc_valid and complete),
            invalid_reason=loc.mentor_tracking.invalid_reason if not loc.source_qc_valid else ("ok" if complete else "window_out_of_bounds"),
            x4_px=g.x_center_px, x_left_edge_px=g.x_left_edge_px, x_right_edge_px=g.x_right_edge_px,
            apparent_width_px=g.x_right_edge_px-g.x_left_edge_px, apparent_width_um=g.lateral_width_um,
            z_top_edge_px=g.z_top_edge_px,z_bottom_edge_px=g.z_bottom_edge_px)
        return details,g
    except ValueError as error:
        details["invalid_reason"] = str(error)
        return details,None


def run_volumes(intake_local):
    prepare()
    manifest = pd.read_csv(ROOT/"data/diameter_stage_v1/volume_manifest.csv")
    selected = manifest[manifest.diameter_um.ne(128)].sort_values(["diameter_um","flow_speed_mm_s","scan_id"])
    expected = pd.read_csv(INTAKE/"processing_summary.csv").set_index("scan_id")
    assert len(selected)==18 and int(selected.bscan_count.sum())==9000
    assert set(selected.scan_id)==set(expected.index) and int(expected.geometry_valid_frames.sum())==8509
    jobs = {j["scan_id"]:j for j in read_json(intake_local/"outputs/local_run_inputs/diameter_stage_v1/export_jobs.json")}
    local_manifest = pd.read_csv(intake_local/"outputs/local_run_inputs/diameter_stage_v1/volume_manifest_local.csv").set_index("scan_id")
    all_frames,all_depths,tracking_identity = [],[],[]
    for row in selected.to_dict("records"):
        scan=row["scan_id"]
        print("QUANTIFY "+scan,flush=True)
        # Local scientific fields are checked only to reject stale locators.
        local=local_manifest.loc[scan]
        for field in ("oct_sha256","flow_dicom_sha256","oct_relative_path","flow_dicom_relative_path"):
            assert local[field]==row[field], (scan,field)
        source=Path(local.oct_absolute_path)
        assert source.as_posix().endswith(row["oct_relative_path"])
        job=jobs[scan]
        assert Path(job["source_file"]).resolve()==source.resolve()
        directory=Path(job["interim_dir"])
        audit=pd.read_csv(INTAKE/"exports"/f"{scan}_frames.csv")
        assert audit.scan_id.eq(scan).all() and np.array_equal(audit.frame_index,np.arange(500))
        tracking_dir=intake_local/"outputs/diameter_stage_v1/tracking"/scan
        meta=read_json(tracking_dir/f"{scan}_mentor_tracking.metadata.json")
        summary=read_json(INTAKE/"tracking"/f"{scan}_summary.json")
        tp=tracking_dir/meta["primary_table"]
        assert sha(tp)==summary["tracking_table_sha256"]==meta["primary_table_sha256"]
        assert summary["tracking_config_sha256"]==sha(TRACK_CONFIG)
        assert meta["effective_tracking_config"]==read_json(TRACK_CONFIG)["tracking"]
        assert meta["flow_dicom_sha256"]==row["flow_dicom_sha256"] and meta["diameter_um"]==row["diameter_um"]
        track=pd.read_csv(tp)  # identical parser to intake geometry computation
        assert track.scan_id.eq(scan).all() and np.array_equal(track.frame_index,np.arange(500))
        assert track.diameter_um.eq(row["diameter_um"]).all()
        frozen_geom=pd.read_csv(tracking_dir/"geometry_validation.csv")
        assert np.array_equal(frozen_geom.frame_index,np.arange(500))
        tracking_identity.append(dict(scan_id=scan,tracking_table_sha256=sha(tp),
                                      intake_geometry_table_sha256=sha(tracking_dir/"geometry_validation.csv")))
        count=0
        for a in audit.itertuples(index=False):
            fi=int(a.frame_index)
            rec=track.iloc[fi]
            detail,g=geometry_record(rec,row)
            assert detail["geometry_qc_valid"]==bool(frozen_geom.iloc[fi].geometry_qc_valid),(scan,fi,"QC mismatch")
            sv,digest=validate_mat(directory/a.file_name,a,row,source)
            metrics,anchors=quantify_raw_sv(sv,g,geometry_qc_valid=detail["geometry_qc_valid"])
            key=dict(scan_id=scan,diameter_um=float(row["diameter_um"]),flow_mm_s=float(row["flow_speed_mm_s"]),frame_index_0based=fi)
            item=dict(**key,**detail,valid=detail["geometry_qc_valid"],
                tracking_status=rec.new_tracking_class,assessability=rec.vessel_presence_prediction,
                localization_source="short_gap_fill" if rec.z_short_gap_filled else ("direct_candidate" if rec.z_candidate_accepted else "unavailable"),
                candidate_status=rec.z_continuity_status, input_mat_sha256=digest,
                input_identity_valid=True,reconstructed_after_intake=False,**metrics)
            all_frames.append(item)
            all_depths.extend(dict(**key,**anchor) for anchor in anchors)
            count+=detail["geometry_qc_valid"]
        assert count==int(expected.loc[scan,"geometry_valid_frames"])==summary["geometry_qc_valid"],scan
        print(f"VERIFIED {scan}: 500 MAT, {count} valid",flush=True)
    fw=pd.DataFrame(all_frames)
    dep=pd.DataFrame(all_depths)
    validate_tables(fw,dep,expected)
    write_csv(OUT/"framewise_all.csv.gz",fw)
    primary=fw[fw.valid].copy()
    write_csv(OUT/"framewise_primary.csv",primary)
    write_csv(OUT/"depth_anchor_framewise.csv.gz",dep)
    summaries=[]
    metrics=("source_mean_raw","tail_mean_raw_500um","ri_tail",*[f"ri_tail_{n}um" for n in (100,200,300,500)])
    for scan,group in fw.groupby("scan_id",sort=True):
        v=group[group.valid]
        summary=dict(scan_id=scan,diameter_um=float(group.diameter_um.iloc[0]),flow_mm_s=float(group.flow_mm_s.iloc[0]),
                     nominal_frames=len(group),geometry_valid_frames=len(v),invalid_frames=len(group)-len(v),geometry_valid_fraction=len(v)/len(group))
        for metric in metrics:
            values=v[metric].to_numpy()
            q1,med,q3=np.quantile(values,[.25,.5,.75])
            summary.update({f"{metric}_{k}":float(x) for k,x in zip(("q1","median","q3","mean","sd"),(q1,med,q3,values.mean(),values.std(ddof=1)))})
        summaries.append(summary)
    write_csv(OUT/"volume_summary.csv",pd.DataFrame(summaries))
    depth_summary=[]
    for (scan,target),group in dep.groupby(["scan_id","target_r_um"],sort=True):
        q1,med,q3=np.quantile(group.RI_r,[.25,.5,.75])
        depth_summary.append(dict(scan_id=scan,diameter_um=float(group.diameter_um.iloc[0]),flow_mm_s=float(group.flow_mm_s.iloc[0]),
            target_r_um=float(target),n_frames=len(group),selected_r_um_median=float(group.selected_r_um.median()),
            max_abs_sampling_error_um=float(group.abs_sampling_error_um.max()),RI_r_q1=q1,RI_r_median=med,RI_r_q3=q3,
            raw_row_signal_median=float(group.S_tail_raw_r.median())))
    ds=pd.DataFrame(depth_summary)
    write_csv(OUT/"depth_summary_10um.csv",ds)
    write_csv(OUT/"depth_selected.csv",ds[ds.target_r_um.isin([0,50,100,200,300,400,500])])
    write_json(OUT/"tracking_identity.json",tracking_identity)
    write_json(LOCAL/"quantification_validation.json",dict(nominal_frames=len(fw),validated_input_mat_frames=len(fw),
        volumes_processed=18,geometry_valid_frames=len(primary),invalid_frames=len(fw)-len(primary),
        per_volume_valid_counts={r["scan_id"]:r["geometry_valid_frames"] for r in summaries},depth_anchor_rows=len(dep),
        input_mat_sha_checks="9000/9000 matched frozen intake hashes",mat_sha_mismatches=0,
        direct_vs_profile_max_relative_error=float(primary.direct_profile_max_relative_error.max()),
        tail_area_identity_max_relative_error=float(primary.tail_area_max_relative_error.max()),
        maximum_depth_sampling_error_um=float(dep.abs_sampling_error_um.max()),
        background_subtracted=False,normalization=False,log=False,clipping=False,interpolation=False,zero_fill=False,
        p_values_computed=False,bscans_treated_as_independent_replicates=False))


def validate_tables(fw,dep,expected):
    assert len(fw)==9000 and not fw.duplicated(["scan_id","frame_index_0based"]).any()
    for scan,group in fw.groupby("scan_id"):
        assert np.array_equal(group.frame_index_0based,np.arange(500))
        assert int(group.valid.sum())==int(expected.loc[scan,"geometry_valid_frames"])
    v=fw[fw.valid]
    assert len(v)==int(expected.geometry_valid_frames.sum())==8509
    assert fw.loc[~fw.valid,list(METRICS)].isna().all().all()
    assert np.isfinite(v[list(METRICS)].to_numpy()).all() and v.source_mean_raw.gt(0).all()
    assert np.array_equal(v.ri_tail,v.ri_tail_500um)
    assert np.array_equal(v.ri_tail,v.tail_mean_raw_500um/v.source_mean_raw)
    assert len(dep)==len(v)*51 and dep.groupby(["scan_id","frame_index_0based"]).size().eq(51).all()
    assert dep.abs_sampling_error_um.max()<=6.7/2+1e-9


def regression(packages_dir, wait_seconds=0):
    prepare()
    refs=pd.read_csv(D128/"no_background_relative_tail_v1_full2422/framewise_primary.csv",float_precision="round_trip").set_index(["scan_id","frame_index_0based"])
    observed=pd.read_csv(D128/"observed_tail_intensity_full2422/observed_tail_intensity_framewise.csv",float_precision="round_trip").set_index(["scan_id","frame_index_0based"])
    geometry=pd.read_csv(FORMAL128/"frame_results.csv",float_precision="round_trip").set_index(["scan_id","frame_index_0based"])
    array_hashes=pd.read_csv(FORMAL128/"arrays_sha256.csv").set_index(["scan_id","frame_index_0based"])
    packages=pd.read_csv(FORMAL128/"download_packages.csv")
    maxima=dict(source=0.,tail=0.,ri=0.)
    seen=set()
    package_audit=[]
    for pkg in packages.itertuples(index=False):
        path=packages_dir/pkg.file
        deadline=time.monotonic()+wait_seconds
        if not path.exists() and wait_seconds:
            print("WAIT FOR FROZEN PACKAGE "+pkg.file,flush=True)
        while not path.exists() and time.monotonic()<deadline:
            time.sleep(5)
        if not path.exists():
            raise FileNotFoundError("Missing frozen D128 package: "+pkg.file)
        assert sha(path)==pkg.sha256 and path.stat().st_size==pkg.bytes,pkg.file
        scan,lo,hi=re.search(r"_(flow\d+)_(\d{3})_(\d{3})\.zip$",pkg.file).groups()
        count=0
        with zipfile.ZipFile(path) as archive:
            for fi in range(int(lo),int(hi)+1):
                key=(scan,fi)
                if key not in refs.index:
                    continue
                assert key not in seen
                row=geometry.loc[key]
                data=archive.read(f"arrays/{scan}/frame_{fi:03d}.npz")
                assert hashlib.sha256(data).hexdigest()==array_hashes.loc[key,"sha256"]
                with np.load(io.BytesIO(data),allow_pickle=False) as arrays:
                    sv=arrays["sv_raw"]
                g=VesselGeometry(float(row.x_left_edge_px),float(row.x_right_edge_px),float(row.z_top_edge_px),128.,12.7,6.7)
                values,_=quantify_raw_sv(sv,g)
                for name,metric,refmetric,obsmetric in (
                    ("source","source_mean_raw","source_mean_raw","source_mean_observed"),
                    ("tail","tail_mean_raw_500um","tail_mean_raw","tail_mean_observed"),
                    ("ri","ri_tail","ri_tail","ri_tail_observed")):
                    maxima[name]=max(maxima[name],relative_error(values[metric],float(refs.loc[key,refmetric])),
                                     relative_error(values[metric],float(observed.loc[key,obsmetric])))
                assert max(maxima.values())<=1e-10,(key,maxima)
                seen.add(key)
                count+=1
        package_audit.append(dict(file=pkg.file,sha256=pkg.sha256,frames_replayed=count))
        print(f"REPLAY {pkg.file}: {count} passed",flush=True)
    assert len(seen)==len(refs)==2422
    write_json(OUT/"d128_nobg_regression.json",dict(status="passed",frames_replayed=len(seen),
        source_max_relative_error=maxima["source"],tail_max_relative_error=maxima["tail"],RI_max_relative_error=maxima["ri"],
        tolerance=1e-10,error_definition="abs(new-reference)/abs(reference)",
        reference_sets=["observed_tail_intensity_full2422","no_background_relative_tail_v1_full2422"],
        packages=package_audit,kernel_sha256=sha(KERNEL)))


def has_absolute_windows_path(content):
    # A URL scheme ending in 's:/' is not a drive; accept only a standalone
    # drive letter, and also reject UNC paths (including JSON-escaped forms).
    return bool(re.search(r"(?<![A-Za-z0-9])[A-Za-z]:[\\/]|\\\\[A-Za-z0-9_.-]+\\",content))


def public_path_audit():
    for path in OUT.iterdir():
        if path.is_file():
            content=gzip.decompress(path.read_bytes()).decode("utf-8") if path.suffix==".gz" else path.read_text(encoding="utf-8")
            if has_absolute_windows_path(content):
                raise ValueError("Absolute Windows path in "+path.name)


def finalize():
    prepare()
    val=read_json(LOCAL/"quantification_validation.json")
    reg=read_json(OUT/"d128_nobg_regression.json")
    assert reg["status"]=="passed" and reg["frames_replayed"]==2422 and reg["kernel_sha256"]==sha(KERNEL)
    before=read_json(LOCAL/"protected_before.json")
    after=protected_snapshot()
    assert before==after,"Protected inputs changed"
    write_json(OUT/"protected_sha256.json",dict(modified_files=0,before=before,after=after))
    val.update(status="COMPLETE_LOCAL_QUANTIFICATION",d128_regression_status="passed",protected_files_unchanged=len(before))
    write_json(OUT/"validation.json",val)
    refs=[D128/"no_background_relative_tail_v1_full2422/framewise_primary.csv",
          D128/"observed_tail_intensity_full2422/observed_tail_intensity_framewise.csv"]
    registry=pd.read_csv(ROOT/"data/diameter_stage_v1/volume_manifest.csv")
    write_json(OUT/"d128_reference.json",dict(included_in_new_outputs=False,
        mappings=[dict(legacy_scan_id=r.legacy_scan_id,registry_scan_id=r.scan_id,oct_sha256=r.oct_sha256,flow_dicom_sha256=r.flow_dicom_sha256)
                  for r in registry[registry.diameter_um.eq(128)].itertuples(index=False)],
        authoritative_results=[dict(path=p.relative_to(ROOT).as_posix(),sha256=sha(p)) for p in refs]))
    inputs=[ROOT/"data/diameter_stage_v1/volume_manifest.csv",INTAKE/"processing_summary.csv",INTAKE/"provenance.json",
            ROOT/"src/svrecttail/geometry.py",ROOT/"src/svrecttail/mentor_tracking.py",TRACK_CONFIG,
            FORMAL128/"frame_results.csv",FORMAL128/"arrays_sha256.csv",FORMAL128/"download_packages.csv",FORMAL128/"run_config.json",
            OUT/"formal_config.json",*refs,*sorted((INTAKE/"exports").glob("*_frames.csv")),*sorted((INTAKE/"tracking").glob("*_summary.json"))]
    write_csv(OUT/"input_sha256.csv",pd.DataFrame([dict(path=p.relative_to(ROOT).as_posix(),sha256=sha(p)) for p in inputs]))
    intake_prov=read_json(INTAKE/"provenance.json")
    write_json(OUT/"provenance.json",dict(repository="bulbel-magnolia/OCTA-",base_branch="codex/sv-diameter-intake-v1",base_sha=BASE,
        working_branch="codex/sv-diameter-formal-quant-v1",output_commit_sha=None,
        output_commit_note="Filled with the metrics commit SHA in a subsequent provenance-only commit; a commit cannot embed its own SHA.",
        python_version=platform.python_version(),packages=dict(numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__),
        matlab_exporter_provenance=dict(intake_provenance_sha256=sha(INTAKE/"provenance.json"),matlab_version=intake_prov["matlab_version"],
            source_sha256={k:v for k,v in intake_prov["source_sha256"].items() if k.startswith("matlab/")},exporter_modified=False,reconstruction_repeated=False),
        volume_manifest_sha256=sha(inputs[0]),tracking_config_sha256=sha(TRACK_CONFIG),geometry_py_sha256=sha(ROOT/"src/svrecttail/geometry.py"),
        new_quantifier_sha256=sha(KERNEL),driver_sha256=sha(Path(__file__)),formal_config_sha256=sha(OUT/"formal_config.json"),
        per_volume_input_reference="data/diameter_stage_v1/volume_manifest.csv; OCT/DICOM SHA and calibration are authoritative",
        calibration="inherited project acquisition calibration; not independently measured",
        per_frame_mat_hash_source="results/diameter_intake_v1/exports/<scan_id>_frames.csv; actual hashes also in framewise_all.csv.gz",
        local_tracking_identity="tracking_identity.json; table SHA anchored in committed intake summary",
        signal_definition=CONFIG["signal"],geometry_definition=CONFIG["geometry"],depth_sampling_rule=CONFIG["depth"],
        experimental_unit=CONFIG["experimental_unit"],bscan=CONFIG["bscan"],background_subtraction=False))
    summary=pd.read_csv(OUT/"volume_summary.csv")
    lines=["# SV Diameter Formal Quantification v1","",
        "RI_tail = mean(raw SV in 500um tail ROI) / mean(raw SV in vessel ROI).",
        "RI(r) = raw SV lateral mean at nearest native row / raw vessel mean.","",
        "Input: linear sv_raw = var(abs(E), 1, 3), denominator N. No background subtraction, log, normalization, gain, clipping or positive truncation.","",
        "Geometry: frozen continuity-first v2.1 X4 centre, X1 apparent lateral span, actual manifest diameter for source axial height; dx=12.7 and dz=6.7 um. Calibration is inherited from the project acquisition protocol.",
        "Source: fractional ellipse (supersample 16). Tail: fractional rectangle with the same lateral span, unrounded physical bottom, guard 0, windows 100/200/300/500 um; primary alias is exactly the 500um value.",
        "RI(r): 51 targets at 0..500 um in 10um steps, nearest native pixel with positive tail overlap (including fractional boundary pixels). No interpolation; maximum error <= dz/2 + 1e-9 um.","",
        "Entry: geometry_qc_valid only. All 9000 input rows retained; 8509 valid and 491 invalid. Invalid metric fields are NA, with reasons preserved.","",
        "| scan_id | nominal | valid | invalid |","|---|---:|---:|---:|"]
    lines += [f"| {r.scan_id} | {r.nominal_frames} | {r.geometry_valid_frames} | {r.invalid_frames} |" for r in summary.itertuples(index=False)]
    lines += ["","Validation: every MAT SHA and embedded identity matched intake; all per-volume QC counts matched; full 2422-frame D128 replay passed. Direct/profile integrals and analytical rectangle areas agree within the frozen numerical tolerances. See validation.json and d128_nobg_regression.json.","",
        "Files: framewise_all.csv.gz (9000 rows); framewise_primary.csv (8509); depth_anchor_framewise.csv.gz (433959); volume_summary.csv (18); depth_summary_10um.csv (918); depth_selected.csv (126).",
        "Framewise fields retain source/tail integrals, fractional areas, raw means, four ratios, geometry and input MAT hashes. Depth fields retain target/native coordinates, fractional overlap and raw row signal. Summary tables use within-volume Q1/median/Q3 and, where specified, mean/sample SD (ddof=1).",
        "formal_config.json freezes definitions; input_sha256.csv, tracking_identity.json, protected_sha256.json and provenance.json record identities. d128_reference.json links historical results without adding D128 rows to these outputs.","",
        "Experimental unit: scan volume. A B-scan is a slow-axis spatial sample, not an independent experiment. Aggregation is descriptive within each volume; no condition comparisons, regression, significance testing, fitting or scientific interpretation were performed."]
    (OUT/"README.md").write_text("\n".join(lines)+"\n",encoding="utf-8",newline="\n")
    public_path_audit()


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phase",choices=["quantify","regression","finalize"])
    parser.add_argument("--intake-local",type=Path)
    parser.add_argument("--packages-dir",type=Path,default=LOCAL/"d128_packages")
    parser.add_argument("--wait-for-packages-seconds",type=int,default=0,
                        help="Bounded wait per package for a separate authenticated download")
    args=parser.parse_args()
    if args.phase=="quantify":
        if args.intake_local is None:
            parser.error("--intake-local is required to locate existing inputs")
        run_volumes(args.intake_local.resolve())
    elif args.phase=="regression":
        regression(args.packages_dir.resolve(),args.wait_for_packages_seconds)
    else:
        finalize()
