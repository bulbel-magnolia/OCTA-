"""Frozen localization and raw SV export orchestration. No tail intensities are computed."""
from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import loadmat

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from svrecttail.mentor_tracking import load_flow_dicom, build_localization_from_tracking
from svrecttail.geometry import ellipse_is_complete, interval_is_complete
from prepare_diameter_stage_inputs import sha256, write_json, write_csv

REPO = Path(__file__).resolve().parents[1]
LOCAL = REPO / "outputs/local_run_inputs/diameter_stage_v1"
RESULTS = REPO / "results/diameter_intake_v1"
LOGS = REPO / "outputs/diameter_stage_v1/logs"
CONFIG = REPO / "config/tracking_config.continuity_first_v2_1.a020_n4.json"


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def load_inputs():
    rows = list(csv.DictReader((LOCAL / "volume_manifest_local.csv").open(encoding="utf-8")))
    return {r["scan_id"]: r for r in rows if r["scan_id"]}, {j["scan_id"]:j for j in read_json(LOCAL / "export_jobs.json")}


def run(command, log):
    LOGS.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["PYTHONPATH"] = str(REPO / "src")
    env["PYTHONIOENCODING"] = "utf-8"
    with Path(log).open("w", encoding="utf-8") as stream:
        completed = subprocess.run(command, cwd=REPO, env=env, stdout=stream, stderr=subprocess.STDOUT)
    if completed.returncode:
        raise RuntimeError(f"Subprocess exit {completed.returncode}; see local log {Path(log).name}")


def distribution(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    return dict(n=int(len(values)), **{k:float(v) for k,v in zip(("min","p25","median","p75","max"), np.quantile(values,[0,.25,.5,.75,1]))}) if len(values) else dict(n=0)


def segments(mask):
    edges = np.diff(np.r_[False, np.asarray(mask,bool), False].astype(int))
    return [dict(first_frame=int(a), last_frame=int(b-1), frame_count=int(b-a)) for a,b in zip(np.flatnonzero(edges == 1), np.flatnonzero(edges == -1))]


def tracking(row, *, make_qc=False):
    scan = row["scan_id"]
    output = REPO / "outputs/diameter_stage_v1/tracking" / scan
    metadata_path = output / f"{scan}_mentor_tracking.metadata.json"
    volume = load_flow_dicom(row["flow_dicom_absolute_path"])
    config = read_json(CONFIG)["tracking"]
    expected = (int(config["expected_frames"]), int(config["expected_z_px"]), int(config["expected_x_px"]))
    if volume.shape != expected or volume.shape[0] != int(row["bscan_count"]):
        raise ValueError("DICOM shape incompatible with frozen tracker; no reshape permitted")
    if not metadata_path.exists():
        run([sys.executable,"-m","svrecttail","mentor-track","--flow-dicom",row["flow_dicom_absolute_path"],"--scan-id",scan,"--diameter-um",row["diameter_um"],"--tracking-config",str(CONFIG),"--output",str(output)], LOGS / f"{scan}_tracking.log")
    meta = read_json(metadata_path)
    assert meta["scan_id"] == scan and meta["diameter_um"] == float(row["diameter_um"])
    assert meta["flow_dicom_sha256"] == row["flow_dicom_sha256"]
    assert all(meta["effective_tracking_config"][k] == v for k,v in config.items()), "Frozen config mismatch"
    path = output / meta["primary_table"]
    assert sha256(path) == meta["primary_table_sha256"]
    table = pd.read_csv(path)
    assert np.array_equal(table.frame_index.to_numpy(), np.arange(volume.shape[0]))
    assert (table.scan_id == scan).all() and (table.diameter_um == float(row["diameter_um"])).all()
    valid = []
    geometry_rows = []
    for _, record in table.iterrows():
        detail = dict(scan_id=scan, frame_index=int(record.frame_index), geometry_qc_valid=False, source_qc_valid=False, full_tail_window_fits=False)
        try:
            loc = build_localization_from_tracking(record, diameter_um=float(row["diameter_um"]), dx_um=float(row["dx_um"]), dz_um=float(row["dz_um"]))
            g = loc.geometry
            source_fits = ellipse_is_complete(volume.shape[1:], g)
            tail_fits = interval_is_complete(volume.shape[1],g.z_bottom_edge_px,g.z_bottom_edge_px+500/g.dz_um) and interval_is_complete(volume.shape[2],g.x_left_edge_px,g.x_right_edge_px)
            detail.update(source_qc_valid=loc.source_qc_valid, full_tail_window_fits=tail_fits,
                          geometry_qc_valid=loc.source_qc_valid and source_fits and tail_fits,
                          physical_diameter_um=g.diameter_um, apparent_width_um=g.lateral_width_um,
                          x4_px=g.x_center_px, x_left_edge_px=g.x_left_edge_px, x_right_edge_px=g.x_right_edge_px,
                          z_top_edge_px=g.z_top_edge_px,z_bottom_edge_px=g.z_bottom_edge_px,
                          tail_start_edge_px=g.z_bottom_edge_px,tail_end_edge_px=g.z_bottom_edge_px+500/g.dz_um,
                          invalid_reason=loc.mentor_tracking.invalid_reason if not loc.source_qc_valid else ("window_out_of_bounds" if not source_fits or not tail_fits else "ok"))
        except ValueError as error:
            detail["invalid_reason"] = str(error)
        valid.append(detail["geometry_qc_valid"])
        geometry_rows.append(detail)
    zvalid = table.z_upper_px.notna().to_numpy()
    summary = dict(scan_id=scan, total_frames=len(table), z_localization_valid=int(zvalid.sum()),
                   geometry_qc_valid=sum(valid), failed=int((table.new_tracking_class == "failed").sum()),
                   not_assessable=int((table.vessel_presence_prediction != "assessable").sum()),
                   direct_candidate=int(table.z_candidate_accepted.fillna(False).sum()),
                   short_gap_fill=int(table.z_short_gap_filled.fillna(False).sum()),
                   missing_segment_count=len(segments(~zvalid)), relock_events=int(table.z_relock_start.fillna(False).sum()),
                   dicom_dimensions=list(volume.shape), finite_fraction=float(np.isfinite(volume).mean()),
                   tracking_config_sha256=sha256(CONFIG), tracking_table_sha256=sha256(path),
                   X1_width_px=distribution(table.local_body_run_width_px), X4_px=distribution(table.x4_centroid_isolated_jump_corrected_px), z_top_px=distribution(table.z_upper_px),
                   missing_segments=segments(~zvalid), relock_frames=table.loc[table.z_relock_start.fillna(False),"frame_index"].astype(int).tolist())
    write_json(RESULTS / "tracking" / f"{scan}_summary.json", summary)
    write_csv(output / "geometry_validation.csv", geometry_rows)
    if make_qc and any(valid):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.patches import Ellipse, Rectangle
        indices = np.flatnonzero(valid)
        selected = sorted(set([int(indices[0]),int(indices[len(indices)//2]),int(indices[-1])]))
        qc_paths = []
        for index in selected:
            d = geometry_rows[index]
            fig, ax = plt.subplots(figsize=(8,6))
            ax.imshow(volume[index],cmap="gray",origin="upper",aspect="auto")
            width = d["x_right_edge_px"]-d["x_left_edge_px"]
            height = float(row["diameter_um"])/float(row["dz_um"])
            ax.add_patch(Ellipse((d["x4_px"],d["z_top_edge_px"]+height/2),width,height,fill=False,edgecolor="#00d8ff",linewidth=1.6))
            ax.add_patch(Rectangle((d["x_left_edge_px"],d["tail_start_edge_px"]),width,500/float(row["dz_um"]),fill=False,edgecolor="#ffbd3d",linewidth=1.6))
            ax.axvline(d["x4_px"],color="#00d8ff",linestyle="--",linewidth=.8,label="X4")
            ax.hlines([d["z_top_edge_px"],d["z_bottom_edge_px"]],d["x_left_edge_px"],d["x_right_edge_px"],colors="#48ee8b",label="physical top / bottom")
            ax.set(xlabel="A-line (0-based pixel centre)",ylabel="Depth pixel (0-based)",title=f"{scan} / B-scan {index}\nPhysical diameter {row['diameter_um']} um; X1 span {width:.0f} px; tail 500 um; guard 0")
            ax.legend(loc="upper right",fontsize=8)
            fig.tight_layout()
            target = RESULTS / "qc" / f"{scan}_frame_{index:03d}.png"
            target.parent.mkdir(exist_ok=True,parents=True)
            fig.savefig(target,dpi=120)
            plt.close(fig)
            qc_paths.append(target.relative_to(REPO).as_posix())
        summary["qc_frames"] = selected
        summary["qc_files"] = qc_paths
    return summary


def validate_mat(path, row, frame):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError("Missing or empty frame " + path.name)
    mat = loadmat(path, simplify_cells=True)
    sv = mat["sv_raw"]
    meta = mat["metadata"]
    expected_shape = (351, int(row["alines_per_bscan"]))
    if sv.shape != expected_shape:
        raise ValueError("sv_raw shape mismatch")
    source = str(Path(row["oct_absolute_path"]).resolve())
    if str(Path(meta["source_file"]).resolve()) != source:
        raise ValueError("Export source identity mismatch")
    if meta["bscan_index_matlab_1based"] != frame + 1 or path.name != f"frame_{frame:03d}.mat":
        raise ValueError("Export B-scan identity mismatch")
    if meta["formal_signal_name"] != "sv_raw" or meta["formal_signal_definition"] != "var(abs(E), 1, 3)" or meta["variance_denominator"] != "N":
        raise ValueError("Formal signal definition mismatch")
    if meta["dimension_order"] != "depth x A-line":
        raise ValueError("Dimension-order mismatch")
    reconstruction = meta["reconstruction"]
    if str(Path(reconstruction["source_file"]).resolve()) != source or reconstruction["bscan_index"] != frame+1:
        raise ValueError("Reconstruction identity mismatch")
    if reconstruction["n_bscan_positions"] != int(row["bscan_count"]):
        raise ValueError("Reconstruction volume count mismatch")
    header = reconstruction["header"]
    if (header["SPL"],header["nX"],header["frame_per_pos"]) != (int(row["spectral_samples"]),int(row["alines_per_bscan"]),int(row["repeat_count"])):
        raise ValueError("Reconstruction acquisition metadata mismatch")
    if not np.isfinite(sv).all():
        raise ValueError("Non-finite raw SV")
    return dict(scan_id=row["scan_id"], frame_index=frame, file_name=path.name, size_bytes=path.stat().st_size,
                sha256=sha256(path), shape=list(sv.shape),finite_fraction=float(np.isfinite(sv).mean()),min=float(sv.min()),max=float(sv.max()),
                source_identity=True,bscan_identity=True,formal_signal_definition=meta["formal_signal_definition"],variance_denominator="N",validation_status="passed")


def validate_volume(job, row):
    directory = Path(job["interim_dir"])
    expected = job["frame_indices_0based"]
    actual = sorted(p.name for p in directory.glob("frame_*.mat"))
    if actual != sorted(f"frame_{f:03d}.mat" for f in expected):
        raise ValueError("Missing or unexpected frame files")
    records = [validate_mat(directory / f"frame_{frame:03d}.mat",row,frame) for frame in expected]
    write_csv(RESULTS / "exports" / f"{row['scan_id']}_frames.csv", records)
    result = dict(scan_id=row["scan_id"],expected_frames=len(expected),validated_frames=len(records),status="passed",total_size_bytes=sum(r["size_bytes"] for r in records))
    write_json(RESULTS / "exports" / f"{row['scan_id']}_validation.json",result)
    return result


def export_volume(job, row, matlab):
    directory = Path(job["interim_dir"])
    directory.mkdir(parents=True,exist_ok=True)
    # Validate every resumed file before calling the frozen wrapper's identity check.
    for frame in job["frame_indices_0based"]:
        path = directory / f"frame_{frame:03d}.mat"
        if path.exists():
            validate_mat(path,row,frame)
    missing = sum(not (directory / f"frame_{f:03d}.mat").exists() for f in job["frame_indices_0based"])
    if shutil.disk_usage(directory).free < missing * 6_000_000 + 1_000_000_000:
        raise OSError("Insufficient disk space for all remaining frames and reserve")
    job_file = LOCAL / f"job_{row['scan_id']}.json"
    write_json(job_file,[job])
    quote = lambda p: str(p).replace("'","''")
    command = f"addpath('matlab'); export_sv_collection('{quote(job_file)}','{row['scan_id']}');"
    run([matlab,"-batch",command],LOGS / f"{row['scan_id']}_export.log")
    return validate_volume(job,row)


def pilot(matlab):
    rows,jobs = load_inputs()
    scan = read_json(LOCAL / "pilot_selection.json")["scan_id"]
    row,job = rows[scan],jobs[scan]
    assert row["pairing_status"] == "complete" and row["acquisition_metadata_status"] == "complete"
    summary = tracking(row,make_qc=True)
    frame = 249 if int(row["bscan_count"]) >= 250 else (int(row["bscan_count"])-1)//2
    path = Path(job["interim_dir"]) / f"frame_{frame:03d}.mat"
    if not path.exists():
        quote = lambda p: str(p).replace("'","''")
        run([matlab,"-batch",f"addpath('matlab'); export_sv_omag_frame('{quote(job['source_file'])}',{frame+1},'{quote(path)}',2);"],LOGS / "pilot_single_frame.log")
    single = validate_mat(path,row,frame)
    run([sys.executable,"-m","svrecttail","inspect-maps",str(path)],LOGS / "pilot_inspect_maps.json")
    result = dict(scan_id=scan,diameter_um=float(row["diameter_um"]),flow_speed_mm_s=float(row["flow_speed_mm_s"]),
                  hard_gates_passed=True,hard_gates=dict(unique_pair=True,metadata_complete=True,dicom_readable=True,frozen_tracking_ran=True,single_frame_reconstruction=True,export_identity_verified=True),
                  tracking=summary,single_frame=single,frame_selection_rule="249 if N>=250 else floor((N-1)/2)",full_volume_status="pending",
                  signal_definition="Var_t(abs(E)), population denominator N",background_subtraction=False,tail_length_um=500,guard_um=0,
                  validated_at_utc=datetime.now(timezone.utc).isoformat())
    write_json(RESULTS / "pilot_validation.json",result)
    write_csv(RESULTS / "pilot_tracking_summary.csv",[{k:v for k,v in summary.items() if not isinstance(v,(dict,list))}])
    return result


def batch(matlab, interim_root=None, pilot_only=False):
    rows,jobs = load_inputs()
    validation = read_json(RESULTS / "pilot_validation.json")
    assert validation["hard_gates_passed"] and all(validation["hard_gates"].values())
    scan = validation["scan_id"]
    if not pilot_only:
        assert validation["full_volume_status"] == "passed", "Pilot whole-volume validation required"
    selected = [scan] if pilot_only else sorted((s for s in jobs if s != scan),key=lambda s:(float(rows[s]["diameter_um"]),float(rows[s]["flow_speed_mm_s"]),rows[s]["oct_relative_path"]))
    status_path = RESULTS / "processing_summary.csv"
    status = {r["scan_id"]:r for r in csv.DictReader(status_path.open(encoding="utf-8"))} if status_path.exists() else {}
    for current in selected:
        row,job = rows[current],jobs[current]
        print("PROCESSING " + current,flush=True)
        record = dict(scan_id=current,diameter_um=row["diameter_um"],flow_speed_mm_s=row["flow_speed_mm_s"],tracking_status="pending",sv_export_status="pending",processing_status="running",exported_frames=0,error="")
        try:
            if interim_root and current != scan:
                candidate = str(Path(interim_root).resolve() / current)
                previous = Path(job["interim_dir"])
                if previous.exists() and any(previous.iterdir()) and str(previous.resolve()) != candidate:
                    raise ValueError("Existing output cannot be relocated implicitly")
                job["interim_dir"] = candidate
            summary = tracking(row)
            record["tracking_status"] = "passed"
            record["geometry_valid_frames"] = summary["geometry_qc_valid"]
            record["z_valid_frames"] = summary["z_localization_valid"]
            write_json(LOCAL / "export_jobs.json",list(jobs.values()))
            exported = export_volume(job,row,matlab)
            record.update(sv_export_status="passed",processing_status="complete",exported_frames=exported["validated_frames"])
            if current == scan:
                validation.update(full_volume_status="passed",full_volume_export=exported)
                write_json(RESULTS / "pilot_validation.json",validation)
        except Exception as error:
            record.update(processing_status="failed",sv_export_status="failed",error=f"{type(error).__name__}: {error}")
            # Full exception/paths remain local; public error messages are sanitized below.
            write_json(LOGS / f"{current}_error.json",record)
            record["error"] = record["error"].replace(str(REPO),"REPOSITORY").replace(row["oct_absolute_path"],row["oct_relative_path"]).replace(row["flow_dicom_absolute_path"],row["flow_dicom_relative_path"])
            if interim_root:
                record["error"] = record["error"].replace(str(Path(interim_root).resolve()),"INTERIM_ROOT")
        record["completed_at_utc"] = datetime.now(timezone.utc).isoformat()
        status[current] = record
        write_csv(status_path,list(status.values()))
        write_csv(LOGS / "batch_status.csv",list(status.values()))
        write_csv(LOGS / "errors.csv",[r for r in status.values() if r["processing_status"] == "failed"],fields=["scan_id","error"] if not any(r["processing_status"] == "failed" for r in status.values()) else None)
        print(json.dumps(record),flush=True)
    summaries = [read_json(p) for p in sorted((RESULTS / "tracking").glob("*_summary.json"))]
    write_csv(RESULTS / "per_volume_tracking_summary.csv",[{k:v for k,v in s.items() if not isinstance(v,(list,dict))} for s in summaries])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage",choices=["pilot","pilot-full","batch"])
    parser.add_argument("--matlab",required=True)
    parser.add_argument("--interim-root",type=Path)
    args = parser.parse_args()
    if args.stage == "pilot":
        print(json.dumps(pilot(args.matlab)))
    else:
        batch(args.matlab,args.interim_root,pilot_only=args.stage == "pilot-full")
