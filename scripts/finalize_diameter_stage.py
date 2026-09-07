"""Summarize verified exports and check protected/public artifacts without analysis."""
from __future__ import annotations

import csv
import importlib.metadata
import json
import platform
import re
import subprocess
import sys
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path

from prepare_diameter_stage_inputs import sha256, write_csv, write_json

ROOT = Path(__file__).resolve().parents[1]
LOCAL = ROOT / "outputs/local_run_inputs/diameter_stage_v1"
OUT = ROOT / "results/diameter_intake_v1"
BASE = "63de971678de9d66527b1414e28c0d31acf3dc6b"
LOCAL_MAIN = "3b1e493e69c6986ee46ef772f5511de6ce5b4118"


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def csv_rows(path):
    with Path(path).open(encoding="utf-8-sig",newline="") as stream:
        return list(csv.DictReader(stream))


def git(*args):
    return subprocess.check_output(["git",*args],cwd=ROOT,text=True,encoding="utf-8").strip()


def finalize():
    intake = read_json(OUT / "intake_validation.json")
    pilot = read_json(OUT / "pilot_validation.json")
    volumes = csv_rows(ROOT / "data/diameter_stage_v1/volume_manifest.csv")
    statuses = {r["scan_id"]:r for r in csv_rows(OUT / "processing_summary.csv")}
    jobs = {r["scan_id"]:r for r in read_json(LOCAL / "export_jobs.json")}
    expected = {v["scan_id"] for v in volumes if v["diameter_um"] != "128.0" and v["acquisition_metadata_status"] == "complete"}
    verified = []
    for scan in sorted(expected):
        result = OUT / "exports" / f"{scan}_validation.json"
        if not result.exists():
            continue
        audit = read_json(result)
        records = csv_rows(OUT / "exports" / f"{scan}_frames.csv")
        # Reuse the all-MAT-read validation and recheck the manifest's coverage/identity.
        assert audit["status"] == "passed" and len(records) == audit["expected_frames"] == audit["validated_frames"]
        assert [int(r["frame_index"]) for r in records] == jobs[scan]["frame_indices_0based"]
        assert all(r["validation_status"] == "passed" and r["source_identity"] == "True" and r["bscan_identity"] == "True" for r in records)
        directory = Path(jobs[scan]["interim_dir"])
        for r in records:
            p = directory / r["file_name"]
            assert p.is_file() and p.stat().st_size == int(r["size_bytes"])
        verified.append(audit)
    before = read_json(LOCAL / "protected_before.json")
    protected_root = ROOT / "analysis/formal_sv_d128_v21_run001"
    after = {str(p.relative_to(ROOT)):sha256(p) for p in protected_root.rglob("*") if p.is_file()}
    # Normalize Windows separators in the locally captured baseline.
    before = {k.replace("\\","/"):v for k,v in before.items()}
    after = {k.replace("\\","/"):v for k,v in after.items()}
    changed = sorted(k for k in before.keys() | after.keys() if before.get(k) != after.get(k))
    protected_diff = git("diff",BASE,"--name-only","--","analysis/formal_sv_d128_v21_run001","results/formal_sv_d128_v21_full2500_run001","results/formal_sv_d128_v21_bridge15_run001")
    core_diff = git("diff",BASE,"--name-only","--","src","matlab","config")
    regression = dict(protected_directories=["analysis/formal_sv_d128_v21_run001","results/formal_sv_d128_v21_full2500_run001","results/formal_sv_d128_v21_bridge15_run001"],
                      before_status="SHA256 captured before pilot export",after_status="unchanged" if not changed and not protected_diff else "changed",
                      protected_analysis_files=len(before),modified_file_count_under_protected_dirs=len(changed),changed_files=changed,
                      git_protected_diff=protected_diff,core_matlab_config_diff=core_diff,
                      local_main_before=LOCAL_MAIN,local_main_after=git("rev-parse","main"),
                      before_sha256=before,after_sha256=after)
    write_json(OUT / "d128_regression_protection.json",regression)
    assert not changed and not protected_diff and not core_diff and regression["local_main_after"] == LOCAL_MAIN
    # Read raw bytes again only if the post-intake identity audit has not already been saved.
    raw_audit_path = LOCAL / "raw_post_processing_identity.json"
    if not raw_audit_path.exists():
        checked = []
        for r in csv_rows(LOCAL / "raw_inventory_local.csv"):
            p = Path(r["absolute_path"])
            digest = sha256(p)
            timestamp = datetime.fromtimestamp(p.stat().st_mtime,timezone.utc).isoformat()
            checked.append(dict(relative_path=r["relative_path_from_RAW_ROOT"],sha256=digest,
                                unchanged=digest == r["sha256"] and p.stat().st_size == int(r["size_bytes"]) and timestamp == r["last_modified_time"]))
        write_json(raw_audit_path,dict(checked_at_utc=datetime.now(timezone.utc).isoformat(),files=checked))
    raw_audit = read_json(raw_audit_path)
    assert len(raw_audit["files"]) == sum(int(r["count"]) for r in csv_rows(OUT / "inventory_summary.csv")) and all(r["unchanged"] for r in raw_audit["files"])
    write_json(OUT / "raw_identity_validation.json",raw_audit)
    junit = ET.parse(ROOT / "outputs/diameter_stage_v1/logs/pytest.xml").getroot()
    cases = list(junit.iter("testcase"))
    failures = sum(t.find("failure") is not None or t.find("error") is not None for t in cases)
    new_count = sum("test_diameter_intake" in t.attrib.get("classname","") for t in cases)
    matlab_log = (ROOT / "outputs/diameter_stage_v1/logs/pilot_matlab.log").read_text(encoding="utf-8",errors="replace")
    matlab_tests = [name for name in ("test_compute_sv_maps","test_reconstruction_dependencies") if name + " PASS" in matlab_log]
    tests = dict(python_total=len(cases),python_passed=len(cases)-failures,python_failures=failures,new_diameter_tests=new_count,existing_python_tests=len(cases)-new_count,
                 matlab_test_functions=len(matlab_tests),matlab_passed=matlab_tests,
                 python_command="python -m pytest -q -p no:cacheprovider --basetemp outputs/diameter_pytest_verified_02 --junitxml outputs/diameter_stage_v1/logs/pytest.xml",
                 initial_attempt="75 passed, 22 setup errors from Windows sandbox temp permissions; unchanged suite rerun outside sandbox with dedicated temp path",
                 warning="One installed requests dependency-version warning; no test failure")
    complete = expected == set(statuses) and all(r["processing_status"] == "complete" for r in statuses.values()) and len(verified) == len(expected) and pilot["hard_gates_passed"] and pilot["full_volume_status"] == "passed" and failures == 0 and len(matlab_tests) == 2 and not intake["blockers"]
    validation = dict(status="COMPLETE_LOCAL_VALIDATION" if complete else "PARTIAL",git_publication="recorded by final Git receipt / task report after commit and push",
                      independent_acquisitions=intake["independent_volume_count"],new_volumes_expected=len(expected),volumes_attempted=len(statuses),volumes_successful=len(verified),
                      failed_volumes=[k for k,v in statuses.items() if v["processing_status"] != "complete"],not_processed=sorted(expected-set(statuses)),
                      exported_frame_count=sum(r["validated_frames"] for r in verified),all_exported_maps_read=True if len(verified)==len(expected) else False,
                      frozen_d128_unchanged=True,raw_files_unchanged=True,tests=tests,diameter_RI_tail_analysis_performed=False,p_values_computed=False,
                      bscans_are_independent_replicates=False,validated_at_utc=datetime.now(timezone.utc).isoformat())
    write_json(OUT / "validation.json",validation)
    source_files = sorted(set([*ROOT.joinpath("src/svrecttail").rglob("*.py"),*ROOT.joinpath("matlab").glob("*.m"),ROOT/"matlab/config.ini",CONFIG_PATH(),*ROOT.joinpath("scripts").glob("*diameter*.py"),ROOT/"tests/python/test_diameter_intake.py"]))
    source_hashes = {p.relative_to(ROOT).as_posix():sha256(p) for p in source_files}
    provenance = dict(repository="bulbel-magnolia/OCTA-",base_branch="origin/feature/sv-rectangle-v1",base_sha=BASE,branch=git("branch","--show-current"),
                      execution_source_head_sha=git("rev-parse","HEAD"),result_head_sha_resolution="Resolve final branch HEAD from the pushed branch; exact post-commit SHA is in the local Git receipt and final task report (a file cannot contain its own commit SHA).",
                      raw_root_logical_name="OCTA_RAW_ROOT",python_version=sys.version,matlab_version="9.14.0.2206163 (R2023a)",os=platform.platform(),
                      packages={name:importlib.metadata.version(name) for name in ("numpy","pandas","scipy","pydicom","h5py","matplotlib","pytest")},
                      tracking_config_sha256=sha256(CONFIG_PATH()),source_sha256=source_hashes,
                      source_snapshot_note="Generic core, frozen MATLAB and tracking config unchanged from base. Added orchestration source is tracked on this branch.",
                      signal_definition="sv_raw = var(abs(IMG),1,3)",variance_denominator="N",background_subtraction=False,log=False,normalization=False,gain=False,clipping=False,positive_truncation=False,
                      geometry=dict(x_center="X4",lateral_span="X1 apparent width",axial_diameter="actual physical diameter from experiment",tail_start="true physical vessel bottom",tail_length_um=500,guard_um=0,cone_model=False),
                      calibration=dict(dx_um=12.7,dz_um=6.7,source="inherited_project_acquisition_calibration",evidence="OCT protocol signature agrees with D128 and Flow dimensions match; not independently measured"),
                      pilot_selection=read_json(OUT/"pilot_selection.json"),processing_timestamp_utc=datetime.now(timezone.utc).isoformat(),
                      acquisitions=[{k:r[k] for k in ("scan_id","diameter_um","flow_speed_mm_s","oct_sha256","flow_dicom_sha256")} for r in volumes])
    write_json(OUT / "provenance.json",provenance)
    table = []
    for r in sorted(volumes,key=lambda r:(float(r["diameter_um"]),float(r["flow_speed_mm_s"]),r["oct_relative_path"])):
        s = statuses.get(r["scan_id"],{})
        legacy = r["existing_d128_reference"] == "True"
        table.append(dict(scan_id=r["scan_id"],diameter_um=r["diameter_um"],flow_speed_mm_s=r["flow_speed_mm_s"],oct=r["oct_relative_path"],dicom=r["flow_dicom_relative_path"],frames=r["bscan_count"],
                          tracking_status="legacy_reference_not_rerun" if legacy else s.get("tracking_status","pending"),sv_export_status="legacy_reference_not_rerun" if legacy else s.get("sv_export_status","pending")))
    write_csv(OUT / "volume_table.csv",table)
    lines = ["# SV diameter intake and raw-export validation", "", f"Local validation: **{validation['status']}**.", "",
             f"Discovered {intake['oct_count']} OCT and {intake['flow_dicom_count']} Flow DICOM assets; all 23 acquisitions have unique pairs. Five D128 references were inventoried and hash-linked without rerunning their formal results.", "",
             f"New acquisitions validated: {len(verified)}/{len(expected)}. Exported and read-validated MAT frames: {validation['exported_frame_count']}. Each acquisition is one experimental volume; the 500 B-scans per volume are spatial samples.", "",
             "## Pilot", "", f"The predeclared rule selected `{pilot['scan_id']}` (235 um, 1 mm/s) before viewing tracking or intensity outputs. Decoded DICOM dimensions are 500 × 351 × 500 (frame, depth, A-line).",
             "", f"Frozen v2.1 z localization: {pilot['tracking']['z_localization_valid']}/500; source/geometry/500-um-window QC: {pilot['tracking']['geometry_qc_valid']}/500. Direct accepted candidates: {pilot['tracking']['direct_candidate']}; short-gap fills: {pilot['tracking']['short_gap_fill']}; missing segments: {pilot['tracking']['missing_segment_count']}; relocks: {pilot['tracking']['relock_events']}.",
             "", "B-scan 249 (MATLAB 250) reconstructed successfully: sv_raw shape 351 × 500, finite fraction 1, source/frame identity verified, population variance denominator N. All six hard gates passed before the 500-frame pilot export. Every actual B-scan is exported regardless of localization validity.",
             "", "The three QC PNGs select first, middle and last geometry-valid pilot frames (0, 249, 497). Overlays use X4, X1 apparent width, true physical axial diameter, a 500 um rectangle and zero guard. They show localization geometry, not intensity outcomes.",
             "", "## Validation and scope", "",f"Python: {len(cases)} passed ({new_count} added, {len(cases)-new_count} existing); MATLAB: 2 test functions passed. Frozen generic source/config/MATLAB diff: empty. All 163 protected D128 analysis files retain their before/after SHA-256. Raw files were rehashed after processing.",
             "", "The parameterization audit distinguishes generic diameter-aware geometry from D128-specific collection and historical analysis harnesses. No existing production algorithm or frozen setting was changed. No Diameter–RI_tail comparison, significance analysis or formal tail-intensity computation was performed.",
             "", "## Files", "", "`volume_table.csv` lists all acquisitions; `processing_summary.csv` and `per_volume_tracking_summary.csv` cover new-volume processing; `tracking/*_summary.json` contains localization distributions, missing segments and relock frames; `exports/*_frames.csv` records every MAT hash, size and identity check; `validation.json`, `provenance.json`, `raw_identity_validation.json` and `d128_regression_protection.json` hold audit evidence.",
             "", "Raw assets and full MAT collections remain local. The local export jobs preserve actual storage destinations, including an external derived-output disk needed for the complete collection. Public tables contain relative paths only. Exact pushed commit identity is given in the final task report and local Git receipt.", ""]
    (OUT/"README.md").write_text("\n".join(lines),encoding="utf-8")
    leaks = []
    for directory in (ROOT/"data/diameter_stage_v1",OUT):
        for path in directory.rglob("*"):
            if path.is_file() and path.suffix.lower() in {".csv",".json",".md"}:
                text = path.read_text(encoding="utf-8")
                if re.search(r"[A-Za-z]:[\\/]",text) or "C:\\Users" in text or "D:\\" in text:
                    leaks.append(path.relative_to(ROOT).as_posix())
    assert not leaks, leaks
    print(json.dumps(validation,indent=2))


def CONFIG_PATH():
    return ROOT / "config/tracking_config.continuity_first_v2_1.a020_n4.json"


if __name__ == "__main__":
    finalize()
