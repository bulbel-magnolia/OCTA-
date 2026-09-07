"""Identity, conservative pairing and physical-diameter regression protection."""
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from scipy.io import savemat

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from prepare_diameter_stage_inputs import assign_ids, pair_oct, parse_condition, relative, sha256, write_csv, preserve_job_destinations
from process_diameter_stage import validate_mat, validate_volume
from svrecttail.geometry import VesselGeometry, ellipse_weights
from svrecttail.mentor.xroi import local_body_features
from svrecttail.mentor_tracking import build_localization_from_tracking
from svrecttail.pipeline import _check_calibration


def test_exact_companion_wins():
    assert pair_oct("64/5.oct", ["64/5.oct-Flow_ed.dcm", "other/5.oct-Flow_ed.dcm"]) == ("64/5.oct-Flow_ed.dcm", "exact_companion", "complete")


def test_missing_companion():
    assert pair_oct("64/5.oct", ["64/5.oct-Stru.dcm"])[2] == "missing"


def test_ambiguous_companion():
    assert pair_oct("64/5.oct", ["64/5.oct-Flow_ed(1).dcm", "64/5.oct-Flow_ed(2).dcm"])[2] == "ambiguous"


def test_unique_suffix_companion():
    assert pair_oct("285/3.oct", ["285/3.oct-Flow_ed(1).dcm"])[2] == "complete"


def test_conflicting_exact_names():
    assert pair_oct("64/5.oct", ["64/5.oct-Flow_ed.dcm", "64/5.OCT-FLOW_ED.DCM"])[2] == "ambiguous"


def test_windows_relative_path():
    assert pair_oct(r"64\5.oct", [r"64\5.oct-Flow_ed.dcm"])[0] == "64/5.oct-Flow_ed.dcm"


@pytest.mark.parametrize("value", [r"D:\data\5.oct", "/data/5.oct", "../5.oct", r"\\host\share\5.oct"])
def test_absolute_paths_rejected(value):
    with pytest.raises(ValueError):
        relative(value)


@pytest.mark.parametrize("path,diameter,flow", [("64/5.oct",64,5),("0.5mm/12.oct",500,12),("256/3.oct",256,3)])
def test_documented_conditions(path,diameter,flow):
    result = parse_condition(path,naming_documented=True)
    assert result["diameter_um"] == diameter and result["flow_speed_mm_s"] == flow


def test_unknown_diameter():
    assert parse_condition("unknown/5.oct",naming_documented=True)["diameter_um"] is None


def test_unknown_flow():
    assert parse_condition("64/scan5.oct",naming_documented=True)["flow_speed_mm_s"] is None


def test_no_undocumented_bare_numeric_inference():
    result = parse_condition("64/5.oct")
    assert result["diameter_um"] is None and result["flow_speed_mm_s"] is None


def test_conflicting_diameter_ancestors():
    assert parse_condition("64/128/5.oct",naming_documented=True)["diameter_um"] is None


def rows():
    return [dict(oct_relative_path=f"64/{p}/5.oct",diameter_um=64.,flow_speed_mm_s=5.,pairing_status="complete",oct_sha256=p,flow_dicom_sha256=p) for p in ("b","a")]


def test_duplicate_condition_volumes_lexical_and_stable():
    first = assign_ids(rows())
    assert {r["oct_relative_path"]:r["scan_id"] for r in first} == {"64/a/5.oct":"D064_F05_V01","64/b/5.oct":"D064_F05_V02"}
    assert sorted(r["scan_id"] for r in assign_ids(list(reversed(rows())),first)) == ["D064_F05_V01","D064_F05_V02"]


def test_prior_identity_never_silently_changes():
    first = assign_ids(rows())
    changed = rows()
    changed[0]["oct_sha256"] = "changed"
    with pytest.raises(ValueError):
        assign_ids(changed,first)


def test_unknown_condition_has_no_invented_scan_id():
    data = rows()
    data[0]["diameter_um"] = None
    assert assign_ids(data)[0]["scan_id"] == ""


def test_intake_rerun_retains_external_output_destination(tmp_path):
    job = dict(scan_id="D064_F05_V01",source_file=str(tmp_path/"scan.oct"),frame_indices_0based=[0,1],interim_dir="default")
    previous = dict(job,interim_dir="previous-derived-disk")
    assert preserve_job_destinations([job],[previous])[0]["interim_dir"] == "previous-derived-disk"


def test_changed_export_job_refused(tmp_path):
    job = dict(scan_id="D064_F05_V01",source_file=str(tmp_path/"scan.oct"),frame_indices_0based=[0,1],interim_dir="default")
    with pytest.raises(ValueError):
        preserve_job_destinations([job],[dict(job,frame_indices_0based=[0])])


def test_sha256_reads_bytes(tmp_path):
    path = tmp_path / "raw.oct"
    path.write_bytes(b"abc")
    assert sha256(path) == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"


def test_public_manifest_no_absolute_path(tmp_path):
    row = dict(scan_id="D064_F05_V01",oct_relative_path="64/5.oct",oct_absolute_path=r"D:\private\64\5.oct")
    target = tmp_path / "manifest.csv"
    write_csv(target,[{k:v for k,v in row.items() if "absolute_path" not in k}])
    content = target.read_text()
    assert "D:" not in content and "absolute_path" not in content and "64/5.oct" in content


def record(diameter):
    return pd.Series(dict(scan_id="synthetic",frame_index=0,alpha=.2,diameter_um=diameter,x_center_px=50,z_upper_px=30,
        x1_local_geometry_px=49,x2_robust_centroid_px=50,x4_centroid_isolated_jump_corrected_px=50,
        local_body_run_width_px=11,expected_lateral_width_px=max(3,round(diameter/12.7)),valid_local_body=True,
        x1_fallback=False,x4_jump_corrected=False,new_tracking_class="high_confidence",vessel_presence_prediction="assessable",
        local_body_background=2,local_body_sigma=1,local_body_peak_cnr=5,local_body_axial_completeness=.8,assessability_score=.9))


@pytest.mark.parametrize("diameter", [64,128,235,256,285,500])
def test_physical_geometry_distinct_from_x1(diameter):
    loc = build_localization_from_tracking(record(diameter),diameter_um=diameter,dx_um=12.7,dz_um=6.7)
    g = loc.geometry
    assert g.diameter_um == diameter
    assert g.lateral_width_um == pytest.approx(11*12.7)
    assert g.lateral_width_um != diameter
    assert g.z_bottom_edge_px == pytest.approx(29.5 + diameter/6.7)
    assert g.x_center_px == 50
    assert loc.source_qc_valid
    assert ellipse_weights((200,100),g).sum() > 0


@pytest.mark.parametrize("diameter", [64,128,256])
def test_local_body_expected_width_and_axial_slab_follow_diameter(diameter):
    image = np.zeros((200,100))
    image[30:100,45:56] = 10
    features = local_body_features(image,pd.Series(dict(x_center_px=50,z_upper_px=30,diameter_um=diameter)))
    assert features["expected_lateral_width_px"] == max(3,round(diameter/12.7))
    assert features["diameter_axial_px"] == round(diameter/6.7)


def test_d128_geometry_matches_preexisting_equations():
    g = build_localization_from_tracking(record(128),diameter_um=128,dx_um=12.7,dz_um=6.7).geometry
    assert g == VesselGeometry(44.5,55.5,29.5,128,12.7,6.7)
    assert g.z_bottom_edge_px == 29.5+128/6.7


@pytest.mark.parametrize("diameter", [64,128,256])
def test_generic_run_accepts_matching_non128_calibration(diameter):
    config = SimpleNamespace(calibration=SimpleNamespace(diameter_um=diameter,dx_um=12.7,dz_um=6.7))
    _check_calibration(pd.Series(dict(scan_id="synthetic",diameter_um=diameter,dx_um=12.7,dz_um=6.7)),config)


def synthetic_mat(tmp_path):
    source = tmp_path / "scan.oct"
    row = dict(scan_id="D064_F05_V01",oct_absolute_path=str(source),alines_per_bscan="3",bscan_count="2",spectral_samples="1024",repeat_count="3")
    meta = dict(source_file=str(source),bscan_index_matlab_1based=1,formal_signal_name="sv_raw",formal_signal_definition="var(abs(E), 1, 3)",variance_denominator="N",dimension_order="depth x A-line",
                reconstruction=dict(source_file=str(source),bscan_index=1,n_bscan_positions=2,header=dict(SPL=1024,nX=3,frame_per_pos=3)))
    path = tmp_path / "frame_000.mat"
    savemat(path,dict(sv_raw=np.ones((351,3)),metadata=meta))
    return path,row,meta


def test_validated_mat_identity(tmp_path):
    path,row,_ = synthetic_mat(tmp_path)
    assert validate_mat(path,row,0)["validation_status"] == "passed"


@pytest.mark.parametrize("key,value", [("source_file","other.oct"),("variance_denominator","N-1"),("bscan_index_matlab_1based",2),("formal_signal_name","sv_cv2")])
def test_invalid_resume_metadata_refused(tmp_path,key,value):
    path,row,meta = synthetic_mat(tmp_path)
    meta[key] = value
    savemat(path,dict(sv_raw=np.ones((351,3)),metadata=meta))
    with pytest.raises(ValueError):
        validate_mat(path,row,0)


def test_silent_missing_frame_refused(tmp_path):
    path,row,_ = synthetic_mat(tmp_path)
    with pytest.raises(ValueError,match="Missing or unexpected"):
        validate_volume(dict(interim_dir=str(tmp_path),frame_indices_0based=[0,1]),row)
