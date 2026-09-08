"""Analytical kernel checks plus frozen-output audit gates."""
import gzip
import importlib.util
import json
import re
from types import SimpleNamespace
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.io import savemat

from svrecttail.geometry import VesselGeometry, ellipse_weights, interval_overlap_weights
from svrecttail.raw_sv_quantification import quantify_raw_sv, WINDOWS_UM, METRICS

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis/formal_sv_diameter_v1"
SPEC = importlib.util.spec_from_file_location("formal_driver", ROOT / "scripts/quantify_sv_diameter_formal_v1.py")
DRIVER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(DRIVER)


@pytest.fixture
def geometry():
    return VesselGeometry(10.15, 23.85, 16.27, 128., 12.7, 6.7)


def test_constant_raw_source_ellipse(geometry):
    m, _ = quantify_raw_sv(np.full((160,40), 17.), geometry)
    assert m["source_mean_raw"] == pytest.approx(17.)
    w = ellipse_weights((160,40),geometry,supersample=16)
    assert m["source_area_um2"] == pytest.approx(w.sum()*12.7*6.7)


@pytest.mark.parametrize("length", WINDOWS_UM)
def test_rectangular_windows_and_area(geometry,length):
    m, _ = quantify_raw_sv(np.full((160,40), 17.), geometry)
    assert m[f"tail_mean_raw_{length}um"] == pytest.approx(17.)
    assert m[f"tail_area_um2_{length}um"] == pytest.approx(geometry.lateral_width_um*length)
    assert m[f"tail_q_raw_{length}um"] == pytest.approx(17.*geometry.lateral_width_um*length)


def test_no_background_subtraction_additive_offset(geometry):
    image = np.repeat((np.arange(160)+1.)[:,None],40,axis=1)
    before,_=quantify_raw_sv(image,geometry)
    after,_=quantify_raw_sv(image+123.,geometry)
    for key in ("source_mean_raw","tail_mean_raw_500um"):
        assert after[key]-before[key] == pytest.approx(123.)
    assert after["ri_tail"] != pytest.approx(before["ri_tail"])


def test_no_gain_normalization_or_log(geometry):
    image=np.repeat((np.arange(160)+1.)[:,None],40,axis=1)
    a,_=quantify_raw_sv(image,geometry)
    b,_=quantify_raw_sv(image*7.,geometry)
    assert b["source_mean_raw"]==pytest.approx(a["source_mean_raw"]*7)
    assert b["tail_mean_raw_500um"]==pytest.approx(a["tail_mean_raw_500um"]*7)
    assert b["ri_tail"]==pytest.approx(a["ri_tail"])


def test_physical_diameter_axial_height_independent_of_x1(geometry):
    wider=VesselGeometry(geometry.x_left_edge_px,geometry.x_right_edge_px,
                         geometry.z_top_edge_px,235.,12.7,6.7)
    assert wider.lateral_width_um==geometry.lateral_width_um
    assert wider.z_bottom_edge_px-geometry.z_bottom_edge_px==pytest.approx((235-128)/6.7)
    assert ellipse_weights((160,40),wider).sum()>ellipse_weights((160,40),geometry).sum()


def test_fractional_bottom_is_not_rounded(geometry):
    assert geometry.z_bottom_edge_px==geometry.z_top_edge_px+128/6.7
    assert geometry.z_bottom_edge_px != round(geometry.z_bottom_edge_px)


def test_fractional_overlap_against_explicit_piecewise_integral(geometry):
    image=np.repeat(np.arange(160,dtype=float)[:,None],40,axis=1)+10.
    m,_=quantify_raw_sv(image,geometry)
    top=geometry.z_bottom_edge_px
    end=top+100/6.7
    integral=sum((z+10.)*max(0.,min(z+.5,end)-max(z-.5,top)) for z in range(160))
    assert m["tail_mean_raw_100um"]==pytest.approx(integral/(end-top))


def test_ratio_and_primary_alias_exact(geometry):
    m,_=quantify_raw_sv(np.arange(6400,dtype=float).reshape(160,40)+1,geometry)
    assert m["ri_tail"]==m["ri_tail_500um"]
    assert m["ri_tail"]==m["tail_mean_raw_500um"]/m["source_mean_raw"]


def test_native_anchors_no_interpolation(geometry):
    image=np.repeat((np.arange(160,dtype=float)**2+1)[:,None],40,axis=1)
    m,anchors=quantify_raw_sv(image,geometry)
    assert len(anchors)==51
    assert [a["target_r_um"] for a in anchors]==list(range(0,501,10))
    zw=interval_overlap_weights(160,geometry.z_bottom_edge_px,geometry.z_bottom_edge_px+500/6.7)
    candidates=np.flatnonzero(zw>0)
    for a in anchors:
        z=a["selected_z_index_0based"]
        distances={i:abs((i-geometry.z_bottom_edge_px)*6.7-a["target_r_um"]) for i in candidates}
        assert z==min(distances,key=distances.get)
        assert a["S_tail_raw_r"]==pytest.approx(z*z+1)
        assert a["RI_r"]==a["S_tail_raw_r"]/m["source_mean_raw"]
        assert a["tail_z_fraction"]==zw[z]
        assert a["abs_sampling_error_um"]<=6.7/2+1e-9


@pytest.mark.parametrize("top", [16.,16.1,16.5,16.99])
def test_depth_half_pixel_tolerance_for_fractional_tops(top):
    g=VesselGeometry(10.,20.,top,128.,12.7,6.7)
    _,anchors=quantify_raw_sv(np.ones((160,40)),g)
    assert max(a["abs_sampling_error_um"] for a in anchors)<=6.7/2+1e-9


def test_invalid_geometry_na_no_zero():
    result,anchors=quantify_raw_sv(None,None,geometry_qc_valid=False)
    assert set(result)==set(METRICS)
    assert all(np.isnan(v) for v in result.values())
    assert anchors==[]


def test_marked_valid_incomplete_geometry_hard_failure():
    g=VesselGeometry(10.,20.,130.,128.,12.7,6.7)
    with pytest.raises(ValueError,match="full source and tail"):
        quantify_raw_sv(np.ones((160,40)),g)


def test_no_positive_truncation_or_ratio_clipping(geometry):
    image=np.full((160,40),20.)
    image[40:,:]=-7.
    result,_=quantify_raw_sv(image,geometry)
    assert result["tail_mean_raw_500um"]<0 and result["ri_tail"]<0


@pytest.mark.parametrize("value",[np.nan,np.inf])
def test_nonfinite_input_hard_failure(geometry,value):
    image=np.ones((160,40)); image[0,0]=value
    with pytest.raises(ValueError,match="finite"):
        quantify_raw_sv(image,geometry)


def test_nonpositive_source_hard_failure(geometry):
    with pytest.raises(ValueError,match="positive"):
        quantify_raw_sv(np.zeros((160,40)),geometry)


def test_direct_vs_profile_nonseparable_image(geometry):
    rng=np.random.default_rng(512)
    image=rng.lognormal(5,2,(160,40))
    result,_=quantify_raw_sv(image,geometry)
    assert result["direct_profile_max_relative_error"]<1e-12
    xw=interval_overlap_weights(40,geometry.x_left_edge_px,geometry.x_right_edge_px)
    zw=interval_overlap_weights(160,geometry.z_bottom_edge_px,geometry.z_bottom_edge_px+500/6.7)
    independent=float(zw @ image @ xw)*12.7*6.7
    assert result["tail_q_raw_500um"]==pytest.approx(independent,rel=1e-12)


def test_all_9000_coverage_and_invalid_na():
    if not (OUT/"framewise_all.csv.gz").exists():
        pytest.skip("Generated integration outputs not yet available")
    fw=pd.read_csv(OUT/"framewise_all.csv.gz",float_precision="round_trip")
    depth=pd.read_csv(OUT/"depth_anchor_framewise.csv.gz",float_precision="round_trip")
    expected=pd.read_csv(ROOT/"results/diameter_intake_v1/processing_summary.csv").set_index("scan_id")
    DRIVER.validate_tables(fw,depth,expected)
    assert fw.input_identity_valid.all() and fw.input_mat_sha256.str.fullmatch(r"[a-f0-9]{64}").all()
    audit=pd.concat([pd.read_csv(p).assign(frame_index_0based=lambda x:x.frame_index) for p in (ROOT/"results/diameter_intake_v1/exports").glob("*_frames.csv")])
    joined=fw.merge(audit[["scan_id","frame_index_0based","sha256"]],on=["scan_id","frame_index_0based"],validate="one_to_one")
    assert len(joined)==9000 and joined.input_mat_sha256.eq(joined.sha256).all()


def test_public_outputs_no_absolute_paths():
    DRIVER.public_path_audit()


@pytest.mark.parametrize("text,expected", [
    ("https://github.com/example/repo",False),
    ("C:/Users/example/file.csv",True),
    (r"D:\raw\scan.mat",True),
    (r"\\server\share\file.csv",True),
    ("relative/scan.csv",False),
])
def test_path_audit_distinguishes_urls_from_windows_paths(text,expected):
    assert DRIVER.has_absolute_windows_path(text)==expected


def test_committed_conditions_and_native_depth_geometry():
    if not (OUT/"framewise_primary.csv").exists():
        pytest.skip("Generated integration outputs not yet available")
    fw=pd.read_csv(OUT/"framewise_primary.csv",float_precision="round_trip")
    manifest=pd.read_csv(ROOT/"data/diameter_stage_v1/volume_manifest.csv")
    conditions=fw.merge(manifest[["scan_id","diameter_um","flow_speed_mm_s"]],on="scan_id",validate="many_to_one",suffixes=("","_manifest"))
    assert conditions.diameter_um.eq(conditions.diameter_um_manifest).all()
    assert conditions.flow_mm_s.eq(conditions.flow_speed_mm_s).all()
    assert not fw.scan_id.str.startswith("D128_").any()
    depth=pd.read_csv(OUT/"depth_anchor_framewise.csv.gz",float_precision="round_trip")
    assert not depth.duplicated(["scan_id","frame_index_0based","target_r_um"]).any()
    joined=depth.merge(fw[["scan_id","frame_index_0based","source_mean_raw","z_bottom_edge_px"]],on=["scan_id","frame_index_0based"],validate="many_to_one")
    assert np.array_equal(joined.S_vessel_raw,joined.source_mean_raw)
    assert np.array_equal(joined.RI_r,joined.S_tail_raw_r/joined.S_vessel_raw)
    expected_r=(joined.selected_z_index_0based-joined.z_bottom_edge_px)*6.7
    np.testing.assert_allclose(joined.selected_r_um,expected_r,rtol=0,atol=1e-12)
    np.testing.assert_allclose(joined.abs_sampling_error_um,abs(expected_r-joined.target_r_um),rtol=0,atol=1e-12)
    assert joined.tail_z_fraction.between(0,1,inclusive="right").all()


def test_volume_descriptive_aggregation():
    if not (OUT/"volume_summary.csv").exists():
        pytest.skip("Generated integration outputs not yet available")
    fw=pd.read_csv(OUT/"framewise_primary.csv",float_precision="round_trip")
    summary=pd.read_csv(OUT/"volume_summary.csv",float_precision="round_trip").set_index("scan_id")
    assert len(summary)==18
    for scan,group in fw.groupby("scan_id"):
        for metric in ("source_mean_raw","tail_mean_raw_500um","ri_tail"):
            row=summary.loc[scan]
            np.testing.assert_allclose([row[f"{metric}_{k}"] for k in ("q1","median","q3")],np.quantile(group[metric],[.25,.5,.75]),rtol=1e-14)
            assert row[f"{metric}_mean"]==pytest.approx(group[metric].mean(),rel=1e-14)
            assert row[f"{metric}_sd"]==pytest.approx(group[metric].std(ddof=1),rel=1e-14)


@pytest.mark.parametrize("corruption", ["sha","source","frame","signal","definition","denominator","shape","finite"])
def test_mat_identity_failures_are_hard(tmp_path,corruption):
    path=tmp_path/"frame_000.mat"
    source=tmp_path/"source.oct"
    meta=dict(source_file=str(source),bscan_index_matlab_1based=1,
        formal_signal_name="sv_raw",formal_signal_definition="var(abs(E), 1, 3)",variance_denominator="N",
        dimension_order="depth x A-line",reconstruction=dict(source_file=str(source),bscan_index=1,
        n_bscan_positions=500,header=dict(SPL=1024,nX=500,frame_per_pos=3)))
    image=np.ones((351,500))
    if corruption=="source": meta["source_file"]=str(tmp_path/"other.oct")
    if corruption=="frame": meta["bscan_index_matlab_1based"]=2
    if corruption=="signal": meta["formal_signal_name"]="corrected"
    if corruption=="definition": meta["formal_signal_definition"]="var(abs(E), 0, 3)"
    if corruption=="denominator": meta["variance_denominator"]="N-1"
    if corruption=="shape": image=np.ones((350,500))
    if corruption=="finite": image[0,0]=np.nan
    savemat(path,dict(sv_raw=image,metadata=meta))
    audit=SimpleNamespace(sha256="0"*64 if corruption=="sha" else DRIVER.sha(path),
                          size_bytes=path.stat().st_size,frame_index=0)
    row=dict(scan_id="synthetic",bscan_count=500,spectral_samples=1024,alines_per_bscan=500,repeat_count=3)
    with pytest.raises(ValueError):
        DRIVER.validate_mat(path,audit,row,source)


def test_d128_full_numerical_regression():
    path=OUT/"d128_nobg_regression.json"
    if not path.exists():
        pytest.skip("D128 full replay has not yet run")
    result=json.loads(path.read_text())
    assert result["status"]=="passed" and result["frames_replayed"]==2422
    assert result["kernel_sha256"]==DRIVER.sha(ROOT/"src/svrecttail/raw_sv_quantification.py")
    for field in ("source_max_relative_error","tail_max_relative_error","RI_max_relative_error"):
        assert result[field]<=1e-10
