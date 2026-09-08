"""Tests for the committed axial partition and strict replay diagnostics."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

from svrecttail.geometry import VesselGeometry, ellipse_weights

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("axial_audit", ROOT / "analysis/sv_diameter_stage2_scientific_v1/d500_source_axial_local_array_audit.py")
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


@pytest.mark.parametrize("diameter,width,top", [(285.,406.4,180.5),(500.,609.6,169.5)])
def test_frozen_subpixel_partition(diameter,width,top):
    g=VesselGeometry(250-width/12.7/2,250+width/12.7/2,top,diameter,12.7,6.7)
    full=ellipse_weights((351,500),g,supersample=16)
    bands=[AUDIT.band_weights((351,500),g,lo,hi) for _,lo,hi in AUDIT.BANDS]
    np.testing.assert_array_equal(sum(bands),full)
    image=np.full((351,500),7.)
    area,q,mean=AUDIT.weighted_stats(image,full,12.7*6.7)
    parts=[AUDIT.weighted_stats(image,w,12.7*6.7) for w in bands]
    assert mean==pytest.approx(7.)
    assert all(p[2]==pytest.approx(7.) for p in parts)
    assert sum(p[0] for p in parts)==pytest.approx(area)
    assert AUDIT.require_band_reconstruction(sum(p[1] for p in parts),q)<1e-12


@pytest.mark.parametrize("qsum",[100.000001,float("nan"),float("inf")])
def test_reconstruction_failure_stops(qsum):
    with pytest.raises(RuntimeError,match="does not reconstruct"):
        AUDIT.require_band_reconstruction(qsum,100.)


def test_descriptive_spearman_ties():
    assert AUDIT.descriptive_spearman([1,1,2,3],[2,2,3,1])==pytest.approx(-1/3)
    assert AUDIT.descriptive_spearman([1,2,3],[3,2,1])==pytest.approx(-1.)


def test_constant_z_has_no_defined_correlation():
    assert np.isnan(AUDIT.descriptive_spearman([1,1,1],[1,2,3]))
