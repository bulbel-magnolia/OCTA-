from pathlib import Path

import numpy as np
import pandas as pd
from pydicom.dataset import FileDataset, FileMetaDataset
from pydicom.uid import ExplicitVRLittleEndian, SecondaryCaptureImageStorage, generate_uid

from svrecttail.mentor_tracking import (
    build_localization_from_tracking,
    load_flow_dicom,
    select_tracking_record,
)


def _tracking_record(*, assessment: str = "assessable") -> pd.Series:
    return pd.Series(
        {
            "scan_id": "flow01",
            "frame_index": 249,
            "alpha": 0.15,
            "x_center_px": 240.0,
            "z_upper_px": 196.0,
            "new_tracking_class": "model_assisted",
            "valid_local_body": True,
            "x1_local_geometry_px": 239.0,
            "x2_robust_centroid_px": 238.1,
            "x4_centroid_isolated_jump_corrected_px": 238.25,
            "x4_jump_corrected": False,
            "x1_fallback": False,
            "local_body_run_width_px": 13,
            "expected_lateral_width_px": 10,
            "local_body_background": 2.0,
            "local_body_sigma": 0.5,
            "local_body_peak_cnr": 8.0,
            "local_body_axial_completeness": 0.8,
            "assessability_score": 0.72,
            "vessel_presence_prediction": assessment,
            "z_edge_local_robust_snr": 4.0,
        }
    )


def test_tracking_geometry_uses_x4_x1_width_and_z_upper() -> None:
    result = build_localization_from_tracking(
        _tracking_record(), diameter_um=128.0, dx_um=12.7, dz_um=6.7
    )
    assert result.geometry.x_center_px == 238.25
    assert result.geometry.x_left_edge_px == 238.25 - 6.5
    assert result.geometry.x_right_edge_px == 238.25 + 6.5
    assert result.geometry.z_top_edge_px == 195.5
    assert np.isclose(result.geometry.z_bottom_edge_px, 195.5 + 128.0 / 6.7)
    assert result.mentor_tracking is not None
    assert result.source_qc_valid


def test_tracking_qc_rejects_uncertain_frame_without_moving_geometry() -> None:
    result = build_localization_from_tracking(
        _tracking_record(assessment="uncertain"),
        diameter_um=128.0,
        dx_um=12.7,
        dz_um=6.7,
    )
    assert result.geometry.x_center_px == 238.25
    assert not result.source_qc_valid


def test_select_tracking_record_requires_unique_scan_and_frame(tmp_path: Path) -> None:
    table = pd.DataFrame([_tracking_record(), _tracking_record()])
    path = tmp_path / "tracking.csv"
    table.to_csv(path, index=False)
    try:
        select_tracking_record(path, scan_id="flow01", frame_index=249)
    except ValueError as error:
        assert "exactly one" in str(error)
    else:
        raise AssertionError("duplicate tracking rows were accepted")


def test_flow_dicom_loader_preserves_frame_z_x_order(tmp_path: Path) -> None:
    array = np.arange(3 * 4 * 5, dtype=np.uint16).reshape(3, 4, 5)
    file_meta = FileMetaDataset()
    file_meta.MediaStorageSOPClassUID = SecondaryCaptureImageStorage
    file_meta.MediaStorageSOPInstanceUID = generate_uid()
    file_meta.TransferSyntaxUID = ExplicitVRLittleEndian
    dataset = FileDataset(str(tmp_path / "flow.dcm"), {}, file_meta=file_meta)
    dataset.SOPClassUID = file_meta.MediaStorageSOPClassUID
    dataset.SOPInstanceUID = file_meta.MediaStorageSOPInstanceUID
    dataset.Rows = 4
    dataset.Columns = 5
    dataset.NumberOfFrames = 3
    dataset.SamplesPerPixel = 1
    dataset.PhotometricInterpretation = "MONOCHROME2"
    dataset.BitsAllocated = 16
    dataset.BitsStored = 16
    dataset.HighBit = 15
    dataset.PixelRepresentation = 0
    dataset.PixelData = array.tobytes()
    path = tmp_path / "flow.dcm"
    dataset.save_as(path, enforce_file_format=True)

    loaded = load_flow_dicom(path)
    assert loaded.shape == (3, 4, 5)
    np.testing.assert_array_equal(loaded, array.astype(np.float32))
