"""Generate a deterministic, non-experimental numeric and QC example."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from svrecttail.detection import detect_tail_extent
from svrecttail.geometry import VesselGeometry, ellipse_weights
from svrecttail.io import FrameMaps
from svrecttail.localization import localize_geometry
from svrecttail.qc import save_detection_qc, save_qc_figure
from svrecttail.quantification import quantify_frame


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    return parser


def main() -> int:
    output = _parser().parse_args().output.resolve()
    existing_names = {path.name for path in output.iterdir()} if output.exists() else set()
    unexpected_existing = existing_names - {"README.md"}
    if unexpected_existing:
        raise FileExistsError(
            f"output directory contains generated files: {output}"
        )
    output.mkdir(parents=True, exist_ok=True)

    shape = (200, 100)
    z = np.arange(shape[0], dtype=np.float64)
    structure = np.repeat((25.0 + 0.01 * z)[:, None], shape[1], axis=1)
    omag = np.ones(shape, dtype=np.float64)
    omag[30:49, 45:55] = 11.0
    sv_raw = np.repeat((5.0 + 0.001 * z)[:, None], shape[1], axis=1)
    expected_geometry = VesselGeometry.from_inclusive_centres(
        x_left_center_px=45,
        x_right_center_px=54,
        z_top_center_px=30,
        diameter_um=128.0,
        dx_um=12.7,
        dz_um=6.7,
    )
    source = ellipse_weights(shape, expected_geometry, supersample=16)
    sv_raw[source > 0] += 3.0
    tail_stop = expected_geometry.z_bottom_edge_px + 500.0 / 6.7
    tail_rows = (z > expected_geometry.z_bottom_edge_px) & (z <= tail_stop)
    sv_raw[np.ix_(tail_rows, np.arange(45, 55))] += 1.0
    maps = FrameMaps(
        sv_raw=sv_raw,
        omag_raw=omag,
        stru_amp=structure,
        sv_cv2=sv_raw / (structure**2 + np.finfo(np.float64).eps),
        metadata={"source": "deterministic synthetic demonstration"},
    )
    localization = localize_geometry(
        maps.omag_raw,
        x_anchor_center_px=50.0,
        z_anchor_center_px=30.0,
        diameter_um=128.0,
        dx_um=12.7,
        dz_um=6.7,
    )
    result = quantify_frame(
        maps.sv_raw,
        localization.geometry,
        tail_length_um=500.0,
        source_qc_valid=localization.source_qc_valid,
    )
    scan_id = "synthetic_demo"
    frame_record = {
        "scan_id": scan_id,
        "data_status": "synthetic_non_experimental",
        "flow_speed_mm_s": np.nan,
        "geometry_qc_valid": localization.source_qc_valid,
        **result.summary(),
    }
    pd.DataFrame([frame_record]).to_csv(
        output / "frame_results.csv", index=False, lineterminator="\n"
    )
    pd.DataFrame(result.profile_records(scan_id)).to_csv(
        output / "profiles.csv", index=False, lineterminator="\n"
    )

    centres = np.arange(shape[0], dtype=np.float64)
    tail_start = result.geometry.z_bottom_edge_px
    signal = result.tail_contrast_profile[
        (centres > tail_start) & (centres <= tail_start + 500.0 / 6.7)
    ]
    signal = signal[: int(np.floor(500.0 / 6.7))]
    offsets = np.array([-0.2, -0.1, 0.0, 0.1, 0.2])[:, None]
    blanks = np.repeat(offsets, signal.size, axis=1)
    detection = detect_tail_extent(signal, blanks, dz_um=6.7)
    pd.DataFrame(
        [
            {
                "scan_id": scan_id,
                "data_status": "synthetic_non_experimental",
                "status": detection.invalid_reason,
                "detected": detection.detected,
                "detectable_length_um": detection.detectable_length_um,
                "right_censored": detection.right_censored,
            }
        ]
    ).to_csv(output / "detection_results.csv", index=False, lineterminator="\n")
    pd.DataFrame(detection.bin_records(scan_id)).to_csv(
        output / "detection_bins.csv", index=False, lineterminator="\n"
    )
    save_qc_figure(
        output / "synthetic_demo_QC01-03.png",
        scan_id=scan_id,
        maps=maps,
        localization=localization,
        result=result,
    )
    save_detection_qc(
        output / "synthetic_demo_QC03_detection.png",
        scan_id=scan_id,
        detection=detection,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
