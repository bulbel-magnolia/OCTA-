"""Run both versioned upper-edge methods over the 15-frame threshold matrix."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC = REPO_ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from svrecttail.mentor import tracking_core  # noqa: E402
from svrecttail.mentor_tracking import load_flow_dicom  # noqa: E402


SPEEDS = (1, 3, 5, 7, 10)
FRAMES = (0, 249, 499)
ALPHAS = (0.15, 0.20, 0.25)
NOISE_MULTIPLIERS = (3.0, 4.0, 5.0)


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--method",
        choices=("legacy", "persistent"),
        required=True,
    )
    return parser.parse_args()


def main() -> int:
    args = arguments()
    args.output.mkdir(parents=True, exist_ok=True)
    method = (
        tracking_core.UPPER_EDGE_LEGACY
        if args.method == "legacy"
        else tracking_core.UPPER_EDGE_PERSISTENT
    )
    rows: list[dict[str, object]] = []
    for speed_rank, speed in enumerate(SPEEDS, start=1):
        path = args.raw_root / f"{speed}.oct-Flow_ed.dcm"
        volume = load_flow_dicom(path)
        for noise_multiplier in NOISE_MULTIPLIERS:
            config = {
                "upper_edge_method": method,
                "upper_edge_alphas": list(ALPHAS),
                "primary_alpha": ALPHAS[0],
                "upper_edge_noise_multiplier": noise_multiplier,
            }
            outputs, _, _ = tracking_core.track_volume(
                volume,
                scan_id=f"flow{speed:02d}",
                diameter_um=128.0,
                config=config,
            )
            for alpha in ALPHAS:
                table = outputs[alpha].set_index("frame_index")
                for position_rank, frame in enumerate(FRAMES, start=1):
                    record = table.loc[frame]
                    rows.append(
                        {
                            "panel": f"{speed_rank}x{position_rank}",
                            "scan_id": f"flow{speed:02d}",
                            "flow_speed_mm_s": speed,
                            "frame_index": frame,
                            "alpha": alpha,
                            "noise_multiplier": noise_multiplier,
                            "z_peak_px": record["z_peak_px"],
                            "seed_z_upper_px": record["seed_z_upper_px"],
                            "z_upper_px": record["z_upper_px"],
                            "upper_component_width_px": record[
                                "upper_component_width_px"
                            ],
                            "peak_snr": record["peak_snr"],
                            "z_candidate_accepted": record[
                                "z_candidate_accepted"
                            ],
                            "tracking_class": record["tracking_class"],
                            "qc_valid": record["qc_valid"],
                            "qc_flags": record["qc_flags"],
                        }
                    )
        del volume
    result = pd.DataFrame.from_records(rows)
    destination = args.output / f"{args.method}_threshold_matrix_15.csv"
    result.to_csv(destination, index=False)
    print(f"wrote {len(result)} rows to {destination}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
