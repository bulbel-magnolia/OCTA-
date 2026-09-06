"""Command-line entry point."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence

import numpy as np

from .config import load_config
from .io import load_frame_maps, load_manifest
from .mentor_tracking import write_mentor_tracking_bundle
from .pipeline import run_batch


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="svrecttail",
        description="Auditable OCTA SV rectangular-tail quantification",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    run = subparsers.add_parser("run", help="run a frozen manifest batch")
    run.add_argument("--config", required=True, type=Path)
    run.add_argument("--manifest", required=True, type=Path)
    run.add_argument("--output", required=True, type=Path)
    run.add_argument(
        "--blank-profiles",
        type=Path,
        help="optional NPZ: each key is scan_id (or all), each value is sample-by-row",
    )

    validate = subparsers.add_parser("validate", help="validate config and manifest")
    validate.add_argument("--config", required=True, type=Path)
    validate.add_argument("--manifest", required=True, type=Path)
    validate.add_argument(
        "--allow-template-gaps",
        action="store_true",
        help="validate schema without requiring source files and anchors",
    )

    inspect = subparsers.add_parser("inspect-maps", help="inspect one interim map file")
    inspect.add_argument("path", type=Path)

    tracking = subparsers.add_parser(
        "mentor-track",
        help="run mentor full-volume localization on one Flow DICOM",
    )
    tracking.add_argument("--flow-dicom", required=True, type=Path)
    tracking.add_argument("--scan-id", required=True)
    tracking.add_argument("--diameter-um", required=True, type=float)
    tracking.add_argument("--output", required=True, type=Path)
    tracking.add_argument(
        "--tracking-config",
        type=Path,
        help="optional JSON object, or a JSON file containing a tracking object",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    if args.command == "run":
        result = run_batch(
            config_path=args.config,
            manifest_path=args.manifest,
            output_dir=args.output,
            blank_profiles_path=args.blank_profiles,
        )
        print(
            json.dumps(
                {
                    "status": "complete",
                    "output_dir": str(result.output_dir),
                    "frame_count": result.frame_count,
                    "valid_frame_count": result.valid_frame_count,
                },
                ensure_ascii=False,
            )
        )
        return 0
    if args.command == "validate":
        config = load_config(args.config)
        manifest = load_manifest(
            args.manifest,
            require_complete=not args.allow_template_gaps,
            require_localization_anchors=(
                args.allow_template_gaps
                or config.localization.mode == "manifest_anchor"
            ),
            require_tracking_files=(
                not args.allow_template_gaps
                and config.localization.mode == "mentor_tracking"
            ),
        )
        print(
            json.dumps(
                {
                    "status": "valid",
                    "schema_version": config.schema_version,
                    "row_count": len(manifest),
                },
                ensure_ascii=False,
            )
        )
        return 0
    if args.command == "inspect-maps":
        maps = load_frame_maps(args.path)
        finite_sv = maps.sv_raw[np.isfinite(maps.sv_raw)]
        print(
            json.dumps(
                {
                    "shape": list(maps.sv_raw.shape),
                    "has_stru_amp": maps.stru_amp is not None,
                    "has_sv_cv2": maps.sv_cv2 is not None,
                    "sv_min": None if finite_sv.size == 0 else float(finite_sv.min()),
                    "sv_max": None if finite_sv.size == 0 else float(finite_sv.max()),
                    "metadata": {key: str(value) for key, value in maps.metadata.items()},
                },
                ensure_ascii=False,
            )
        )
        return 0
    if args.command == "mentor-track":
        tracking_config = None
        if args.tracking_config is not None:
            with args.tracking_config.open("r", encoding="utf-8") as handle:
                tracking_config = json.load(handle)
            if "tracking" in tracking_config:
                tracking_config = tracking_config["tracking"]
        bundle = write_mentor_tracking_bundle(
            args.flow_dicom,
            scan_id=args.scan_id,
            diameter_um=args.diameter_um,
            output_dir=args.output,
            tracking_config=tracking_config,
        )
        print(
            json.dumps(
                {
                    "status": "complete",
                    "primary_table": str(bundle.primary_table_path),
                    "metadata": str(bundle.metadata_path),
                    "frame_count": bundle.frame_count,
                },
                ensure_ascii=False,
            )
        )
        return 0
    raise AssertionError(f"unhandled command: {args.command}")
