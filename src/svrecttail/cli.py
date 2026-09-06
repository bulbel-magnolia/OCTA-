"""Command-line entry point."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence

import numpy as np

from .config import load_config
from .io import load_frame_maps, load_manifest
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
            args.manifest, require_complete=not args.allow_template_gaps
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
    raise AssertionError(f"unhandled command: {args.command}")
