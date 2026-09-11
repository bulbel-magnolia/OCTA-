"""Configuration, manifest, and MATLAB v7.3 Flow I/O helpers.

All arrays returned by :func:`load_flow` use the explicit Python order
``[frame, z, x]``.  The exporter writes MATLAB ``[z, x, frame]``; MATLAB v7.3
therefore appears to h5py as ``[frame, x, z]`` and is transposed here.
"""

from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path
from typing import Any

import h5py
import numpy as np
import pandas as pd


LOG = logging.getLogger(__name__)


def load_config(path: str | Path) -> tuple[dict[str, Any], Path]:
    path = Path(path).resolve()
    with path.open("r", encoding="utf-8") as handle:
        cfg = json.load(handle)
    required = {"study_id", "source", "output", "scales", "metrics", "aggregation"}
    missing = required.difference(cfg)
    if missing:
        raise ValueError(f"Configuration is missing sections: {sorted(missing)}")
    return cfg, path


def config_digest(config_path: str | Path) -> str:
    return hashlib.sha256(Path(config_path).read_bytes()).hexdigest()


def load_manifest(cfg: dict[str, Any], config_path: Path) -> pd.DataFrame:
    manifest_path = Path(cfg["source"]["manifest_csv"])
    if not manifest_path.is_absolute():
        manifest_path = config_path.parent / manifest_path
    table = pd.read_csv(manifest_path)
    expected = {
        "dataset_id",
        "diameter_um",
        "velocity_mm_s",
        "technical_replicate",
        "scan_number",
        "raw_path",
    }
    missing = expected.difference(table.columns)
    if missing:
        raise ValueError(f"Manifest is missing columns: {sorted(missing)}")
    if table["dataset_id"].duplicated().any():
        duplicated = table.loc[table["dataset_id"].duplicated(), "dataset_id"].tolist()
        raise ValueError(f"Duplicate dataset_id values: {duplicated}")
    return table


def output_root(cfg: dict[str, Any]) -> Path:
    return Path(cfg["output"]["root"]).resolve()


def export_path(cfg: dict[str, Any], dataset_id: str) -> Path:
    return output_root(cfg) / cfg["output"]["export_subdir"] / f"{dataset_id}.mat"


def read_export_metadata(flow_path: str | Path) -> dict[str, Any]:
    sidecar = Path(f"{Path(flow_path)}.metadata.json")
    if not sidecar.exists():
        return {}
    with sidecar.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _choose_dataset(handle: h5py.File) -> str:
    for candidate in ("p_bld_ed", "p_bld_ed_sparse"):
        if candidate in handle:
            return candidate
    available = [key for key, value in handle.items() if isinstance(value, h5py.Dataset)]
    raise KeyError(
        "No p_bld_ed or p_bld_ed_sparse dataset found; "
        f"available datasets are {available}"
    )


def load_flow(flow_path: str | Path) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    """Load float Flow and return ``(flow[frame,z,x], source_frames, metadata)``.

    The metadata sidecar is authoritative.  Shape inference is allowed only for
    legacy sparse validation files and requires an unambiguous 351-pixel axial
    dimension; otherwise the function fails loudly.
    """

    flow_path = Path(flow_path)
    if not flow_path.exists():
        raise FileNotFoundError(flow_path)
    metadata = read_export_metadata(flow_path)
    with h5py.File(flow_path, "r") as handle:
        name = _choose_dataset(handle)
        raw = np.asarray(handle[name], dtype=np.float32)

    if raw.ndim != 3:
        raise ValueError(f"{flow_path}: expected a 3-D Flow volume, got {raw.shape}")
    axis_order = str(metadata.get("axis_order_h5py", "")).replace(" ", "").lower()
    if axis_order == "frame,x,z":
        flow = raw.transpose(0, 2, 1)
    elif axis_order == "frame,z,x":
        flow = raw
    elif not axis_order:
        if raw.shape[-1] == 351 and raw.shape[1] != 351:
            LOG.warning("%s has no sidecar; inferring legacy h5py [frame,x,z]", flow_path)
            flow = raw.transpose(0, 2, 1)
        elif raw.shape[1] == 351 and raw.shape[-1] != 351:
            LOG.warning("%s has no sidecar; inferring [frame,z,x]", flow_path)
            flow = raw
        else:
            raise ValueError(
                f"{flow_path}: ambiguous axis order {raw.shape}; a metadata sidecar is required"
            )
    else:
        raise ValueError(f"{flow_path}: unsupported axis_order_h5py={axis_order!r}")

    if not np.issubdtype(flow.dtype, np.floating):
        raise TypeError(f"{flow_path}: Flow must be floating point, got {flow.dtype}")
    if not np.isfinite(flow).all():
        raise ValueError(f"{flow_path}: Flow contains NaN or infinite values")

    source_frames = metadata.get("source_frame_indices")
    source_frames_are_matlab = source_frames is not None
    if source_frames is None and "sparse" in name:
        # Legacy validation exports store the MATLAB source positions as a MAT
        # variable rather than in the JSON sidecar.
        try:
            with h5py.File(flow_path, "r") as handle:
                if "frame_indices" in handle:
                    source_frames = np.asarray(handle["frame_indices"]).squeeze().astype(int)
                    source_frames_are_matlab = True
        except OSError:
            source_frames = None
    if source_frames is None:
        source_frames = np.arange(flow.shape[0], dtype=int)
    else:
        source_frames = np.asarray(source_frames, dtype=int)
        # MATLAB exporter stores source positions as 1..500.  Internally all
        # Python frame indices are 0-based.
        if source_frames_are_matlab and len(source_frames) and source_frames.min() >= 1:
            source_frames = source_frames - 1
    if len(source_frames) != flow.shape[0]:
        raise ValueError(
            f"{flow_path}: {len(source_frames)} source frame indices for {flow.shape[0]} frames"
        )
    return np.ascontiguousarray(flow), source_frames, metadata


def alpha_token(alpha: float) -> str:
    return f"{alpha:.2f}".replace(".", "p")


def tracking_csv_path(tracking_root: Path, dataset_id: str, alpha: float) -> Path:
    scan_dir = tracking_root / dataset_id
    preferred = [
        scan_dir / f"frame_tracking_alpha{int(round(alpha * 100)):03d}.csv",
        scan_dir / f"frame_tracking_alpha{alpha_token(alpha)}.csv",
        scan_dir / f"frame_tracking_alpha_{alpha_token(alpha)}.csv",
        scan_dir / f"tracking_alpha{alpha_token(alpha)}.csv",
    ]
    for path in preferred:
        if path.exists():
            return path
    candidates = sorted(scan_dir.glob("*.csv"))
    for path in candidates:
        try:
            header = pd.read_csv(path, nrows=4)
        except Exception:
            continue
        if "alpha" in header and np.isclose(header["alpha"].astype(float), alpha).all():
            return path
    raise FileNotFoundError(
        f"No tracking CSV for {dataset_id}, alpha={alpha:.2f} under {scan_dir}"
    )


def load_tracking(path: str | Path, alpha: float | None = None) -> pd.DataFrame:
    table = pd.read_csv(path)
    aliases = {
        "x_c_px": "x_center_px",
        "z_u_px": "z_upper_px",
        "confidence_class": "tracking_class",
    }
    table = table.rename(columns={key: value for key, value in aliases.items() if key in table})
    required = {"frame_index", "x_center_px", "z_upper_px", "tracking_class"}
    missing = required.difference(table.columns)
    if missing:
        raise ValueError(f"{path}: tracking CSV is missing {sorted(missing)}")
    if alpha is not None:
        if "alpha" not in table:
            table["alpha"] = alpha
        table = table.loc[np.isclose(table["alpha"].astype(float), alpha)].copy()
    table["frame_index"] = table["frame_index"].astype(int)
    if "qc_valid" in table:
        raw = table["qc_valid"]
        if raw.dtype == object:
            mapped = raw.astype(str).str.strip().str.lower().map(
                {"true": True, "false": False, "1": True, "0": False}
            )
            if mapped.isna().any():
                raise ValueError(f"{path}: qc_valid contains non-boolean values")
            table["qc_valid"] = mapped.astype(bool)
        else:
            table["qc_valid"] = raw.astype(bool)
    if table["frame_index"].duplicated().any():
        raise ValueError(f"{path}: duplicate frame_index values")
    return table.sort_values("frame_index").reset_index(drop=True)


def write_run_metadata(
    path: str | Path,
    *,
    config_path: Path,
    extra: dict[str, Any] | None = None,
) -> None:
    from datetime import datetime, timezone

    from . import __version__

    payload: dict[str, Any] = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "tailq_version": __version__,
        "config_path": str(config_path),
        "config_sha256": config_digest(config_path),
        "random_seed": 20260414,
    }
    if extra:
        payload.update(extra)
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
