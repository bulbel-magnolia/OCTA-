"""Interim map loading and manifest validation."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import h5py
import numpy as np
import pandas as pd
from numpy.typing import NDArray
from scipy.io import loadmat


@dataclass(frozen=True)
class FrameMaps:
    sv_raw: NDArray[np.float64]
    omag_raw: NDArray[np.float64]
    stru_amp: NDArray[np.float64] | None = None
    sv_cv2: NDArray[np.float64] | None = None
    metadata: dict[str, Any] = field(default_factory=dict)


REQUIRED_MANIFEST_COLUMNS = {
    "scan_id",
    "source_file",
    "vessel_id",
    "phantom_id",
    "session_id",
    "diameter_um",
    "flow_speed_mm_s",
    "position_label",
    "dx_um",
    "dz_um",
    "bscan_index",
    "slow_axis_position_um",
    "temporal_repeat_id",
    "temporal_repeat_count",
    "scan_time_interval_s",
    "acquisition_order",
    "reconstruction_version",
    "x_anchor_center_px",
    "z_anchor_center_px",
    "geometry_source",
    "background_excluded_side",
    "background_exclusion_reason",
}


def load_manifest(
    path: str | Path,
    *,
    require_complete: bool = True,
    require_localization_anchors: bool = True,
) -> pd.DataFrame:
    """Load a one-row-per-frame manifest and validate formal inputs."""

    manifest_path = Path(path)
    table = pd.read_csv(
        manifest_path,
        dtype={
            "scan_id": "string",
            "source_file": "string",
            "vessel_id": "string",
            "phantom_id": "string",
            "session_id": "string",
            "temporal_repeat_id": "string",
            "background_excluded_side": "string",
            "background_exclusion_reason": "string",
        },
    )
    missing = sorted(REQUIRED_MANIFEST_COLUMNS - set(table.columns))
    if missing:
        raise ValueError(f"manifest missing columns: {', '.join(missing)}")
    if table.empty:
        raise ValueError("manifest must contain at least one frame")
    if table["scan_id"].isna().any() or table["scan_id"].str.strip().eq("").any():
        raise ValueError("scan_id must be non-empty")
    required_numeric = [
        "diameter_um",
        "flow_speed_mm_s",
        "dx_um",
        "dz_um",
        "bscan_index",
    ]
    optional_numeric = [
        "slow_axis_position_um",
        "temporal_repeat_count",
        "scan_time_interval_s",
        "acquisition_order",
        "x_anchor_center_px",
        "z_anchor_center_px",
    ]
    numeric = required_numeric + optional_numeric
    for column in numeric:
        original = table[column]
        converted = pd.to_numeric(original, errors="coerce")
        supplied = original.notna() & original.astype("string").str.strip().ne("")
        invalid_numeric = supplied & converted.isna()
        if invalid_numeric.any():
            ids = ", ".join(table.loc[invalid_numeric, "scan_id"].astype(str))
            raise ValueError(f"{column} contains non-numeric values for: {ids}")
        table[column] = converted
    excluded_side = table["background_excluded_side"].fillna("").str.strip().str.lower()
    exclusion_reason = table["background_exclusion_reason"].fillna("").str.strip()
    if (~excluded_side.isin(["", "left", "right"])).any():
        raise ValueError("background_excluded_side must be blank, left, or right")
    if (excluded_side.eq("") != exclusion_reason.eq("")).any():
        raise ValueError(
            "background_excluded_side and background_exclusion_reason must be supplied together"
        )
    table["background_excluded_side"] = excluded_side.replace("", pd.NA)
    scan_level_columns = (
        "vessel_id",
        "phantom_id",
        "session_id",
        "diameter_um",
        "flow_speed_mm_s",
        "dx_um",
        "dz_um",
        "temporal_repeat_count",
        "scan_time_interval_s",
        "acquisition_order",
        "reconstruction_version",
    )
    inconsistent = [
        column
        for column in scan_level_columns
        if table.groupby("scan_id", dropna=False)[column].nunique(dropna=True).gt(1).any()
    ]
    if inconsistent:
        raise ValueError(
            "scan-level fields vary within scan_id: " + ", ".join(inconsistent)
        )
    if require_complete:
        required_values = [
            "source_file",
            "position_label",
            "reconstruction_version",
            *required_numeric,
            "geometry_source",
        ]
        if require_localization_anchors:
            required_values.extend(["x_anchor_center_px", "z_anchor_center_px"])
        incomplete = table[required_values].isna().any(axis=1)
        for column in (
            "source_file",
            "position_label",
            "reconstruction_version",
            "geometry_source",
        ):
            incomplete |= table[column].astype("string").str.strip().eq("").fillna(True)
        if incomplete.any():
            ids = ", ".join(table.loc[incomplete, "scan_id"].astype(str))
            raise ValueError(f"incomplete formal manifest rows: {ids}")
        missing_sources = [
            f"{row.scan_id}: {row.source_file}"
            for row in table.itertuples(index=False)
            if not resolve_source_path(manifest_path, str(row.source_file)).is_file()
        ]
        if missing_sources:
            raise FileNotFoundError(
                "manifest source files not found: " + "; ".join(missing_sources)
            )
        if (table[["diameter_um", "dx_um", "dz_um"]] <= 0).any().any():
            raise ValueError("diameter and pixel pitches must be positive")
        if (table["flow_speed_mm_s"] <= 0).any():
            raise ValueError("formal flow speed must be positive")
        if (table["bscan_index"] < 0).any():
            raise ValueError("bscan_index must be non-negative")
        if (table[["x_anchor_center_px", "z_anchor_center_px"]] < 0).any().any():
            raise ValueError("localization anchors must be non-negative")
        if (table["acquisition_order"].dropna() < 0).any():
            raise ValueError("acquisition_order must be non-negative when supplied")
        if (table["temporal_repeat_count"].dropna() <= 0).any():
            raise ValueError("temporal_repeat_count must be positive when supplied")
        if (table["scan_time_interval_s"].dropna() <= 0).any():
            raise ValueError("scan_time_interval_s must be positive when supplied")
        for column in ("bscan_index", "acquisition_order", "temporal_repeat_count"):
            values = table[column].dropna().to_numpy(dtype=float)
            if not np.equal(values, np.floor(values)).all():
                raise ValueError(f"{column} must contain integer values")
        repeat_key = table.get(
            "temporal_repeat_id", pd.Series("", index=table.index)
        ).fillna("").astype(str)
        keys = pd.DataFrame(
            {
                "scan_id": table["scan_id"].astype(str),
                "bscan_index": table["bscan_index"],
                "temporal_repeat_id": repeat_key,
            }
        )
        if keys.duplicated().any():
            raise ValueError(
                "scan_id, bscan_index, and temporal_repeat_id must identify one frame"
            )
    return table


def _field(container: Any, name: str) -> Any:
    if hasattr(container, name):
        return getattr(container, name)
    if isinstance(container, dict) and name in container:
        return container[name]
    return None


def _matlab_metadata_value(value: Any) -> Any:
    """Convert scipy MATLAB structs into JSON-compatible Python containers."""

    if hasattr(value, "_fieldnames"):
        return {
            name: _matlab_metadata_value(getattr(value, name))
            for name in value._fieldnames
        }
    if isinstance(value, np.ndarray):
        if value.size == 1:
            return _matlab_metadata_value(value.item())
        return [_matlab_metadata_value(item) for item in value.tolist()]
    if isinstance(value, (list, tuple)):
        return [_matlab_metadata_value(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def _validate_maps(
    sv_raw: Any,
    omag_raw: Any,
    stru_amp: Any | None,
    sv_cv2: Any | None,
    metadata: dict[str, Any],
) -> FrameMaps:
    sv = np.asarray(sv_raw, dtype=np.float64).squeeze()
    omag = np.asarray(omag_raw, dtype=np.float64).squeeze()
    structure = None if stru_amp is None else np.asarray(stru_amp, dtype=np.float64).squeeze()
    cv2 = None if sv_cv2 is None else np.asarray(sv_cv2, dtype=np.float64).squeeze()
    if sv.ndim != 2 or omag.ndim != 2 or sv.shape != omag.shape:
        raise ValueError("sv_raw and omag_raw must be matching 2-D arrays")
    if structure is not None and structure.shape != sv.shape:
        raise ValueError("stru_amp must match sv_raw shape")
    if cv2 is not None and cv2.shape != sv.shape:
        raise ValueError("sv_cv2 must match sv_raw shape")
    if (
        np.isinf(sv).any()
        or np.isinf(omag).any()
        or (structure is not None and np.isinf(structure).any())
        or (cv2 is not None and np.isinf(cv2).any())
    ):
        raise ValueError("input maps cannot contain infinity")
    finite_sv = sv[np.isfinite(sv)]
    if finite_sv.size and float(finite_sv.min()) < -1e-12 * max(1.0, float(finite_sv.max())):
        raise ValueError("sv_raw contains negative values; supply linear uncorrected variance")
    return FrameMaps(
        sv_raw=sv,
        omag_raw=omag,
        stru_amp=structure,
        sv_cv2=cv2,
        metadata=metadata,
    )


def _load_classic_mat(path: Path) -> FrameMaps:
    content = loadmat(path, squeeze_me=True, struct_as_record=False)
    root = content.get("result", content)
    sv_raw = _field(root, "sv_raw")
    omag_raw = _field(root, "omag_raw")
    stru_amp = _field(root, "stru_amp")
    sv_cv2 = _field(root, "sv_cv2")
    if sv_cv2 is None:
        sv_cv2 = _field(root, "cv2")
    if sv_cv2 is None:
        sv_cv2 = _field(root, "sv_norm")
    if sv_raw is None or omag_raw is None:
        raise ValueError("MAT file must contain sv_raw and omag_raw")
    metadata: dict[str, Any] = {"source_path": str(path), "format": "mat-v5-v7"}
    embedded_metadata = _field(root, "metadata")
    if embedded_metadata is not None:
        converted_metadata = _matlab_metadata_value(embedded_metadata)
        metadata["export_metadata"] = converted_metadata
        if isinstance(converted_metadata, dict):
            for name in (
                "source_file",
                "bscan_index_matlab_1based",
                "formal_signal_definition",
                "dimension_order",
            ):
                if name in converted_metadata:
                    metadata[name] = converted_metadata[name]
    for name in ("source_file", "bscan_index", "x_idx", "z_top", "SV_definition"):
        value = _field(root, name)
        if value is not None:
            if isinstance(value, np.ndarray) and value.size == 1:
                value = value.item()
            metadata[name] = value
    if "x_idx" in metadata or "z_top" in metadata:
        metadata["stored_anchor_index_base"] = 1
    return _validate_maps(sv_raw, omag_raw, stru_amp, sv_cv2, metadata)


def _read_hdf5_dataset(handle: h5py.File, names: tuple[str, ...]) -> NDArray[Any] | None:
    for name in names:
        if name in handle and isinstance(handle[name], h5py.Dataset):
            return np.asarray(handle[name])
    return None


def _load_hdf5(path: Path) -> FrameMaps:
    with h5py.File(path, "r") as handle:
        sv_raw = _read_hdf5_dataset(handle, ("sv_raw", "result/sv_raw"))
        omag_raw = _read_hdf5_dataset(handle, ("omag_raw", "result/omag_raw"))
        stru_amp = _read_hdf5_dataset(handle, ("stru_amp", "result/stru_amp"))
        sv_cv2 = _read_hdf5_dataset(
            handle,
            ("sv_cv2", "cv2", "result/sv_cv2", "result/cv2", "result/sv_norm"),
        )
    if sv_raw is None or omag_raw is None:
        raise ValueError("HDF5 file must contain sv_raw and omag_raw datasets")
    return _validate_maps(
        sv_raw,
        omag_raw,
        stru_amp,
        sv_cv2,
        {"source_path": str(path), "format": "hdf5"},
    )


def load_frame_maps(path: str | Path) -> FrameMaps:
    """Load linear maps from MAT, HDF5, or NumPy NPZ interim storage."""

    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(source)
    suffix = source.suffix.lower()
    if suffix == ".npz":
        with np.load(source, allow_pickle=False) as content:
            if "sv_raw" not in content or "omag_raw" not in content:
                raise ValueError("NPZ file must contain sv_raw and omag_raw")
            stru_amp = content["stru_amp"] if "stru_amp" in content else None
            sv_cv2 = content["sv_cv2"] if "sv_cv2" in content else None
            return _validate_maps(
                content["sv_raw"],
                content["omag_raw"],
                stru_amp,
                sv_cv2,
                {"source_path": str(source), "format": "npz"},
            )
    if suffix in {".h5", ".hdf5"}:
        return _load_hdf5(source)
    if suffix == ".mat":
        try:
            return _load_classic_mat(source)
        except NotImplementedError:
            return _load_hdf5(source)
        except ValueError as error:
            if "Unknown mat file type" in str(error):
                return _load_hdf5(source)
            raise
    raise ValueError(f"unsupported interim map format: {suffix}")


def resolve_source_path(manifest_path: str | Path, source_file: str) -> Path:
    """Resolve a source path relative to its manifest without changing CWD."""

    source = Path(source_file)
    if source.is_absolute():
        return source
    return Path(manifest_path).resolve().parent / source
