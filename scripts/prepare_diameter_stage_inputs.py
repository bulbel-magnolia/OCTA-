"""Read-only raw intake; explicit naming convention, content hashes, no SV analysis."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import struct
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path, PureWindowsPath

import pydicom

RULE = "complete pair and complete metadata; prefer non-128; sort diameter, flow, lexical OCT relative path; select first; never reselect based on images"


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(path, rows, fields=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = fields or list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")


def relative(path):
    text = str(path).replace("\\", "/")
    if PureWindowsPath(text).is_absolute() or text.startswith("/") or ".." in text.split("/"):
        raise ValueError("Expected relative path")
    return text


def parse_condition(path, *, naming_documented=False):
    parts = relative(path).split("/")
    diameter = flow = None
    dsource = fsource = "unknown"
    # Bare numbers are allowed ONLY with the experiment's documented naming rule.
    values = set()
    for part in parts[:-1]:
        match = re.fullmatch(r"(?:D)?(\d+(?:\.\d+)?)\s*(um|µm|μm|mm)?", part, re.I)
        if match and (match[2] or naming_documented):
            value = float(match[1]) * (1000 if (match[2] or "").lower() == "mm" else 1)
            if value > 0:
                values.add(value)
    if len(values) == 1:
        diameter = values.pop()
        dsource = "parent_directory_name;experiment_naming_note" if naming_documented else "directory_explicit_unit"
    match = re.fullmatch(r"(\d+(?:\.\d+)?)(?:mm_s)?\.oct", parts[-1], re.I)
    if match and (naming_documented or "mm_s" in parts[-1]):
        flow = float(match[1])
        fsource = "oct_filename;experiment_naming_note" if naming_documented else "filename_explicit_unit"
    return dict(diameter_um=diameter, flow_speed_mm_s=flow, diameter_source=dsource, flow_source=fsource)


def pair_oct(oct_path, dicom_paths):
    oct_path = relative(oct_path)
    paths = [relative(p) for p in dicom_paths]
    exact = [p for p in paths if p.casefold() == (oct_path + "-Flow_ed.dcm").casefold()]
    if exact:
        return (exact[0], "exact_companion", "complete") if len(exact) == 1 else (None, "exact_companion", "ambiguous")
    name = oct_path.rsplit("/", 1)[-1]
    stem = name[:-4]
    pattern = re.compile(r"^" + re.escape(stem) + r"(?:\.oct)?[-_]Flow[-_]ed(?:\(\d+\))?\.dcm$", re.I)
    candidates = [p for p in paths if pattern.fullmatch(p.rsplit("/", 1)[-1])]
    parent = oct_path.rpartition("/")[0]
    same_dir = [p for p in candidates if p.rpartition("/")[0] == parent]
    candidates = same_dir or candidates
    return (candidates[0], "unique_stem_same_directory" if same_dir else "unique_stem_recursive", "complete") if len(candidates) == 1 else (None, "stem_search", "ambiguous" if candidates else "missing")


def assign_ids(rows, prior=None):
    groups = defaultdict(list)
    for row in rows:
        row.update(scan_id="", volume_index="", legacy_scan_id="", existing_d128_reference=False)
        if row["pairing_status"] == "complete" and row["diameter_um"] is not None and row["flow_speed_mm_s"] is not None:
            groups[(row["diameter_um"], row["flow_speed_mm_s"])].append(row)
    def token(value, width):
        return f"{int(value):0{width}d}" if value == int(value) else str(value).replace(".", "p")
    for (diameter, flow), group in groups.items():
        for index, row in enumerate(sorted(group, key=lambda r: r["oct_relative_path"]), 1):
            row.update(scan_id=f"D{token(diameter,3)}_F{token(flow,2)}_V{index:02d}", volume_index=index,
                       existing_d128_reference=diameter == 128)
    if prior:
        old = {r["oct_relative_path"]: r for r in prior}
        for row in rows:
            previous = old.get(row["oct_relative_path"])
            if previous and any(str(previous[k]) != str(row[k]) for k in ("scan_id", "oct_sha256", "flow_dicom_sha256")):
                raise ValueError("Existing acquisition identity would change; create a new intake revision")
    return rows


def preserve_job_destinations(jobs, previous):
    """Keep previously chosen derived storage locations after an intake rerun."""
    old = {job["scan_id"]: job for job in previous}
    for job in jobs:
        prior = old.get(job["scan_id"])
        if prior:
            if Path(prior["source_file"]).resolve() != Path(job["source_file"]).resolve() or prior["frame_indices_0based"] != job["frame_indices_0based"]:
                raise ValueError("Existing export job identity would change")
            job["interim_dir"] = prior["interim_dir"]
    return jobs


def oct_header(path):
    names = ["bob", "SPL", "nX", "nY_raw", "Boffset", "Blength_minus1", "Xcenter", "Xspan", "Ycenter", "Yspan", "frame_per_pos", "reserved", "sizeBck"]
    fmt = "<IdIIIIddddIII"
    with Path(path).open("rb") as stream:
        data = stream.read(struct.calcsize(fmt))
    h = dict(zip(names, struct.unpack(fmt, data)))
    for field in ("SPL", "nX", "nY_raw", "frame_per_pos"):
        if h[field] <= 0 or h[field] != int(h[field]):
            raise ValueError("Invalid OCT header " + field)
    if h["nY_raw"] % h["frame_per_pos"]:
        raise ValueError("Incomplete repeat group")
    expected = h["bob"] + (int(h["SPL"]) * h["nX"] + 2) * 2 * h["nY_raw"]
    if Path(path).stat().st_size != expected:
        raise ValueError("OCT size does not match header and all repeated frames")
    if h["Boffset"] + h["Blength_minus1"] + 1 > h["SPL"] or h["SPL"] < 400:
        raise ValueError("Incompatible spectral range")
    if h["bob"] < struct.calcsize(fmt) + 2 * h["sizeBck"]:
        raise ValueError("Invalid data offset")
    return h


def protocol(h):
    return tuple(h[k] for k in ("SPL", "nX", "nY_raw", "Boffset", "Blength_minus1", "Xspan", "Yspan", "frame_per_pos"))


def prepare(raw_root, repo, naming_note):
    raw_root, repo = Path(raw_root).resolve(), Path(repo).resolve()
    local = repo / "outputs/local_run_inputs/diameter_stage_v1"
    public = repo / "data/diameter_stage_v1"
    results = repo / "results/diameter_intake_v1"
    note = Path(naming_note)
    note_text = note.read_text(encoding="utf-8-sig")
    if "文件名即流速" not in note_text or "文件夹名即血管直径" not in note_text:
        raise ValueError("Naming note does not establish the required convention")
    files = sorted((p for p in raw_root.rglob("*") if p.is_file() and p.suffix.lower() in {".oct", ".dcm"}), key=lambda p: p.relative_to(raw_root).as_posix())
    inventory, headers, dicom_meta = [], {}, {}
    for index, path in enumerate(files):
        rel = path.relative_to(raw_root).as_posix()
        stat = path.stat()
        asset = "raw_oct" if path.suffix.lower() == ".oct" else "other_dicom"
        info = {}
        if path.suffix.lower() == ".dcm":
            try:
                ds = pydicom.dcmread(path, stop_before_pixels=True)
                info = dict(NumberOfFrames=int(getattr(ds, "NumberOfFrames", 1)), Rows=int(ds.Rows), Columns=int(ds.Columns),
                            pixel_spacing=list(map(float, ds.PixelSpacing)) if hasattr(ds, "PixelSpacing") else None,
                            dicom_header_status="readable")
                if re.search(r"[-_]Flow[-_]ed(?:\(\d+\))?\.dcm$", path.name, re.I):
                    asset = "localization_flow_dicom"
                dicom_meta[rel] = info
            except Exception as error:
                info = dict(dicom_header_status=type(error).__name__)
        else:
            try:
                headers[rel] = oct_header(path)
                info["oct_header_status"] = "readable"
            except Exception as error:
                info["oct_header_status"] = str(error)
        digest = sha256(path)
        after = path.stat()
        if (stat.st_size, stat.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
            raise ValueError("Raw file changed during hashing: " + rel)
        inventory.append(dict(file_name=path.name, absolute_path=str(path), relative_path_from_RAW_ROOT=rel,
                              extension=path.suffix.lower(), size_bytes=stat.st_size, sha256=digest,
                              last_modified_time=datetime.fromtimestamp(stat.st_mtime, timezone.utc).isoformat(),
                              inferred_asset_type=asset, **info))
        print(f"HASHED {index+1}/{len(files)} {rel}", flush=True)
    by_path = {r["relative_path_from_RAW_ROOT"]: r for r in inventory}
    flow_paths = [r["relative_path_from_RAW_ROOT"] for r in inventory if r["inferred_asset_type"] == "localization_flow_dicom"]
    rows, pairs = [], []
    reference_protocols = {protocol(h) for p, h in headers.items() if parse_condition(p, naming_documented=True)["diameter_um"] == 128}
    for item in inventory:
        if item["inferred_asset_type"] != "raw_oct":
            continue
        rel = item["relative_path_from_RAW_ROOT"]
        companion, method, status = pair_oct(rel, flow_paths)
        dcm = by_path.get(companion, {})
        condition = parse_condition(rel, naming_documented=True)
        h = headers.get(rel)
        dm = dicom_meta.get(companion, {})
        compatible = bool(h and protocol(h) in reference_protocols and dm.get("Rows") == 351 and dm.get("Columns") == h["nX"] and dm.get("NumberOfFrames") == h["nY_raw"] // h["frame_per_pos"])
        row = dict(**condition, independent_volume=True, oct_relative_path=rel, oct_absolute_path=item["absolute_path"], oct_size_bytes=item["size_bytes"], oct_sha256=item["sha256"],
                   flow_dicom_relative_path=companion or "", flow_dicom_absolute_path=dcm.get("absolute_path", ""), flow_dicom_size_bytes=dcm.get("size_bytes", ""), flow_dicom_sha256=dcm.get("sha256", ""),
                   dx_um=12.7 if compatible else None, dz_um=6.7 if compatible else None,
                   bscan_count=h["nY_raw"] // h["frame_per_pos"] if h else None, alines_per_bscan=h["nX"] if h else None,
                   repeat_count=h["frame_per_pos"] if h else None, spectral_samples=int(h["SPL"]) if h else None,
                   acquisition_metadata_status="complete" if compatible and all(condition[k] is not None for k in ("diameter_um", "flow_speed_mm_s")) and status == "complete" else "blocked",
                   calibration_source="inherited_project_acquisition_calibration" if compatible else "unknown",
                   calibration_evidence="OCT spectral, repeat, scan-span and dimension signature matches D128; Flow dimensions agree; calibration not independently measured" if compatible else "protocol or header mismatch",
                   pairing_status=status, pairing_method=method, metadata_inference_notes="Experiment naming note; OCT header; B-scans are spatial samples nested within acquisition")
        rows.append(row)
        pairs.append(dict(oct_relative_path=rel, oct_sha256=item["sha256"], dicom_relative_path=companion or "", dicom_sha256=dcm.get("sha256", ""), pairing_method=method, pairing_status=status))
    # Reject accidental many-to-one assignment even when an individual stem search was unique.
    usage = Counter(r["flow_dicom_relative_path"] for r in rows if r["pairing_status"] == "complete")
    for row, pair in zip(rows, pairs):
        if usage[row["flow_dicom_relative_path"]] > 1:
            row["pairing_status"] = pair["pairing_status"] = "ambiguous"
            row["acquisition_metadata_status"] = "blocked"
    prior_path = public / "volume_manifest.csv"
    prior = list(csv.DictReader(prior_path.open(encoding="utf-8"))) if prior_path.exists() else None
    assign_ids(rows, prior)
    # Legacy linkage requires content identity with the frozen manifest, not merely a filename.
    legacy_index = {}
    for path in (repo / "results/formal_sv_d128_v21_full2500_run001").glob("input_sha256.csv"):
        for entry in csv.DictReader(path.open(encoding="utf-8-sig")):
            legacy_index.setdefault(entry.get("scan_id", ""), []).append(entry)
    for row in rows:
        if row["existing_d128_reference"]:
            legacy = f"flow{int(row['flow_speed_mm_s']):02d}"
            text = json.dumps(legacy_index.get(legacy, []))
            if row["oct_sha256"] in text and row["flow_dicom_sha256"] in text:
                row["legacy_scan_id"] = legacy
    eligible = sorted((r for r in rows if r["acquisition_metadata_status"] == "complete"), key=lambda r: (r["diameter_um"], r["flow_speed_mm_s"], r["oct_relative_path"]))
    non128 = [r for r in eligible if r["diameter_um"] != 128]
    pilot = (non128 or eligible or [None])[0]
    selection = dict(rule=RULE, scan_id=pilot["scan_id"] if pilot else None,
                     selected_before_tracking=True, selected_at_utc=datetime.now(timezone.utc).isoformat())
    if (local / "pilot_selection.json").exists():
        existing = json.loads((local / "pilot_selection.json").read_text(encoding="utf-8"))
        if existing["scan_id"] != selection["scan_id"]:
            raise ValueError("Pilot selection is frozen; intake change requires explicit revision")
        selection = existing
    jobs = [dict(scan_id=r["scan_id"], source_file=r["oct_absolute_path"], interim_dir=str(repo / "outputs/local_interim/diameter_stage_v1" / r["scan_id"]), frame_indices_0based=list(range(r["bscan_count"]))) for r in non128]
    if (local / "export_jobs.json").exists():
        preserve_job_destinations(jobs, json.loads((local / "export_jobs.json").read_text(encoding="utf-8")))
    write_csv(local / "raw_inventory_local.csv", inventory)
    write_csv(local / "pairing_report_local.csv", [dict(p, oct_absolute_path=r["oct_absolute_path"], dicom_absolute_path=r["flow_dicom_absolute_path"]) for p,r in zip(pairs,rows)])
    write_csv(local / "volume_manifest_local.csv", rows)
    write_json(local / "export_jobs.json", jobs)
    write_json(local / "pilot_selection.json", selection)
    write_json(local / "acquisition_headers.json", dict(oct=headers, dicom=dicom_meta))
    write_csv(public / "raw_inventory.csv", [{k:v for k,v in r.items() if k != "absolute_path"} for r in inventory])
    write_csv(public / "volume_manifest.csv", [{k:v for k,v in r.items() if "absolute_path" not in k} for r in rows])
    write_csv(public / "pairing_report.csv", pairs)
    write_json(results / "pilot_selection.json", selection)
    write_json(results / "naming_provenance.json", dict(relative_path=note.relative_to(raw_root).as_posix(), sha256=sha256(note), rules=["directory is diameter in um unless explicit mm", "numeric OCT filename is flow speed in mm/s"], calibration_reference_protocols=[list(p) for p in sorted(reference_protocols)]))
    write_csv(results / "input_sha256.csv", [dict(relative_path=r["relative_path_from_RAW_ROOT"], sha256=r["sha256"], size_bytes=r["size_bytes"]) for r in inventory])
    write_csv(results / "inventory_summary.csv", [dict(asset_type=k, count=v) for k,v in Counter(r["inferred_asset_type"] for r in inventory).items()])
    write_csv(results / "pairing_summary.csv", [dict(pairing_status=k, count=sum(r["pairing_status"] == k for r in rows)) for k in ("complete", "missing", "ambiguous")])
    paired = {r["flow_dicom_relative_path"] for r in rows if r["pairing_status"] == "complete"}
    write_json(results / "intake_validation.json", dict(oct_count=len(rows), flow_dicom_count=len(flow_paths), orphan_flow_dicoms=sorted(set(flow_paths)-paired), independent_volume_count=sum(r["pairing_status"] == "complete" for r in rows), eligible_new_volumes=len(non128), diameters=sorted({r["diameter_um"] for r in rows if r["diameter_um"] is not None}), flows=sorted({r["flow_speed_mm_s"] for r in rows if r["flow_speed_mm_s"] is not None}), blockers=[r["oct_relative_path"] for r in rows if r["acquisition_metadata_status"] != "complete"]))
    return selection


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-root", type=Path, required=True)
    parser.add_argument("--naming-note", type=Path, required=True)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    print(json.dumps(prepare(args.raw_root, args.repo, args.naming_note)))
