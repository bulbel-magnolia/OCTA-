# Diameter-stage input registry

This registry contains one row per paired acquisition. B-scans are slow-axis spatial positions nested within scan volumes; they are not independent experimental repeats.

- `raw_inventory.csv`: all recursively discovered OCT and DICOM assets, with actual byte-read SHA-256, sizes, timestamps and header checks.
- `pairing_report.csv`: unique exact companion or unique stem pairing, including explicit missing/ambiguous states.
- `volume_manifest.csv`: physical diameter, flow speed, acquisition metadata, stable scan ID and calibration provenance.

The raw-root experiment note documents numeric folder names as diameter in um and numeric OCT filenames as flow in mm/s. Explicit `0.5mm` converts to 500 um. Its digest and parsing rules are in `results/diameter_intake_v1/naming_provenance.json`. No image appearance is used to infer conditions.

Within a diameter × flow condition, lexical OCT relative-path order assigns V01, V02, etc. An existing registry is checked on rerun: changed hashes or changed IDs are refused. Inserting an earlier lexical acquisition requires a new reviewed intake revision. The pilot selection is likewise frozen independently of tracking quality.

The D128 reference rows retain their raw inventory. `legacy_scan_id` is populated only when OCT and Flow DICOM content hashes match the existing formal input registry. D128 historical outputs are not regenerated.

All local absolute paths and MATLAB jobs reside under ignored `outputs/local_run_inputs/diameter_stage_v1/`. Local per-frame MAT files are never committed. `export_jobs.json` records the actual local interim destination; it may use a separate disk when the repository disk cannot hold the full export. No raw OCT or DICOM is copied into this repository.

Reproduce intake from the repository root:

```powershell
python scripts/prepare_diameter_stage_inputs.py --raw-root $rawRoot --naming-note $experimentNamingNote
python scripts/process_diameter_stage.py pilot --matlab $matlabExecutable
python scripts/process_diameter_stage.py pilot-full --matlab $matlabExecutable
python scripts/process_diameter_stage.py batch --matlab $matlabExecutable --interim-root $derivedOutputRoot
```

Pilot hard gates precede whole-volume export; validated pilot whole-volume export precedes the remaining batch. Frozen exporters refuse source/frame identity mismatches, and the wrapper validates resumed files before reuse. Every expected exported MAT is read and checked. A single acquisition failure is recorded without selecting a different pilot or discarding a poor-looking volume.

These inputs prepare later raw-SV, no-background rectangular-tail analysis. This stage does not calculate Diameter–RI_tail results or statistical comparisons.
