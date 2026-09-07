# SV diameter intake and raw-export validation

Local validation: **COMPLETE_LOCAL_VALIDATION**.

Discovered 23 OCT and 23 Flow DICOM assets; all 23 acquisitions have unique pairs. Five D128 references were inventoried and hash-linked without rerunning their formal results.

New acquisitions validated: 18/18. Exported and read-validated MAT frames: 9000. Each acquisition is one experimental volume; the 500 B-scans per volume are spatial samples.

## Pilot

The predeclared rule selected `D235_F01_V01` (235 um, 1 mm/s) before viewing tracking or intensity outputs. Decoded DICOM dimensions are 500 × 351 × 500 (frame, depth, A-line).

Frozen v2.1 z localization: 492/500; source/geometry/500-um-window QC: 448/500. Direct accepted candidates: 319; short-gap fills: 173; missing segments: 1; relocks: 1.

B-scan 249 (MATLAB 250) reconstructed successfully: sv_raw shape 351 × 500, finite fraction 1, source/frame identity verified, population variance denominator N. All six hard gates passed before the 500-frame pilot export. Every actual B-scan is exported regardless of localization validity.

The three QC PNGs select first, middle and last geometry-valid pilot frames (0, 249, 497). Overlays use X4, X1 apparent width, true physical axial diameter, a 500 um rectangle and zero guard. They show localization geometry, not intensity outcomes.

## Validation and scope

Python: 102 passed (43 added, 59 existing); MATLAB: 2 test functions passed. Frozen generic source/config/MATLAB diff: empty. All 163 protected D128 analysis files retain their before/after SHA-256. Raw files were rehashed after processing.

The parameterization audit distinguishes generic diameter-aware geometry from D128-specific collection and historical analysis harnesses. No existing production algorithm or frozen setting was changed. No Diameter–RI_tail comparison, significance analysis or formal tail-intensity computation was performed.

## Files

`volume_table.csv` lists all acquisitions; `processing_summary.csv` and `per_volume_tracking_summary.csv` cover new-volume processing; `tracking/*_summary.json` contains localization distributions, missing segments and relock frames; `exports/*_frames.csv` records every MAT hash, size and identity check; `validation.json`, `provenance.json`, `raw_identity_validation.json` and `d128_regression_protection.json` hold audit evidence.

Raw assets and full MAT collections remain local. The local export jobs preserve actual storage destinations, including an external derived-output disk needed for the complete collection. Public tables contain relative paths only. Exact pushed commit identity is given in the final task report and local Git receipt.
