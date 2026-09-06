# Fixed-surface 15-frame pilot review bundle

This directory contains the GitHub review package for the 15-frame pilot run generated on 2026-09-06. It is a pilot localization and rectangular-tail quantification run, not a final detection dataset.

## Run summary

- Input: 15 previously reconstructed OMAG/SV frames from five flow-speed scans; raw `.oct` files remain outside Git.
- Fixed surface reference: zero-based axial coordinate `z = 176`.
- Automatic coarse localization: search the OMAG image near the expected vessel depth and find the global lateral peak.
- Local refinement: apply the current X1 local-geometry rule to determine the left edge, right edge, central A-line, and vessel top edge.
- Valid frames: 13/15.
- Frames requiring review: `flow07_Bscan000` and `flow10_Bscan499`. Their horizontal localization succeeded, while top-edge refinement used the physical-depth fallback because no sustained top edge was found.
- Manual adjustments: 0.
- Tail detection: not evaluated because matched blank profiles were not supplied. Thirteen valid frames are marked `not_run_missing_matched_blank_profiles`; the two review frames are marked `not_run_upstream_qc_failed`.

`run_complete.json` records base commit `f44f9c6` with a dirty worktree because the run was executed immediately before the localization changes were committed. The exact implementation used by this run was subsequently committed as `7e4ae29`.

## Physical calibration

The computation does not treat pixels as square:

- lateral scale: `dx = 12.7 um/A-line`;
- axial scale: `dz = 6.7 um/pixel`;
- pixel area used for physical integration: `dx * dz = 85.09 um^2`;
- nominal vessel diameter: `128 um`, corresponding to about `10.08` lateral columns and `19.10` axial rows;
- surface-to-vessel-top prior: `200 um * 1.12 / 6.7 um/pixel = 33.43` axial rows, so the predicted top is near `z = 209.43`.

The PNG panels use display-oriented axis scaling, so the ellipse can look stretched on screen. Geometry and area calculations use the separate physical scales above.

## Included artifacts

- `frame_results.csv`, `localization.csv`, and `scan_summary.csv`: frame-level coordinates, QC, ratios, and scan summaries.
- `profiles.csv` and `profiles/`: merged and per-frame vessel/background/tail profiles.
- `sensitivity_results.csv`: one-pixel localization and background-spacing sensitivity checks.
- `detection_results.csv` and `detection_bins.csv`: explicit not-evaluable detection status for this run.
- `qc/`: one six-panel QC image per frame.
- `flow03_Bscan000_surface_guided_localization.png`: one representative source image with the complete localized region.
- `flow03_Bscan000_explicit_boundaries.png`: representative image with central A-line and all four boundaries explicitly labelled.
- `surface_guided_localization_overview_15.png`: overview of all 15 frames.
- `run_config.json`, `manifest.csv`, `manifest_original.csv.bin`, `run_complete.json`, and `logs/`: frozen settings and execution metadata. The `.bin` file preserves the original manifest bytes used for the recorded run hash; `manifest.csv` is its line-ending-normalized, browser-readable copy.

The 15 per-frame intermediate MAT arrays total 100,758,235 bytes and duplicate reconstructable image-level data, so they are kept in the ignored local run directory instead of this review package. `arrays_sha256.csv` records each omitted array's size and SHA-256 digest for exact traceability.

Local absolute paths in `frame_results.csv` and `localization.csv` were replaced with manifest-relative source paths before publication. Numerical results and run metadata were otherwise preserved. The SHA-256 of `manifest_original.csv.bin` matches the `manifest_sha256` field in `run_complete.json`.
