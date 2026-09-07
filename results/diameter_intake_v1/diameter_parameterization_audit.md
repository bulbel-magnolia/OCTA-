# Diameter parameterization audit

Base: `63de971678de9d66527b1414e28c0d31acf3dc6b`, branch `feature/sv-rectangle-v1`.

Scope: read `src/svrecttail/`, `src/svrecttail/mentor/`, `matlab/`, frozen v2.1 tracking config, D128 collection harnesses and the D128 no-background analysis script. Search literal `128` and trace physical diameter through localization, source geometry and tail placement. No scientific intensity comparison was run.

## Generic diameter-aware path

| Component | Verified behavior |
|---|---|
| `cli.py`, `mentor_tracking.py` | `mentor-track` requires `--diameter-um`; it reaches `track_volume`, row metadata and the tracking bundle. The new runner supplies the acquisition's physical diameter. |
| `pipeline.py:_check_calibration` | A quantitative run requires manifest diameter to match that run's configured calibration diameter. This is an arbitrary-diameter consistency gate, not a fixed-128 assumption. Historical D128 run configs cannot be reused for mixed-diameter formal quantification; future formal runs need appropriately frozen diameter-specific configs. The new intake runner does not invoke formal quantification. |
| `mentor/tracking_core.py` | Lateral tracking, central/support widths, confidence windows and paired-edge depth use diameter divided by the frozen lateral/axial pixel pitch. No fixed 128 um vessel is imposed. |
| `mentor/xroi.py:local_body_features` | Reads row `diameter_um`; axial body slab is rounded diameter/6.7, expected lateral width is max(3, round(diameter/12.7)). X1 is a measured continuous image run. |
| `mentor_tracking.py:build_localization_from_tracking` | X4 sets lateral centre; X1 run width sets apparent lateral span; z-upper minus 0.5 sets the upper pixel edge; the passed physical diameter is stored in `VesselGeometry`. |
| `geometry.py:VesselGeometry` | Physical bottom is `z_top_edge_px + diameter_um/dz_um`, without rounding. The source ellipse uses physical axial diameter and measured lateral width. |
| `quantification.py` | Generic default tail gap is 0, length 500 um; rectangular tail starts at physical vessel bottom and shares the vessel's lateral span. Reviewed only; no quantification invoked in this task. |
| `matlab/load_and_reconstruct_common_oct.m` | Reconstruction reads spectral, A-line and repeat metadata from the actual OCT header; physical vessel diameter is not a reconstruction parameter. |
| `matlab/compute_sv_maps.m` | Formal raw signal remains `var(abs(IMG),1,3)`, population denominator N. No log, normalization, gain, clipping, positive truncation or background subtraction is introduced in raw SV. |

X1 width and physical diameter remain different quantities. Tests hold X1 at 11 pixels while changing physical diameter over 64, 128, 235, 256, 285 and 500 um. Separate tests check the diameter-dependent local-body prior at 64, 128 and 256 um.

## Intentional D128-specific code retained

`src/svrecttail/collection.py` contains literal `diameter_um=128` in its D128 formal collection and frozen-coordinate verification harness. Its run names, bridge gate, flow list and counts explicitly target `formal_sv_d128_v21_*_run001`. It is not the generic localization/geometry core and is not used for new-diameter intake or export.

Historical pilot run configs, threshold comparison and plotting scripts, synthetic example settings, and D128 analysis scripts also contain 128. They remain unchanged. In particular, `analysis/formal_sv_d128_v21_run001/no_background_relative_tail_v1_full2422/analyze_no_background_relative_tail.py` has a D128 release, expected counts, flow identifiers and fixed128 sensitivity; it was reviewed and never executed for new inputs.

The `128` occurrences in `matlab/config.ini` are display-range values, not physical-diameter logic in the common reconstruction path. Reconstruction reads the frozen phase threshold from this file; display ranges do not define `sv_raw`.

## Frozen configuration and geometry boundary

`config/tracking_config.continuity_first_v2_1.a020_n4.json` is byte-identical to base. Alpha 0.20, noise multiplier 4, continuity gates, jump thresholds and assessability rules are unchanged.

The tracker retains legacy diagnostic fields for its 2-pixel guard and 300 um window (`z_tail_start_px`, `window_300_fits`). These are not the formal 500 um/guard-0 geometry. The new compatibility runner constructs `VesselGeometry` through the existing adapter, then checks the exact physical bottom and 500 um rectangular tail without changing the tracking config. Public QC uses these formal edges. Localization-only background/excess calculations inside the frozen tracker do not alter the separately reconstructed raw SV signal.

Calibration is inherited at dx=12.7 um/A-line and dz=6.7 um/pixel only after each OCT acquisition signature and paired DICOM dimensions match D128. This is a documented project-calibration inheritance, not an independent calibration measurement; OCT scan-span fields are not silently interpreted as micrometres.

## Changes and regression protection

No implicit fixed-128 error was found in the generic core. No production core, MATLAB exporter or frozen config was changed. New scripts implement intake, stable identities, compatibility validation and export orchestration outside the D128 collection harness. New tests protect multi-diameter geometry and input/output identities.

`d128_regression_protection.json` records before/after SHA-256 and unchanged protected files. Python and MATLAB test results are recorded in `validation.json`. The D128 geometry test verifies the exact pre-existing mathematical construction at diameter 128. No old test expectations were edited.
