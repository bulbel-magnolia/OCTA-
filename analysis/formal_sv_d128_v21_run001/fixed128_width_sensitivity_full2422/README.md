# Fixed 128 µm lateral-width geometry sensitivity — Task 2

This directory is a **geometry QC / sensitivity analysis**. Frozen run001 remains the primary analysis and is not modified.

## Geometry contrast

- Baseline: X4 centre + X1 apparent lateral width + frozen v2.1 z geometry.
- Variant: X4 centre + **fixed 128 µm lateral width** + the same frozen z geometry.
- Source vertical diameter remains 128 µm.
- Tail width follows the corresponding source width; tail length remains 500 µm; guard remains 0.
- Background rule remains skip 3 / take 5 per side and is recomputed from the **variant vessel edges**.
- `sv_raw`, signed background subtraction, 16×16 ellipse fractional area, calibration, and all other quantification rules are unchanged.

No parameter is chosen according to flow monotonicity or statistical significance. B-scans are spatial positions within a scan volume, not independent biological replicates; no p-values are reported.

## Input and replay validation

- Frozen-valid frames requested: **2422**.
- Fixed-128 valid frames: **2422**.
- Release packages verified: **25/25**.
- Valid-frame NPZ SHA-256 checks passed: **2422**.
- Baseline geometry replay max relative error across audited scalar quantities: **1.750e-14**.
- X4 vs baseline geometry-centre max difference: **2.842e-14 px**.

## Overall sensitivity

| metric | baseline median | fixed-128 median | median |Δ| (%) | P95 |Δ| (%) |
|---|---:|---:|---:|---:|
| Q_V | 2.18226e+12 | 2.00718e+12 | 7.439% | 14.991% |
| Q_T | 1.33454e+12 | 1.10822e+12 | 16.316% | 33.058% |
| R_B (relative tail burden) | 0.609904 | 0.547097 | 9.354% | 24.819% |
| RI_tail (Relative Tail Intensity) | 0.122646 | 0.110018 | 9.355% | 24.821% |

## Apparent-width context

| scan | flow (mm/s) | valid frames | X1 apparent width median [Q1,Q3] µm | fixed width (µm) |
|---|---:|---:|---:|---:|
| flow01 | 1 | 486 | 177.800 [165.100, 190.500] | 128.0 |
| flow03 | 3 | 468 | 177.800 [165.100, 190.500] | 128.0 |
| flow05 | 5 | 491 | 177.800 [165.100, 190.500] | 128.0 |
| flow07 | 7 | 493 | 165.100 [152.400, 177.800] | 128.0 |
| flow10 | 10 | 484 | 177.800 [165.100, 190.500] | 128.0 |

## Files

- `fixed128_framewise.csv`: paired baseline/fixed-128 values for every frozen-valid frame.
- `fixed128_sensitivity_summary.csv`: overall and per-scan descriptive sensitivity statistics.
- `apparent_width_summary.csv`: frozen X1 apparent-width distribution by scan.
- `release_package_audit.csv`: Release package and per-package array-verification accounting.
- `input_sha256.csv`: repository input hashes used by the audit.
- `validation.json`: inclusion, replay, array-hash, and variant-validity checks.
- `provenance.json`: frozen/input/workflow provenance and fixed geometry definition.
- `audit_fixed128_width.py`: reproducible script.

## Interpretation boundary

This comparison quantifies how much the frozen endpoints change when X1-dependent apparent width is replaced by the known 128 µm phantom diameter. It is not used to choose the version with a more favorable flow pattern. Scientific flow interpretation follows only after this QC is complete.

## Provenance

- Frozen formal handoff SHA: `acef5eb9d5ac356f1acf11aee885895da79a22e3`.
- Task 1 result SHA: `f71e8df806ffca56eacb2df5167310108536a026`.
- Workflow-trigger SHA: `598383ce873c4d250b7d99b8d6f74411a1e71867`.
- 2-D array Release: `formal-sv-d128-v21-run001`.
