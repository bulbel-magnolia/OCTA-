# X1 midpoint vs X4 full-volume stability audit

This is a localization-only descriptive audit. It does **not** modify `results/formal_sv_d128_v21_full2500_run001`, the frozen formal geometry, or any Q_vessel/Q_tail/R endpoint.

## Definitions

- X1 midpoint: `x1_local_geometry_px`, the midpoint of the selected contiguous local-body run.
- X4: `x4_centroid_isolated_jump_corrected_px`, i.e. robust X2 centroid after isolated 3-frame jump correction.
- Signed difference: `delta_px = X4 - X1`; positive values mean X4 is to the larger-x side of X1.
- Lateral scale: 12.7 um/A-line; nominal vessel diameter: 128.0 um.
- No B-scan is treated as an independent biological replicate; results are descriptive across spatial frames.

## Whole 2500-frame result

- Paired X1/X4 coordinates: 2453/2500 (98.12%).
- Signed X4-X1: mean 0.0515 px; median 0.0311 px.
- Absolute difference: median 0.5608 px (7.12 um); p95 1.6389 px (20.81 um); maximum 3.6333 px (46.14 um).
- Within 0.5 px: 45.17%; within 1.0 px: 78.31%; >2 px: 1.47%.
- X1 fallback frames: 0; X4 isolated-jump corrections: 9.

## Assessable subset

- Frames: 2422.
- Signed median X4-X1: 0.0312 px.
- Absolute difference: median 0.5577 px; p95 1.5991 px; maximum 3.3342 px.
- Within 0.5 px: 45.38%; within 1.0 px: 78.78%.

## Adjacent-frame trajectory stability

| Scan | X1 median step (px) | X1 p95 step | X4 median step (px) | X4 p95 step |
|---|---:|---:|---:|---:|
| flow01 | 1.0000 | 2.1000 | 0.5566 | 1.6366 |
| flow03 | 1.0000 | 2.2000 | 0.5212 | 1.6780 |
| flow05 | 1.0000 | 2.5000 | 0.5374 | 1.6924 |
| flow07 | 1.0000 | 2.5000 | 0.5767 | 1.8110 |
| flow10 | 1.0000 | 2.5000 | 0.5465 | 1.6833 |

## Largest discrepancy

Largest |X4-X1| occurs at `flow03` frame 431: X1=246.5000 px, X4=242.8667 px, delta=-3.6333 px (-46.14 um), assessability=`uncertain`, X1_fallback=False, X4_jump_corrected=False.

See `x1_x4_summary.csv` for flow-stratified values, `x1_x4_extremes_top50.csv` for the largest discrepancies, and `x1_x4_framewise.csv` for all 2500 frames.

## Provenance

Workflow source SHA before generated commit: `4ce4d334332674b54e2368a96c035605beb95a3f`.
The five input tracking CSV SHA-256 values are recorded in `audit_metadata.csv`.
