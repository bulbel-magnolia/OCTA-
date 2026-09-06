# 128 µm SV-OCTA: five-flow descriptive analysis

This is the first scientific flow-pattern analysis after completion of the frozen localization/quantification work, Relative Tail Intensity derivation, and fixed-128 geometry sensitivity QC.

## Statistical unit and scope

Each flow condition currently corresponds to one complete scan volume. The included B-scans are spatial positions within that volume, not independent vessels or independent experimental replicates. Results below describe volume-level spatial distributions and effect patterns. No ANOVA, t-test, p-value, or population-level flow-effect claim is made.

The primary endpoint values are the frozen run001 geometry (X4 centre + X1 apparent width). Fixed-128 values are shown only as geometry-sensitivity context and do not replace the primary analysis.

## Primary scan-level results

| flow (mm/s) | valid frames | Q_V median [Q1,Q3] ×10^12 | Sbar_V median [Q1,Q3] ×10^8 | Q_T median [Q1,Q3] ×10^12 | RI_tail median [Q1,Q3] | Sbar_T median [Q1,Q3] ×10^7 |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 486 | 2.254 [1.930, 2.625] | 1.265 [1.095, 1.484] | 1.363 [1.173, 1.601] | 0.1200 [0.1034, 0.1393] | 1.539 [1.328, 1.769] |
| 3 | 468 | 2.290 [1.940, 2.708] | 1.318 [1.085, 1.508] | 1.322 [1.122, 1.548] | 0.1147 [0.0950, 0.1386] | 1.495 [1.308, 1.695] |
| 5 | 491 | 2.338 [1.893, 2.762] | 1.314 [1.073, 1.567] | 1.321 [1.112, 1.539] | 0.1140 [0.0955, 0.1377] | 1.495 [1.270, 1.743] |
| 7 | 493 | 1.674 [1.401, 1.998] | 1.008 [0.842, 1.208] | 1.273 [1.068, 1.500] | 0.1533 [0.1253, 0.1820] | 1.544 [1.285, 1.847] |
| 10 | 484 | 2.336 [1.991, 2.694] | 1.319 [1.125, 1.572] | 1.372 [1.163, 1.591] | 0.1170 [0.1006, 0.1400] | 1.579 [1.364, 1.841] |

## Descriptive flow pattern

Using the 1 mm/s scan median only as a descriptive reference:

| flow (mm/s) | Δ median Q_V | Δ median Sbar_V | Δ median Q_T | Δ median RI_tail |
|---:|---:|---:|---:|---:|
| 1 | +0.00% | +0.00% | +0.00% | +0.00% |
| 3 | +1.60% | +4.15% | -3.04% | -4.46% |
| 5 | +3.72% | +3.85% | -3.10% | -5.05% |
| 7 | -25.75% | -20.31% | -6.63% | +27.71% |
| 10 | +3.64% | +4.23% | +0.62% | -2.50% |

### Geometry-QC context

The fixed-128 sensitivity changed absolute magnitudes, especially Q_T, but did not remove the distinct flow07 relative-intensity pattern. RI_tail medians under fixed-128 geometry are:

| flow (mm/s) | primary RI_tail median | fixed-128 RI_tail median | fixed128 vs primary |
|---:|---:|---:|---:|
| 1 | 0.1200 | 0.1078 | -10.18% |
| 3 | 0.1147 | 0.1044 | -9.01% |
| 5 | 0.1140 | 0.1000 | -12.24% |
| 7 | 0.1533 | 0.1431 | -6.66% |
| 10 | 0.1170 | 0.1059 | -9.53% |

## Relative intensity with depth RI(r)

For cross-volume description, every valid frame is sampled at fixed physical depth targets from 0 to 500 µm in 10 µm steps. For each target, the nearest original axial profile sample is selected and retained only within half one axial pixel (3.35 µm). This is nearest-sample reporting: **no interpolation and no curve fitting**. Negative background-corrected RI(r) values are retained.

Observed maximum nearest-sample depth mismatch: **3.350 µm**.

Selected depths are shown below; `ri_depth_summary_10um.csv` contains all 51 depth targets.

| flow (mm/s) | depth (µm) | RI(r) median [Q1,Q3] | negative-frame fraction |
|---:|---:|---:|---:|
| 1 | 0 | 0.5122 [0.3334, 0.7121] | 0.010 |
| 1 | 50 | 0.3400 [0.2104, 0.4857] | 0.027 |
| 1 | 100 | 0.1865 [0.1035, 0.2909] | 0.076 |
| 1 | 200 | 0.0870 [0.0170, 0.1679] | 0.200 |
| 1 | 300 | 0.0345 [-0.0249, 0.0933] | 0.360 |
| 1 | 400 | 0.0109 [-0.0514, 0.0677] | 0.442 |
| 1 | 500 | 0.0058 [-0.0437, 0.0559] | 0.463 |
| 3 | 0 | 0.4788 [0.2968, 0.7391] | 0.013 |
| 3 | 50 | 0.2953 [0.1781, 0.4614] | 0.019 |
| 3 | 100 | 0.1845 [0.0886, 0.2870] | 0.079 |
| 3 | 200 | 0.0891 [0.0190, 0.1575] | 0.184 |
| 3 | 300 | 0.0373 [-0.0141, 0.1005] | 0.318 |
| 3 | 400 | 0.0240 [-0.0271, 0.0678] | 0.357 |
| 3 | 500 | 0.0103 [-0.0341, 0.0571] | 0.447 |
| 5 | 0 | 0.4631 [0.3014, 0.6902] | 0.010 |
| 5 | 50 | 0.2985 [0.1811, 0.4359] | 0.029 |
| 5 | 100 | 0.1802 [0.0915, 0.2751] | 0.073 |
| 5 | 200 | 0.0706 [0.0096, 0.1497] | 0.222 |
| 5 | 300 | 0.0391 [-0.0149, 0.0987] | 0.326 |
| 5 | 400 | 0.0142 [-0.0389, 0.0693] | 0.436 |
| 5 | 500 | 0.0064 [-0.0414, 0.0570] | 0.473 |
| 7 | 0 | 0.6290 [0.4184, 0.9310] | 0.002 |
| 7 | 50 | 0.3921 [0.2248, 0.6350] | 0.032 |
| 7 | 100 | 0.2442 [0.1127, 0.3844] | 0.083 |
| 7 | 200 | 0.1022 [0.0277, 0.2033] | 0.197 |
| 7 | 300 | 0.0408 [-0.0260, 0.1171] | 0.337 |
| 7 | 400 | 0.0204 [-0.0515, 0.0871] | 0.422 |
| 7 | 500 | 0.0130 [-0.0561, 0.0836] | 0.448 |
| 10 | 0 | 0.4852 [0.3305, 0.7144] | 0.010 |
| 10 | 50 | 0.2856 [0.1909, 0.4157] | 0.019 |
| 10 | 100 | 0.1800 [0.0837, 0.2849] | 0.079 |
| 10 | 200 | 0.0950 [0.0243, 0.1582] | 0.174 |
| 10 | 300 | 0.0402 [-0.0106, 0.0999] | 0.300 |
| 10 | 400 | 0.0146 [-0.0328, 0.0678] | 0.411 |
| 10 | 500 | 0.0122 [-0.0392, 0.0521] | 0.428 |

## Interpretation

Across the current five volumes, Q_V and Sbar_V are similar for 1, 3, 5 and 10 mm/s and lower in the 7 mm/s volume. Q_T varies less across the five volumes. Sbar_T is also comparatively stable, so the elevated RI_tail in flow07 is driven primarily by the lower vessel-reference intensity Sbar_V rather than by a proportional rise in absolute tail signal. The fixed-128 geometry QC preserves this RI_tail ordering, so the pattern is not explained solely by the narrower X1 apparent width in flow07.

These statements describe the present scan volumes. Independent acquisitions are required before assigning the observed volume-to-volume pattern to a general flow-speed effect.

RI(r) is reported as a signed relative signal at fixed physical depths. No detection depth is inferred because a matched blank is unavailable; detection depth remains `not_evaluated / NA`.

## Outputs

- `scan_metrics_summary.csv`: spatial descriptive distributions for Q_V, Sbar_V, Q_T, RI_tail and supporting Sbar_T/R_B.
- `flow_metric_pattern.csv`: scan medians and descriptive changes relative to the 1 mm/s scan.
- `geometry_qc_context.csv`: primary run001 vs fixed-128 scan medians; QC only.
- `ri_depth_anchor_framewise.csv.gz`: nearest original RI(r) sample for every valid frame × 51 physical depth targets.
- `ri_depth_summary_10um.csv`: scan-wise RI(r) median/IQR and signed-value counts at all depth targets.
- `ri_depth_selected.csv`: 0/50/100/200/300/400/500 µm subset for compact review.
- `input_sha256.csv`, `validation.json`, `provenance.json`, and this reproducible script.

## Provenance

- Task 1 / primary derived metric source SHA: `f71e8df806ffca56eacb2df5167310108536a026`.
- Task 2 geometry-QC result SHA: `02cc44a2bb2f5821e444d631f9a96bcaeb78e59e`.
- Workflow-trigger SHA: `fcaf2c944f1b46f0fe40ceaedb1c6d0a3f2f8a80`.
