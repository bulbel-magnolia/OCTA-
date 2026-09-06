# Relative Tail Intensity — formal SV d128 v2.1 run001

This directory is a post-processing derivative of the frozen formal run. The frozen directory `results/formal_sv_d128_v21_full2500_run001` is not modified, and no `.oct` reconstruction is performed here.

## Formal inclusion

- Included frames: **2422** with frozen `valid == True`.
- Excluded frames remain excluded: **78**; no zero filling or interpolation is used.
- Linear background-corrected SV values are used as stored; negative values are retained.
- Flow/OMAG localization images are not used in the numerical Q_V/Q_T/RI integration.
- Detection depth remains outside this analysis (`not_evaluated / NA` without a matched blank).

## Definitions

- `q_vessel` = Q_V: background-corrected source-ellipse integral.
- `source_area_um2` = A_V: frozen fractional source area.
- `source_mean` = Sbar_V = Q_V / A_V.
- `q_tail` = Q_T: background-corrected 0–500 µm rectangular-tail integral.
- `tail_area_um2` = A_T: **frozen fractional tail area** from the formal geometry; no ideal-width approximation is substituted.
- `tail_mean` = Sbar_T = Q_T / A_T.
- `ratio_tail_to_vessel` = R_B = Q_T / Q_V: **relative tail burden**, not tail intensity.
- `ri_tail` = RI_tail = (Q_T/A_T)/(Q_V/A_V): **Relative Tail Intensity**.
- `T` = T(r) = V(r) − B(r), retaining signed background-corrected depth signal.
- `ri_r` = RI(r) = T(r)/Sbar_V.

The depth-profile outputs preserve every available formal profile row belonging to a valid frame and keep `tail_z_fraction`. No cross-frame depth interpolation or curve fitting is introduced in Task 1.

## Scan-level descriptive summary

Frames are spatial samples within one scan volume per flow condition. The table is descriptive; the B-scans are **not** treated as independent biological replicates and no p-values are computed.

| scan | flow (mm/s) | valid frames | Q_V median [Q1,Q3] | Sbar_V median [Q1,Q3] | Q_T median [Q1,Q3] | Sbar_T median [Q1,Q3] | R_B median [Q1,Q3] | RI_tail median [Q1,Q3] |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| flow01 | 1 | 486 | 2.25421e+12 [1.92982e+12, 2.62507e+12] | 1.26523e+08 [1.09538e+08, 1.48365e+08] | 1.3631e+12 [1.17279e+12, 1.60143e+12] | 1.53908e+07 [1.32837e+07, 1.76914e+07] | 0.596952 [0.514004, 0.693034] | 0.120037 [0.103363, 0.139329] |
| flow03 | 3 | 468 | 2.29035e+12 [1.93976e+12, 2.70829e+12] | 1.31772e+08 [1.08461e+08, 1.50792e+08] | 1.32172e+12 [1.12172e+12, 1.54778e+12] | 1.49527e+07 [1.30827e+07, 1.69502e+07] | 0.570396 [0.472698, 0.689414] | 0.114687 [0.0950173, 0.138642] |
| flow05 | 5 | 491 | 2.33811e+12 [1.89315e+12, 2.7624e+12] | 1.3139e+08 [1.07311e+08, 1.56733e+08] | 1.32078e+12 [1.11228e+12, 1.53937e+12] | 1.49486e+07 [1.27032e+07, 1.74335e+07] | 0.566833 [0.474892, 0.685031] | 0.113971 [0.0954609, 0.137733] |
| flow07 | 7 | 493 | 1.67368e+12 [1.40145e+12, 1.99816e+12] | 1.00828e+08 [8.41865e+07, 1.2079e+08] | 1.27273e+12 [1.06777e+12, 1.50015e+12] | 1.54424e+07 [1.28539e+07, 1.84713e+07] | 0.762375 [0.623122, 0.905114] | 0.153293 [0.12531, 0.182027] |
| flow10 | 10 | 484 | 2.33626e+12 [1.99129e+12, 2.69409e+12] | 1.31879e+08 [1.12518e+08, 1.57202e+08] | 1.37154e+12 [1.16326e+12, 1.59053e+12] | 1.57947e+07 [1.36412e+07, 1.84083e+07] | 0.582027 [0.500445, 0.69624] | 0.117034 [0.100615, 0.139999] |

## Outputs

- `relative_tail_intensity_framewise.csv`: one row per formally valid frame.
- `relative_tail_intensity_scan_summary.csv`: spatial descriptive statistics by scan/flow.
- `profiles/*_relative_intensity.csv.gz`: valid-frame depth profiles retaining T(r), fractional tail coverage, and adding RI(r).
- `profile_index.csv`: row/frame accounting for every profile chunk.
- `input_sha256.csv`: SHA-256 hashes of frozen scalar/profile inputs used by this derivation.
- `validation.json`: identity checks, inclusion counts, profile coverage, and signed-value checks.
- `provenance.json`: frozen source, workflow source, definitions, and analysis restrictions.
- `derive_relative_tail_intensity.py`: reproducible derivation script.

## Provenance

- Frozen formal handoff SHA: `acef5eb9d5ac356f1acf11aee885895da79a22e3`.
- Workflow-trigger source SHA: `4717cc00a5cd944296c0ad6b26181b4e1c58bb7f`.
- Formal 2D-array release tag: `formal-sv-d128-v21-run001`.
