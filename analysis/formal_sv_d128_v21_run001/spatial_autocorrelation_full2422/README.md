# Spatial autocorrelation — formal SV d128 run001

This analysis quantifies **within-volume spatial dependence** of the four frozen scalar metrics. It does not test a flow effect and does not treat B-scans as independent experimental replicates.

## Method

- Input: all 2422 Task-1 formal valid frames.
- Frame index 0–499 is retained as the spatial coordinate. Invalid positions remain missing; no interpolation or zero fill is used.
- ACF is evaluated from lag 0 to 200 B-scan positions.
- The frozen manifest contains no slow-axis spacing/position in µm, so correlation lengths are reported in **B-scan-position units only**.
- `descriptive_spatial_ess` is a volume-internal information-count approximation. It is **not** an independent biological/experimental sample size.

## Per-volume correlation scales

| scan | metric | ACF lag1 | ACF lag5 | ACF lag10 | 1/e lag | 0.1 lag | first <=0 lag | tau+ | spatial ESS* |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| flow01 | Q_V | 0.581 | 0.117 | 0.010 | 2.27 | 5.36 | 33 | 7.82 | 62.1 |
| flow01 | Sbar_V | 0.417 | 0.082 | -0.028 | 1.32 | 4.66 | 7 | 3.26 | 149.3 |
| flow01 | Q_T | 0.491 | 0.168 | 0.134 | 2.52 | 12.70 | 33 | 10.36 | 46.9 |
| flow01 | RI_tail | 0.395 | 0.253 | 0.126 | 1.30 | 19.91 | 78 | 17.69 | 27.5 |
| flow03 | Q_V | 0.536 | 0.185 | 0.030 | 2.60 | 7.92 | 49 | 14.95 | 31.3 |
| flow03 | Sbar_V | 0.428 | 0.149 | 0.075 | 1.40 | 9.68 | 35 | 10.46 | 44.8 |
| flow03 | Q_T | 0.510 | 0.188 | -0.038 | 1.95 | 6.68 | 10 | 5.01 | 93.4 |
| flow03 | RI_tail | 0.449 | 0.251 | 0.107 | 1.78 | 17.75 | 32 | 10.40 | 45.0 |
| flow05 | Q_V | 0.630 | 0.252 | 0.154 | 2.73 | 29.95 | 61 | 17.91 | 27.4 |
| flow05 | Sbar_V | 0.510 | 0.229 | 0.114 | 2.09 | 6.66 | 53 | 13.43 | 36.6 |
| flow05 | Q_T | 0.515 | 0.187 | 0.146 | 2.70 | 7.73 | 14 | 6.29 | 78.1 |
| flow05 | RI_tail | 0.415 | 0.237 | 0.205 | 1.43 | 35.33 | 77 | 21.20 | 23.2 |
| flow07 | Q_V | 0.533 | 0.155 | 0.074 | 2.57 | 7.21 | 23 | 8.07 | 61.1 |
| flow07 | Sbar_V | 0.462 | 0.077 | 0.014 | 1.68 | 4.82 | 11 | 3.82 | 129.1 |
| flow07 | Q_T | 0.539 | 0.195 | 0.091 | 2.44 | 9.84 | 23 | 7.80 | 63.2 |
| flow07 | RI_tail | 0.416 | 0.191 | 0.093 | 1.79 | 8.88 | 23 | 7.99 | 61.7 |
| flow10 | Q_V | 0.503 | 0.165 | 0.112 | 2.30 | 6.44 | 37 | 9.42 | 51.4 |
| flow10 | Sbar_V | 0.448 | 0.104 | 0.030 | 2.01 | 5.05 | 8 | 3.95 | 122.5 |
| flow10 | Q_T | 0.523 | 0.228 | 0.013 | 2.91 | 6.70 | 11 | 5.43 | 89.1 |
| flow10 | RI_tail | 0.309 | 0.227 | 0.069 | 0.92 | 6.69 | 17 | 5.52 | 87.7 |

`tau+` = initial-positive integrated autocorrelation span. `spatial ESS*` = n_valid / tau+ and is descriptive only.

## Cross-volume overview

| metric | median ACF lag1 (range) | median 1/e lag (range) | median first <=0 lag (range) | median tau+ (range) |
|---|---:|---:|---:|---:|
| Q_V | 0.54 (0.50–0.63) | 2.57 (2.27–2.73) | 37.00 (23.00–61.00) | 9.42 (7.82–17.91) |
| Sbar_V | 0.45 (0.42–0.51) | 1.68 (1.32–2.09) | 11.00 (7.00–53.00) | 3.95 (3.26–13.43) |
| Q_T | 0.51 (0.49–0.54) | 2.52 (1.95–2.91) | 14.00 (10.00–33.00) | 6.29 (5.01–10.36) |
| RI_tail | 0.41 (0.31–0.45) | 1.43 (0.92–1.79) | 32.00 (17.00–78.00) | 10.40 (5.52–21.20) |

## Interpretation boundary

These correlation scales describe how rapidly the measured metric changes along the slow-axis positions **inside one acquisition**. They can later inform spacing/blocking and subsampling design. They do not create independent replicates for flow inference; independent acquisitions remain the experimental unit for that purpose.

## Outputs

- `acf_by_scan_metric.csv`: complete lag-by-lag ACF table.
- `acf_selected_lags.csv`: compact selected-lag table.
- `correlation_length_summary.csv`: per-scan/per-metric correlation scales.
- `cross_scan_autocorrelation_summary.csv`: descriptive five-volume overview.
- `validation.json`, `provenance.json`, `input_sha256.csv`.
- `analyze_spatial_autocorrelation.py`: reproducible script.
