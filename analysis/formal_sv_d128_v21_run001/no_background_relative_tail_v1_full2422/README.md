# No-background Relative Tail Intensity — primary SV reanalysis v1

From this analysis forward, **RI_tail** means the no-background raw-SV ratio:

`RI_tail = mean(raw SV in tail ROI) / mean(raw SV in vessel ROI)`

`RI(r) = mean(raw SV across the vessel lateral span at depth r) / mean(raw SV in vessel ROI)`

No local background is subtracted. Historical background-corrected outputs remain in their old directories only and are not used by this primary analysis.

## Five-flow primary result

| scan | flow | valid | RI_tail median [Q1,Q3] | change vs flow01 | vessel raw mean | tail raw mean |
|---|---:|---:|---:|---:|---:|---:|
| flow01 | 1 | 486 | 0.268943 [0.237327,0.301055] | +0.00% | 1.49475e+08 | 3.99738e+07 |
| flow03 | 3 | 468 | 0.256329 [0.225350,0.294434] | -4.69% | 1.55468e+08 | 3.9259e+07 |
| flow05 | 5 | 491 | 0.256646 [0.223266,0.299072] | -4.57% | 1.54205e+08 | 3.91762e+07 |
| flow07 | 7 | 493 | 0.323180 [0.282615,0.368924] | +20.17% | 1.23277e+08 | 4.0209e+07 |
| flow10 | 10 | 484 | 0.259467 [0.223592,0.293673] | -3.52% | 1.54339e+08 | 3.97747e+07 |

## Fixed-128 geometry stress test

- RI_tail framewise |Δ| median/P95: **12.08% / 22.03%**.
- Source raw mean |Δ| median/P95: **22.10% / 45.16%**.

## Spatial autocorrelation of new RI_tail

Across five volumes, median lag-1 ACF = **0.459**; median 1/e crossing = **1.72 frames**; median initial-positive correlation span = **10.67 frames**. This is spatial dependence, not biological n.

## Systematic subsampling stability of new RI_tail

| stride | ~positions/500 | median |error| | P95 |error| | max |error| |
|---:|---:|---:|---:|---:|
| 2 | 250 | 0.82% | 1.59% | 1.97% |
| 5 | 100 | 1.12% | 3.59% | 4.84% |
| 10 | 50 | 2.45% | 5.31% | 6.71% |
| 20 | 25 | 3.10% | 9.96% | 15.98% |
| 25 | 20 | 3.24% | 10.42% | 15.99% |
| 50 | 10 | 5.06% | 12.51% | 21.33% |

## Selected RI(r) medians

| flow | 0 | 50 | 100 | 200 | 300 | 400 | 500 µm |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0.5977 | 0.4465 | 0.3290 | 0.2373 | 0.1862 | 0.1659 | 0.1575 |
| 3 | 0.5723 | 0.4192 | 0.3185 | 0.2333 | 0.1850 | 0.1594 | 0.1522 |
| 5 | 0.5623 | 0.4261 | 0.3173 | 0.2213 | 0.1793 | 0.1613 | 0.1536 |
| 7 | 0.7158 | 0.5252 | 0.4046 | 0.2814 | 0.2207 | 0.1992 | 0.1927 |
| 10 | 0.5734 | 0.4160 | 0.3226 | 0.2300 | 0.1818 | 0.1570 | 0.1493 |

## Validation

- valid frames: **2422**
- raw-array packages verified: **25/25**
- per-frame NPZ SHA verified: **2422**
- baseline observed raw-SV replay max relative error: **3.013e-15**
- no interpolation, no zero-fill, no background subtraction, no p-values.
- one scan volume per flow: all flow comparisons remain descriptive until independent-volume replication.

## Files

- `framewise_primary.csv`
- `scan_flow_summary.csv`
- `depth_selected.csv` / `depth_full_10um.csv`
- `fixed128_observed_framewise.csv` / `fixed128_geometry_summary.csv`
- `acf_ri_tail.csv` / `acf_summary.csv`
- `subsampling_phase_results.csv` / `subsampling_summary.csv`
- `validation.json`, `provenance.json`, `input_sha256.csv`
