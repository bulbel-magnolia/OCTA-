# Final two audits for no-background SV RI_tail v1

Primary endpoint remains `RI_tail = mean(raw SV in tail ROI) / mean(raw SV in vessel ROI)`. No background subtraction is used.

## 1. Aggregation sensitivity

| aggregation method | flow07 | flow07 vs flow01 | flow07 vs other4 median | rank | max within-scan change vs primary |
|---|---:|---:|---:|---:|---:|
| mean_framewise_RI | 0.327819 | +20.46% | +25.02% | 1/5 | 2.53% |
| median_of_5_spatial_segment_medians | 0.320495 | +19.42% | +23.18% | 1/5 | 3.27% |
| pooled_area_weighted_raw_ratio | 0.314214 | +19.42% | +24.81% | 1/5 | 3.11% |
| primary_median_framewise_RI | 0.323180 | +20.17% | +25.24% | 1/5 | 0.00% |
| ratio_of_mean_raw_means | 0.312682 | +19.89% | +25.23% | 1/5 | 3.95% |
| ratio_of_median_raw_means | 0.326168 | +21.97% | +27.47% | 1/5 | 1.49% |

- Flow07 ranks first under **6/6** pre-specified aggregation methods.

## 2. Residual geometry/localization QC

### Flow07 geometry medians

| metric | flow07 | flow01 | other4 scan median | flow07 vs other4 |
|---|---:|---:|---:|---:|
| lateral_width_um | 165.1 | 177.8 | 177.8 | -7.14% |
| source_area_um2 | 16598.9 | 17875.5 | 17875 | -7.14% |
| z_top_edge_px | 199.5 | 198.5 | 195 | +2.31% |
| x4_px | 238.514 | 239.122 | 241.404 | -1.20% |

### RI_tail vs geometry: Spearman by scan

| scan | width | source area | z_top | X4 |
|---|---:|---:|---:|---:|
| flow01 | 0.239 | 0.244 | -0.210 | 0.190 |
| flow03 | 0.248 | 0.240 | -0.156 | 0.201 |
| flow05 | 0.180 | 0.167 | -0.347 | 0.358 |
| flow07 | 0.166 | 0.163 | -0.242 | 0.214 |
| flow10 | 0.235 | 0.237 | -0.209 | 0.203 |

### Restricted-subset sensitivity

| subset | flow07 RI | vs flow01 | vs other4 median | rank | vessel vs other4 | tail vs other4 |
|---|---:|---:|---:|---:|---:|---:|
| all_valid | 0.323180 | +20.17% | +25.24% | 1/5 | -20.09% | +1.75% |
| central80_width | 0.322607 | +19.69% | +25.15% | 1/5 | -19.52% | +1.92% |
| central80_z | 0.325885 | +21.54% | +26.10% | 1/5 | -20.16% | +1.64% |
| central80_z_and_width | 0.323834 | +20.99% | +25.85% | 1/5 | -20.00% | +1.57% |
| direct_and_central80_z_width | 0.322233 | +20.42% | +23.64% | 1/5 | -18.53% | +1.63% |
| direct_only | 0.321311 | +19.87% | +23.59% | 1/5 | -18.42% | +1.67% |

Localization-source composition and category-specific medians are in `localization_source_summary.csv`. Restrictions are QC sensitivity subsets only; they do not create independent replicates.

## Validation

- valid frames: **2422**
- scan counts: `{'flow01': 486, 'flow03': 468, 'flow05': 491, 'flow07': 493, 'flow10': 484}`
- no background subtraction, interpolation, p-values, or B-scan pseudo-replication.

## Files

- `aggregation_all_methods.csv`, `aggregation_summary.csv`, `aggregation_dispersion_by_scan.csv`
- `geometry_signal_correlations.csv`, `geometry_scan_summary.csv`, `geometry_flow07_context.csv`
- `localization_source_summary.csv`, `restricted_subset_scan_summary.csv`, `restricted_subset_flow07_context.csv`
- `validation.json`, `provenance.json`, `input_sha256.csv`
