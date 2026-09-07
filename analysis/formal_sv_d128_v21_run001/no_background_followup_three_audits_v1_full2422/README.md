# No-background SV follow-up: three descriptive audits v1

Primary definition remains `RI_tail = mean(raw SV in tail ROI) / mean(raw SV in vessel ROI)`. No background subtraction is used.
B-scans are spatial samples inside one scan volume. No p-values or independent-flow inference are reported.

## 1. Source-tail-ratio coupling

| scan | source vs RI Spearman | tail vs RI Spearman | source vs tail Spearman | source IQR/median | tail IQR/median | RI IQR/median |
|---|---:|---:|---:|---:|---:|---:|
| flow01 | -0.908 | -0.148 | 0.507 | 0.268 | 0.112 | 0.237 |
| flow03 | -0.913 | -0.035 | 0.389 | 0.280 | 0.112 | 0.270 |
| flow05 | -0.924 | -0.174 | 0.501 | 0.319 | 0.114 | 0.295 |
| flow07 | -0.887 | -0.089 | 0.494 | 0.313 | 0.134 | 0.267 |
| flow10 | -0.916 | -0.129 | 0.477 | 0.299 | 0.118 | 0.270 |

Flow07 volume-level context:

- vessel raw mean: flow07 vs flow01 **-17.53%**; vs median of other four scan medians **-20.09%**.
- tail raw mean: flow07 vs flow01 **+0.59%**; vs median of other four scan medians **+1.75%**.
- RI_tail: flow07 vs flow01 **+20.17%**; vs median of other four scan medians **+25.24%**.

The negative source-vs-RI association is expected partly from the ratio definition itself. It is reported as a coupling audit, not as a biological association.

## 2. Exact tail-window sensitivity

| window (µm) | flow01 | flow03 | flow05 | flow07 | flow10 | flow07 vs flow01 |
|---:|---:|---:|---:|---:|---:|---:|
| 100 | 0.4821 | 0.4707 | 0.4648 | **0.5733** | 0.4659 | **+18.92%** |
| 200 | 0.3884 | 0.3734 | 0.3703 | **0.4650** | 0.3740 | **+19.72%** |
| 300 | 0.3314 | 0.3171 | 0.3147 | **0.3980** | 0.3193 | **+20.09%** |
| 500 | 0.2689 | 0.2563 | 0.2566 | **0.3232** | 0.2595 | **+20.17%** |

All four windows use exact native-pixel integration from the vessel bottom; no 10-µm anchor approximation is used.

## 3. Five-segment within-volume spatial stability

| scan | RI segment min | RI segment max | range / full median | max |segment-full| |
|---|---:|---:|---:|---:|
| flow01 | 0.2512 | 0.2962 | 16.7% | 10.1% |
| flow03 | 0.2197 | 0.2802 | 23.6% | 14.3% |
| flow05 | 0.2316 | 0.3041 | 28.2% | 18.5% |
| flow07 | 0.2973 | 0.3584 | 18.9% | 10.9% |
| flow10 | 0.2227 | 0.2799 | 22.1% | 14.2% |

Flow07 RI relative to the median of the other four scans in the same 100-frame segment:

- frames 0-99: **+15.98%**, rank 1/5.
- frames 100-199: **+41.09%**, rank 1/5.
- frames 200-299: **+20.58%**, rank 1/5.
- frames 300-399: **+22.05%**, rank 1/5.
- frames 400-499: **+24.40%**, rank 1/5.

## Validation

- valid frames: **2422**
- exact-window raw-array packages verified: **25/25**
- valid-frame NPZ SHA verified: **2422**
- 500-µm RI replay max relative error: **1.304e-15**
- log identity `log(RI)=log(tail)-log(source)` max absolute numerical error: **3.997e-15**
- no background subtraction, interpolation, zero-fill, p-values, or B-scan pseudo-replication.

## Files

- `coupling_summary.csv`, `coupling_source_deciles.csv`, `coupling_source_quintiles.csv`, `flow07_coupling_context.csv`
- `window_framewise.csv.gz`, `window_scan_summary.csv`, `window_flow07_context.csv`
- `segment_summary.csv`, `segment_stability_by_scan.csv`, `segment_flow07_context.csv`
- `validation.json`, `provenance.json`, `input_sha256.csv`, `release_package_audit.csv`
