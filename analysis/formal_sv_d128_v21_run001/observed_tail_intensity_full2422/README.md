# Observed SV tail intensity — no background subtraction

Question: **How bright does the tail appear in the SV image relative to the vessel itself?**

`RI_tail_observed = (Q_T_observed/A_T)/(Q_V_observed/A_V)` and `RI_observed(r)=V(r)/source_mean_observed`. Both source and tail use linear `sv_raw` without background subtraction. The previous background-corrected endpoint is retained separately.

Frozen X4+X1 geometry, v2.1 z geometry, 500 µm tail, guard=0, fractional weights and `valid == True` are unchanged. No p-values are computed; B-scans are spatial positions, not independent biological replicates.

## Validation

- Valid frames: **2422**
- Release packages: **25/25**
- NPZ SHA checks: **2422**
- Frozen corrected-metric replay max relative error: **1.750e-14**
- Area replay max relative error: **2.139e-15**
- Max depth-anchor sampling error: **3.350 µm**

## Scan-level observed RI_tail

| scan | flow (mm/s) | median | Q1 | Q3 |
|---|---:|---:|---:|---:|
| flow01 | 1 | 0.268943 | 0.237327 | 0.301055 |
| flow03 | 3 | 0.256329 | 0.22535 | 0.294434 |
| flow05 | 5 | 0.256646 | 0.223266 | 0.299072 |
| flow07 | 7 | 0.32318 | 0.282615 | 0.368924 |
| flow10 | 10 | 0.259467 | 0.223592 | 0.293673 |

`RI_tail_observed=0.20` means the average signal actually visible in the tail rectangle is 20% of the average signal actually visible in the vessel ROI. It intentionally includes local SV background.

Files: `observed_tail_intensity_framewise.csv`, `observed_depth_anchor_framewise.csv.gz`, `observed_tail_intensity_scan_summary.csv`, `observed_depth_summary.csv`, plus validation/provenance/SHA audits.
