# SV Diameter Formal Quantification v1

RI_tail = mean(raw SV in 500um tail ROI) / mean(raw SV in vessel ROI).
RI(r) = raw SV lateral mean at nearest native row / raw vessel mean.

Input: linear sv_raw = var(abs(E), 1, 3), denominator N. No background subtraction, log, normalization, gain, clipping or positive truncation.

Geometry: frozen continuity-first v2.1 X4 centre, X1 apparent lateral span, actual manifest diameter for source axial height; dx=12.7 and dz=6.7 um. Calibration is inherited from the project acquisition protocol.
Source: fractional ellipse (supersample 16). Tail: fractional rectangle with the same lateral span, unrounded physical bottom, guard 0, windows 100/200/300/500 um; primary alias is exactly the 500um value.
RI(r): 51 targets at 0..500 um in 10um steps, nearest native pixel with positive tail overlap (including fractional boundary pixels). No interpolation; maximum error <= dz/2 + 1e-9 um.

Entry: geometry_qc_valid only. All 9000 input rows retained; 8509 valid and 491 invalid. Invalid metric fields are NA, with reasons preserved.

| scan_id | nominal | valid | invalid |
|---|---:|---:|---:|
| D235_F01_V01 | 500 | 448 | 52 |
| D235_F03_V01 | 500 | 448 | 52 |
| D235_F05_V01 | 500 | 424 | 76 |
| D235_F07_V01 | 500 | 444 | 56 |
| D235_F10_V01 | 500 | 461 | 39 |
| D285_F01_V01 | 500 | 453 | 47 |
| D285_F03_V01 | 500 | 474 | 26 |
| D285_F05_V01 | 500 | 492 | 8 |
| D285_F07_V01 | 500 | 458 | 42 |
| D285_F10_V01 | 500 | 494 | 6 |
| D500_F01_V01 | 500 | 500 | 0 |
| D500_F02_V01 | 500 | 493 | 7 |
| D500_F03_V01 | 500 | 490 | 10 |
| D500_F05_V01 | 500 | 471 | 29 |
| D500_F07_V01 | 500 | 477 | 23 |
| D500_F09_V01 | 500 | 482 | 18 |
| D500_F10_V01 | 500 | 500 | 0 |
| D500_F12_V01 | 500 | 500 | 0 |

Validation: every MAT SHA and embedded identity matched intake; all per-volume QC counts matched; full 2422-frame D128 replay passed. Direct/profile integrals and analytical rectangle areas agree within the frozen numerical tolerances. See validation.json and d128_nobg_regression.json.

Files: framewise_all.csv.gz (9000 rows); framewise_primary.csv (8509); depth_anchor_framewise.csv.gz (433959); volume_summary.csv (18); depth_summary_10um.csv (918); depth_selected.csv (126).
Framewise fields retain source/tail integrals, fractional areas, raw means, four ratios, geometry and input MAT hashes. Depth fields retain target/native coordinates, fractional overlap and raw row signal. Summary tables use within-volume Q1/median/Q3 and, where specified, mean/sample SD (ddof=1).
formal_config.json freezes definitions; input_sha256.csv, tracking_identity.json, protected_sha256.json and provenance.json record identities. d128_reference.json links historical results without adding D128 rows to these outputs.

Experimental unit: scan volume. A B-scan is a slow-axis spatial sample, not an independent experiment. Aggregation is descriptive within each volume; no condition comparisons, regression, significance testing, fitting or scientific interpretation were performed.
