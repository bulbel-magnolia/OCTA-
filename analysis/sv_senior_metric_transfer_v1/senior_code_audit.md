# Senior source audit

Authoritative metric route: run_single_scan → formal_9scan_freeze; primary helper calls were recovered from actual source before analysis. The original source snapshots are unchanged.

| Source | Definition | Lines |
|---|---|---|
| `senior_source/python/tailq/metrics.py` | `_profile_geometry` | 52–83 |
| `senior_source/python/tailq/metrics.py` | `extract_profiles` | 86–130 |
| `senior_source/python/tailq/metrics.py` | `_normalizer` | 133–157 |
| `senior_source/python/tailq/metrics.py` | `_window_pixels` | 160–161 |
| `senior_source/python/tailq/metrics.py` | `_interval_values` | 346–386 |
| `senior_source/python/13_formal_9scan/formal_9scan_freeze.py` | `_interval_result` | 318–341 |
| `senior_source/python/13_formal_9scan/formal_9scan_freeze.py` | `_frame_metrics` | 344–492 |
| `senior_source/python/13_formal_9scan/formal_9scan_freeze.py` | `_scan_row` | 495–538 |
| `senior_source/python/02_track_vessel/tracking_core.py` | `load_flow_volume` | 229–332 |
| `senior_source/python/02_track_vessel/tracking_core.py` | `locate_lateral_track` | 446–643 |
| `senior_source/python/02_track_vessel/tracking_core.py` | `extract_profiles` | 646–713 |
| `senior_source/python/02_track_vessel/tracking_core.py` | `locate_peak_path` | 716–728 |
| `senior_source/python/02_track_vessel/tracking_core.py` | `_robust_trajectory` | 743–785 |
| `senior_source/python/02_track_vessel/tracking_core.py` | `track_one_alpha` | 788–1078 |
| `senior_source/python/10_xroi_refinement/develop_xroi_assessability.py` | `_isolated_jump_correct` | 133–146 |
| `senior_source/python/10_xroi_refinement/develop_xroi_assessability.py` | `_local_body_features` | 149–267 |
| `senior_source/python/10_xroi_refinement/develop_xroi_assessability.py` | `_add_assessability_score` | 270–296 |
| `senior_source/python/run_single_scan.py` | `_make_features` | 164–233 |
| `senior_source/python/run_single_scan.py` | `_make_metrics` | 236–360 |
| `senior_source/python/tailq/io.py` | `load_flow` | 88–154 |

The source-specific differences are recorded in senior_method_contract.json: central third wording vs 0.40 configuration, formal 1-pixel zero endpoint vs generic zero offset, exact formal division vs generic epsilon addition, and schema alias vessel_presence_prediction → assessability_class. The frozen contract is retained verbatim in validation/senior_method_contract_preexecution.json.

The post-execution loader audit found float32 loading for localization/features. This fidelity correction was applied to all 15 Stage B scans with no parameter tuning; all preliminary Stage B tables were retained under validation/stageB_before_loader_dtype_correction_DIAGNOSTIC_ONLY. Stage A and all metrics continue to use the unchanged original retained SV_raw values. The checkpoint is diagnostic only and never pooled into final outputs.

A synthetic zero-tail test initially zeroed [203,232), one pixel deeper than its stated D128 z=180 fixture [202,231). Correcting the fixture, without changing endpoint implementation, passed the test. The 90 fixed real-frame checks and independent scan medians also passed.

Initial script execution stopped before reading arrays due to Windows default text decoding; UTF-8 was made explicit. No scientific outputs existed at that point.

Volume coverage exposes metric-only assessability, full-process primary eligibility, and all-position mathematical validity separately. This column clarification changes no scan medians or support.
