# Systematic subsampling stability — formal SV d128 run001

This analysis tests whether spatially uniform subsets reproduce the **full-volume median** of the frozen valid-frame endpoints. It does not change run001, interpolate invalid frames, or treat B-scans as independent biological replicates.

## Sampling design

For each stride `s`, every phase `0..s-1` is evaluated: a phase keeps formally valid frames whose original `frame_index_0based % s == phase`. Thus no favorable starting offset is selected. Tested strides correspond approximately to 250, 100, 50, 25, 20 and 10 planned positions per 500-frame volume.

## Primary endpoint stability

| stride | approx positions / 500 | Q_T P95 abs error | Q_T max abs error | RI_tail P95 abs error | RI_tail max abs error | joint max |
|---:|---:|---:|---:|---:|---:|---:|
| 2 | 250.0 | 1.75% | 1.88% | 2.40% | 2.49% | 2.49% |
| 5 | 100.0 | 3.97% | 4.35% | 4.75% | 7.19% | 7.19% |
| 10 | 50.0 | 5.13% | 6.57% | 6.13% | 12.31% | 12.31% |
| 20 | 25.0 | 9.99% | 14.27% | 12.12% | 14.88% | 14.88% |
| 25 | 20.0 | 13.05% | 18.93% | 14.97% | 17.82% | 18.93% |
| 50 | 10.0 | 16.26% | 26.89% | 16.80% | 30.13% | 30.13% |

The error is the percentage difference between a systematic subset median and that scan's full-valid-volume median. P95 and maximum summarize **all five scans and all phase offsets**.

## Interpretation boundary

These results quantify spatial subsampling efficiency for one already-acquired volume. They do not estimate the number of independent acquisitions required per flow condition. The latter requires between-volume variance from true independent repeats.

`planning_criteria_candidates.csv` reports which tested grids satisfy 1%, 2%, 5% or 10% absolute-error criteria for **every** Q_T and RI_tail scan-phase combination. Those thresholds are planning aids, not significance or validation thresholds.

## Outputs

- `full_volume_reference.csv`: full-valid-volume medians and quartiles.
- `systematic_phase_results.csv`: every scan × stride × phase × metric result.
- `subsampling_summary_by_scan.csv`: phase robustness within each scan.
- `cross_scan_subsampling_summary.csv`: pooled descriptive phase-error summaries across scans.
- `primary_endpoint_planning_matrix.csv`: compact Q_T / RI_tail planning view.
- `planning_criteria_candidates.csv`: sparsest tested grid meeting each descriptive all-phase error criterion.
- `phase_coverage_validation.csv`, `validation.json`, `provenance.json`, `input_sha256.csv`.
- `analyze_subsampling_stability.py`: reproducible script.
