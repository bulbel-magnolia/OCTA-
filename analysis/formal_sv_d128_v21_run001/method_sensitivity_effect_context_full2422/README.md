# Observed flow contrasts vs method-sensitivity context

This audit puts the current **volume-level** flow contrasts on the same percent scale as three audited perturbations. It does **not** estimate a single measurement-error SD and does not establish a flow-speed effect.

- `background +2`: local method perturbation.
- `fixed-128 width`: geometry **stress test**, not a random uncertainty estimate.
- `systematic subsampling`: sampling-design approximation; 100-position (stride 5) is used as the local planning reference, while 250-position results are also retained.

## Method sensitivity across all valid frames / scan phases

| metric | background +2 median / P95 | fixed-128 median / P95 | ~100 positions median / P95 | ~250 positions median / P95 |
|---|---:|---:|---:|---:|
| Q_V | 1.05% / 3.59% | 7.44% / 14.99% | 1.89% / 4.36% | 0.49% / 1.04% |
| Sbar_V | 1.05% / 3.59% | 24.87% / 51.02% | 1.56% / 3.97% | 0.69% / 1.67% |
| Q_T | 5.35% / 16.54% | 16.32% / 33.06% | 0.88% / 3.97% | 0.58% / 1.75% |
| RI_tail | 4.93% / 15.86% | 9.35% / 24.82% | 1.85% / 4.75% | 1.17% / 2.40% |

## Current volume-level contrasts relative to 1 mm/s

| flow | metric | observed median change | local-method P95 reference | geometry-stress P95 | context |
|---:|---|---:|---:|---:|---|
| 3 | Q_T | -3.04% | 16.54% | 33.06% | within_local_method_sensitivity |
| 3 | Q_V | +1.60% | 4.36% | 14.99% | within_local_method_sensitivity |
| 3 | RI_tail | -4.46% | 15.86% | 24.82% | within_local_method_sensitivity |
| 3 | Sbar_V | +4.15% | 3.97% | 51.02% | above_local_but_within_geometry_stress |
| 5 | Q_T | -3.10% | 16.54% | 33.06% | within_local_method_sensitivity |
| 5 | Q_V | +3.72% | 4.36% | 14.99% | within_local_method_sensitivity |
| 5 | RI_tail | -5.05% | 15.86% | 24.82% | within_local_method_sensitivity |
| 5 | Sbar_V | +3.85% | 3.97% | 51.02% | within_local_method_sensitivity |
| 7 | Q_T | -6.63% | 16.54% | 33.06% | within_local_method_sensitivity |
| 7 | Q_V | -25.75% | 4.36% | 14.99% | above_local_and_geometry_stress |
| 7 | RI_tail | +27.71% | 15.86% | 24.82% | above_local_and_geometry_stress |
| 7 | Sbar_V | -20.31% | 3.97% | 51.02% | above_local_but_within_geometry_stress |
| 10 | Q_T | +0.62% | 16.54% | 33.06% | within_local_method_sensitivity |
| 10 | Q_V | +3.64% | 4.36% | 14.99% | within_local_method_sensitivity |
| 10 | RI_tail | -2.50% | 15.86% | 24.82% | within_local_method_sensitivity |
| 10 | Sbar_V | +4.23% | 3.97% | 51.02% | above_local_but_within_geometry_stress |

## Interpretation

The classification is deliberately descriptive. `within_local_method_sensitivity` means the observed absolute scan-median contrast does not exceed the larger P95 of background +2 and ~100-position subsampling. `above_local_but_within_geometry_stress` means it exceeds those local perturbations but remains within the P95 of the stronger fixed-128 geometry stress test. `above_local_and_geometry_stress` means it exceeds both reference scales.

None of these labels is a p-value or a causal flow-speed claim. Each flow still has one independent scan volume, so independent-volume replication is required.

## Outputs

- `method_sensitivity_summary.csv`
- `observed_effect_vs_method_sensitivity.csv`
- `flow07_context_summary.csv`
- `scan_median_reference.csv`
- `input_sha256.csv`
- `validation.json`
- `provenance.json`
- `analyze_method_sensitivity_effect_context.py`
