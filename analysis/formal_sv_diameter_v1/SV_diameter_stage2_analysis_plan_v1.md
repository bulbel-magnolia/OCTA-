# SV Diameter Stage 2 — scientific analysis plan v1

Date: 2026-09-08

Starting point: formal no-background multi-diameter quantification is complete. This stage analyzes the scientific structure of the SV results without changing the frozen endpoint or geometry.

## 1. Primary question

Determine how physical vessel diameter and flow condition are associated with the observed SV tail metrics, with particular emphasis on separating:

1. vessel/source raw SV;
2. absolute tail raw SV;
3. relative tail intensity `RI_tail = tail mean / vessel mean`;
4. depth-dependent `RI(r)`.

The main objective is to explain the non-monotonic diameter pattern now visible in RI_tail, especially the fall from 128→235→285 um and the return to higher RI at 500 um.

## 2. Experimental-unit boundary

Every diameter × flow cell currently contains one independent scan volume.

Therefore:

- B-scans remain slow-axis spatial samples, not independent experimental replicates;
- no B-scan-level t-test, ANOVA or p-value is permitted;
- all cross-condition conclusions remain scan-volume-level descriptive comparisons;
- inferential flow/diameter effects require independent-volume replication in a future acquisition stage.

## 3. Primary comparison grid

Primary common grid:

- diameters: 128, 235, 285, 500 um
- flows: 1, 3, 5, 7, 10 mm/s

This produces 20 matched condition volumes.

D500 flows 2, 9 and 12 mm/s are retained as secondary within-D500 descriptive extension data and are not used to redefine the common grid.

## 4. Stage 2A — common-grid descriptive matrix

Construct a single 4 × 5 matrix for each of the following volume-level endpoints:

- `source_mean_raw`
- `tail_mean_raw_100um`
- `tail_mean_raw_200um`
- `tail_mean_raw_300um`
- `tail_mean_raw_500um`
- `ri_tail_100um`
- `ri_tail_200um`
- `ri_tail_300um`
- `ri_tail_500um` / primary `ri_tail`

For each fixed flow, calculate descriptive diameter contrasts:

- 128→235
- 235→285
- 285→500
- 128→500

For each fixed diameter, summarize the variation across common flows.

Report both absolute differences and percent differences where meaningful. These are descriptive contrasts, not effect estimates with inferential significance.

## 5. Stage 2B — source–tail–ratio decomposition

For every common-grid volume, place the three quantities side by side:

`source_mean_raw`, `tail_mean_raw_500um`, `RI_tail`.

Primary questions:

- Is a change in RI accompanied by a change in absolute tail raw SV?
- Is it mainly driven by the vessel denominator?
- Do numerator and denominator move in opposite directions?

The D500 pattern receives a targeted audit because its current RI returns to approximately D128 levels while its absolute tail raw SV remains lower than D235/D285 and its vessel raw SV is markedly lower.

No mechanistic claim is made until geometry and source-signal diagnostics below are completed.

## 6. Stage 2C — depth-profile analysis

Use the frozen `RI(r)` data from 0 to 500 um.

Primary summaries:

- full median RI(r) curve by scan volume;
- selected depths: 0, 50, 100, 200, 300, 400, 500 um;
- comparison of 0–100, 100–200, 200–300 and 300–500 um behavior using the already-computed window metrics and depth profiles.

Questions:

- At which depths is the 128→235→285 separation largest?
- Does the difference converge toward a common deep raw-SV floor?
- Is the D500 high RI primarily a near-vessel phenomenon or a broadly elevated normalized floor?

Interpretation rule:

Deep non-zero RI(r) is not a detection depth or tail endpoint because the primary curve is not background-subtracted.

## 7. Stage 2D — slow-axis spatial robustness

For every volume, divide the nominal 500-frame slow axis into the same five fixed spatial segments:

- 0–99
- 100–199
- 200–299
- 300–399
- 400–499

Within each segment, report:

- valid-frame count;
- median source raw SV;
- median 500 um tail raw SV;
- median RI_tail.

Determine whether the diameter ranking is distributed across the volume or driven by one localized region.

Do not use the segments as independent biological replicates.

## 8. Stage 2E — geometry and QC robustness

Because geometry-valid fractions differ among the new diameter volumes, audit whether the main diameter pattern is explained by selection or geometry composition.

At minimum examine:

- geometry-valid fraction by volume;
- X1 apparent width distribution;
- X4 distribution;
- z_top distribution;
- source ROI area distribution;
- localization / candidate source composition where available;
- relation of these geometry variables to `source_mean_raw`, `tail_mean_raw_500um` and `RI_tail` within each volume.

Predefined restricted subsets should include, where data permit:

1. direct-candidate-supported frames only;
2. central 80% X1 apparent width within each volume;
3. central 80% z_top within each volume;
4. central 80% width + z_top;
5. direct-only + central geometry.

The primary question is whether the 128→235→285 decline and the D500 rebound remain visible after excluding geometry extremes and localization-composition differences.

## 9. Stage 2F — D500 denominator audit

The low D500 vessel raw SV is currently the most important unresolved feature.

First use committed framewise data to test whether the D500 source decrease is associated with:

- apparent X1 width;
- source ROI area;
- z_top;
- slow-axis location;
- geometry/QC composition;
- flow condition.

If the D500 source decrease persists across spatial segments and restricted geometry subsets, perform a second targeted local-array audit only if needed.

A targeted local-array audit may examine the axial distribution of raw SV inside the frozen physical source ellipse, for example comparing predefined relative axial subregions of the vessel. Such diagnostics must remain secondary and must not replace the frozen primary source definition.

Do not redesign the source ROI simply because D500 has a lower denominator.

## 10. Stage 2G — D500 extra-flow extension

After the 4 × 5 common grid is understood, add D500 flows 2, 9 and 12 mm/s as within-D500 descriptive extension points.

Use them to assess whether the D500 pattern across flow remains consistent with the common-grid observations.

Do not use these extra flows to create an unbalanced cross-diameter model at this stage.

## 11. Analyses deliberately deferred

Do not yet perform:

- polynomial or nonlinear diameter fitting;
- Flow × Diameter regression surface;
- ANOVA/t-tests/p-values;
- mixed-effects inference;
- detection-depth estimation;
- cone-tail modelling;
- OMAG-vs-SV harmonized comparison.

OMAG remains contextual reference only during this SV-specific stage.

## 12. Planned outputs

Create a new analysis directory under the current analysis branch, for example:

`analysis/sv_diameter_stage2_scientific_v1/`

Recommended outputs:

- `common_grid_volume_metrics.csv`
- `diameter_contrasts_by_flow.csv`
- `flow_variation_by_diameter.csv`
- `source_tail_ratio_decomposition.csv`
- `depth_selected_common_grid.csv`
- `depth_profile_summary.csv`
- `slow_axis_segment_summary.csv`
- `geometry_qc_summary.csv`
- `restricted_subset_summary.csv`
- `d500_denominator_audit.csv`
- figures for the common grid, source/tail/ratio decomposition and RI(r)
- `validation.json`
- `provenance.json`
- `README.md`

## 13. Stage-2 completion criterion

Stage 2 is complete when the following are clear from the existing data:

1. whether the 128→235→285 RI decline is spatially and geometrically robust;
2. whether the large 128→235 drop and smaller 235→285 drop persist across all common flows and restricted subsets;
3. whether the D500 RI rebound is primarily explained by the vessel denominator, the tail numerator, or both;
4. how the diameter differences evolve with depth;
5. whether any apparent flow pattern is consistent across diameters at the descriptive volume level.

Only after these questions are resolved should the project decide whether any additional model fitting, new acquisition replication, or cross-algorithm comparison is scientifically justified.
