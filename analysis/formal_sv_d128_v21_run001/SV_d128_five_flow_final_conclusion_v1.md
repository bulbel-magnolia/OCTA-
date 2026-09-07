# SV 128 µm five-flow analysis — consolidated conclusion v1

## Scope

This document consolidates the current conclusions from the 128 µm artificial-vessel SV analysis before moving to vessel-diameter analysis.

Current primary endpoint:

`RI_tail = mean(raw SV in 500 µm tail ROI) / mean(raw SV in vessel ROI)`

Depth profile:

`RI(r) = mean(raw SV across the vessel lateral span at depth r) / mean(raw SV in vessel ROI)`

No background subtraction is used in the current primary definition. Historical background-corrected outputs remain archived in their original directories and are not used for the current primary interpretation.

## Dataset and statistical scope

- Vessel inner diameter: 128 µm.
- Flow volumes: 1, 3, 5, 7, and 10 mm/s.
- One independently acquired scan volume per flow condition.
- 2500 B-scans planned; 2422 passed the frozen formal validity criteria.
- Valid B-scans by flow: 486, 468, 491, 493, and 484.
- B-scans are slow-axis spatial samples within one volume, not independent experimental replicates.
- Therefore all flow comparisons in this dataset remain descriptive at the scan-volume level; no causal flow-speed effect is claimed and no flow-level p-values are reported.

## Primary five-flow result

Volume-level median `RI_tail` values:

| flow (mm/s) | RI_tail median |
|---:|---:|
| 1 | 0.268943 |
| 3 | 0.256329 |
| 5 | 0.256646 |
| 7 | 0.323180 |
| 10 | 0.259467 |

Relative to 1 mm/s, the 7 mm/s volume is approximately +20.17% in the primary `RI_tail` definition. Relative to the median of the other four scan-volume medians, it is approximately +25.24%.

## Source–tail decomposition

The elevated `RI_tail` in the current 7 mm/s volume is driven mainly by the vessel denominator rather than by a large increase in raw tail signal.

Relative to the median of the other four scan-volume medians:

- vessel raw mean in flow07: approximately -20.09%
- tail raw mean in flow07: approximately +1.75%
- `RI_tail` in flow07: approximately +25.24%

Thus the current descriptive pattern is:

`vessel raw SV lower + tail raw SV broadly similar -> relative tail intensity higher`.

The correct interpretation is therefore not that the 7 mm/s tail signal is absolutely much stronger, but that the tail is stronger relative to the vessel signal in this particular scan volume.

## Depth-window robustness

Exact native-pixel integration was recomputed for 100, 200, 300, and 500 µm tail windows.

| tail window | flow01 | flow03 | flow05 | flow07 | flow10 | flow07 vs flow01 |
|---:|---:|---:|---:|---:|---:|---:|
| 100 µm | 0.4821 | 0.4707 | 0.4648 | 0.5733 | 0.4659 | +18.92% |
| 200 µm | 0.3884 | 0.3734 | 0.3703 | 0.4650 | 0.3740 | +19.72% |
| 300 µm | 0.3314 | 0.3171 | 0.3147 | 0.3980 | 0.3193 | +20.09% |
| 500 µm | 0.2689 | 0.2563 | 0.2566 | 0.3232 | 0.2595 | +20.17% |

The 7 mm/s volume ranks first at all tested tail-window lengths. Therefore the observed ranking is not created by the choice of a 500 µm tail window.

## Within-volume spatial robustness

Each 500-frame volume was divided into five contiguous 100-frame slow-axis segments.

For flow07, `RI_tail` remained rank 1/5 in every spatial segment. Relative to the median of the other four scans within the same segment, flow07 was higher by approximately:

- frames 0–99: +15.98%
- frames 100–199: +41.09%
- frames 200–299: +20.58%
- frames 300–399: +22.05%
- frames 400–499: +24.40%

In the same five segments, the flow07 vessel raw mean was the lowest among the five volumes in every segment, whereas tail raw mean differed from the other scans by only a few percent.

Therefore the flow07 pattern is distributed across the volume rather than being driven by one isolated spatial region.

## Aggregation robustness

Six pre-specified volume-level aggregation methods were audited:

1. primary median of framewise RI
2. mean of framewise RI
3. ratio of median raw tail mean to median raw vessel mean
4. ratio of mean raw tail mean to mean raw vessel mean
5. pooled area-weighted raw ratio
6. median of five spatial-segment RI medians

Flow07 ranked first under all 6/6 methods. Its contrast versus flow01 remained approximately +19% to +22% across these aggregation choices.

This indicates that the flow07 ranking is not an artifact of the primary `median(framewise RI)` aggregation convention.

## Geometry and localization robustness

Residual QC evaluated X1 apparent width, source ROI area, z-top position, X4 position, and localization-source composition.

Flow07 has a somewhat smaller apparent width and source area than the other scans, but the high-RI pattern persists after restricted-subset analyses.

Flow07 remains rank 1/5 under:

- all valid frames
- central 80% of apparent width
- central 80% of z position
- central 80% of both z and width
- direct-localization frames only
- direct-localization plus central-80% z and width

Under direct-only frames, flow07 remains approximately +23.59% above the median of the other four scans, while its vessel raw mean remains approximately -18.42% and its tail raw mean approximately +1.67% relative to the other four.

Localization-source composition is similar across the five volumes and does not explain the flow07 result.

## Additional robustness results retained in the analysis tree

The current primary endpoint has also been audited for:

- fixed 128 µm lateral-width geometry stress
- spatial autocorrelation along slow-axis B-scan index
- systematic subsampling stability
- source–tail–ratio coupling
- exact depth profiles `RI(r)`

The fixed-128 geometry variant preserves the same qualitative flow07 ranking. Systematic sampling of approximately 100 uniformly distributed B-scan positions per volume reproduces the full-volume `RI_tail` median within 5% for all tested scan/phase combinations in the current dataset, although full acquisition remains the reference dataset.

## Consolidated scientific conclusion

In the present 128 µm artificial-vessel dataset, the 1, 3, 5, and 10 mm/s scan volumes show broadly similar no-background raw-SV relative tail intensity, whereas the current 7 mm/s scan volume consistently shows higher `RI_tail`.

Across source–tail decomposition, exact 100–500 µm tail windows, five spatial segments, six volume-level aggregation methods, fixed-width stress testing, restricted geometry subsets, and localization-source restrictions, the same qualitative pattern is retained.

The most consistent internal explanation is that the 7 mm/s volume has lower vessel raw SV signal while its tail raw SV signal remains broadly comparable with the other volumes, thereby increasing the tail-to-vessel ratio.

This establishes a robust characteristic of the current flow07 scan volume, but it does **not** establish that 7 mm/s flow causally produces this characteristic, because each flow condition currently contains only one independent scan volume.

## Transition to diameter analysis

The 128 µm five-flow analysis is considered sufficiently complete for the current dataset. Further work should now proceed to vessel-diameter analysis using the same primary no-background `RI_tail` definition and the same interpretation discipline. Diameter analysis should preserve condition-level scan volumes as the experimental units and should avoid treating individual B-scans as independent replicates.

## Related result directories

- `analysis/formal_sv_d128_v21_run001/no_background_relative_tail_v1_full2422/`
- `analysis/formal_sv_d128_v21_run001/no_background_followup_three_audits_v1_full2422/`
- `analysis/formal_sv_d128_v21_run001/no_background_final_two_audits_v1_full2422/`

