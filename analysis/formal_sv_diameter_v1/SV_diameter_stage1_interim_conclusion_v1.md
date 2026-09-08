# SV Diameter Stage 1 — interim descriptive conclusions

Date: 2026-09-08

Base formal quantification: `codex/sv-diameter-formal-quant-v1` at `6cae8806d2070d606057873d382f37f09883ce5a`.

This file records the first scientific readout after `SV Diameter Formal Quantification v1`. It does not replace the frozen D128 conclusion and does not introduce inferential statistics. Each diameter × flow condition currently contains one independent scan volume; B-scans are slow-axis spatial samples nested within that volume.

## 1. Frozen endpoint

Primary no-background endpoint:

`RI_tail = mean(raw SV in 500 um tail ROI) / mean(raw SV in vessel ROI)`

Depth profile:

`RI(r) = raw SV lateral mean at depth r / mean(raw SV in vessel ROI)`

The source ROI is the frozen ellipse defined by X4 centre, X1 apparent lateral span, continuity-first v2.1 z geometry, and the true physical vessel diameter. The tail starts at the true physical vessel bottom, guard = 0, and uses the same X1 lateral span. No background subtraction, log transform, normalization, clipping, interpolation or zero-fill is used.

## 2. Formal data coverage

New-diameter formal quantification contains 18 independent scan volumes and 9000 nominal B-scans. Of these, 8509 frames satisfy the frozen geometry QC and enter the primary endpoint. The new-diameter set contains:

- 235 um: flows 1, 3, 5, 7, 10 mm/s
- 285 um: flows 1, 3, 5, 7, 10 mm/s
- 500 um: flows 1, 2, 3, 5, 7, 9, 10, 12 mm/s

The existing D128 formal reference contains flows 1, 3, 5, 7 and 10 mm/s.

The common comparison grid for the next stage is therefore:

- diameters: 128, 235, 285, 500 um
- common flows: 1, 3, 5, 7, 10 mm/s

## 3. Main Stage-1 finding: RI_tail decreases from 128 to 235 to 285 um

For the five common flows, the volume-level RI_tail medians are:

| Flow (mm/s) | D128 | D235 | D285 | D128→D235 change | D235→D285 change |
|---:|---:|---:|---:|---:|---:|
| 1 | 0.268943 | 0.177773 | 0.162519 | -33.90% | -8.58% |
| 3 | 0.256329 | 0.176817 | 0.162068 | -31.02% | -8.34% |
| 5 | 0.256646 | 0.173374 | 0.158516 | -32.45% | -8.57% |
| 7 | 0.323180 | 0.175291 | 0.160953 | -45.76% | -8.18% |
| 10 | 0.259467 | 0.169371 | 0.160131 | -34.72% | -5.46% |

Across the five common flows, the median of the five volume-level medians is:

- D128: 0.259467
- D235: 0.175291
- D285: 0.160953

Thus the dominant drop occurs between 128 and 235 um:

- D128→D235: -32.44%
- D235→D285: -8.18%

This pattern is present at every common flow. The 128→235 decrease is therefore substantially larger than the additional 235→285 decrease. The present data support a descriptive staged decline rather than an approximately uniform linear decline over 128–285 um.

## 4. Source–tail–ratio decomposition

To avoid interpreting a ratio as an absolute tail signal, source and tail raw SV are retained separately.

Using the median across the five common-flow volume medians:

| Diameter | vessel raw SV | 500 um tail raw SV | RI_tail |
|---:|---:|---:|---:|
| 128 um | 1.54205e8 | 3.97747e7 | 0.259467 |
| 235 um | 1.83653e8 | 3.17953e7 | 0.175291 |
| 285 um | 1.87972e8 | 3.01960e7 | 0.160953 |
| 500 um | 1.02259e8 | 2.65971e7 | 0.257642 |

The 128→235→285 RI decline is accompanied by both a rise in vessel raw SV and a fall in absolute tail raw SV.

The 500 um condition is different: absolute tail raw SV remains lower than at 235/285 um, while vessel raw SV is also markedly lower. Its RI_tail therefore returns to a level close to D128. This observation is descriptive and is not yet treated as an established diameter effect.

## 5. Preliminary depth-profile observation

The 235 and 285 um RI(r) profiles decline with depth and the 235 um profile is generally above the 285 um profile across the measured 0–500 um range. The separation is visually larger near the vessel and becomes smaller at deeper positions.

Because the primary RI(r) is not background-subtracted, non-zero deep values must not be interpreted as a tail endpoint or detection depth. Deep RI(r) can approach the raw SV background floor divided by the vessel raw SV.

The 500 um profiles show a higher deep relative floor than 235/285 um in several flows. Given the markedly lower D500 vessel denominator, the next stage must explicitly separate near-vessel tail signal, deep raw-SV floor and vessel denominator before mechanistic interpretation.

## 6. Current interpretation boundary

Current evidence supports the following descriptive statement:

> Across the presently available independent scan volumes, SV relative tail intensity decreases strongly from 128 to 235 um and decreases further, more modestly, from 235 to 285 um. The 500 um volumes show a higher RI_tail than 235/285 um despite lower absolute tail raw SV, indicating that source-signal changes are an important component of the ratio.

Current evidence does not support:

- treating B-scans as independent experimental replicates;
- claiming that diameter causally produces these differences;
- fitting a monotonic diameter-response relation before the D500 pattern and geometry/QC composition are audited;
- defining a tail endpoint from RI(r).

## 7. Authoritative source files

New-diameter formal results:

- `analysis/formal_sv_diameter_v1/framewise_primary.csv`
- `analysis/formal_sv_diameter_v1/volume_summary.csv`
- `analysis/formal_sv_diameter_v1/depth_summary_10um.csv`
- `analysis/formal_sv_diameter_v1/depth_selected.csv`
- `analysis/formal_sv_diameter_v1/validation.json`
- `analysis/formal_sv_diameter_v1/provenance.json`

D128 reference:

- `analysis/formal_sv_d128_v21_run001/no_background_relative_tail_v1_full2422/scan_flow_summary.csv`
- `analysis/formal_sv_d128_v21_run001/no_background_relative_tail_v1_full2422/depth_selected.csv`
- `analysis/formal_sv_d128_v21_run001/SV_d128_five_flow_final_conclusion_v1.md`
