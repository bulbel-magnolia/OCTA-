# Background +2 A-line full-volume sensitivity audit

## What was changed

This is a **QC sensitivity audit**, not a replacement of frozen run001. Baseline background skips 3 A-lines from each vessel edge and then takes 5 A-lines per side. The variant skips 5 A-lines and still takes 5 per side, so both background strips move **2 A-lines = 25.4 μm outward**. Source/tail geometry, SV definition, 500 μm tail window, signed residual handling and inclusion rules are unchanged.

Spatial B-scans in one volume are not independent experiments; this audit is descriptive and reports no p-values.

## Headline results

- Baseline-valid frames audited: **2422**.
- +2 variant valid: **2422/2422 (100.00%)**.
- |ΔQ_tail|: median **5.352%**, P95 **16.544%**, max **71.943%**.
- |ΔQ_vessel|: median **1.054%**, P95 **3.593%**.
- |ΔR|: median **0.0299**, P95 **0.0978**.
- Largest per-flow median |ΔQ_tail|: **flow01 = 5.618%**.
- Frames with |ΔQ_tail| >5% / >10% / >20%: **1291 / 524 / 60**.

Large relative ΔQ_tail can occur when signed baseline Q_tail is close to zero. Therefore the framewise file retains absolute values/deltas, and the extremes file is diagnostic rather than an exclusion list.

## Files

- `background_plus2_framewise.csv`: all baseline-valid frames.
- `background_plus2_summary.csv`: overall, per-flow and localization-source robust summaries.
- `background_plus2_extremes.csv`: top 50 by |relative ΔQ_tail|.
- `audit_metadata.csv`: perturbation definition and row counts.
- `audit_background_plus2.py`: reproducible script.

## Interpretation rule

This audit tests whether a reasonable outward shift of local background materially changes Q_vessel, Q_tail or R. It must not be used to select whichever background gives a more monotonic or statistically favorable flow trend. The frozen directory `results/formal_sv_d128_v21_full2500_run001/` is untouched.
