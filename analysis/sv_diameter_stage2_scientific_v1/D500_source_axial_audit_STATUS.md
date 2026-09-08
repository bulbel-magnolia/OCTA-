# D500 frozen-source axial local-array audit — status

Date: 2026-09-08

Status: **IMPLEMENTED / SCIENTIFIC EXECUTION BLOCKED BY NON-PERSISTED LOCAL ARRAYS**

## Purpose

Secondary diagnostic for the Stage-2 finding that D500 has a markedly lower frozen vessel/source raw-SV denominator than D285 while its absolute 500-um tail raw SV is also lower. This audit does not alter the frozen primary source definition.

## Frozen analysis retained

- Signal: `sv_raw = var(abs(E),1,3)`, denominator N.
- No log, normalization, gain, clipping, positive truncation or background subtraction.
- Geometry: frozen X4 centre, X1 apparent lateral span, continuity-first v2.1 z-top, physical vessel diameter.
- Source ROI: unchanged fractional ellipse, supersample 16.
- Experimental unit: scan volume; B-scans remain slow-axis spatial positions.

## Secondary axial subdivision

The unchanged frozen source ellipse is partitioned by relative physical axial coordinate into:

- upper third: 0 to 1/3 of physical diameter from z_top;
- middle third: 1/3 to 2/3;
- lower third: 2/3 to 1.

Subpixel samples are partitioned inside the same 16x16 ellipse supersampling grid. The three band weight arrays must reconstruct the full frozen ellipse. The audit reports band raw-SV mean, band Q contribution, lower/upper ratio, flow summaries, slow-axis fifths and predefined central/direct restrictions. D285 same-flow volumes are the primary comparator for D500.

## Mandatory replay gates

For every analyzed MAT:

1. MAT SHA-256 must equal `analysis/formal_sv_diameter_v1/framewise_primary.csv::input_mat_sha256`.
2. Recomputed full-ellipse source area, source Q and source mean must reproduce the frozen framewise values within relative tolerance `1e-10`.
3. Upper + middle + lower band Q must reconstruct the full source Q.

No scientific output is valid unless all gates pass.

## Input availability audit

The intake record explicitly states that local per-frame MAT files are never committed and that raw OCT/DICOM data are not copied into the repository. The committed repository contains per-frame MAT names, sizes and SHA-256 values, but not the MAT bytes. The only public raw-array release currently found is the frozen D128 release; it does not contain D285/D500 arrays.

The currently connected File Library and Google Drive searches did not locate the D500 raw OCT files or new-diameter MAT export. The conversation-local `/mnt/data` contains the older selected B-scan RAR packages and code packages, not the D500 full-volume MAT export.

Therefore the exact upper/middle/lower scientific numbers cannot be reconstructed from committed summary tables alone. `source_mean_raw` is a weighted mean over the whole ellipse; the internal axial distribution is not identifiable without `sv_raw` pixels.

## Current Stage-2 evidence before local-array execution

The existing frozen/Stage-2 results already establish that:

- D500 source raw SV is lower than D285 across all 25 common-flow × slow-axis-fifth comparisons;
- D500 RI rebound persists in all predefined restricted geometry/QC subsets;
- within D500, source raw SV has strong spatial association with z_top in several volumes, but central-z restrictions do not remove the cross-diameter low denominator;
- deep D500 and D285 absolute raw-SV values can converge while D500 RI remains higher, consistent with a low-denominator normalized floor.

The local-array audit is therefore targeted at locating the low source signal within the frozen vessel ellipse, not at deciding whether the low D500 denominator exists.

## Executable implementation

`analysis/sv_diameter_stage2_scientific_v1/d500_source_axial_local_array_audit.py`

Example once the retained local MAT root is mounted:

```bash
python analysis/sv_diameter_stage2_scientific_v1/d500_source_axial_local_array_audit.py \
  --mat-root <retained-formal-MAT-root>
```

Planned outputs are written under `analysis/sv_diameter_stage2_scientific_v1/d500_source_axial_audit/` and include framewise bands, volume summaries, D500-vs-D285 band contrasts, slow-axis robustness, `validation.json` and `provenance.json`.
