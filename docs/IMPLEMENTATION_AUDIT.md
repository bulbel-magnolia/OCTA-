# Implementation audit

Date: 2026-09-06

## Confirmed local components

- MATLAB R2023a (`9.14.0.2206163`) runs in batch mode.
- `load_and_reconstruct_common_oct.m` exists and returns a complex
  depth-by-A-line-by-repeat array after FFT, crop `50:400`, subpixel
  registration, and phase compensation.
- The formal SV implementation uses `var(abs(IMG), 1, 3)`.
- The mentor delivery package passes its Python structural self-check.
- The existing `SV_Tail_Precise_SelfTest` passes in MATLAB R2023a.

## Reusable localization boundary

The first implementation reuses the mentor method's local-body X1 rule on a
co-registered OMAG frame. The selected continuous component is exposed as left
and right pixel-centre indices, then converted to physical pixel edges before
area weighting. The X4 robust centroid is retained only as a possible tracking
anchor and is not used as the final geometric centre.

## Available data

No original `.oct` file is currently present on the mounted C, D, or F data
locations searched. A historical archive contains one reconstructed frame for
each speed (1, 3, 5, 7, and 10 mm/s), plus two manually selected blank regions.
Those files support a non-formal smoke test only. They cannot supply the planned
front/middle/rear 15-frame pilot or a condition-level result table.

Unknown acquisition identities remain missing in the manifest. Filenames are
used only as flow labels where that mapping was explicitly established; they
are never interpreted as vessel diameter.

## Implemented verification

- The 29-test Python synthetic suite covers the formal N-denominator signal, fractional
  geometry, fixed bilateral background, retained negative residuals, source
  and tail integrals, incomplete-window NA behavior, X1 localization,
  top-edge re-establishment, detection persistence, and right censoring.
- MATLAB R2023a independently verifies `var(abs(E), 1, 3)` on a known array;
  a second test executes the copied subpixel registration, phase compensation,
  and OMAG eigenspace-filter dependency chain on deterministic complex data.
- Five historical vessel frames complete the full 500 um batch path; two
  historical blank selections fail source QC and do not produce a ratio.
- A run freezes configuration and manifest hashes, software versions, scalar
  and scan-level tables, per-frame full profiles, complete MAT arrays, background
  pixels, area weights, sensitivity results, six-panel QC figures, and logs.
- Unknown identity or acquisition metadata stays missing and is counted rather
  than inferred. A declared one-sided background requires an excluded side plus
  a factual reason; undeclared one-sided input fails background QC.
- Recoverable frame-level input/localization errors produce NA records and do
  not prevent the remaining manifest rows from completing.
