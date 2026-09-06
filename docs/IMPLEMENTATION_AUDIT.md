# Implementation audit

## 2026-09-06 upper-edge v2 candidate audit

The versioned `persistent_core_paired_edge_v2` candidate adds lateral row
support, sustained axial evidence, paired inside/outside contrast, and
physical-diameter balance to upper-edge localization.  Its explicit pilot
configuration raises the peak/noise thresholds from 0.15/3 to 0.20/4 without
changing downstream QC.  All 15 selected frames remain valid; the five
user-flagged panels move by +4, +5, +12, +12, and +1 axial pixels.

The 2500-frame audit identified a separate failure mode in flow03 near frames
425–465: sparse valid edge candidates plus deeper false candidates cause the
existing trajectory completion step to descend as much as 46 pixels relative
to v1.  The selected 0/249/499 frames are outside that event.  For this reason
v2 is retained as a review candidate, the default remains
`legacy_connected_component_v1`, and the earlier frozen result is not
overwritten.

## 2026-09-06 current-state correction

本节取代下方早期的 Reusable localization boundary 与 Available data 状态。
5 个原始 OCT 文件及对应的 5 个完整 500 帧 Flow DICOM 已定位并校验；
15 个前/中/后 sv_raw 单帧中间图也已生成。当前实现纳入师兄原始
tracking_core.py（与交付文件 SHA-256 完全一致），并仅抽取 X1/X2/X4 与
可评估性模块。师兄的拖尾 AUC 指标未接入。

当前正式定位由整卷 Viterbi、X4、X1 宽度和 z_upper 组成，不再使用固定
表面或表面距离。5 卷共 2500 帧均形成轨迹；所选 15 帧全部为 assessable，
全部通过定位、背景、窗口和最终指标 QC。完整证据见
MENTOR_TRACKING_INTEGRATION.md 和 results/pilot_mentor_tracking_15。

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
