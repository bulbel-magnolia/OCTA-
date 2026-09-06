# Synthetic demonstration

This directory contains deterministic, fully synthetic output for interface and
QC inspection. It is not experimental evidence and must not be included in a
flow-speed analysis.

The example has a 128 um source, 12.7 um lateral pitch, 6.7 um axial pitch,
and a constant signed tail excess in the 500 um formal window. It demonstrates:

- scalar source/tail integration and the ratio in `frame_results.csv`;
- the complete signed `V`, `B_left`, `B_right`, `B`, `T`, and `P` profiles;
- 5-row detection bins, with the final complete bin at 469.0 um marked as
  right-censored;
- representative QC01-QC03 and detection-threshold figures.

Regenerate into a new empty directory from the repository root:

```powershell
$env:PYTHONPATH = (Resolve-Path 'src').Path
python scripts/generate_synthetic_demo.py --output <empty-output-directory>
```
