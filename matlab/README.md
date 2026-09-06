# MATLAB reconstruction and export layer

`export_sv_omag_frame.m` invokes the verified project reconstruction path and
writes one linear interim MAT file. No display transform enters the export.

The following upstream files were copied from the local OCTA OMAG
implementation so the reconstruction dependency boundary is explicit:

- `load_and_reconstruct_common_oct.m`
- `OCTA_F_SubPixReg.m`
- `SSOCT_F_PhCompV3.m`
- `OCTA_F_ED_Clutter_EigFeed.m`
- `iniset.m`, `ini2struct.m`, `config.ini`

Source directory at integration time:
`算法差异分析/OCTA/OMAG` (the common loader came from `算法差异分析/OCTA`).
The registration, phase-compensation, OMAG, and INI parser files retain their
upstream algorithmic implementation and comments; only text formatting was
normalized. The common loader is the maintained adapter that exposes one
depth-by-A-line-by-repeat group and freezes the bundled phase configuration.
The formal SV wrapper and exporter are maintained in this repository.

`load_and_reconstruct_common_oct.m` resolves the phase threshold only from the
bundled `matlab/config.ini`; it does not search the caller's current directory.
The selected threshold and configuration identifier are stored in reconstruction
metadata.

MATLAB R2023a test:

```matlab
addpath(fullfile(pwd, 'tests', 'matlab'));
test_compute_sv_maps;
test_reconstruction_dependencies;
```
