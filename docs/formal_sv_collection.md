# Formal SV collection runner

The approved collection uses baseline `ee3012c04be0d02289520a34409d0855a3e5ca58` and the existing continuity-first v2.1 primary tables. `collection.py` verifies the frozen algorithms/configs against that Git revision. It does not run the tracker. MATLAB's new batch wrapper calls the existing frame exporter without changing reconstruction.

Commands (replace the local roots; never publish local job files):

```powershell
python scripts/collect_formal_sv.py prepare --stage bridge --raw-root <raw-root> --output results/formal_sv_d128_v21_bridge15_run001 --archive <archive-parent>/bridge15
# MATLAB: addpath('matlab'); export_sv_collection('<archive-parent>/bridge15/local_jobs.json');
python scripts/collect_formal_sv.py process --scan flow01 --output results/formal_sv_d128_v21_bridge15_run001 --archive <archive-parent>/bridge15
```

Process all five scan IDs. Review all 15 mapping images and numerical checks, then record `bridge_validation.json` with `stage_b_authorized_by_checks=true` only when the approved Stage A criteria pass. Full preparation requires this record; use `--stage full`, the full2500 result directory and `<archive-parent>/full2500`. Interim files are shared under `<archive-parent>/interim`. Each frame's archive is compressed losslessly and verified by a round trip before checkpoint completion. Resumption verifies array/config/tracking hashes. Interrupted MATLAB exports have temporary names and cannot masquerade as completed frame exports. `--remove-interim` optionally removes only successfully archived temporary MAT files within that archive parent's interim directory.

Formal SV, geometry, background and sensitivity definitions come from the frozen code and run configuration. Missing geometry retains its own full input arrays and all-NA quantitative profiles. Invalid quantitative values are marked diagnostic and excluded from scan summaries. No blank acquisition exists, so detection remains not evaluated. The auxiliary DICOM array is for localization checks; quantitative SV always comes from the OCT reconstruction.

`package` produces a self-contained bridge ZIP and full-run ZIPs of 100 frames, including complete B-scans, fractional source/tail weights, background column indices, portable reconstruction metadata, tables/configs and matching QC/profile subsets. All ZIP and NPZ checksums accompany the Release downloads. Original OCT/DICOM/MAT files are not committed.

For the approved full run, after all workers finish, execute `summarize`, then `python scripts/verify_formal_sv.py --output <full-result-directory> --archive <archive-parent>/full2500`. `python scripts/finalize_formal_sv_run.py` checks frozen coordinates and bridge/full agreement and assembles the fixed run001 README/QC gallery and output hashes. Finally run `package`. Verification and finalization fail on inconsistent inputs rather than changing scientific parameters. Text checksums refer to the exact bytes in ZIP files; Git may normalize text line endings.
