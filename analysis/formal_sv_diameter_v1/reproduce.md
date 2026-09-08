# Reproduction and validation

Use the repository at the recorded code commit with the package versions in provenance.json. The intake worktree supplies local file locators and frozen tracking tables; committed volume_manifest.csv supplies conditions, calibration and input identities.

1. Obtain the 25 official assets listed in results/formal_sv_d128_v21_full2500_run001/download_packages.csv from the formal-sv-d128-v21-run001 GitHub Release. Keep these local under outputs/formal_quant_v1/d128_packages. The regression phase checks package size/SHA and each used NPZ SHA before calculating metrics.
2. Run `python scripts/quantify_sv_diameter_formal_v1.py quantify --intake-local ../sv-diameter-intake-v1`. Adjust the local locator argument for the host. Existing MAT files are read, hashed and identity-checked; the script does not overwrite or reconstruct them.
3. Run `python scripts/quantify_sv_diameter_formal_v1.py regression` for all 2422 frozen-valid D128 frames. This uses the new kernel and both frozen no-background reference tables, without running historical scientific-analysis scripts.
4. Run `python scripts/quantify_sv_diameter_formal_v1.py finalize` after both phases pass. This checks protection hashes and writes the lightweight evidence records.
5. Run `python -m pytest -q`. Tests include the analytical kernel, MAT identity rejection, all-frame coverage, native-depth coordinates, descriptive aggregation and D128 replay evidence. A host-specific temporary directory may be supplied with --basetemp.

The source and result SHA records identify the executed definitions. A subsequent provenance-only commit records the metrics commit SHA; a commit cannot contain its own hash. Tests and publication receipts are recorded after calculation. No raw OCT/DICOM, MAT, NPZ or release package belongs in this output directory or Git commit.
