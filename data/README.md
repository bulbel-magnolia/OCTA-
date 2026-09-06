# Data layout

This repository is configured to store OCT/MAT/HDF5/NumPy binary data through
Git LFS. Do not bypass the LFS rules for raw scans.

- `raw/`: immutable source scans and acquisition sidecars.
- `interim/`: reconstructed linear maps used as pipeline inputs.
- `processed/<run_id>/`: frozen run configuration, manifest, tables, profiles,
  arrays, QC figures, and logs.

Raw `.oct` files in the current acquisition format are approximately 1.5 GB
per scan. Confirm the repository's available Git LFS storage before adding a
batch. Each committed run must include the exact manifest and run configuration
used to produce it.
