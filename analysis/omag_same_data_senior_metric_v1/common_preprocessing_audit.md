# Common preprocessing identity audit

Gate: PASS. Fixed raw frames: 45. OMAG independent/reuse checks: PASS.

The raw reader follows the frozen senior exporter: three repeats, double FFT, crop 50:400, Nsub=10, Colshift=false, phase arg 1, configuration Thr=85 dB, median shift true. An independent call to the retained-SV common reader gives the same complex IMG. On that IMG, var(abs(IMG),1,3) is compared pixelwise with retained SV. Tolerances were frozen in analysis_plan.json before reconstruction (rtol 1e-10, atol 1e-6). No parameter search occurred. See validation CSVs for full errors, dtypes and array hashes.

The retained SV exporter already computes second-output OMAG NumEV=2 on exactly this IMG, without filtering. All retained files must additionally pass their frozen SHA and metadata checks before reuse. Raw OMAG remains float64, preserving the retained signal with no storage quantization; full-pipeline localization casts to float32 exactly as the senior loader and previous SV-B adapter. Metrics use original raw precision, matching SV-A/SV-B comparison precision.

config.ini contains an unrelated display crop setting; both audited computational readers explicitly use 50:400. Actual dependency paths and hashes, including iniset and config, are recorded in reconstruction_contract.json. The independent reader shares the audited mathematical dependencies and differs in raw-file reader code; it is a computational reproducibility check, not an independent acquisition.

All six senior/retained dependency files have exactly matching nonempty lines after trimming surrounding whitespace; actual file SHA differences are newline/blank-line/trailing-space differences. See validation/dependency_equivalence.csv. Fixed SV and OMAG results are bitwise equal in all 45 frames, not merely within tolerance.
