## SCutils v0.8.1

### 🛠️ Bug Fixes

- Unclutters namespace between dplyr and stats.

### 🏷️ Release tag

- GitHub prerelease tag/title: `v0.8.1`.

#### Previous Release Notes (v0.8.0)

##### 🚀 New Features

- **Improved QC Metrics**: The `CalculateQC()` function now computes log-normalized counts per cell and detected features per cell for the normalized layer, adding them as `nCount_logRNA` and `nFeature_logRNA`.
- **Cell Cycle Scoring**: The function now includes cell cycle scoring capabilities.
- **MALAT1-based QC Thresholding**: The function now supports MALAT1-based QC thresholding, adding `MALAT1.threshold` (numeric) and `MALAT1.pass` (boolean) columns to `meta.data`.

##### 🛠️ Bug Fixes

- Fixed issues with layer-related options in various plotting functions.
- Unclutters namespace between dplyr and stats.