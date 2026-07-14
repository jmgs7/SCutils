## SCutils v0.8.2-beta

### 🚀 New Features

- The `CalculateQC()` function now admits additional arguments to pass to the internal computation of feature threshold by the ComputeFeatureThreshold module (internal, not user-facing). This allows for more flexibility in adjusting the parameters for the computation of the MALAT1 threshold test, enabling users to customize the analysis according to their specific needs.

- We have included a Utils module where we will upload utility functions. The first function of this module, `ExtractFeatureTestResults()`, retrieves and summarizes the per-cell results from the MALAT1 test of `CalculateQC()`. This function is designed to facilitate the extraction of feature threshold test results from a Seurat object, providing users with a convenient way to access and analyze these results.

### 🛠️ Bug Fixes

- Improvement of code readability and maintainability by refactoring the `FeatureDensityPlot()` function. This includes better organization of the code, clearer variable names, and enhanced documentation to make it easier for users to understand and utilize the function effectively.
- Enhanced documentation of the code in several modules.

### 🏷️ Release tag

- GitHub prerelease tag/title: `v0.8.2-beta`.

**NOTE**: This is a beta release, so some features may still be under development or need further testing. We encourage users to test the new features and provide feedback. Please report any issues or bugs encountered during testing to help us improve the package.

---

### Previous Release Notes (v0.8.0 series)

### v0.8.1

### 🛠️ Bug Fixes

- Unclutters namespace between dplyr and stats.

---

#### v0.8.0

##### 🚀 New Features

- **Improved QC Metrics**: The `CalculateQC()` function now computes log-normalized counts per cell and detected features per cell for the normalized layer, adding them as `nCount_logRNA` and `nFeature_logRNA`.
- **Cell Cycle Scoring**: The function now includes cell cycle scoring capabilities.
- **MALAT1-based QC Thresholding**: The function now supports MALAT1-based QC thresholding, adding `MALAT1.threshold` (numeric) and `MALAT1.pass` (boolean) columns to `meta.data`.

##### 🛠️ Bug Fixes

- Fixed issues with layer-related options in various plotting functions.
- Unclutters namespace between dplyr and stats.

