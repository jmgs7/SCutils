## SCutils v0.6.4-beta

This prerelease refactors `FeatureDensityPlot()` to use standard sequential `lapply()` only and removes the unused `mc.cores` option.

### ✨ FeatureDensityPlot updates

- Replaced multicore/future execution with base `lapply()`.
- Chose a conditional nested `lapply()` design:
  - outer loop over features,
  - inner loop over groups only when `split.plot = TRUE`.
- Removed the `mc.cores` argument from the function interface.
- Kept line styling consistent:
  - `vline` remains dashed red,
  - median overlays remain dashed black.

### 🧹 Code clarity improvements

- Refactored internal variable and helper names to use dot-separated style.
- Expanded inline comments to explain validation, grouping, line drawing, and return-shape logic.
- Added clearer sectioning in the function body for easier maintenance.

### 🧾 Documentation and package updates

- Updated roxygen docs and regenerated package help files.
- Updated README to reflect the lapply-only implementation.
- Bumped package version to `0.6.4`.

### 🗑️ Cleanup

- Removed older release note files from the repository.

### 🏷️ Release tag

- GitHub prerelease tag/title: `v0.6.4-beta`.
