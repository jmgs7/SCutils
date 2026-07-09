## SCutils v0.6.2-beta

This prerelease improves `FeatureDensityPlot()` with new line controls, title customization, and clearer output behavior for multi-feature calls.

### ✨ FeatureDensityPlot updates

- Added `plot.median` argument (default: `TRUE`) to draw median line(s) independently from `vline`.
- Enforced line color separation for clarity:
  - `vline` reference line(s): **red**
  - `plot.median` median line(s): **black**
- Changed `split.plot` default to `TRUE`.
- Added `plot.title` argument to support custom titles (single title or one per feature).
- Centered global feature titles in split-plot outputs.
- Updated multi-feature return behavior:
  - single feature -> single plot object
  - multiple features -> named list of plots (names match `features`)

### 🧾 Documentation and maintainability

- Expanded inline code comments in `FeatureDensityPlot.R` to improve readability and maintenance.
- Updated roxygen documentation and regenerated `man/FeatureDensityPlot.Rd`.
- Updated README function section to reflect new arguments, defaults, line-color behavior, and return shape.

### 📦 Versioning note

- R package version updated to `0.6.2` (required valid DESCRIPTION format).
- GitHub prerelease tag/title remains `v0.6.2-beta`.

### 🔎 Compatibility

- No breaking changes to other exported functions.
