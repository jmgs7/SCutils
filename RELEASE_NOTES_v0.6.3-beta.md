## SCutils v0.6.3-beta

This prerelease migrates `FeatureDensityPlot()` to the `future` framework and improves plotting consistency and code style.

### ✨ FeatureDensityPlot parallelization migration

- Replaced nested `parallel::mclapply()` with `future.apply::future_lapply()`.
- Implemented single-level task parallelization across feature/group tasks to avoid nested parallel apply behavior.
- Added support for user-controlled parallel backends via `future::plan(...)`.
- Kept optional per-call worker override through `mc.cores` using a temporary local plan.

### 🎨 Plot behavior updates

- Updated `vline` rendering to be consistently **dashed red** in all branches.
- Preserved independent median overlays in **black** when `plot.median = TRUE`.

### 🧹 Code quality and style

- Refactored internals to be more idiomatic and readable (Google-style oriented naming/structure for internal objects and helpers).
- Expanded sectioned comments to clarify parallelization design, line-layering semantics, and assembly logic.

### 🧾 Documentation and package updates

- Updated roxygen documentation for `FeatureDensityPlot()` to reflect:
  - future-based execution,
  - dashed `vline` behavior,
  - current return object semantics.
- Updated README function section and examples for `future::plan(...)` usage.
- Bumped package version to `0.6.3`.
- Regenerated `NAMESPACE` and Rd documentation with roxygen2.

### 🏷️ Release tag

- GitHub prerelease tag/title: `v0.6.3-beta`.
