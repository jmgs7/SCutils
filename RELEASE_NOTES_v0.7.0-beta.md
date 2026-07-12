## SCutils v0.7.0-beta

This prerelease significantly extends `FeatureScatterGradient()` with grouping, multi-layer, and configurable correlation support; migrates `FeatureDensityPlot()` to `Seurat::FetchData()` with per-feature layer support; and adds per-feature custom title support to `VlnPlotGradient()`.

### ✨ FeatureScatterGradient major update

- **Grouping support**: added `group.by` parameter (default `"ident"`). When non-`NULL`, the function returns a combined patchwork of per-group scatter panels instead of a single plot.
- **Per-group correlations**: in grouped mode, correlation is recomputed independently for each group level and displayed as the subtitle of each panel.
- **Configurable correlation method**: new `corr.method` parameter accepts `"pearson"` (default), `"spearman"`, or `"kendall"`.
- **Per-feature layer selection**: new `layer1`, `layer2`, and `layer.gradient` parameters allow fetching each of the three axes from different Seurat v5 assay layers.
- **Custom grouped title**: new `plot.title` parameter sets a custom main title for grouped patchwork output; ignored in ungrouped mode.
- **Global gradient scale in grouped mode**: gradient limits are computed globally from all data and applied consistently via a single shared patchwork legend.
- **Asymmetric limit clamping**: `lower.limit` and `upper.limit` are now independently optional (`NULL`); when only one is set, the missing bound is inferred numerically to prevent duplicate patchwork guides.
- **FetchData-based data access**: all feature and grouping data are resolved via `Seurat::FetchData()`, consistent with `VlnPlotGradient()` and `FeatureDensityPlot()`.
- **Axis label conventions**: axis labels follow `feature_layer` for assay-backed features when a layer is specified, and bare `feature` for metadata-backed features.
- **Duplicate legend fix**: resolved a bug where multiple redundant gradient legends appeared in grouped mode due to NA limit handling in patchwork.
- **Input hardening**: comprehensive validation added for all new and existing parameters.

### ✨ FeatureDensityPlot major update

- **Migrated to `Seurat::FetchData()`**: features are now resolved via `FetchData()`, enabling support for metadata columns, assay features from any layer, and dimensional reduction variables (e.g. `"PC_1"`).
- **Per-feature layer selection**: new `layer` parameter accepts a scalar (applied to all features) or a per-feature character vector, mirroring `VlnPlotGradient()` semantics.
- **Stable internal key design**: duplicate feature names with different layers are handled via position-based internal keys, preserving correct per-occurrence layer assignment.
- **Default grouping variable changed**: `group.by` now defaults to `"ident"` (was `"active.ident"`) for direct `FetchData()` compatibility.
- **Output name convention**: returned list names use `feature_layer` for assay-backed features and bare `feature` for metadata-backed features; layer suffixes are cleanly stripped for metadata columns.
- **Bug fix**: metadata-backed features no longer receive spurious layer suffixes in returned list names.

### ✨ VlnPlotGradient: per-feature custom titles

- **New `plot.title` parameter**: accepts `NULL` (default, preserves historical naming), a single string (recycled to all panels), or a character vector of length `length(features)` (one title per panel).
- Follows the same title-resolution contract as `FeatureDensityPlot()`.

### 🧹 Code clarity and hardening

- Consistent use of `Seurat::FetchData()` across `FeatureDensityPlot()`, `VlnPlotGradient()`, and `FeatureScatterGradient()` for all data access.
- Expanded input validation with explicit, informative error messages across all updated functions.
- Formatting and comment improvements throughout the codebase.

### 🧾 Documentation and package updates

- Updated roxygen docs and regenerated package help files for all modified functions.
- Updated README to reflect new parameters and function signatures.
- Bumped package version to `0.7.0`.

### 🏷️ Release tag

- GitHub prerelease tag/title: `v0.7.0-beta`.
