## SCutils v0.7.1-beta

This prerelease fixes the all the layer-related options so now `layer`, `layer1`, `layer2`, and `layer.gradient` accept `NULL`, `NA`, or a single non-NA character value. Empty strings and `"null"` (case-insensitive) are also treated as `NULL`. This allows for more flexible per-feature layer selection in `FeatureScatterGradient()`, `FeatureDensityPlot()`, and `VlnPlotGradient()`.

### 🏷️ Release tag

- GitHub prerelease tag/title: `v0.7.1-beta`.
