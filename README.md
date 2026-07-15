# SCutils

An R package providing single-cell RNA-seq utility functions to complement Seurat-based workflows.

## Status and version

- Current package version: **0.8.5-beta**
- Status: **beta**; some features are still under development and may change in future releases. Please report any issues on GitHub.

## Installation

```r
# install.packages("remotes")
remotes::install_github("jmgs7/SCutils")
```

## Functions

### `BatchOpenH5()`

Batch-opens `.h5` single-cell matrices, writes/loads BPCells matrix directories (`*_BP`), supports both 10X and AnnData HDF5 inputs, and can optionally convert ENSEMBL IDs to gene symbols; it can also generate per-cell sample provenance metadata. Processing uses `parallel::mclapply()` (on Windows, cores are forced to 1).

```r
BatchOpenH5(
  files,
  relative = TRUE,
  BP.data.dir = NULL,
  platform = "10X",
  use.names = TRUE,
  ensembl.to.symbol = FALSE,
  species = "human",
  generate.metadata = FALSE,
  mc.cores = length(files)
)
```

| Parameter | Default | Description |
| :-- | :-- | :-- |
| `files` | — | Character vector of input `.h5` file paths |
| `relative` | `TRUE` | Store matrix directory paths as relative paths inside each matrix object |
| `BP.data.dir` | `NULL` | Output directory for `*_BP` folders; defaults to `dirname(files[[1]])` |
| `platform` | `"10X"` | Input format: `"10X"` or `"anndata"` |
| `use.names` | `TRUE` | For 10X only: replace feature IDs with `/matrix/features/name` |
| `ensembl.to.symbol` | `FALSE` | Convert ENSEMBL IDs to symbols via `ConvertEnsembleToSymbol2()` (applied only when `use.names = FALSE`) |
| `species` | `"human"` | Species passed to symbol conversion: `"human"` or `"mouse"` |
| `generate.metadata` | `FALSE` | If `TRUE`, returns a named list with `data.list` plus per-cell metadata (`cell.tag`, `sample.procedence`) |
| `mc.cores` | `length(files)` | Cores used by `mclapply()`; forced to `1` on Windows |

***

### `CalculateCDR()`

Calculates the **Cellular Detection Rate** (fraction of detected features per cell), z-scales it, and adds a `CDR` column to `SeuratObject@meta.data`. Handles both Seurat v3 (`@assays$RNA@counts`) and v5/split-layer (`@assays$RNA@layers$counts`) object formats.

```r
SeuratObject <- CalculateCDR(SeuratObject)
```

| Parameter | Default | Description |
| :-- | :-- | :-- |
| `SeuratObject` | — | A Seurat object |

***

### `CalculateQC()`

Computes common single-cell QC metrics and appends them directly to `SeuratObject@meta.data`.

**Always-computed columns** (from raw counts):

- `percent.mt`: The percentage of mitochondrial gene counts per cell, computed from the `^MT-` gene pattern (human) or `^mt-` (mouse).
- `percent.ribo`: The percentage of ribosomal gene counts per cell, computed from the `^RP[SL]` gene pattern (human) or `^Rp[sl]` (mouse).
- `percent.hb`: The percentage of hemoglobin gene counts per cell, computed from the `^HB[^(P)]"` gene pattern (human) or `^Hb^(p)]"` (mouse).
- `percent.ig`: The percentage of immunoglobulin gene counts per cell, computed from the `^IG` gene pattern (human) or `^Ig` (mouse).
- `percent.plat`: The percentage of platelet gene counts per cell, computed from the `^PPBP|^PF4` gene pattern (human) or `^Ppbp|^Pf4` (mouse).
- `percent.MALAT1`: The percentage of MALAT1 gene counts per cell, computed from the `MALAT1` gene pattern (human) or `Malat1` (mouse).
- `percent.S100A9`: The percentage of S100A9 gene counts per cell, computed from the `S100A9` gene pattern (human) or `S100a9` (mouse).
- `percent.S100A8`: The percentage of S100A8 gene counts per cell, computed from the `S100A8` gene pattern (human) or `S100a8` (mouse).
- `percent.FCGR3B`: The percentage of FCGR3B gene counts per cell, computed from the `FCGR3B` gene pattern (human) or `Fcgr3b` (mouse).
- `log10_nFeature_RNA`: The log10 of the number of features per cell, computed from the `nFeature_RNA` column.
- `log10_nCount_RNA`: The log10 of the number of counts per cell, computed from the `nCount_RNA` column.
- `complexity`: The complexity of the library per cell, computed as `log10(nFeature_RNA) / log10(nCount_RNA)`.

**Conditionally computed** (requires normalized counts in `data` layer):

- `nCount_logRNA`, `nFeature_logRNA`: log-normalized count and detected feature totals per cell.
- Cell cycle scoring (S and G2M scores, `Phase`) — when `perform.cell.cycle.scoring = TRUE`.
- MALAT1 thresholding — when `perform.MALAT1.test = TRUE`; adds `MALAT1.threshold` (numeric) and `MALAT1.pass` (logical) columns. See the [feature_threshold](https://github.com/jmgs7/feature_threshold) repository for details.

The function supports split Seurat v5 objects: if normalized counts are stored in standard `data.*` layers or a concrete layer vector is passed via `layers`, metrics are computed per-layer.

```r
SeuratObject <- CalculateQC(
  SeuratObject,
  assay = "RNA",
  layers = NULL,
  species = "human" , # Also 'mouse' available
  perform.cell.cycle.scoring = TRUE,
  perform.MALAT1.test = TRUE
)
```

| Parameter | Default | Description |
| :-- | :-- | :-- |
| `SeuratObject` | — | Seurat object to annotate |
| `assay` | `"RNA"` | Assay to read counts from |
| `layers` | `NULL` | Normalized layer(s) to use; `NULL` uses the default `data` layer; a character vector processes each layer independently |
| `species` | `"human"` | Species of the dataset, used to select the appropriate gene sets for cell cycle scoring. Human and Mouse currently supported. |
| `perform.cell.cycle.scoring` | `TRUE` | If `TRUE`, performs cell cycle scoring |
| `perform.MALAT1.test` | `TRUE` | If `TRUE`, performs MALAT1-based QC thresholding |
| `...` | | Additional arguments forwarded to the internal `.CalculateFeatureThresholdSeurat()` |

***

### `CellsHistoGradient()`

Creates a bar plot showing the number of cells per group, with bars filled by a viridis gradient based on cell count.

```r
CellsHistoGradient(
  SeuratObject,
  group.by = NULL,
  scale.colors = "viridis",
  breaks = scales::extended_breaks()
)
```

| Parameter | Default | Description |
| :-- | :-- | :-- |
| `SeuratObject` | — | A Seurat object |
| `group.by` | `NULL` | Metadata column for grouping; if `NULL`, uses active identity |
| `scale.colors` | `"viridis"` | Viridis palette option: `"magma"` (`"A"`), `"inferno"` (`"B"`), `"plasma"` (`"C"`), `"viridis"` (`"D"`), `"cividis"` (`"E"`), `"rocket"` (`"F"`), `"mako"` (`"G"`), `"turbo"` (`"H"`) |
| `breaks` | `scales::extended_breaks()` | Y-axis break specification |


***

### `FeatureScatterGradient()`

Extends `Seurat::FeatureScatter()` by coloring points with a continuous viridis gradient defined by a third feature. Supports optional per-group panels with independent per-group correlation subtitles, configurable correlation methods, per-feature layer selection for all three axes, and a single shared legend across grouped panels. The gradient scale is computed globally across all panels.

```r
FeatureScatterGradient(
  SeuratObject,
  feature1,
  feature2,
  gradient,
  group.by = "ident",
  scale.colors = "viridis",
  lower.limit = NULL,
  upper.limit = NULL,
  corr.method = "pearson",
  layer1 = NULL,
  layer2 = NULL,
  layer.gradient = NULL,
  plot.title = NULL,
  pt.size = 0.5,
  common.scales = TRUE,
  collect.axes = FALSE
)
```

| Parameter | Default | Description |
| :-- | :-- | :-- |
| `feature1` | — | X-axis feature (metadata column, assay feature, or reduction variable) |
| `feature2` | — | Y-axis feature; same resolution rules as `feature1` |
| `gradient` | — | Feature used for the continuous color scale, resolved via `FetchData()` |
| `group.by` | `"ident"` | Metadata column for per-group panels; `NULL` for a single combined plot |
| `scale.colors` | `"viridis"` | Viridis palette option or letter code |
| `lower.limit` | `NULL` | Lower gradient clamp; `NULL` infers from data |
| `upper.limit` | `NULL` | Upper gradient clamp; `NULL` infers from data |
| `corr.method` | `"pearson"` | Correlation method: `"pearson"`, `"spearman"`, or `"kendall"` |
| `layer1` | `NULL` | Assay layer for `feature1`; `NULL` uses Seurat default; `""` and `"null"` treated as `NULL` |
| `layer2` | `NULL` | Assay layer for `feature2`; semantics mirror `layer1` |
| `layer.gradient` | `NULL` | Assay layer for `gradient`; metadata-backed features ignore this |
| `plot.title` | `NULL` | Custom main title; ungrouped plots default to the correlation string |
| `pt.size` | `0.5` | Point size |
| `common.scales` | `TRUE` | If `TRUE`, all sub-plots share the same x and y axis limits |
| `collect.axes` | `FALSE` | If `TRUE`, collects axes across sub-plots using patchwork |


***

### `VlnPlotGradient()`

Creates violin plots colored by a per-identity aggregated gradient (`"nCells"` or the mean of a chosen feature), ordered by gradient value. Supports per-feature assay layer selection via `layer` and optional per-feature custom titles via `plot.title`. When `layer` is specified for assay-backed features, the default panel title is `feature_layer`; metadata-backed features always use bare names.

```r
VlnPlotGradient(
  SeuratObject,
  features,
  gradient = "nCells",
  group.by = "ident",
  scale.colors = "viridis",
  lower.limit = 0,
  upper.limit = NULL,
  pt.size = 0,
  ncol = NULL,
  layer = NULL,
  layer.gradient = NULL,
  plot.title = NULL
)
```

| Parameter | Default | Description |
| :-- | :-- | :-- |
| `SeuratObject` | — | A Seurat object |
| `features` | — | Character vector of features (genes or metadata columns) to plot |
| `gradient` | `"nCells"` | `"nCells"` or a feature name whose per-group mean colors the violins |
| `group.by` | `"ident"` | Metadata grouping column; `NULL` uses active identity |
| `scale.colors` | `"viridis"` | Viridis palette option |
| `lower.limit` | `0` | Lower gradient clamp |
| `upper.limit` | `NULL` | Upper gradient clamp |
| `pt.size` | `0` | Size of jittered points overlaid on violins; `0` hides points |
| `ncol` | `NULL` | Number of columns in the combined plot |
| `layer` | `NULL` | Assay layer(s): scalar applies to all features; vector of length equal to `features` applies per-feature; metadata-backed features ignore this |
| `layer.gradient` | `NULL` | Assay layer for the `gradient` feature; metadata-backed features ignore this |
| `plot.title` | `NULL` | Custom title(s): `NULL`, one string for all, or one string per feature |


***

### `FeatureDensityPlot()`

Creates kernel-density plots for one or more features extracted from a Seurat object via `FetchData()`, supporting metadata columns, assay features from any layer, and reduction variables. Features grouped by `group.by` can be shown as overlaid densities or split into per-group panels. Supports dashed red reference lines (`vline`), independent black median lines (`plot.median`), custom titles, and per-feature layer selection. When multiple features are requested, returns a named list of plots; names follow `feature_layer` for assay-backed features and bare `feature` for metadata-backed features (suffixes `_NULL` and `_NA` are stripped).

```r
FeatureDensityPlot(
  SeuratObject,
  features,
  group.by = "ident",
  split.plot = TRUE,
  scale.colors = "viridis",
  ncol = NULL,
  vline = NULL,
  layer = NULL,
  plot.median = TRUE,
  plot.title = NULL,
  nmad = 2,
  alpha = 0.3,
  pt.size = 0,
  common.scales = TRUE,
  collect.axes = FALSE
)
```

| Parameter | Default | Description |
| :-- | :-- | :-- |
| `SeuratObject` | — | A Seurat object |
| `features` | — | Character vector of metadata columns, assay features, or reduction variables |
| `group.by` | `"ident"` | Grouping variable: metadata column, `"ident"` (active identity), or `NULL` for no grouping |
| `split.plot` | `TRUE` | If `TRUE`, one panel per group level; if `FALSE`, overlay all groups in one panel per feature |
| `scale.colors` | `"viridis"` | Viridis palette option for group color mapping |
| `ncol` | `NULL` | Number of columns for split panels within each feature plot |
| `vline` | `NULL` | Reference line: `NULL`, `"mean"`, `"median"`, `"upper"` (median + `nmad` MADs), `"lower"` (median − `nmad` MADs), `"both"`, a numeric value, or a per-feature/per-group character/numeric vector |
| `layer` | `NULL` | Assay layer(s): scalar (all features) or vector (per-feature); metadata features ignore this |
| `plot.median` | `TRUE` | If `TRUE`, draws an independent black dashed median line for each group |
| `plot.title` | `NULL` | Custom title(s): `NULL` uses feature names, one string for all, or one per feature |
| `nmad` | `2` | Number of MADs used when `vline` is `"upper"`, `"lower"`, or `"both"` |
| `alpha` | `0.3` | Fill transparency for density geometries |
| `pt.size` | `0` | Rug line width; `0` disables the rug |
| `common.scales` | `TRUE` | If `TRUE`, all sub-plots share the same x and y axis limits |
| `collect.axes` | `FALSE` | If `TRUE`, collects axes across sub-plots using patchwork |


***

### `QCMetricsBoxplot()`

Creates boxplots to visualize the distribution of QC metrics (default: `nFeature_RNA`, `nCount_RNA`, `percent.mt`) grouped by an entity metadata column. Individual cells are overlaid as jitter points, optionally colored by a gradient feature. Multiple metrics are combined into a single patchwork figure.

```r
QCMetricsBoxplot(
  SeuratObject,
  entity_type,
  entity_name = NULL,
  qc_metrics = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  gradient_col = NULL,
  scale.colors = "viridis",
  lower.limit = 0,
  upper.limit = NULL,
  pt.size = 1,
  pt.alpha = 0.6,
  fill_color = "lightblue",
  outlier.size = 1,
  ncol = NULL,
  gradient_legend = FALSE
)
```

| Parameter | Default | Description |
| :-- | :-- | :-- |
| `SeuratObject` | — | A Seurat object |
| `entity_type` | — | Metadata column name used to group cells (e.g. `"orig.ident"`) |
| `entity_name` | `NULL` | If specified, filters the plot to one value of `entity_type` |
| `qc_metrics` | `c("nFeature_RNA", "nCount_RNA", "percent.mt")` | QC metadata columns to plot |
| `gradient_col` | `NULL` | Metadata column or feature used to color jitter points; if `NULL`, the current metric is used |
| `scale.colors` | `"viridis"` | Viridis palette option |
| `lower.limit` | `0` | Lower gradient scale clamp |
| `upper.limit` | `NULL` | Upper gradient scale clamp |
| `pt.size` | `1` | Size of individual cell points |
| `pt.alpha` | `0.6` | Transparency of individual cell points |
| `fill_color` | `"lightblue"` | Fill color of the boxplot body |
| `outlier.size` | `1` | Size of boxplot outlier points |
| `ncol` | `NULL` | Number of columns in the combined plot |
| `gradient_legend` | `FALSE` | Whether to show the gradient legend |


***

### `scGSEAmarkers()`

Runs **fgsea** per cluster using the output of `Seurat::FindAllMarkers()` and returns a list of enriched pathways for each cluster. Genes are ranked by `avg_log2FC`; only genes passing the `padj.threshold` are included in the ranked list.

```r
results <- scGSEAmarkers(
  cluster_markers,
  reference_markers,
  padj.threshold = 1e-6,
  only.pos = TRUE,
  workers = 4
)
```

| Parameter | Default | Description |
| :-- | :-- | :-- |
| `cluster_markers` | — | `data.frame` from `FindAllMarkers()` containing `cluster`, `gene`, `avg_log2FC`, and `p_val_adj` columns |
| `reference_markers` | — | GSEA gene-set database in `fgsea`-compatible list format |
| `padj.threshold` | `1e-6` | Adjusted p-value cutoff for marker inclusion |
| `only.pos` | `TRUE` | If `TRUE`, returns only positively enriched pathways (NES > 0) |
| `workers` | `4` | Number of cores for parallel fgsea computation |


***

### `ExtractFeatureTestResults()`

Extracts the per-cell results of `.CalculateFeatureThresholdSeurat()` for a given feature from a Seurat object's metadata and returns a tidy `data.frame` with one row per unique batch level, including the computed threshold and pass/fail statistics.

```r
ExtractFeatureTestResults(
  SeuratObject,
  batch.col,
  feature = "MALAT1"
)
```

| Parameter | Default | Description |
| :-- | :-- | :-- |
| `SeuratObject` | — | Seurat object containing `<feature>.threshold` and `<feature>.pass` columns in `meta.data` |
| `batch.col` | — | Metadata column name identifying the batch variable |
| `feature` | `"MALAT1"` | Feature name whose threshold results to extract |


***

## Internal functions

The following functions are internal helpers not exported from the namespace. They are documented here for reference.

### `.ComputeFeatureThreshold()`

Computes an expression threshold for a single feature from a normalized counts vector using kernel density estimation and smoothing splines. The threshold corresponds to the left quadratic intercept of the density curve, with fallback to a conservative default value.

### `.CalculateFeatureThresholdSeurat()`

Applies `.ComputeFeatureThreshold()` across all layers of a Seurat assay and stores per-cell results as `<feature>.threshold` (numeric) and `<feature>.pass` (logical) in `meta.data`.

### `.ConvertEnsembleToSymbol2()`

Converts ENSEMBL gene IDs to gene symbols for a raw count matrix using `biomaRt`. Returns the filtered matrix with symbol rownames. Requires `biomaRt` and `dplyr`; not exported.

### `.SetCommonScales()`

Sets common x and y axis limits across a `patchwork` joined ggplot object. Used internally by `FeatureScatterGradient()` and `FeatureDensityPlot()` when `common.scales = TRUE`.

***

## Dependencies

| Package | Source |
| :-- | :-- |
| Seurat | CRAN  |
| SeuratObject | CRAN |
| ggplot2 | CRAN |
| dplyr | CRAN |
| patchwork | CRAN |
| fgsea | Bioconductor |
| BiocParallel | Bioconductor |
| tibble | CRAN |
| scales | CRAN |
| BPCells | GitHub (`bnprks/BPCells`) |
| stats | base R |

***

## Changelog

### 0.8.4-beta

- The `FeatureDensityPlot()` and `FeatureScatterGradient()`functions now set a common x and y axes scale by default. You can override this behavior by setting the `common.scales` argument to `FALSE`. This enhancement allows for better comparison of feature distributions across different groups or conditions, as it ensures that the scales are consistent and comparable. They also include a new argument `collect.axes` (default `FALSE`) to draw a unique x and y axis per plot when `split.plot = TRUE`. This provides users with more flexibility in visualizing their data, allowing for clearer interpretation of feature distributions and relationships across different groups or conditions.
- Added support for mouse genome in `CalculateQC()` function. The function now accepts a `species` argument, which can be set to either `"human"` (default) or `"mouse"`. This allows users to perform quality control calculations on datasets from both human and mouse samples, ensuring that the appropriate gene sets are used for cell cycle scoring and other QC metrics.

### 0.8.3 (beta)

- The `FeatureDensityPlot()` and `FeatureScatterGradient()`functions now set a common x and y axes scale by default. You can override this behavior by setting the `common.scales` argument to `FALSE`. This enhancement allows for better comparison of feature distributions across different groups or conditions, as it ensures that the scales are consistent and comparable. They also include a new argument `collect.axes` (default `FALSE`) to draw a unique x and y axis per plot when `split.plot = TRUE`. This provides users with more flexibility in visualizing their data, allowing for clearer interpretation of feature distributions and relationships across different groups or conditions.
- Added support for mouse genome in `CalculateQC()` function. The function now accepts a `species` argument, which can be set to either `"human"` (default) or `"mouse"`. This allows users to perform quality control calculations on datasets from both human and mouse samples, ensuring that the appropriate gene sets are used for cell cycle scoring and other QC metrics.

### 0.8.2 (beta)

- Introduced `QCMetricsBoxplot()` for grouped QC boxplots with gradient-colored jitter overlays.
- Added `...` forwarding in `CalculateQC()` to pass custom arguments to the internal `.CalculateFeatureThresholdSeurat()`.
- Added `ExtractFeatureTestResults()` (Utils module) to retrieve and summarize per-cell MALAT1 threshold test results.
- Enhanced `FeatureDensityPlot()` to accept a vector of `vline` values per group when plotting a single feature with `split.plot = TRUE`.
- Refactored `.ComputeFeatureThreshold()` and `.CalculateFeatureThresholdSeurat()` as documented internal helpers.

### 0.8.1

- Uncluttered namespace conflicts between `dplyr` and `stats`.

### 0.8.0

- `CalculateQC()` now computes `nCount_logRNA` and `nFeature_logRNA` from normalized counts.
- Added cell cycle scoring support to `CalculateQC()`.
- Added MALAT1-based QC thresholding to `CalculateQC()`, producing `MALAT1.threshold` and `MALAT1.pass` metadata columns.
- Fixed layer-related bugs in plotting functions.

### 0.7.1

- Hardened `FeatureScatterGradient()` argument validation and single-layer handling, treating `""` and `"null"` as omitted layers and enforcing numeric gradient limit ordering when both bounds are specified.
- Improved grouped `FeatureScatterGradient()` behaviour: NA grouping values are excluded, per-group correlations appear as subtitles, and a single global viridis gradient scale with shared legend is assembled via `patchwork::plot_layout(guides = "collect")`.
- Updated `FeatureDensityPlot()` to use `group.by = "ident"` by default for direct `FetchData()` compatibility; output names cleaned up to use `feature_layer` for assay-backed features and bare feature names for metadata-backed features.
- Enhanced `VlnPlotGradient()` to support per-feature custom plot titles via `plot.title` while maintaining `feature_layer` naming for assay-backed features when a layer is specified.
- Fixed several plotting edge cases: handling NAs and empty strings in all layer options as `NULL`, avoiding errors when using default grouping by identity, and eliminating duplicate legends in grouped gradient plots.

### 0.7.0 (beta)

- Roxygenized the package and aligned documentation for the main plotting utilities and QC helpers.
- Introduced the gradient-based visualization suite (`FeatureScatterGradient`, `FeatureDensityPlot`, `VlnPlotGradient`, `QCMetricsBoxplot`) with consistent viridis-based palettes and Seurat-compatible semantics.


***

## Credits

- `CalculateCDR` adapted from [conquer_comparison](https://github.com/csoneson/conquer_comparison/blob/master/scripts/apply_MASTcpmDetRate.R) by Charlotte Soneson.

 - MALAT1 thresholding logic adapted from [BaderLab's MALAT1_threshold] (https://github.com/BaderLab/MALAT1_threshold)

## License

GPLv3 © José Manuel Gómez Silva

## Contact and contributions

If you find bugs or have suggestions, please open an issue or pull request on GitHub.

When contributing code:

- Follow the existing style (based partially on Google's R style guide: CamelCase for exported functions, dot.case for arguments/variables, explicit `package::function` calls).
- Include clear roxygen documentation and tests where appropriate.
- Respect licensing and attribution, particularly for any upstream code you incorporate.
