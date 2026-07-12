# SCutils

An R package providing single-cell RNA-seq utility functions to complement Seurat-based workflows.

## Installation

```r
# install.packages("remotes")
remotes::install_github("jmgs7/SCutils")
```

## Functions

### `BatchOpenH5(files, relative, BP.data.dir, platform, use.names, ensembl.to.symbol, species, generate.metadata, mc.cores)`
Batch-opens `.h5` single-cell matrices, writes/loads BPCells matrix directories (`*_BP`), supports both 10X and AnnData HDF5 inputs, and can optionally convert ENSEMBL IDs to gene symbols. It can also generate per-cell sample provenance metadata. Processing uses `parallel::mclapply()` (on Windows, cores are forced to 1).

```r
BatchOpenH5(
  files = c("/data/sample1.h5", "/data/sample2.h5"),
  relative = TRUE,
  BP.data.dir = "/data/BP",
  platform = "10X",
  use.names = TRUE,
  ensembl.to.symbol = FALSE,
  species = "human",
  generate.metadata = TRUE,
  mc.cores = 4
)
```

| Parameter | Default | Description |
|---|---|---|
| `files` | — | Character vector of input `.h5` file paths |
| `relative` | `TRUE` | Store matrix directories as relative paths in each matrix object (`./<basename(BP.data.dir)>/<matrix_dir>`) |
| `BP.data.dir` | `NULL` | Output directory for `*_BP` folders; defaults to `dirname(files[[1]])` |
| `platform` | `"10X"` | Input format: `"10X"` or `"anndata"` |
| `use.names` | `TRUE` | For 10X only: replace feature IDs with `/matrix/features/name` |
| `ensembl.to.symbol` | `FALSE` | Convert ENSEMBL IDs to symbols with `ConvertEnsembleToSymbol2()` (applied only when `use.names = FALSE`) |
| `species` | `"human"` | Species passed to symbol conversion |
| `generate.metadata` | `FALSE` | If `TRUE`, returns `data.list` plus metadata (`cell.tag`, `sample.procedence`) |
| `mc.cores` | `length(files)` | Cores used by `mclapply()` |

---

### `CalculateCDR(SeuratObject)`
Calculates **Cellular Detection Rate** (fraction of detected features per cell), scales it, and adds `CDR` to `SeuratObject@meta.data`.

```r
SeuratObject <- CalculateCDR(SeuratObject)
```

---

### `CalculateQC(SeuratObject)`
Computes common single-cell QC metrics and appends them to `SeuratObject@meta.data`.

Added columns include:
- `percent.mt`, `percent.ribo`, `percent.hb`, `percent.ig`, `percent.plat`
- `percent.MALAT1`, `percent.S100A9`, `percent.S100A8`, `percent.FCGR3B`
- `log10_nFeature_RNA`, `log10_nCount_RNA`, `complexity`

```r
SeuratObject <- CalculateQC(SeuratObject)
```

| Parameter | Default | Description |
|---|---|---|
| `SeuratObject` | — | Seurat object to annotate with QC metadata |

---

### `CellsHistoGradient(SeuratObject, group.by, scale.colors, breaks)`
Creates a bar plot with number of cells per group, filled by a viridis gradient.

```r
CellsHistoGradient(
  SeuratObject,
  group.by = "seurat_clusters",
  scale.colors = "plasma"
)
```

| Parameter | Default | Description |
|---|---|---|
| `group.by` | `NULL` | Metadata grouping column; if `NULL`, uses active identity |
| `scale.colors` | `"viridis"` | Viridis palette option |
| `breaks` | `scales::extended_breaks()` | Y-axis breaks |

---

### `FeatureScatterGradient(SeuratObject, feature1, feature2, gradient, group.by, scale.colors, lower.limit, upper.limit, corr.method, layer1, layer2, layer.gradient, plot.title, pt.size)`
Extends feature scatter plots by coloring points with a continuous gradient defined by a third feature. Supports optional grouping into per-group panels with independent per-group correlation display, configurable correlation method (`pearson`, `spearman`, `kendall`), and per-feature layer selection for all three axes (`layer1`, `layer2`, `layer.gradient`). In grouped mode, the gradient color scale is computed globally across all panels and a single shared legend is collected by patchwork.

```r
FeatureScatterGradient(
  SeuratObject,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA",
  gradient = "percent.mt",
  upper.limit = 100,
  scale.colors = "viridis"
)
```

| Parameter | Default | Description |
|---|---|---|
| `feature1` | — | X-axis feature |
| `feature2` | — | Y-axis feature |
| `gradient` | — | Feature used for color scale |
| `group.by` | `"ident"` | Metadata column or `"ident"` for per-group panels; `NULL` for a single plot |
| `scale.colors` | `"viridis"` | Viridis palette option |
| `lower.limit` | `NULL` | Lower gradient clamp; `NULL` infers from data |
| `upper.limit` | `NULL` | Upper gradient clamp; `NULL` infers from data |
| `corr.method` | `"pearson"` | Correlation method: `"pearson"`, `"spearman"`, or `"kendall"` |
| `layer1` | `NULL` | Assay layer for `feature1`; `NULL` uses Seurat default |
| `layer2` | `NULL` | Assay layer for `feature2`; `NULL` uses Seurat default |
| `layer.gradient` | `NULL` | Assay layer for `gradient`; `NULL` uses Seurat default |
| `plot.title` | `NULL` | Custom main title for grouped plots; ignored in ungrouped mode |
| `pt.size` | `0.5` | Point size |

---

### `VlnPlotGradient(SeuratObject, features, gradient, group.by, scale.colors, lower.limit, upper.limit, pt.size, ncol, layer, plot.title)`
Creates violin plots colored by a per-identity aggregated gradient (`nCells` or mean of selected feature), ordered by gradient value. Supports per-feature assay layer selection via `layer` (mirroring `FeatureDensityPlot()` semantics) and optional per-feature custom titles via `plot.title`. When `layer` is provided for an assay-backed feature, the default panel title is `feature_layer`; metadata-backed features always use bare names.

```r
VlnPlotGradient(
  SeuratObject,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  gradient = "nCells",
  group.by = "seurat_clusters",
  scale.colors = "mako"
)
```

| Parameter | Default | Description |
|---|---|---|
| `features` | — | Features (genes/metadata) to plot |
| `gradient` | — | `"nCells"` or feature used to color violins |
| `group.by` | `"ident"` | Metadata grouping column; if `NULL`, active identity |
| `scale.colors` | `"viridis"` | Viridis palette option |
| `lower.limit` | `0` | Lower gradient clamp |
| `upper.limit` | `NULL` | Upper gradient clamp |
| `pt.size` | `0` | Size of jittered points on violins (`0` hides points) |
| `ncol` | `NULL` | Number of columns in combined plot |
| `layer` | `NULL` | Assay layer(s): scalar (all features) or vector (per-feature) |
| `plot.title` | `NULL` | Custom title(s): `NULL`, one string, or one per feature |

---

### `FeatureDensityPlot(SeuratObject, features, group.by, split.plot, scale.colors, ncol, vline, layer, plot.median, plot.title, nmad, alpha, pt.size)`
Creates density plots for one or more features directly from a Seurat object via `Seurat::FetchData()`, supporting metadata columns, assay features from any layer, and dimensional reduction variables. Supports grouped overlays or split-by-group panels, dashed red reference lines (`vline`), independent median overlays (`plot.median`, black), custom titles, and per-feature layer selection via `layer`. The default grouping variable has changed from `"active.ident"` to `"ident"` for direct compatibility with `FetchData()`. When multiple features are requested, the function returns a named list of plots (one plot per feature); names use `feature_layer` for assay-backed features and bare `feature` for metadata-backed features.

```r
FeatureDensityPlot(
  SeuratObject,
  features = c("percent.mt", "nCount_RNA", "nFeature_RNA"),
  group.by = "batch",
  split.plot = TRUE,
  scale.colors = "viridis",
  vline = "upper",
  plot.median = TRUE,
  plot.title = c("Mitochondrial %", "UMI Counts", "Detected Features"),
  nmad = 2.5,
  alpha = 0.3
)
```

| Parameter | Default | Description |
|---|---|---|
| `features` | — | Metadata columns, assay features, or reduction variables to plot |
| `group.by` | `"ident"` | Grouping variable: metadata column, `"ident"` (active identity), or `NULL` for no grouping |
| `split.plot` | `TRUE` | If `TRUE`, one panel per group level; if `FALSE`, overlay groups in one panel per feature |
| `scale.colors` | `"viridis"` | Viridis palette option for grouped color mapping |
| `ncol` | `NULL` | Number of columns for split panels within each feature plot |
| `vline` | `NULL` | Reference line mode: `NULL`, `"mean"`, `"median"`, `"upper"`, `"lower"`, `"both"`, numeric, or per-feature character vector |
| `layer` | `NULL` | Assay layer(s): scalar (all features) or vector (per-feature); `NULL` uses Seurat default |
| `plot.median` | `TRUE` | Draw independent median line(s) in black |
| `plot.title` | `NULL` | Custom title(s): `NULL`, one title for all features, or one title per feature |
| `nmad` | `2` | Number of MADs used when `vline` is `"upper"`, `"lower"`, or `"both"` |
| `alpha` | `0.3` | Density fill transparency |
| `pt.size` | `0` | Rug line width (`0` disables rugs) |

---

### `QCMetricsBoxplot(SeuratObject, entity_type, entity_name, qc_metrics, gradient_col, scale.colors, lower.limit, upper.limit, pt.size, pt.alpha, fill_color, outlier.size, ncol, gradient_legend)`
Creates QC boxplots (default: `nFeature_RNA`, `nCount_RNA`, `percent.mt`) grouped by an entity column, with outliers overlaid as gradient-colored jitter points.

```r
QCMetricsBoxplot(
  SeuratObject,
  entity_type = "orig.ident",
  entity_name = "Sample_1",
  qc_metrics = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  gradient_col = NULL,
  scale.colors = "plasma",
  pt.size = 1.5,
  pt.alpha = 0.7,
  ncol = 3
)
```

| Parameter | Default | Description |
|---|---|---|
| `entity_type` | — | Metadata column to group by |
| `entity_name` | `NULL` | Optional value in `entity_type` to filter one entity |
| `qc_metrics` | `c("nFeature_RNA", "nCount_RNA", "percent.mt")` | QC metadata columns to plot |
| `gradient_col` | `NULL` | Feature used for point color; if `NULL`, uses current metric |
| `scale.colors` | `"viridis"` | Viridis palette option |
| `lower.limit` | `0` | Lower gradient clamp |
| `upper.limit` | `NULL` | Upper gradient clamp |
| `pt.size` | `1` | Outlier point size |
| `pt.alpha` | `0.6` | Outlier point alpha |
| `fill_color` | `"lightblue"` | Box fill color |
| `outlier.size` | `1` | Boxplot outlier size setting |
| `ncol` | `NULL` | Number of columns in combined plot |
| `gradient_legend` | `FALSE` | Toggle legend positioning behavior |

---

### `scGSEAmarkers(cluster_markers, reference_markers, padj.threshold, only.pos, workers)`
Runs **fgsea** per cluster using `FindAllMarkers()` output and returns enriched pathways for each cluster.

```r
results <- scGSEAmarkers(
  cluster_markers = markers_df,
  reference_markers = msigdb_c8,
  padj.threshold = 1e-6,
  only.pos = TRUE,
  workers = 4
)
```

| Parameter | Default | Description |
|---|---|---|
| `padj.threshold` | `1e-6` | Adjusted p-value cutoff |
| `only.pos` | `TRUE` | Keep only positive NES pathways |
| `workers` | `4` | Cores for parallel fgsea |

---

## Dependencies

| Package | Source |
|---|---|
| Seurat | CRAN / Bioconductor |
| ggplot2 | CRAN |
| dplyr | CRAN |
| patchwork | CRAN |
| fgsea | Bioconductor |
| BiocParallel | Bioconductor |
| tibble | CRAN |
| scales | CRAN |
| BPCells | CRAN / GitHub |
| Azimuth | CRAN / GitHub |
| parallel | base R |

---

## Credits

- `CalculateCDR` adapted from [conquer_comparison](https://github.com/csoneson/conquer_comparison/blob/master/scripts/apply_MASTcpmDetRate.R) by Charlotte Soneson.

## License

GPLv3 © José Manuel Gómez Silva
