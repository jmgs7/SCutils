# SCutils

An R package providing single-cell RNA-seq utility functions to complement Seurat-based workflows.

## Installation

```r
# Install from GitHub using devtools
# install.packages("devtools")
devtools::install_github("jmgs7/SCutils")
```

## Functions

### `CalculateCDR(SeuratObject)`
Calculates the **Cellular Detection Rate** for each cell in a Seurat object. The CDR is the fraction of features with count > 0, scaled to 0–1. Useful as a covariate in differential expression analyses (e.g., MAST).

```r
SeuratObject <- CalculateCDR(SeuratObject)
# Adds a "CDR" column to SeuratObject@meta.data
```

---

### `CellsHistoGradient(SeuratObject, group.by, scale.colors, breaks)`
Generates a **bar plot** showing the number of cells per identity/group, with bars colored by a viridis gradient based on cell count.

```r
CellsHistoGradient(seurat_obj, group.by = "seurat_clusters", scale.colors = "plasma")
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `group.by` | `NULL` (active ident) | Metadata column to group by |
| `scale.colors` | `"viridis"` | Viridis palette (`"magma"`, `"inferno"`, `"plasma"`, etc.) |
| `breaks` | `scales::extended_breaks()` | Y-axis breaks |

---

### `FeatureScatterGradient(SeuratObject, feature1, feature2, gradient, ...)`
Extends `Seurat::FeatureScatter()` by **coloring dots with a gradient** based on a third feature. Ideal for QC visualization (e.g., nCount vs nFeature colored by % mitochondrial genes).

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
|-----------|---------|-------------|
| `gradient` | — | Feature whose values set the color |
| `upper.limit` | `NULL` | Upper clamp for the gradient scale |
| `lower.limit` | `0` | Lower clamp (only when `upper.limit` set) |
| `scale.colors` | `"viridis"` | Viridis palette |

---

### `VlnPlotGradient(SeuratObject, features, gradient, ...)`
Extends `Seurat::VlnPlot()` by **coloring each violin** with a gradient based on a per-identity aggregated value (e.g., number of cells, mean expression).

```r
VlnPlotGradient(
  SeuratObject,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  gradient = "nCells",
  scale.colors = "viridis"
)
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `gradient` | — | Feature or `"nCells"` to color violins by |
| `group.by` | `NULL` (active ident) | Metadata column to group by |
| `scale.colors` | `"viridis"` | Viridis palette |
| `upper.limit` | `NULL` | Upper clamp for gradient |
| `pt.size` | `0.1` | Size of overlaid jitter points |
| `ncol` | `NULL` | Number of columns in combined plot |

---

### `QCMetricsBoxplot(SeuratObject, entity_name, entity_type, ...)`
Creates **boxplots** to visualize QC metrics (nFeature_RNA, nCount_RNA, percent.mt) for cells within specific entities. Individual cells are displayed as gradient-colored points, allowing rapid visual assessment of QC metric distributions by sample, batch, cluster, or other grouping variables.

```r
QCMetricsBoxplot(
  SeuratObject,
  entity_name = "orig.ident",
  qc_metrics = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  scale.colors = "plasma",
  pt.size = 1.5,
  pt.alpha = 0.7
)
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `entity_name` | — | Metadata column to group by (e.g., "orig.ident", "batch", "seurat_clusters") |
| `entity_type` | `NULL` | If specified, filters to show only cells from this entity; if `NULL`, shows all |
| `qc_metrics` | `c("nFeature_RNA", "nCount_RNA", "percent.mt")` | QC metrics to visualize |
| `gradient_col` | `NULL` | Feature for point gradient coloring; if `NULL`, each metric colors by its own values |
| `scale.colors` | `"viridis"` | Viridis palette (`"magma"`, `"inferno"`, `"plasma"`, etc.) |
| `pt.size` | `1` | Size of points representing individual cells |
| `pt.alpha` | `0.6` | Transparency of points (0-1) |
| `fill_color` | `"lightblue"` | Boxplot fill color |
| `outlier.size` | `1` | Size of boxplot outliers |
| `upper.limit` | `NULL` | Upper clamp for gradient scale |
| `lower.limit` | `0` | Lower clamp (only when `upper.limit` set) |
| `ncol` | `NULL` | Number of columns in combined plot (default: 3) |

---

### `scGSEAmarkers(cluster_markers, reference_markers, ...)`
Runs **fgsea** for every cluster in a Seurat `FindAllMarkers()` output against a reference gene-set database. Useful for automated cluster annotation (e.g., using MSigDB C8).

```r
results <- scGSEAmarkers(
  cluster_markers = markers_df,  # output of Seurat::FindAllMarkers()
  reference_markers = msigdb_c8, # fgsea-formatted gene set list
  padj.threshold = 1e-6,
  only.pos = TRUE,
  workers = 4
)
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `padj.threshold` | `1e-6` | Adjusted p-value cutoff |
| `only.pos` | `TRUE` | Return only positively enriched sets (NES > 0) |
| `workers` | `4` | Cores for parallel computation |

---

## Dependencies

| Package | Source |
|---------|--------|
| Seurat | CRAN / Bioconductor |
| ggplot2 | CRAN |
| dplyr | CRAN |
| patchwork | CRAN |
| fgsea | Bioconductor |
| BiocParallel | Bioconductor |
| tibble | CRAN |
| scales | CRAN |

---

## Credits

- `CalculateCDR` adapted from [conquer_comparison](https://github.com/csoneson/conquer_comparison/blob/master/scripts/apply_MASTcpmDetRate.R) by Charlotte Soneson.

## License

MIT © José Manuel Gómez Silva
