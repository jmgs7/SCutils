# SCutils

An R package providing single-cell RNA-seq utility functions to complement Seurat-based workflows.

## Installation

```r
# install.packages("devtools")
devtools::install_github("jmgs7/SCutils")
```

## Functions

### `BatchOpenH5(files, BP.data.dir, platform, ensembl.to.symbol, species, generate.metadata, mc.cores)`
Batch-opens `.h5*` files, writes/loads BPCells matrix directories, optionally converts ENSEMBL IDs to gene symbols or uses genes symbols if available in the h5 matrix data. Can generate sample provenance metadata. File processing is parallelized with `parallel::mclapply()` (Windows falls back to 1 core).

```r
BatchOpenH5(
  files = c("/data/sample1.h5ad", "/data/sample2.h5ad"),
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
| `files` | — | Character vector of `.h5*` file paths |
| `BP.data.dir` | `NULL` | Output directory for `*_BP` folders; defaults to first file directory |
| `platform` | `"10X"` | Input format: `"10X"` or `"anndata"` |
| `use.names` | `TRUE` | Use gene symbols instead of IDs. Requires 10X object and names to be in the default /matrix/features/name directory |
| `ensembl.to.symbol` | `FALSE` | Convert ENSEMBL IDs to symbols with Azimuth |
| `species` | `"human"` | Species for ID conversion |
| `generate.metadata` | `FALSE` | Return sample provenance metadata (`cell.tag`, `sample.procedence`) |
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

### `FeatureScatterGradient(SeuratObject, feature1, feature2, gradient, upper.limit, lower.limit, scale.colors)`
Extends feature scatter plots by coloring points with a third feature gradient and reports Pearson correlation.

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
| `upper.limit` | `NULL` | Optional upper gradient clamp |
| `lower.limit` | `0` | Lower gradient clamp (used with `upper.limit`) |
| `scale.colors` | `"viridis"` | Viridis palette option |

---

### `VlnPlotGradient(SeuratObject, features, gradient, group.by, scale.colors, lower.limit, upper.limit, pt.size, ncol)`
Creates violin plots colored by a per-identity aggregated gradient (`nCells` or mean of selected feature), ordered by gradient value.

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
| `group.by` | `NULL` | Metadata grouping column; if `NULL`, active identity |
| `scale.colors` | `"viridis"` | Viridis palette option |
| `lower.limit` | `0` | Lower gradient clamp |
| `upper.limit` | `NULL` | Upper gradient clamp |
| `pt.size` | `0.1` | Size of jittered points on violins (`0` hides points) |
| `ncol` | `NULL` | Number of columns in combined plot |

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

MIT © José Manuel Gómez Silva
