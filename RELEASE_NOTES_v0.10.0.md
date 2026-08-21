## SCutils v0.10.0

### 🚀 New features.

**New functions:**

- `VariableFeaturePlot2()`: An enhanced version of Seurat's `VariableFeaturePlot()` that allows for labeling the top variable features on the plot. It also includes additional customization options for the plot aesthetics.

**New internal functions:**

- `FilterVariableFeatures()`: Filters the high variable features of a Seurat object after applying `Seurat::FindVariableFeatures()` to delete low-informative genes such as mitochondrial, ribosomal, antisense and non-coding genes.
- `SelectPCs()`: Estimates the intrinsic dimensionality of a Seurat object based on the PCA embeddings and returns a vector with the most informative PCs.

### 🛠️ Bug fixes and improvements.

- `CalculateQC()` now automatically uses the default assay if `assay = NULL`.