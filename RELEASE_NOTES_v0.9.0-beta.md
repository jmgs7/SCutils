## SCutils v0.9.0

### 🚀 New features.

- `FeatureTest()` function: This function allows users to perform feature testing on a Seurat object, enabling the identification of cells that meet specific criteria based on user-defined features and thresholds. For a given feature, users can specify a threshold and an operator (e.g., >, <, <=, ...>) to define the testing criteria. The a new column named `<feature>.pass` is created to indicate whether each cell passes the test. The function also supports multiple features and thresholds, allowing for complex filtering criteria. The results of the feature test are stored in the Seurat object's metadata, and users can extract and summarize these results using the `ExtractFeatureTestResults()` function. This function is very useful in combination with `CalculateQC()` to perform quality control on single-cell RNA-seq data, allows to visualize in a PCA plot or UMAP/tSNE plot the cells that pass or fail the feature test, and can be used to filter out low-quality cells from downstream analyses in a more controlled manner than directly applying filters with Seurat's `subset()` function.

### 🛠️ Bug fixes and improvements.

- `ExtractFeatureTestResults()` now returns the thresholds optionaly to improve compatibility with `FeatureTest()`.
- Improved documentation of multiple functions.