## SCutils v0.8.5-beta

### 🛠️ Bug fixes and improvements

- `FeatureDensityPlot()` now sorts plots alphabetically by group as `FeatureScatterGradient()` does, ensuring consistent ordering of plots across different visualizations.
- Now we apply coord_cartesian() instead of scale_x_continuous() and scale_y_continuous() in these functions to avoid issues with axis limits when using `common.scales = TRUE`. This change ensures that the axes are consistent across all plots without losing any data points or changing the density plot kernels.
- We have updated the `percent.plat` metric calculation in `CalculateQC()`. Now it calculates the percentage of counts matching the following platelet-specific genes: `PPBP`, `PF4`, `GP9`, `ITGA2B`, `TUBB1`, `ITGB3`, and `GP1BA`. This addition allows for better identification of platelet contamination in single-cell RNA-seq data.

### 🏷️ Release tag

- GitHub prerelease tag/title: `v0.8.4-beta`.

**NOTE**: This is a beta release, so some features may still be under development or need further testing. We encourage users to test the new features and provide feedback. Please report any issues or bugs encountered during testing to help us improve the package.