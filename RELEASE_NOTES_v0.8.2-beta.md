## SCutils v0.8.3-beta

### 🚀 New Features.

- The `FeatureDensityPlot()` and `FeatureScatterGradient()`functions now set a common x and y axes scale by default. You can override this behavior by setting the `common.scales` argument to `FALSE`. This enhancement allows for better comparison of feature distributions across different groups or conditions, as it ensures that the scales are consistent and comparable. They also include a new argument `collect.axes` (default `FALSE`) to draw a unique x and y axis per plot when `split.plot = TRUE`. This provides users with more flexibility in visualizing their data, allowing for clearer interpretation of feature distributions and relationships across different groups or conditions.

### 🏷️ Release tag

- GitHub prerelease tag/title: `v0.8.3-beta`.

**NOTE**: This is a beta release, so some features may still be under development or need further testing. We encourage users to test the new features and provide feedback. Please report any issues or bugs encountered during testing to help us improve the package.