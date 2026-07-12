#' @title FeatureScatterGradient
#'
#' @description
#' `FeatureScatterGradient()` extends `Seurat::FeatureScatter()` by coloring
#' each point according to a continuous gradient defined by a third feature,
#' while optionally computing correlations per group and faceting by group.
#'
#' The function is designed to:
#' - Preserve the core behaviour of `FeatureScatter`: scatter of two features
#'   (`feature1` on the x-axis and `feature2` on the y-axis), correlation
#'   annotation, and Seurat-like aesthetics.
#' - Use `Seurat::FetchData()` for all data access, with optional per-feature
#'   layer selection (`layer1`, `layer2`) to support Seurat v5 assay layers.
#' - Provide grouping via `group.by`, consistent with `VlnPlotGradient()` and
#'   `FeatureDensityPlot()` (grouping resolved from metadata or identity, without
#'   mutating `Seurat::Idents()`).
#' - Compute and display correlation values either globally or independently per
#'   group, using a configurable correlation method (`corr.method`).
#'
#' When `group.by` is `NULL`, a single scatter plot is returned, with the title
#' set to the global correlation string (e.g. `"Pearson = 0.85"`). When
#' `group.by` is non-`NULL`, the function returns a combined patchwork of
#' per-group scatter panels, with a main title summarising the feature pairing
#' (or an optional `plot.title`) and per-panel titles showing the group-specific
#' correlation values.
#'
#' Axis labels follow `feature_layer` naming for assay-backed features when a
#' layer is specified, and plain `feature` naming for metadata-backed features,
#' matching the naming conventions used in `VlnPlotGradient()` and
#' `FeatureDensityPlot()`.
#'
#' @param SeuratObject A Seurat object.
#' @param feature1 Character scalar. Name of the feature to plot on the x-axis.
#'   Can be a metadata column, assay feature, dimensional reduction variable, or
#'   any variable resolvable by `Seurat::FetchData()`.
#' @param feature2 Character scalar. Name of the feature to plot on the y-axis.
#'   Same resolution rules as for `feature1`.
#' @param gradient Character scalar. Name of the feature whose values will be
#'   used to set the color gradient of the points. Resolved via `Seurat::FetchData()`.
#' @param group.by Character scalar or `NULL`. Grouping variable used to create
#'   per-group panels when non-`NULL`. If `NULL`, no grouping is applied and a
#'   single scatter plot is returned. When non-`NULL`, must be either a metadata
#'   column name or `"ident"`, which refers to the active identity classes and
#'   is resolved via `Seurat::FetchData()` using the special `"ident"` keyword.
#' @param scale.colors Character scalar. Viridis palette option used for the
#'   gradient color scale (`"magma"`, `"inferno"`, `"plasma"`, `"viridis"`,
#'   `"cividis"`, `"rocket"`, `"mako"`, `"turbo"` or their letter codes). Default
#'   is `"viridis"`.
#' @param lower.limit Numeric scalar or `NULL`. Lower limit of the gradient
#'   color scale. Only applied when `upper.limit` is non-`NULL`. When `NULL`,
#'   the lower limit is inferred automatically by ggplot when `upper.limit` is
#'   also `NULL`.
#' @param upper.limit Numeric scalar or `NULL`. Upper limit of the gradient
#'   color scale. When `NULL`, the gradient limits are set automatically by
#'   ggplot. When non-`NULL`, `lower.limit` must also be non-`NULL` and strictly
#'   smaller than `upper.limit`.
#' @param corr.method Character scalar. Method used to compute the correlation
#'   between `feature1` and `feature2`. Accepted values are `"pearson"`
#'   (default), `"spearman"`, and `"kendall"`. Passed to `stats::cor()` as the
#'   `method` argument.
#' @param layer1 Character scalar or `NULL`. Assay layer from which `feature1`
#'   should be obtained (e.g. `"counts"`, `"data"`, `"scale.data"`). When
#'   `NULL`, Seurat's default layer is used. Empty strings and `"null"`
#'   (case-insensitive) are treated as `NULL`. Metadata-backed features ignore
#'   `layer1` at the data access level.
#' @param layer2 Character scalar or `NULL`. Assay layer from which `feature2`
#'   should be obtained. Semantics mirror `layer1`. Metadata-backed features
#'   ignore `layer2`.
#' @param plot.title Character scalar or `NULL`. Optional custom main title for
#'   grouped plots. When `group.by` is `NULL`, this argument is ignored and the
#'   plot title is always set to the global correlation string. When `group.by`
#'   is non-`NULL`, this string is used as the combined patchwork title if
#'   provided; otherwise, the default main title is the concatenation
#'   `"feature1(_layer1) VS feature2(_layer2)"`.
#' @param pt.size Numeric scalar. Point size for the scatter plot. Default is
#'   `0.5`.
#'
#' @return
#' If `group.by` is `NULL`, returns a single `ggplot2` object representing the
#' scatter of `feature1` vs `feature2` with points colored by `gradient` and a
#' title displaying the global correlation value.
#'
#' If `group.by` is non-`NULL`, returns a combined `ggplot2`/`patchwork` object
#' with one panel per group level. Each panel shows the scatter of `feature1`
#' vs `feature2` for that group, colored by `gradient`, and a panel title
#' displaying the group-specific correlation value. The combined plot has a
#' main title and shared legend.
#'
#' @details
#' **Data access**: All variables (`feature1`, `feature2`, `gradient`, and
#' `group.by`) are retrieved via `Seurat::FetchData()` using a common cell
#' order (`colnames(SeuratObject)`), ensuring alignment across features. Layers
#' for `feature1` and `feature2` are controlled by `layer1` and `layer2`.
#' When grouping is requested, cells with missing `group.by` values are excluded
#' before correlations and panels are computed.
#'
#' **Metadata vs assay naming**: Axis labels and the default grouped main title
#' follow the convention `feature_layer` for assay-backed features when a layer
#' is specified, and plain `feature` for metadata-backed features. Trailing
#' suffixes such as `"_"`, `"_NA"`, or `"_NULL"` are stripped for robustness,
#' matching naming safeguards used in other SCutils plotting functions.
#'
#' **Grouping and correlation**: When `group.by` is `NULL`, correlation is
#' computed across all cells using the selected `corr.method` and displayed on
#' the plot. When `group.by` is non-`NULL`, correlation is recomputed separately
#' for each group level and displayed in the title of each panel. Groups with
#' insufficient non-NA data yield an `"NA"` correlation value.
#'
#' **Aesthetics**: Ungrouped plots mimic Seurat's `FeatureScatter` aesthetics
#' (single panel, Seurat-like theme, point coloring), while grouped plots are
#' styled to be cohesive with `FeatureDensityPlot()`
#' (consistent fonts, legend placement, and patchwork layout).
#'
#' @examples
#' \dontrun{
#' # Global correlation of QC metrics colored by percent.mt
#' FeatureScatterGradient(
#'   SeuratObject = SeuratObject,
#'   feature1     = "nCount_RNA",
#'   feature2     = "nFeature_RNA",
#'   gradient     = "percent.mt"
#' )
#'
#' # Grouped scatter by identity, with per-group correlations
#' FeatureScatterGradient(
#'   SeuratObject = SeuratObject,
#'   feature1     = "nCount_RNA",
#'   feature2     = "nFeature_RNA",
#'   gradient     = "percent.mt",
#'   group.by     = "ident",
#'   plot.title   = "QC metrics across identities"
#' )
#'
#' # Assay-backed features from different layers, colored by metadata
#' FeatureScatterGradient(
#'   SeuratObject = SeuratObject,
#'   feature1     = "CD3D",
#'   feature2     = "CD3E",
#'   gradient     = "percent.mt",
#'   layer1       = "counts",
#'   layer2       = "data",
#'   corr.method  = "spearman"
#' )
#' }
#'
#' @import Seurat
#' @import ggplot2
#' @import patchwork
#' @export
FeatureScatterGradient <- function(
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
  plot.title = NULL,
  pt.size = 0.5
) {
  # Validate Seurat object: required for FetchData semantics and structure.
  if (!inherits(SeuratObject, "Seurat")) {
    stop("'SeuratObject' must be a Seurat object.")
  }

  # Validate feature1: must be a single non-NA character.
  if (!is.character(feature1) || length(feature1) != 1L || is.na(feature1)) {
    stop("'feature1' must be a single non-NA character value.")
  }

  # Validate feature2: must be a single non-NA character.
  if (!is.character(feature2) || length(feature2) != 1L || is.na(feature2)) {
    stop("'feature2' must be a single non-NA character value.")
  }

  # Validate gradient: must be a single non-NA character.
  if (!is.character(gradient) || length(gradient) != 1L || is.na(gradient)) {
    stop("'gradient' must be a single non-NA character value.")
  }

  # Validate group.by: NULL or a single non-NA character.
  if (!is.null(group.by)) {
    if (!is.character(group.by) || length(group.by) != 1L || is.na(group.by)) {
      stop("'group.by' must be NULL or a single non-NA character value.")
    }
  }

  # Validate pt.size: single non-negative numeric.
  if (
    !is.numeric(pt.size) ||
      length(pt.size) != 1L ||
      is.na(pt.size) ||
      pt.size < 0
  ) {
    stop("'pt.size' must be a single non-negative numeric value.")
  }

  # Validate lower.limit: NULL or single numeric.
  if (!is.null(lower.limit)) {
    if (
      !is.numeric(lower.limit) ||
        length(lower.limit) != 1L ||
        is.na(lower.limit)
    ) {
      stop("'lower.limit' must be NULL or a single non-NA numeric value.")
    }
  }

  # Validate upper.limit: NULL or single numeric, and greater than lower.limit when both are set.
  if (!is.null(upper.limit)) {
    if (
      !is.numeric(upper.limit) ||
        length(upper.limit) != 1L ||
        is.na(upper.limit)
    ) {
      stop("'upper.limit' must be NULL or a single non-NA numeric value.")
    }
    if (is.null(lower.limit)) {
      stop("'lower.limit' must be specified when 'upper.limit' is non-NULL.")
    }
    if (upper.limit <= lower.limit) {
      stop("'upper.limit' must be strictly greater than 'lower.limit'.")
    }
  }

  # Validate corr.method: must be one of the allowed correlation methods.
  valid.methods <- c("pearson", "spearman", "kendall")
  if (
    !is.character(corr.method) ||
      length(corr.method) != 1L ||
      is.na(corr.method)
  ) {
    stop("'corr.method' must be a single non-NA character value.")
  }
  corr.method.lower <- tolower(corr.method)
  if (!corr.method.lower %in% valid.methods) {
    stop(sprintf(
      "'corr.method' must be one of '%s'.",
      paste(valid.methods, collapse = "', '")
    ))
  }

  # Validate scale.colors: must be a known viridis option or letter code.
  valid.viridis <- c(
    "magma",
    "inferno",
    "plasma",
    "viridis",
    "cividis",
    "rocket",
    "mako",
    "turbo",
    "A",
    "B",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H"
  )
  if (
    !is.character(scale.colors) ||
      length(scale.colors) != 1L ||
      is.na(scale.colors)
  ) {
    stop("'scale.colors' must be a single non-NA character value.")
  }
  if (!tolower(scale.colors) %in% tolower(valid.viridis)) {
    stop(sprintf(
      "'scale.colors' must be one of %s.",
      paste(valid.viridis, collapse = ", ")
    ))
  }

  # Validate layer1: NULL or single character.
  if (!is.null(layer1)) {
    if (!is.character(layer1) || length(layer1) != 1L || is.na(layer1)) {
      stop("'layer1' must be NULL or a single non-NA character value.")
    }
  }

  # Validate layer2: NULL or single character.
  if (!is.null(layer2)) {
    if (!is.character(layer2) || length(layer2) != 1L || is.na(layer2)) {
      stop("'layer2' must be NULL or a single non-NA character value.")
    }
  }

  # Validate plot.title: NULL or single character.
  if (!is.null(plot.title)) {
    if (
      !is.character(plot.title) || length(plot.title) != 1L || is.na(plot.title)
    ) {
      stop("'plot.title' must be NULL or a single non-NA character value.")
    }
  }

  # Normalize single-layer inputs so callers can pass the same sentinel values
  # already accepted elsewhere in SCutils (`NA` is rejected above because these
  # arguments are scalar-only, but empty strings and "null" should behave like
  # an omitted layer rather than being forwarded to FetchData()).
  normalizeSingleLayer <- function(layer.value) {
    if (is.null(layer.value)) {
      return(NULL)
    }

    layer.value.lower <- tolower(layer.value)
    if (layer.value.lower %in% c("", "null")) {
      return(NULL)
    }

    return(layer.value)
  }

  # Apply layer normalization once so both FetchData access and axis-label
  # construction use the same resolved layer semantics.
  layer1 <- normalizeSingleLayer(layer1)
  layer2 <- normalizeSingleLayer(layer2)

  # Helper to format correlation method label for titles (capitalised).
  formatMethodLabel <- function(method.name) {
    if (method.name == "pearson") {
      return("Pearson")
    }
    if (method.name == "spearman") {
      return("Spearman")
    }
    if (method.name == "kendall") {
      return("Kendall")
    }
    return(method.name)
  }

  # Helper to construct axis labels and title tokens based on metadata vs assay
  # origin and layer specification.
  makeAxisLabel <- function(feature.name, layer.value, metadata.names) {
    # Determine whether the feature is backed by metadata.
    is.metadata <- feature.name %in% metadata.names

    # Metadata-backed features ignore layer suffixes and use bare feature names.
    if (is.metadata) {
      return(feature.name)
    }

    # Assay-backed features: append layer when non-NULL; otherwise use bare name.
    if (is.null(layer.value)) {
      label <- feature.name
    } else {
      label <- paste0(feature.name, "_", layer.value)
    }

    # Strip undesired trailing suffixes such as "_", "_NA", "_NULL".
    label <- gsub("(_NA|_NULL)$", "", label)
    label <- gsub("_$", "", label)

    return(label)
  }

  # Helper to fetch feature values via Seurat::FetchData(), with optional layer
  # and robust error messaging.
  fetchFeatureValues <- function(
    object,
    feature.name,
    layer.value = NULL,
    cells.use
  ) {
    # Call FetchData with appropriate layer and cell order.
    fetched <- tryCatch(
      Seurat::FetchData(
        object = object,
        vars = feature.name,
        cells = cells.use,
        layer = layer.value,
        clean = FALSE
      ),
      error = function(e) {
        stop(sprintf(
          "Cannot fetch feature '%s' (layer '%s'): %s",
          feature.name,
          if (is.null(layer.value)) "NULL" else layer.value,
          e$message
        ))
      }
    )

    # Validate that the requested feature is present in the fetched data.
    if (!feature.name %in% colnames(fetched)) {
      stop(sprintf(
        "The feature '%s' could not be fetched from the object. Check the feature name and layer.",
        feature.name
      ))
    }

    # Extract values aligned to cells.use.
    values <- fetched[cells.use, feature.name]

    # Coerce to numeric explicitly; NA-producing coercions are tolerated.
    if (!is.numeric(values)) {
      suppressWarnings(
        values <- as.numeric(values)
      )
    }

    return(values)
  }

  # Prepare common cell order for all FetchData calls.
  cells.use <- colnames(SeuratObject)

  # Cache metadata names for axis and title naming logic.
  metadata.names <- colnames(SeuratObject@meta.data)

  # Resolve grouping vector when group.by is non-NULL; default to "ident" when
  # grouping is requested but group.by equals "ident".
  group.values <- NULL
  group.levels <- NULL
  if (!is.null(group.by)) {
    # Grouping variable name for FetchData.
    group.fetch.var <- group.by

    # Fetch grouping values via FetchData.
    group.data <- tryCatch(
      Seurat::FetchData(
        object = SeuratObject,
        vars = group.fetch.var,
        cells = cells.use,
        clean = FALSE
      ),
      error = function(e) {
        stop(sprintf(
          "Cannot fetch grouping variable '%s': %s",
          group.fetch.var,
          e$message
        ))
      }
    )

    # Validate presence of grouping column.
    if (!group.fetch.var %in% colnames(group.data)) {
      stop(sprintf(
        "Grouping variable '%s' was not found in the fetched data.",
        group.fetch.var
      ))
    }

    # Extract grouping values as character.
    group.values <- as.character(group.data[cells.use, group.fetch.var])

    # Remove missing grouping assignments before turning the values into a factor.
    # This avoids creating a literal "NA" group panel, which would be misleading
    # and inconsistent with FeatureDensityPlot(), where missing grouping values
    # are excluded before plotting.
    valid.group.cells <- !is.na(group.values)
    cells.use <- cells.use[valid.group.cells]
    group.values <- group.values[valid.group.cells]

    # Determine group levels, ordered as in FeatureDensityPlot (sorted unique).
    group.levels <- sort(unique(group.values))
    if (length(group.levels) == 0L) {
      stop(
        "No groups available to plot after filtering missing grouping values."
      )
    }
    group.values <- factor(group.values, levels = group.levels)
  }

  # Fetch feature1 values via FetchData with optional layer1.
  x.values <- fetchFeatureValues(
    object = SeuratObject,
    feature.name = feature1,
    layer.value = layer1,
    cells.use = cells.use
  )

  # Fetch feature2 values via FetchData with optional layer2.
  y.values <- fetchFeatureValues(
    object = SeuratObject,
    feature.name = feature2,
    layer.value = layer2,
    cells.use = cells.use
  )

  # Fetch gradient values (always from default layer; metadata vs assay handled
  # automatically by FetchData).
  gradient.values <- fetchFeatureValues(
    object = SeuratObject,
    feature.name = gradient,
    layer.value = NULL,
    cells.use = cells.use
  )

  # Construct axis labels for x and y based on metadata vs assay origin and
  # layer specification.
  axis.x.label <- makeAxisLabel(
    feature.name = feature1,
    layer.value = layer1,
    metadata.names = metadata.names
  )
  axis.y.label <- makeAxisLabel(
    feature.name = feature2,
    layer.value = layer2,
    metadata.names = metadata.names
  )

  # Prepare gradient scale limits: only used when upper.limit is non-NULL.
  gradient.limits <- NULL
  if (!is.null(upper.limit)) {
    gradient.limits <- c(lower.limit, upper.limit)
  }

  # Label for the gradient color scale (use raw gradient name).
  gradient.label <- gradient

  # Format correlation method label for titles.
  method.label <- formatMethodLabel(corr.method.lower)

  # Ungrouped case: group.by is NULL; compute global correlation and build a
  # single scatter plot mimicking Seurat::FeatureScatter aesthetics.
  if (is.null(group.by)) {
    # Build data frame for plotting: x, y, gradient.
    df <- data.frame(
      x = x.values,
      y = y.values,
      gradient = gradient.values,
      stringsAsFactors = FALSE
    )

    # Compute global correlation using selected method and pairwise complete obs.
    corr.global <- suppressWarnings(
      stats::cor(
        x.values,
        y.values,
        method = corr.method.lower,
        use = "pairwise.complete.obs"
      )
    )

    # Round correlation to two decimal places.
    corr.global.round <- round(corr.global, 2)

    # Build title string: e.g. "Pearson = 0.85".
    corr.title <- paste0(method.label, " = ", corr.global.round)

    # Build ggplot: scatter of x vs y, colored by gradient.
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(x = x, y = y, color = gradient)
    ) +
      ggplot2::geom_point(size = pt.size, alpha = 0.7) +
      ggplot2::scale_color_viridis_c(
        name = gradient.label,
        option = scale.colors,
        limits = gradient.limits,
        na.value = "grey70"
      ) +
      ggplot2::labs(
        title = corr.title,
        x = axis.x.label,
        y = axis.y.label
      ) +
      ggplot2::theme_classic(base_size = 11) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          face = "bold",
          size = 11,
          hjust = 0.5
        ),
        axis.text = ggplot2::element_text(size = 9),
        axis.title = ggplot2::element_text(size = 10)
      )

    # Return single scatter plot.
    return(p)
  }

  # Grouped case: group.by is non-NULL; compute per-group correlations and
  # build one scatter panel per group level, then combine via patchwork.
  # This matches the grouped outputs and aesthetic cohesion of
  # FeatureDensityPlot() and VlnPlotGradient().

  # Build full data frame with x, y, gradient, and group factor.
  df.full <- data.frame(
    x = x.values,
    y = y.values,
    gradient = gradient.values,
    group = group.values,
    stringsAsFactors = FALSE
  )

  # Prepare list to hold per-group ggplot objects.
  plot.list <- vector(mode = "list", length = length(group.levels))

  # Iterate over group levels to build per-group panels.
  for (i in seq_along(group.levels)) {
    # Current group label.
    current.group <- group.levels[[i]]

    # Subset data frame to current group.
    df.group <- df.full[df.full$group == current.group, , drop = FALSE]

    # Extract x and y values for correlation within this group.
    x.group <- df.group$x
    y.group <- df.group$y

    # Compute group-specific correlation using selected method; handle NA cases.
    corr.group <- suppressWarnings(
      stats::cor(
        x.group,
        y.group,
        method = corr.method.lower,
        use = "pairwise.complete.obs"
      )
    )

    # Round correlation to two decimal places when available.
    corr.group.round <- if (!is.na(corr.group)) {
      round(corr.group, 2)
    } else {
      NA_real_
    }

    # Build per-group title: e.g. "Pearson = 0.85" or "Pearson = NA".
    corr.title.group <- paste0(
      method.label,
      " = ",
      if (is.na(corr.group.round)) "NA" else corr.group.round
    )

    # Build per-group ggplot: scatter of x vs y, colored by gradient.
    p.group <- ggplot2::ggplot(
      df.group,
      ggplot2::aes(x = x, y = y, color = gradient)
    ) +
      ggplot2::geom_point(size = pt.size, alpha = 0.7) +
      ggplot2::scale_color_viridis_c(
        name = gradient.label,
        option = scale.colors,
        limits = gradient.limits,
        na.value = "grey70"
      ) +
      ggplot2::labs(
        title = corr.title.group,
        x = axis.x.label,
        y = axis.y.label
      ) +
      ggplot2::theme_classic(base_size = 9) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          face = "bold",
          size = 9,
          hjust = 0.5
        ),
        axis.text = ggplot2::element_text(size = 7),
        axis.title = ggplot2::element_text(size = 8)
      )

    # Store panel in list.
    plot.list[[i]] <- p.group
  }

  # Determine heuristic number of columns for patchwork layout: use square-root
  # rule, as in FeatureDensityPlot() and related grouped plots.
  ncol <- ceiling(sqrt(length(group.levels)))

  # Build default main title combining axis labels.
  default.main.title <- paste(axis.x.label, "VS", axis.y.label)

  # Use custom plot.title when provided; otherwise use default.
  main.title <- if (is.null(plot.title)) default.main.title else plot.title

  # Combine panels into a single patchwork object with shared legend.
  combined <- patchwork::wrap_plots(plot.list, ncol = ncol) +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "right")

  # Add the grouped main title through patchwork annotation rather than
  # ggplot labels. Attaching it as a ggplot label to the combined patchwork
  # object causes downstream ggplot build/print failures because patchwork's
  # root plot expects title styling to live in the annotation theme.
  combined <- combined +
    patchwork::plot_annotation(
      title = main.title,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          face = "bold",
          size = 11,
          hjust = 0.5
        )
      )
    )

  # Return combined grouped scatter plot.
  return(combined)
}
