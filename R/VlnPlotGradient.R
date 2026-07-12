#' @title VlnPlotGradient
#'
#' @description
#' `VlnPlotGradient()` extends `Seurat::VlnPlot()` functionality by coloring each
#' single-cell violin plot with a continuous gradient that encodes a per-group
#' summary of a feature (for example, number of cells per identity or the mean
#' expression of a marker gene).
#'
#' The function is designed to match the visual aesthetics and combined-output
#' behaviour of `Seurat::VlnPlot()` while providing a more flexible backend:
#' - Data access is based entirely on `Seurat::FetchData()`, which can retrieve
#'   metadata columns, dimensional reduction variables (e.g. `PC_1`), and assay
#'   features from arbitrary layers (such as "counts", "data", or
#'   "scale.data") in Seurat v5 objects.
#' - Per-feature assay layer selection is supported via the `layer` argument,
#'   closely mirroring the behaviour of `FeatureDensityPlot()`.
#' - Optional custom per-feature titles are supported via `plot.title`.
#' - Grouping is controlled via `FetchData()` (using the special keyword
#'   "ident" for active identities), avoiding any temporary changes to
#'   `Seurat::Idents()` and keeping identity classes stable.
#'
#' Titles for individual violins are constructed from the feature name and the
#' per-feature layer choice. Metadata-backed features (columns in
#' `SeuratObject@meta.data`) never receive a layer suffix, while assay-backed
#' features appended with an underscore and the chosen layer (e.g. `CD3D_counts`
#' or `CD3D_data`). This behaviour matches the `FeatureDensityPlot()` semantics
#' while preserving the single combined plot output of `Seurat::VlnPlot()`.
#'
#' @param SeuratObject A Seurat object.
#' @param features Character vector of features to plot. Each entry can be a
#'   metadata column name, an assay feature name, or any variable retrievable via
#'   `Seurat::FetchData()` (including dimensional reduction variables such as
#'   "PC_1" or the special keyword "ident").
#' @param gradient Character scalar. The feature used to compute the per-group
#'   gradient value. Use "nCells" to color by the number of cells per group
#'   (identity). Any metadata column or feature name resolvable by
#'   `Seurat::FetchData()` is also accepted; its mean per group is used.
#' @param group.by Character scalar or `NULL`. Name of the grouping variable
#'   used for the x-axis identities. When `group.by = "ident"` (default), the
#'   active identity classes are pulled via `Seurat::FetchData()` using the
#'   special "ident" keyword. When `group.by` is a metadata column name, that
#'   column is used instead. If `group.by` is `NULL`, identities are set from
#'   the active identity via "ident".
#' @param scale.colors Character scalar. Viridis palette option used for the
#'   gradient scales. Accepted values include "magma" / "A", "inferno" / "B",
#'   "plasma" / "C", "viridis" / "D", "cividis" / "E", "rocket" / "F",
#'   "mako" / "G", "turbo" / "H". Default is "viridis".
#' @param lower.limit Numeric scalar. Lower limit of the gradient scale. Only
#'   applied when `upper.limit` is non-`NULL`. Default is `0`.
#' @param upper.limit Numeric scalar or `NULL`. Upper limit of the gradient
#'   scale. When `NULL` (default), limits are set automatically by ggplot.
#' @param pt.size Numeric scalar. Size of jittered points overlaid on the
#'   violins. When `0` (default), no points are drawn. Values greater than
#'   `0` enable `geom_jitter()` overlays.
#' @param ncol Integer scalar or `NULL`. Number of columns in the combined
#'   patchwork layout. When `NULL`, the number of columns equals
#'   `length(features)`, but if `length(features) > 5`, the number of columns #'   is set to `ceiling(sqrt(length(features)))` to avoid excessive horizontal
#'   stretching.
#' @param layer Character scalar or vector, `NA` or `NULL`. Assay layer(s) used
#'   to retrieve feature expression values. When `NA`or `NULL`, Seurat's default 
#'   layer resolution is used. When a length-1 character vector is supplied, the 
#'   same layer is applied to all features. When a character vector of length 
#'   equal to `length(features)` is supplied, each feature uses its corresponding
#'   layer entry. Entries that are "null" (case-insensitive), or empty
#'   strings are treated as `NULL` for that feature. Metadata-backed features
#'   ignore `layer` at the data access level, but assay-backed features obey it.
#' @param plot.title Optional custom title(s). `NULL` uses internally generated
#'   per-feature names (metadata features use `feature`; assay features use
#'   `feature_layer` when layer is provided). A length-1 string is recycled to
#'   all features. A character vector with length equal to `length(features)`
#'   applies one title per feature.
#'
#' @return A single `ggplot2` / `patchwork` object with one violin panel per
#'   requested feature. The panels share a common gradient legend, and the x-axis
#'   identities are ordered from highest to lowest gradient value.
#'
#' @details
#' **Data access**: All feature and grouping values are retrieved via
#' `Seurat::FetchData()`. Grouping uses the `group.by` variable, which defaults
#' to "ident" (the active identities). Gradient values are computed either as
#' cell counts per group (for `gradient = "nCells"`) or as mean feature values
#' per group (for any other gradient feature). Feature expression values for
#' the violins are fetched per feature using the optional `layer` argument.
#'
#' **Layer handling**: The `layer` parameter mirrors `FeatureDensityPlot()`.
#' A length-1 character value applies to all features; a vector matching
#' `length(features)` specifies per-feature layers. Metadata columns ignore
#' `layer` when fetching data. For assay-backed features, the default title of
#' each violin is constructed as `feature_layer` when a non-`NULL` layer is
#' used; metadata-backed features retain bare names without any suffix.
#'
#' **Title handling**: `plot.title = NULL` preserves the historical naming
#' behaviour (`feature_layer` for assay-backed features with layer, `feature`
#' otherwise). A length-1 `plot.title` is recycled to all features, and a
#' length-`length(features)` vector maps one custom title per feature.
#'
#' **Efficiency**: A top-level `lapply()` iterates over features to build the
#' per-feature violin panels. Data extraction is performed per feature via
#' `FetchData()`, and grouping / gradient values are computed once and reused
#' across features. There are no explicit cell-wise loops; all operations work
#' on vectors or small data frames.
#'
#' **Output**: Unlike `FeatureDensityPlot()`, `VlnPlotGradient()` always returns
#' a single combined `ggplot2` / `patchwork` object, matching the default
#' behaviour of `Seurat::VlnPlot()` when `combine = TRUE`.
#'
#' @examples
#'
#' \dontrun{
#' # QC metrics colored by number of cells per identity
#' VlnPlotGradient(
#'   SeuratObject,
#'   features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
#'   gradient = "nCells"
#' )
#'
#' # Mixed metadata and assay features, with per-feature layers
#' VlnPlotGradient(
#'   SeuratObject,
#'   features = c("percent.mt", "CD3D"),
#'   gradient = "nCount_RNA",
#'   group.by = "ident",
#'   layer = c(NA, "counts")
#' )
#'
#' # Duplicate gene with different layers, both plotted side-by-side
#' VlnPlotGradient(
#'   SeuratObject,
#'   features = c("CD3D", "CD3D"),
#'   gradient = "CD3D",
#'   group.by = "ident",
#'   layer = c("counts", "data")
#' )
#' }
#'
#' @import ggplot2
#' @import dplyr
#' @import patchwork
#' @import Seurat
#'
#' @export
VlnPlotGradient <- function(
  SeuratObject,
  features,
  gradient,
  group.by = "ident",
  scale.colors = "viridis",
  lower.limit = 0,
  upper.limit = NULL,
  pt.size = 0,
  ncol = NULL,
  layer = NULL,
  plot.title = NULL
) {
  # Validate the Seurat object early so downstream code can assume a consistent
  # structure and type.
  if (!inherits(SeuratObject, "Seurat")) {
    stop("'SeuratObject' must be a Seurat object.")
  }

  # Validate feature names: must be a non-empty character vector.
  if (!is.character(features) || length(features) == 0L) {
    stop("'features' must be a non-empty character vector.")
  }

  # Validate gradient: must be a single character string.
  if (!is.character(gradient) || length(gradient) != 1L || is.na(gradient)) {
    stop("'gradient' must be a single non-NA character value.")
  }

  # Validate group.by: either NULL or a single character string.
  if (!is.null(group.by)) {
    if (!is.character(group.by) || length(group.by) != 1L || is.na(group.by)) {
      stop("'group.by' must be NULL or a single non-NA character value.")
    }
  }

  # Normalize group.by: if NULL, default to "ident" so that the active identities
  # are used as grouping variable. This keeps behaviour explicit and avoids
  # temporary changes to Idents(SeuratObject).
  group.fetch.var <- if (is.null(group.by)) "ident" else group.by

  # Validate numeric scalars: lower.limit, upper.limit, pt.size.
  if (
    !is.numeric(lower.limit) || length(lower.limit) != 1L || is.na(lower.limit)
  ) {
    stop("'lower.limit' must be a single non-NA numeric value.")
  }
  if (!is.null(upper.limit)) {
    if (
      !is.numeric(upper.limit) ||
        length(upper.limit) != 1L ||
        is.na(upper.limit)
    ) {
      stop("'upper.limit' must be NULL or a single non-NA numeric value.")
    }
    if (upper.limit <= lower.limit) {
      stop(
        "'upper.limit' must be strictly greater than 'lower.limit' when both are specified."
      )
    }
  }
  if (
    !is.numeric(pt.size) ||
      length(pt.size) != 1L ||
      is.na(pt.size) ||
      pt.size < 0
  ) {
    stop("'pt.size' must be a single non-negative numeric value.")
  }

  # Validate ncol: NULL or a positive integer.
  if (!is.null(ncol)) {
    if (!is.numeric(ncol) || length(ncol) != 1L || is.na(ncol) || ncol <= 0) {
      stop("'ncol' must be NULL or a single positive numeric value.")
    }
  }

  # Validate plot.title so custom title behavior matches FeatureDensityPlot()
  # semantics while preserving backward-compatible defaults.
  if (!is.null(plot.title)) {
    # plot.title must be character when provided because titles are strings.
    if (!is.character(plot.title)) {
      stop("'plot.title' must be NULL or a character vector.")
    }

    # Allow a single title to be recycled to every panel, or one title per
    # requested feature position.
    if (!(length(plot.title) == 1L || length(plot.title) == length(features))) {
      stop("'plot.title' must have length 1 or length(features).")
    }
  }

  # Validate scale.colors: must be a known viridis option.
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

  # Helper to normalize the layer argument into a per-feature list of character
  # values or NULLs. This mirrors FeatureDensityPlot() semantics.
  normalizeLayerPerFeature <- function(layer.value, n.features) {
    # When no layer is supplied, return a list of NULLs (one per feature).
    if (is.null(layer.value)) {
      return(rep(list(NULL), n.features))
    }

    # Layers must be character when provided.
    if (!is.character(layer.value)) {
      stop("'layer' must be NULL or a character vector.")
    }

    # Length-1 layer is recycled to all features.
    if (length(layer.value) == 1L) {
      return(rep(list(layer.value), n.features))
    }

    # Per-feature layer specification must match the number of features.
    if (length(layer.value) != n.features) {
      stop(sprintf(
        "'layer' has length %d but 'features' has length %d; it must be length 1 or length(features).",
        length(layer.value),
        n.features
      ))
    }

    # Map each entry to either a character string or NULL (for NA / "null" / empty).
    return(lapply(layer.value, function(entry) {
      if (is.na(entry)) {
        return(NULL)
      }
      entry.lower <- tolower(entry)
      if (entry.lower %in% c("null", "")) {
        return(NULL)
      }
      return(entry)
    }))
  }

  # Prepare per-feature layer specification.
  layer.per.feature <- normalizeLayerPerFeature(layer, length(features))

  # Cache metadata column names once for efficient metadata detection when
  # constructing plot titles.
  metadata.names <- colnames(SeuratObject@meta.data)

  # Resolve grouping values via FetchData so that grouping always comes from the
  # requested variable (active identity via "ident" or a metadata column).
  cells.use <- colnames(SeuratObject)
  grouping.data <- tryCatch(
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

  # Extract grouping values as character; this vector is used both for gradient
  # computation and for x-axis identities.
  if (!group.fetch.var %in% colnames(grouping.data)) {
    stop(sprintf(
      "Grouping variable '%s' was not found in the fetched data.",
      group.fetch.var
    ))
  }
  group.values <- as.character(grouping.data[cells.use, group.fetch.var])

  # Gradient computation -------------------------------------------------------
  # Special case: gradient = "nCells" uses simple per-group cell counts.
  if (gradient == "nCells") {
    # Build a data frame with one row per cell, storing its group.
    gradient.ident.data <- data.frame(
      group = group.values,
      stringsAsFactors = FALSE
    )

    # Count cells per group using dplyr::count().
    gradient.values <- gradient.ident.data %>%
      dplyr::count(group, name = "nCells") %>%
      dplyr::rename(gradient_val = nCells)

    # Legend label for this special case.
    gradient.label <- "nCells"
  } else {
    # General case: compute mean gradient feature per group.
    gradient.data <- tryCatch(
      Seurat::FetchData(
        object = SeuratObject,
        vars = gradient,
        cells = cells.use,
        clean = FALSE
      ),
      error = function(e) {
        stop(sprintf(
          "Cannot fetch gradient feature '%s': %s",
          gradient,
          e$message
        ))
      }
    )

    # Validate that the gradient feature was successfully fetched.
    if (!gradient %in% colnames(gradient.data)) {
      stop(sprintf(
        "The gradient feature '%s' could not be fetched from the object. Check the feature name.",
        gradient
      ))
    }

    # Extract raw gradient values per cell and align them to cells.use.
    gradient.raw <- gradient.data[cells.use, gradient]

    # Coerce to numeric explicitly; any non-numeric entries become NA.
    if (!is.numeric(gradient.raw)) {
      suppressWarnings(
        gradient.raw <- as.numeric(gradient.raw)
      )
    }

    # Assemble per-cell data with group and raw gradient values.
    gradient.ident.data <- data.frame(
      group = group.values,
      gradient_raw = gradient.raw,
      stringsAsFactors = FALSE
    )

    # Summarise mean gradient per group, removing NA values.
    gradient.values <- gradient.ident.data %>%
      dplyr::group_by(group) %>%
      dplyr::summarise(
        gradient_val = mean(gradient_raw, na.rm = TRUE),
        .groups = "drop"
      )

    # Legend label for the general mean case.
    gradient.label <- paste0("mean(", gradient, ")")
  }

  # Order groups by descending gradient value so violins appear from highest to
  # lowest gradient on the x-axis.
  gradient.values <- gradient.values %>%
    dplyr::arrange(dplyr::desc(gradient_val))
  ordered.levels <- gradient.values$group

  # Feature data extraction ----------------------------------------------------
  # For each feature, fetch its values for all cells, optionally using a layer.
  feature.data <- lapply(seq_along(features), function(feature.id) {
    # Current feature name.
    feature.name <- features[[feature.id]]
    # Per-feature layer (NULL or character).
    feature.layer <- layer.per.feature[[feature.id]]

    # Fetch feature values via Seurat::FetchData().
    fetched <- tryCatch(
      Seurat::FetchData(
        object = SeuratObject,
        vars = feature.name,
        cells = cells.use,
        layer = feature.layer,
        clean = FALSE
      ),
      error = function(e) {
        stop(sprintf(
          "Cannot fetch feature '%s' (layer '%s'): %s",
          feature.name,
          if (is.null(feature.layer)) "NULL" else feature.layer,
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

    # Coerce to numeric explicitly; NA-producing coercions are tolerated and
    # handled later by ggplot.
    if (!is.numeric(values)) {
      suppressWarnings(
        values <- as.numeric(values)
      )
    }

    return(values)
  })

  # Construct per-feature default names with metadata-aware suffix stripping.
  # These names preserve the historical behavior and are used whenever
  # 'plot.title' is not provided by the caller.
  plot.names <- vapply(
    seq_along(features),
    function(feature.id) {
      feature.name <- features[[feature.id]]
      feature.layer <- layer.per.feature[[feature.id]]

      # Determine whether the feature is backed by metadata.
      is.metadata <- feature.name %in% metadata.names

      # Metadata-backed features ignore layer in titles.
      if (is.metadata) {
        return(feature.name)
      }

      # Assay-backed features: append layer when non-NULL; otherwise use bare name.
      if (is.null(feature.layer)) {
        return(feature.name)
      }

      return(paste0(feature.name, "_", feature.layer))
    },
    FUN.VALUE = character(1L)
  )

  # Resolve final panel titles.
  # - If plot.title is NULL, keep legacy behavior by using computed plot.names.
  # - If plot.title has length 1, recycle it to all feature panels.
  # - If plot.title has length length(features), use per-feature custom titles.
  plot.titles <- if (is.null(plot.title)) {
    plot.names
  } else if (length(plot.title) == 1L) {
    rep(plot.title, length(features))
  } else {
    plot.title
  }

  # Global gradient scale limits shared across all panels.
  scale.limits <- if (!is.null(upper.limit)) {
    c(lower.limit, upper.limit)
  } else {
    NULL
  }

  # Build one ggplot panel per feature using lapply() for efficiency.
  plot.list <- lapply(seq_along(features), function(feature.id) {
    feature.name <- features[[feature.id]]
    feature.values <- feature.data[[feature.id]]

    # Assemble per-cell data frame: group, value, and gradient_val.
    cell.data <- data.frame(
      group = group.values,
      value = feature.values,
      stringsAsFactors = FALSE
    ) %>%
      dplyr::left_join(gradient.values, by = "group") %>%
      dplyr::mutate(
        group = factor(group, levels = ordered.levels)
      )

    # Construct the ggplot object: violin plus optional jitter, viridis
    # gradients, Seurat-like theme.
    p <- ggplot2::ggplot(
      cell.data,
      ggplot2::aes(x = group, y = value, fill = gradient_val)
    ) +
      ggplot2::geom_violin(
        scale = "width",
        trim = TRUE,
        adjust = 1,
        linewidth = 0.3,
        color = "black"
      ) +
      {
        if (pt.size > 0) {
          ggplot2::geom_jitter(
            ggplot2::aes(color = gradient_val),
            width = 0.3,
            size = pt.size,
            alpha = 0.6,
            show.legend = FALSE
          )
        }
      } +
      ggplot2::scale_fill_viridis_c(
        name = gradient.label,
        option = scale.colors,
        limits = scale.limits,
        na.value = "grey70"
      ) +
      ggplot2::scale_color_viridis_c(
        option = scale.colors,
        limits = scale.limits,
        guide = "none"
      ) +
      ggplot2::labs(
        title = plot.titles[[feature.id]],
        x = NULL,
        y = NULL
      ) +
      ggplot2::theme_classic(base_size = 11) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          face = "bold",
          size = 11,
          hjust = 0.5
        ),
        axis.text.x = ggplot2::element_text(
          angle = 45,
          hjust = 1,
          vjust = 1,
          size = 9
        ),
        axis.text.y = ggplot2::element_text(size = 9),
        axis.title.y = ggplot2::element_text(size = 9),
        legend.title = ggplot2::element_text(size = 9),
        legend.text = ggplot2::element_text(size = 8),
        legend.key.height = ggplot2::unit(0.8, "cm"),
        legend.key.width = ggplot2::unit(0.25, "cm"),
        panel.grid.major.y = ggplot2::element_line(
          color = "grey90",
          linewidth = 0.25
        ),
        panel.border = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(linewidth = 0.4)
      )

    return(p)
  })

  # Combine panels using patchwork, sharing a single gradient legend.
  ncol <- if (!is.null(ncol)) {
    ncol
  } else if (length(features) > 5) {
    ceiling(sqrt(length(features)))
  } else {
    ncol <- length(features)
  }

  combined <- patchwork::wrap_plots(plot.list, ncol = ncol) +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "right")

  return(combined)
}
