#' @title FeatureDensityPlot
#' @description
#' `FeatureDensityPlot()` draws density plots for features from a Seurat object
#' without splitting the object. It supports metadata columns, assay features
#' (from any layer), and dimensional reduction variables via `Seurat::FetchData()`.
#'
#' Grouping is resolved from the specified `group.by` variable, avoiding expensive
#' object partitioning for large datasets.
#'
#' The function supports:
#' - overlayed grouped densities,
#' - split (faceted-by-group) density panels,
#' - optional vertical reference lines (`vline`) in dashed red,
#' - optional independent median overlays (`plot.median`) in dashed black,
#' - custom plot titles via `plot.title`,
#' - multi-feature output as a named list (one plot per feature),
#' - per-feature layer selection via `layer`.
#'
#' The implementation uses sequential base `lapply()` for efficiency.
#'
#' @param SeuratObject A Seurat object.
#' @param features Character vector of features to plot (metadata columns, assay
#'   feature names, or reduction variable names like `"PC_1"`).
#' @param group.by Character scalar. Grouping variable: a metadata column name,
#'   `"active.ident"` (default), or `NULL` for no grouping.
#' @param split.plot Logical. If `TRUE` (default), creates one panel per group level
#'   for each feature. If `FALSE`, groups are overlaid in one panel per feature.
#' @param scale.colors Character scalar. Viridis palette option used for grouped
#'   density colors (`"viridis"`, `"magma"`, `"plasma"`, `"inferno"`,
#'   `"cividis"`, `"rocket"`, `"mako"`, `"turbo"`). Default is `"viridis"`.
#' @param ncol Integer or `NULL`. Number of columns for split panels within each
#'   feature plot when `split.plot = TRUE`. If `NULL`, inferred automatically.
#' @param vline Optional reference-line specification (drawn in dashed red).
#'   Accepted values: `NULL`, `"mean"`, `"median"`, `"upper"`, `"lower"`,
#'   `"both"`, a numeric value, or a character vector with one entry per feature.
#'   For vector input, each element follows the same rules as the scalar form.
#'   `"upper"`, `"lower"`, and `"both"` use median +/- `nmad`*MAD.
#'   A length-1 vector is recycled to all features.
#' @param layer Character or `NULL`. Specifies which assay layer to extract feature
#'   values from (e.g., `"data"`, `"counts"`, `"scale.data"`). `NULL` uses the
#'   default layer. May be a single string (applied to all features) or a character
#'   vector of length equal to `features` (per-feature layer specification).
#' @param plot.median Logical. If `TRUE` (default), draws median line(s) in dashed
#'   black, independently of `vline`.
#' @param plot.title Optional custom title(s). `NULL` uses feature names as titles.
#'   A length-1 string applies to all features. A character vector with length equal
#'   to `length(features)` applies one title per feature.
#' @param nmad Numeric. Number of MADs used when `vline` is `"upper"`, `"lower"`,
#'   or `"both"`. Default is `2`.
#' @param alpha Numeric in `[0, 1]`. Fill alpha for density geometries.
#' @param pt.size Numeric. If `0` (default), no rug is drawn. If `> 0`, adds a rug
#'   (`geom_rug()`) with this line width.
#'
#' @return If `length(features) == 1`, returns a `ggplot2`/`patchwork` plot object.
#'   If `length(features) > 1`, returns a named list of plot objects (one per feature;
#'   list names equal `features`).
#'
#' @details
#' **Data access**: Features are resolved via `Seurat::FetchData()`, which tries
#' metadata columns first, then keyed variables, then assay features, then the
#' special `ident` keyword.
#'
#' **Layer handling**: When `layer` is a vector of different values, each feature
#' is fetched with its assigned layer. Metadata columns ignore the layer parameter;
#' assay features use it to select the desired matrix.
#'
#' **Efficiency**: A top-level `lapply()` iterates over features. A second nested
#' `lapply()` is used only when `split.plot = TRUE` to iterate over group levels.
#' This structure keeps per-feature data local and reduces indexing overhead.
#'
#' @examples
#' \dontrun{
#' # Metadata-only (backward compatible)
#' plot.list <- FeatureDensityPlot(
#'   SeuratObject,
#'   features = c("percent.mt", "nCount_RNA"),
#'   group.by = "batch",
#'   split.plot = TRUE,
#'   vline = "upper",
#'   plot.median = TRUE
#' )
#'
#' # Assay features from log-normalized data
#' FeatureDensityPlot(
#'   SeuratObject,
#'   features = "CD3D",
#'   layer = "data"
#' )
#'
#' # Mixed per-feature layers
#' FeatureDensityPlot(
#'   SeuratObject,
#'   features = c("percent.mt", "CD3D", "nCount_RNA"),
#'   layer = c(NA, "counts", NA)  # CD3D from raw counts, others from metadata
#' )
#' }
#'
#' @import ggplot2
#' @import dplyr
#' @import patchwork
#' @import Seurat
#' @export
FeatureDensityPlot <- function(
  SeuratObject,
  features,
  group.by = "active.ident",
  split.plot = TRUE,
  scale.colors = "viridis",
  ncol = NULL,
  vline = NULL,
  layer = NULL,
  plot.median = TRUE,
  plot.title = NULL,
  nmad = 2,
  alpha = 0.3,
  pt.size = 0
) {
  # ─────────────────────────────────────────────────────────────────────────────
  # 1) Input validation
  #
  # Validate early so the plotting code can stay compact and deterministic.
  # ─────────────────────────────────────────────────────────────────────────────
  if (!inherits(SeuratObject, "Seurat")) {
    stop("'SeuratObject' must be a Seurat object.")
  }

  if (!is.character(features) || length(features) == 0L) {
    stop("'features' must be a non-empty character vector.")
  }

  if (
    !is.logical(split.plot) || length(split.plot) != 1L || is.na(split.plot)
  ) {
    stop("'split.plot' must be a single logical value.")
  }

  if (
    !is.logical(plot.median) || length(plot.median) != 1L || is.na(plot.median)
  ) {
    stop("'plot.median' must be a single logical value.")
  }

  if (
    !is.numeric(alpha) ||
      length(alpha) != 1L ||
      is.na(alpha) ||
      alpha < 0 ||
      alpha > 1
  ) {
    stop("'alpha' must be a single numeric value between 0 and 1.")
  }

  if (
    !is.numeric(pt.size) ||
      length(pt.size) != 1L ||
      is.na(pt.size) ||
      pt.size < 0
  ) {
    stop("'pt.size' must be a single non-negative numeric value.")
  }

  if (!is.numeric(nmad) || length(nmad) != 1L || is.na(nmad) || nmad < 0) {
    stop("'nmad' must be a single non-negative numeric value.")
  }

  if (!is.null(plot.title)) {
    if (!is.character(plot.title)) {
      stop("'plot.title' must be NULL or a character vector.")
    }
    if (!(length(plot.title) == 1L || length(plot.title) == length(features))) {
      stop("'plot.title' must have length 1 or length(features).")
    }
  }

  normalize.vline.entry <- function(
    x,
    valid.keywords = c("mean", "median", "upper", "lower", "both")
  ) {
    if (is.null(x) || is.na(x) || (is.character(x) && tolower(x) == "null")) {
      return(NULL)
    }
    if (is.character(x)) {
      x.lower <- tolower(x)
      if (x.lower %in% valid.keywords) {
        return(x.lower)
      }
      x.numeric <- suppressWarnings(as.numeric(x))
      if (!is.na(x.numeric)) {
        return(x.numeric)
      }
      stop(sprintf(
        "vline entry '%s' is not a valid keyword, 'NULL', or a numeric value.",
        x
      ))
    }
    if (is.numeric(x)) {
      if (length(x) == 1L && !is.na(x)) {
        return(x)
      }
      stop("numeric vline must be a single value, not a vector.")
    }
    stop(sprintf(
      "vline entry has unexpected type: %s",
      paste(class(x), collapse = ", ")
    ))
  }

  valid.vline.values <- c("mean", "median", "upper", "lower", "both")
  if (is.null(vline)) {
    vline.per.feature <- rep(list(NULL), length(features))
  } else if (is.numeric(vline)) {
    if (length(vline) != 1L || is.na(vline)) {
      stop("'vline' numeric input must be a single non-NA value.")
    }
    vline.per.feature <- rep(list(vline), length(features))
  } else if (is.character(vline)) {
    if (length(vline) == 1L) {
      vline.per.feature <- rep(
        list(normalize.vline.entry(vline, valid.vline.values)),
        length(features)
      )
    } else if (length(vline) == length(features)) {
      vline.per.feature <- lapply(
        vline,
        normalize.vline.entry,
        valid.keywords = valid.vline.values
      )
    } else {
      stop(sprintf(
        "'vline' has length %d but 'features' has length %d; it must be length 1 or length(features).",
        length(vline),
        length(features)
      ))
    }
  } else {
    stop("'vline' must be NULL, numeric, or character.")
  }

  if (is.null(layer)) {
    layer.per.feature <- rep(list(NULL), length(features))
  } else if (is.character(layer)) {
    if (length(layer) == 1L) {
      layer.per.feature <- rep(list(layer), length(features))
    } else if (length(layer) == length(features)) {
      layer.per.feature <- lapply(layer, function(x) {
        if (is.na(x) || tolower(x) == "null") {
          return(NULL)
        }
        x
      })
    } else {
      stop(sprintf(
        "'layer' has length %d but 'features' has length %d; it must be length 1 or length(features).",
        length(layer),
        length(features)
      ))
    }
  } else {
    stop("'layer' must be NULL or a character vector.")
  }

  # ─────────────────────────────────────────────────────────────────────────────
  # 2) Data extraction via Seurat::FetchData
  #
  # We construct per-feature internal keys so duplicated feature names remain
  # distinct internally. This is essential for cases like features = c("MALAT1",
  # "MALAT1") with different layer entries; name-based indexing would collapse
  # those entries and silently reuse one layer.
  # ─────────────────────────────────────────────────────────────────────────────
  cells.use <- colnames(SeuratObject)
  has.grouping <- !is.null(group.by)

  # Internal stable keys (one key per requested feature occurrence).
  # These keys are used only for internal storage and lookup.
  feature.keys <- paste0(".feature_", seq_along(features))

  # Resolve grouping values once so all feature fetches align to the same cells.
  if (has.grouping) {
    if (!is.character(group.by) || length(group.by) != 1L) {
      stop("'group.by' must be NULL or a single character value.")
    }
    group.fetch.var <- if (group.by == "active.ident") "ident" else group.by
    group.data <- Seurat::FetchData(
      object = SeuratObject,
      vars = group.fetch.var,
      cells = cells.use,
      clean = FALSE
    )
    group.values <- as.character(group.data[[group.fetch.var]])
    group.label <- group.by
  } else {
    group.fetch.var <- NULL
    group.values <- rep("all", length(cells.use))
    group.label <- "all"
  }

  # Fetch each feature by position (not by name). Position-based alignment is the
  # contract that links features, layer.per.feature, vline.per.feature, and titles.
  feature.data <- lapply(seq_along(features), function(feature.id) {
    feature.name <- features[[feature.id]]
    feature.layer <- layer.per.feature[[feature.id]]

    fetched <- Seurat::FetchData(
      object = SeuratObject,
      vars = feature.name,
      cells = cells.use,
      layer = feature.layer,
      clean = FALSE
    )

    if (!feature.name %in% colnames(fetched)) {
      stop(sprintf(
        "The feature '%s' could not be fetched from the object. Check the feature name and layer.",
        feature.name
      ))
    }

    # Explicitly reindex by cells.use to preserve row alignment even if any future
    # Seurat internals alter FetchData row ordering.
    fetched[cells.use, feature.name, drop = TRUE]
  })
  names(feature.data) <- feature.keys

  # Assemble internal plotting data with unique keys so duplicate feature names do
  # not overwrite or alias one another.
  plot.data <- data.frame(row.names = cells.use)
  for (feature.id in seq_along(features)) {
    feature.key <- feature.keys[[feature.id]]
    plot.data[[feature.key]] <- feature.data[[feature.key]]
  }

  # Append grouping column and remove rows without grouping values.
  plot.data$.group <- group.values
  plot.data <- plot.data[!is.na(plot.data$.group), , drop = FALSE]

  # Density requires numeric x values. Coercion is explicit and NA-producing
  # coercions are handled downstream by removing missing values per feature.
  for (feature.id in seq_along(features)) {
    feature.key <- feature.keys[[feature.id]]
    if (!is.numeric(plot.data[[feature.key]])) {
      suppressWarnings(
        plot.data[[feature.key]] <- as.numeric(plot.data[[feature.key]])
      )
    }
  }

  # We will the define the plot nlames based on the features and the layer used for each feature. This is important for cases where the same feature name is used with different layers, as it ensures that each plot has a unique name.
  # This variable will be also used to name the plots in the returned list, ensuring that each plot can be easily identified based on the feature and layer used.
  plot.names <- paste0(features, "_", layer.per.feature)
  plot.names <- gsub("_+$", "", plot.names) # Remove trailing underscores if layer is NULL
  plot.names <- gsub("_NULL$", "", plot.names) # Remove trailing _NULL if layer is "NULL" (string)
  plot.names <- gsub("_NA$", "", plot.names) # Remove trailing _NA if layer is NA

  # Titles stay user-facing and index-aligned to features.
  feature.titles <- if (is.null(plot.title)) {
    plot.names
  } else if (length(plot.title) == 1L) {
    rep(plot.title, length(features))
  } else {
    plot.title
  }

  # Cache group levels once for split-plot branch.
  group.levels <- unique(plot.data$.group)
  if (length(group.levels) == 0L) {
    stop("No groups available to plot after filtering missing grouping values.")
  }

  # ─────────────────────────────────────────────────────────────────────────────
  # 3) Helper functions
  # ─────────────────────────────────────────────────────────────────────────────
  compute.vline.positions <- function(values, vline.spec, nmad.value) {
    values <- values[!is.na(values)]
    if (length(values) == 0L || is.null(vline.spec)) {
      return(data.frame(xintercept = numeric(0), stringsAsFactors = FALSE))
    }

    if (is.numeric(vline.spec)) {
      return(data.frame(xintercept = vline.spec, stringsAsFactors = FALSE))
    }

    median.value <- stats::median(values)

    if (vline.spec == "mean") {
      return(data.frame(xintercept = mean(values), stringsAsFactors = FALSE))
    }

    if (vline.spec == "median") {
      return(data.frame(xintercept = median.value, stringsAsFactors = FALSE))
    }

    mad.value <- stats::mad(values, na.rm = TRUE)
    if (is.na(mad.value)) {
      return(data.frame(xintercept = numeric(0), stringsAsFactors = FALSE))
    }

    if (vline.spec == "upper") {
      return(data.frame(
        xintercept = median.value + nmad.value * mad.value,
        stringsAsFactors = FALSE
      ))
    }

    if (vline.spec == "lower") {
      return(data.frame(
        xintercept = median.value - nmad.value * mad.value,
        stringsAsFactors = FALSE
      ))
    }

    if (vline.spec == "both") {
      return(data.frame(
        xintercept = c(
          median.value - nmad.value * mad.value,
          median.value + nmad.value * mad.value
        ),
        stringsAsFactors = FALSE
      ))
    }

    return(data.frame(xintercept = numeric(0), stringsAsFactors = FALSE))
  }

  compute.median.positions <- function(values) {
    values <- values[!is.na(values)]
    if (length(values) == 0L) {
      return(data.frame(xintercept = numeric(0), stringsAsFactors = FALSE))
    }
    return(data.frame(
      xintercept = stats::median(values),
      stringsAsFactors = FALSE
    ))
  }

  build.single.density.panel <- function(
    feature.df,
    x.label,
    panel.title,
    vline.spec
  ) {
    current.plot <- ggplot2::ggplot(feature.df, ggplot2::aes(x = value)) +
      ggplot2::geom_density(fill = "lightblue", alpha = alpha, color = "black")

    if (pt.size > 0) {
      current.plot <- current.plot +
        ggplot2::geom_rug(sides = "b", linewidth = pt.size, alpha = 0.35)
    }

    vline.df <- compute.vline.positions(feature.df$value, vline.spec, nmad)
    if (nrow(vline.df) > 0L) {
      current.plot <- current.plot +
        ggplot2::geom_vline(
          data = vline.df,
          ggplot2::aes(xintercept = xintercept),
          color = "red",
          linetype = "dashed",
          linewidth = 0.6,
          show.legend = FALSE
        )
    }

    if (plot.median) {
      median.df <- compute.median.positions(feature.df$value)
      if (nrow(median.df) > 0L) {
        current.plot <- current.plot +
          ggplot2::geom_vline(
            data = median.df,
            ggplot2::aes(xintercept = xintercept),
            color = "black",
            linetype = "dashed",
            linewidth = 0.7,
            show.legend = FALSE
          )
      }
    }

    current.plot +
      ggplot2::theme_bw() +
      ggplot2::labs(title = panel.title, x = x.label, y = "Density") +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  }

  # ─────────────────────────────────────────────────────────────────────────────
  # 4) Feature-level plotting via top-level lapply
  # ─────────────────────────────────────────────────────────────────────────────
  feature.plots <- lapply(seq_along(features), function(feature.id) {
    # feature.name/feature.title are user-facing labels; feature.key is internal
    # and unique per position, which preserves correct behavior for duplicates.
    feature.name <- features[[feature.id]]
    feature.key <- feature.keys[[feature.id]]
    feature.title <- feature.titles[[feature.id]]
    feature.vline <- vline.per.feature[[feature.id]]

    feature.df <- plot.data[, c(feature.key, ".group"), drop = FALSE]
    names(feature.df)[1] <- "value"
    feature.df <- feature.df[!is.na(feature.df$value), , drop = FALSE]

    if (nrow(feature.df) == 0L) {
      stop(sprintf(
        "Feature '%s' has no non-missing numeric values to plot.",
        feature.name
      ))
    }

    if (!has.grouping) {
      return(build.single.density.panel(
        feature.df = feature.df,
        x.label = feature.name,
        panel.title = feature.title,
        vline.spec = feature.vline
      ))
    }

    if (!split.plot) {
      current.plot <- ggplot2::ggplot(
        feature.df,
        ggplot2::aes(x = value, color = .group, fill = .group)
      ) +
        ggplot2::geom_density(alpha = alpha)

      if (pt.size > 0) {
        current.plot <- current.plot +
          ggplot2::geom_rug(
            ggplot2::aes(color = .group),
            sides = "b",
            linewidth = pt.size,
            alpha = 0.35,
            show.legend = FALSE
          )
      }

      if (!is.null(feature.vline)) {
        vline.df <- feature.df |>
          dplyr::group_by(.group) |>
          dplyr::group_modify(
            ~ compute.vline.positions(.x$value, feature.vline, nmad)
          ) |>
          dplyr::ungroup()

        if (nrow(vline.df) > 0L) {
          current.plot <- current.plot +
            ggplot2::geom_vline(
              data = vline.df,
              ggplot2::aes(xintercept = xintercept, group = .group),
              color = "red",
              linetype = "dashed",
              linewidth = 0.55,
              show.legend = FALSE
            )
        }
      }

      if (plot.median) {
        median.df <- feature.df |>
          dplyr::group_by(.group) |>
          dplyr::summarise(
            xintercept = stats::median(value, na.rm = TRUE),
            .groups = "drop"
          )

        if (nrow(median.df) > 0L) {
          current.plot <- current.plot +
            ggplot2::geom_vline(
              data = median.df,
              ggplot2::aes(xintercept = xintercept, group = .group),
              color = "black",
              linetype = "dashed",
              linewidth = 0.65,
              show.legend = FALSE
            )
        }
      }

      return(
        current.plot +
          ggplot2::scale_color_viridis_d(
            option = scale.colors,
            name = group.label
          ) +
          ggplot2::scale_fill_viridis_d(
            option = scale.colors,
            name = group.label
          ) +
          ggplot2::theme_bw() +
          ggplot2::labs(
            title = feature.title,
            x = feature.name,
            y = "Density"
          ) +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
      )
    }

    group.plots <- lapply(group.levels, function(group.name) {
      group.df <- feature.df[feature.df$.group == group.name, , drop = FALSE]
      build.single.density.panel(
        feature.df = group.df,
        x.label = feature.name,
        panel.title = group.name,
        vline.spec = feature.vline
      )
    })

    ncol.groups <- if (!is.null(ncol)) {
      ncol
    } else {
      ceiling(sqrt(length(group.plots)))
    }

    return(
      patchwork::wrap_plots(group.plots, ncol = ncol.groups) +
        patchwork::plot_annotation(
          title = feature.title,
          theme = ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
          )
        )
    )
  })

  # ─────────────────────────────────────────────────────────────────────────────
  # 5) Return-shape contract
  # ─────────────────────────────────────────────────────────────────────────────

  names(feature.plots) <- plot.names

  if (length(feature.plots) == 1L) {
    return(feature.plots[[1L]])
  }

  return(feature.plots)
}
