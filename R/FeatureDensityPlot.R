#' @title FeatureDensityPlot
#' @description
#' `FeatureDensityPlot()` draws density plots for one or more metadata features from a
#' Seurat object without splitting the object. Grouping is resolved directly from
#' `SeuratObject@meta.data` (or from `active.ident`), which avoids expensive object
#' partitioning for large datasets.
#'
#' The function supports:
#' - overlayed grouped densities,
#' - split (faceted-by-group) density panels,
#' - optional vertical reference lines (`vline`) in dashed red,
#' - optional independent median overlays (`plot.median`) in dashed black,
#' - custom plot titles via `plot.title`,
#' - multi-feature output as a named list (one plot per feature).
#'
#' The implementation is sequential and uses base `lapply()`.
#'
#' @param SeuratObject A Seurat object.
#' @param features Character vector of metadata columns to plot on the x-axis.
#' @param group.by Character scalar. Grouping variable. Use a metadata column name,
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
#'   `"both"`, or a numeric value. `"upper"`, `"lower"`, and `"both"` use
#'   median +/- `nmad`*MAD.
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
#' This function intentionally uses sequential base `lapply()`.
#'
#' Efficiency design:
#' - A top-level `lapply()` iterates over `features`.
#' - A second nested `lapply()` is used only when `split.plot = TRUE` to iterate over
#'   groups within each feature.
#'
#' This conditional nested structure avoids building a global feature-by-group task
#' table and avoids costly post-hoc reassembly. In sequential execution, it keeps
#' per-feature data local and reduces indexing overhead while preserving readability.
#'
#' @examples
#' \dontrun{
#' plot.list <- FeatureDensityPlot(
#'   SeuratObject,
#'   features = c("percent.mt", "nCount_RNA"),
#'   group.by = "batch",
#'   split.plot = TRUE,
#'   vline = "upper",
#'   plot.median = TRUE
#' )
#'
#' FeatureDensityPlot(
#'   SeuratObject,
#'   features = "percent.mt",
#'   plot.title = "Mitochondrial Percentage Density"
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
  plot.median = TRUE,
  plot.title = NULL,
  nmad = 2,
  alpha = 0.3,
  pt.size = 0
) {
  # ─────────────────────────────────────────────────────────────────────────────
  # 1) Input validation
  #
  # We keep strict validation because plotting functions often fail later with less
  # informative messages if inputs are malformed. Failing early makes debugging much
  # faster for users.
  # ─────────────────────────────────────────────────────────────────────────────
  if (!inherits(SeuratObject, "Seurat")) {
    stop("'SeuratObject' must be a Seurat object.")
  }

  if (!is.character(features) || length(features) == 0) {
    stop("'features' must be a non-empty character vector.")
  }

  if (!is.logical(split.plot) || length(split.plot) != 1 || is.na(split.plot)) {
    stop("'split.plot' must be a single logical value.")
  }

  if (!is.logical(plot.median) || length(plot.median) != 1 || is.na(plot.median)) {
    stop("'plot.median' must be a single logical value.")
  }

  if (!is.numeric(alpha) || length(alpha) != 1 || is.na(alpha) || alpha < 0 || alpha > 1) {
    stop("'alpha' must be a single numeric value between 0 and 1.")
  }

  if (!is.numeric(pt.size) || length(pt.size) != 1 || is.na(pt.size) || pt.size < 0) {
    stop("'pt.size' must be a single non-negative numeric value.")
  }

  if (!is.numeric(nmad) || length(nmad) != 1 || is.na(nmad) || nmad < 0) {
    stop("'nmad' must be a single non-negative numeric value.")
  }

  if (!is.null(plot.title)) {
    if (!is.character(plot.title)) {
      stop("'plot.title' must be NULL or a character vector.")
    }
    if (!(length(plot.title) == 1 || length(plot.title) == length(features))) {
      stop("'plot.title' must have length 1 or length(features).")
    }
  }

  valid.vline.values <- c("mean", "median", "upper", "lower", "both")
  if (!is.null(vline)) {
    if (is.character(vline)) {
      if (length(vline) != 1 || !tolower(vline) %in% valid.vline.values) {
        stop("'vline' must be one of: NULL, 'mean', 'median', 'upper', 'lower', 'both', or a numeric value.")
      }
      vline <- tolower(vline)
    } else if (!is.numeric(vline) || length(vline) != 1 || is.na(vline)) {
      stop("'vline' must be one of: NULL, 'mean', 'median', 'upper', 'lower', 'both', or a numeric value.")
    }
  }

  # ─────────────────────────────────────────────────────────────────────────────
  # 2) Data extraction and grouping resolution
  #
  # We only touch metadata and never split the Seurat object itself. This preserves
  # memory and keeps runtime predictable for large datasets.
  # ─────────────────────────────────────────────────────────────────────────────
  metadata.df <- SeuratObject@meta.data

  missing.features <- setdiff(features, colnames(metadata.df))
  if (length(missing.features) > 0) {
    stop(
      paste0(
        "These features are not metadata columns: ",
        paste(missing.features, collapse = ", ")
      )
    )
  }

  has.grouping <- !is.null(group.by)
  if (has.grouping) {
    if (!is.character(group.by) || length(group.by) != 1) {
      stop("'group.by' must be NULL or a single character value.")
    }

    if (group.by == "active.ident") {
      group.values <- as.character(Seurat::Idents(SeuratObject))
      group.label <- "active.ident"
    } else {
      if (!group.by %in% colnames(metadata.df)) {
        stop(paste0("'group.by' column '", group.by, "' not found in SeuratObject@meta.data."))
      }
      group.values <- as.character(metadata.df[[group.by]])
      group.label <- group.by
    }
  } else {
    group.values <- rep("all", nrow(metadata.df))
    group.label <- "all"
  }

  plot.data <- metadata.df[, features, drop = FALSE]
  plot.data$.group <- group.values
  plot.data <- plot.data[!is.na(plot.data$.group), , drop = FALSE]

  # Density requires numeric x values. Coercion is explicit and NA-producing coercions
  # are handled downstream by removing missing values per feature.
  for (feature.name in features) {
    if (!is.numeric(plot.data[[feature.name]])) {
      suppressWarnings(plot.data[[feature.name]] <- as.numeric(plot.data[[feature.name]]))
    }
  }

  group.levels <- unique(plot.data$.group)
  if (length(group.levels) == 0) {
    stop("No groups available to plot after filtering missing grouping values.")
  }

  feature.titles <- if (is.null(plot.title)) {
    features
  } else if (length(plot.title) == 1) {
    rep(plot.title, length(features))
  } else {
    plot.title
  }

  # ─────────────────────────────────────────────────────────────────────────────
  # 3) Helper functions
  #
  # The helpers isolate statistical line logic so the main plotting branches stay
  # readable. This separation also makes future maintenance safer.
  # ─────────────────────────────────────────────────────────────────────────────
  compute.vline.positions <- function(values, vline.spec, nmad.value) {
    values <- values[!is.na(values)]
    if (length(values) == 0 || is.null(vline.spec)) {
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
      return(data.frame(xintercept = median.value + nmad.value * mad.value, stringsAsFactors = FALSE))
    }

    if (vline.spec == "lower") {
      return(data.frame(xintercept = median.value - nmad.value * mad.value, stringsAsFactors = FALSE))
    }

    if (vline.spec == "both") {
      return(data.frame(
        xintercept = c(median.value - nmad.value * mad.value, median.value + nmad.value * mad.value),
        stringsAsFactors = FALSE
      ))
    }

    return(data.frame(xintercept = numeric(0), stringsAsFactors = FALSE))
  }

  compute.median.positions <- function(values) {
    values <- values[!is.na(values)]
    if (length(values) == 0) {
      return(data.frame(xintercept = numeric(0), stringsAsFactors = FALSE))
    }
    return(data.frame(xintercept = stats::median(values), stringsAsFactors = FALSE))
  }

  build.single.density.panel <- function(feature.df, x.label, panel.title) {
    current.plot <- ggplot2::ggplot(feature.df, ggplot2::aes(x = value)) +
      ggplot2::geom_density(fill = "lightblue", alpha = alpha, color = "black")

    if (pt.size > 0) {
      current.plot <- current.plot +
        ggplot2::geom_rug(sides = "b", linewidth = pt.size, alpha = 0.35)
    }

    # Color/linetype contract is intentionally fixed for interpretability:
    # red dashed for vline references, black dashed for medians.
    vline.df <- compute.vline.positions(feature.df$value, vline, nmad)
    if (nrow(vline.df) > 0) {
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
      if (nrow(median.df) > 0) {
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

    return(
      current.plot +
        ggplot2::theme_bw() +
        ggplot2::labs(title = panel.title, x = x.label, y = "Density") +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    )
  }

  # ─────────────────────────────────────────────────────────────────────────────
  # 4) Feature-level plotting via top-level lapply
  #
  # Efficiency choice:
  # - Always iterate features with lapply.
  # - Only in split mode, use a second nested lapply over groups.
  #
  # This avoids global task-index construction and keeps per-feature data localized.
  # ─────────────────────────────────────────────────────────────────────────────
  feature.plots <- lapply(seq_along(features), function(feature.id) {
    feature.name <- features[[feature.id]]
    feature.title <- feature.titles[[feature.id]]

    feature.df <- plot.data[, c(feature.name, ".group"), drop = FALSE]
    names(feature.df)[1] <- "value"
    feature.df <- feature.df[!is.na(feature.df$value), , drop = FALSE]

    if (nrow(feature.df) == 0) {
      stop(paste0("Feature '", feature.name, "' has no non-missing numeric values to plot."))
    }

    if (!has.grouping) {
      return(build.single.density.panel(feature.df, x.label = feature.name, panel.title = feature.title))
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

      if (!is.null(vline)) {
        vline.df <- feature.df |>
          dplyr::group_by(.group) |>
          dplyr::group_modify(~compute.vline.positions(.x$value, vline, nmad)) |>
          dplyr::ungroup()

        if (nrow(vline.df) > 0) {
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
          dplyr::summarise(xintercept = stats::median(value, na.rm = TRUE), .groups = "drop")

        if (nrow(median.df) > 0) {
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
          ggplot2::scale_color_viridis_d(option = scale.colors, name = group.label) +
          ggplot2::scale_fill_viridis_d(option = scale.colors, name = group.label) +
          ggplot2::theme_bw() +
          ggplot2::labs(title = feature.title, x = feature.name, y = "Density") +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
      )
    }

    # Split branch: nested lapply across groups (only when needed).
    group.plots <- lapply(group.levels, function(group.name) {
      group.df <- feature.df[feature.df$.group == group.name, , drop = FALSE]
      return(
        build.single.density.panel(
          feature.df = group.df,
          x.label = feature.name,
          panel.title = group.name
        )
      )
    })

    ncol.groups <- if (!is.null(ncol)) ncol else ceiling(sqrt(length(group.plots)))

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
  #
  # Preserve existing behavior:
  # - one feature -> single plot object
  # - multiple features -> named list of plots
  # ─────────────────────────────────────────────────────────────────────────────
  names(feature.plots) <- features

  if (length(feature.plots) == 1) {
    return(feature.plots[[1]])
  }

  return(feature.plots)
}
