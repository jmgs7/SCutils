#' @title FeatureDensityPlot
#' @description
#' `FeatureDensityPlot()` draws density plots for one or more metadata features from a
#' Seurat object without splitting the object. Grouping is handled directly from
#' `SeuratObject@meta.data` (or `active.ident`), which is substantially faster than
#' `SplitObject()` workflows for large datasets.
#'
#' The function supports:
#' - overlayed grouped densities,
#' - split (faceted-by-group) density panels,
#' - optional vertical reference lines (`vline`) in dashed red,
#' - optional independent median overlays (`plot.median`) in black,
#' - custom plot titles via `plot.title`,
#' - multi-feature output as a named list (one plot per feature),
#' - future-based parallel plotting with `future.apply::future_lapply()`.
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
#' @param plot.median Logical. If `TRUE` (default), draws median line(s) in black,
#'   independently of `vline`.
#' @param plot.title Optional custom title(s). `NULL` uses feature names as titles.
#'   A length-1 string applies to all features. A character vector with length equal
#'   to `length(features)` applies one title per feature.
#' @param nmad Numeric. Number of MADs used when `vline` is `"upper"`, `"lower"`,
#'   or `"both"`. Default is `2`.
#' @param alpha Numeric in `[0, 1]`. Fill alpha for density geometries.
#' @param pt.size Numeric. If `0` (default), no rug is drawn. If `> 0`, adds a rug
#'   (`geom_rug()`) with this line width.
#' @param mc.cores Integer or `NULL`. Optional worker override for a temporary local
#'   `future::plan(future::multisession, workers = mc.cores)` used during plotting.
#'   If `NULL`, the currently active future plan is respected.
#'
#' @return If `length(features) == 1`, returns a `ggplot2`/`patchwork` plot object.
#'   If `length(features) > 1`, returns a named list of plot objects (one per feature;
#'   list names equal `features`).
#'
#' @details
#' Parallelization uses `future.apply::future_lapply()`. For split mode
#' (`split.plot = TRUE`), work units are parallelized over feature/group pairs so
#' plotting can scale across both dimensions without nested parallel apply calls.
#'
#' @examples
#' \dontrun{
#' # Respect an externally configured plan.
#' future::plan(future::multisession, workers = 4)
#' plot_list <- FeatureDensityPlot(
#'   SeuratObject,
#'   features = c("percent.mt", "nCount_RNA"),
#'   group.by = "batch",
#'   split.plot = TRUE,
#'   vline = "upper",
#'   plot.median = TRUE
#' )
#' future::plan(future::sequential)
#'
#' # Or override workers only for this call.
#' FeatureDensityPlot(
#'   SeuratObject,
#'   features = "percent.mt",
#'   plot.title = "Mitochondrial Percentage Density",
#'   mc.cores = 2
#' )
#' }
#'
#' @import ggplot2
#' @import dplyr
#' @import patchwork
#' @import Seurat
#' @import future
#' @importFrom future.apply future_lapply
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
    pt.size = 0,
    mc.cores = NULL
) {
  # ─────────────────────────────────────────────────────────────────────────────
  # 1) Input validation
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

  if (!is.null(mc.cores)) {
    if (!is.numeric(mc.cores) || length(mc.cores) != 1 || is.na(mc.cores) || mc.cores < 1) {
      stop("'mc.cores' must be NULL or a positive integer.")
    }
    mc.cores <- as.integer(mc.cores)
  }

  if (!is.null(plot.title)) {
    if (!is.character(plot.title)) {
      stop("'plot.title' must be NULL or a character vector.")
    }
    if (!(length(plot.title) == 1 || length(plot.title) == length(features))) {
      stop("'plot.title' must have length 1 or length(features).")
    }
  }

  valid_vline_values <- c("mean", "median", "upper", "lower", "both")
  if (!is.null(vline)) {
    if (is.character(vline)) {
      if (length(vline) != 1 || !tolower(vline) %in% valid_vline_values) {
        stop("'vline' must be one of: NULL, 'mean', 'median', 'upper', 'lower', 'both', or a numeric value.")
      }
      vline <- tolower(vline)
    } else if (!is.numeric(vline) || length(vline) != 1 || is.na(vline)) {
      stop("'vline' must be one of: NULL, 'mean', 'median', 'upper', 'lower', 'both', or a numeric value.")
    }
  }

  # ─────────────────────────────────────────────────────────────────────────────
  # 2) Data extraction and grouping resolution
  # ─────────────────────────────────────────────────────────────────────────────
  metadata_df <- SeuratObject@meta.data

  missing_features <- setdiff(features, colnames(metadata_df))
  if (length(missing_features) > 0) {
    stop(
      paste0(
        "These features are not metadata columns: ",
        paste(missing_features, collapse = ", ")
      )
    )
  }

  has_grouping <- !is.null(group.by)
  if (has_grouping) {
    if (!is.character(group.by) || length(group.by) != 1) {
      stop("'group.by' must be NULL or a single character value.")
    }

    if (group.by == "active.ident") {
      group_values <- as.character(Seurat::Idents(SeuratObject))
      group_label <- "active.ident"
    } else {
      if (!group.by %in% colnames(metadata_df)) {
        stop(paste0("'group.by' column '", group.by, "' not found in SeuratObject@meta.data."))
      }
      group_values <- as.character(metadata_df[[group.by]])
      group_label <- group.by
    }
  } else {
    group_values <- rep("all", nrow(metadata_df))
    group_label <- "all"
  }

  plot_data <- metadata_df[, features, drop = FALSE]
  plot_data$.group <- group_values
  plot_data <- plot_data[!is.na(plot_data$.group), , drop = FALSE]

  for (feature_name in features) {
    if (!is.numeric(plot_data[[feature_name]])) {
      suppressWarnings(plot_data[[feature_name]] <- as.numeric(plot_data[[feature_name]]))
    }
  }

  group_levels <- unique(plot_data$.group)
  if (length(group_levels) == 0) {
    stop("No groups available to plot after filtering missing grouping values.")
  }

  feature_titles <- if (is.null(plot.title)) {
    features
  } else if (length(plot.title) == 1) {
    rep(plot.title, length(features))
  } else {
    plot.title
  }

  # ─────────────────────────────────────────────────────────────────────────────
  # 3) Future plan handling
  # If `mc.cores` is provided, apply a temporary local plan and restore on exit.
  # ─────────────────────────────────────────────────────────────────────────────
  if (!is.null(mc.cores)) {
    previous_plan <- future::plan()
    on.exit(future::plan(previous_plan), add = TRUE)
    future::plan(future::multisession, workers = mc.cores)
  }

  # ─────────────────────────────────────────────────────────────────────────────
  # 4) Helper functions for statistical lines and panel generation
  # ─────────────────────────────────────────────────────────────────────────────
  ComputeVlinePositions <- function(values, vline_spec, nmad_value) {
    values <- values[!is.na(values)]
    if (length(values) == 0 || is.null(vline_spec)) {
      return(data.frame(xintercept = numeric(0), stringsAsFactors = FALSE))
    }

    if (is.numeric(vline_spec)) {
      return(data.frame(xintercept = vline_spec, stringsAsFactors = FALSE))
    }

    median_value <- stats::median(values)

    if (vline_spec == "mean") {
      return(data.frame(xintercept = mean(values), stringsAsFactors = FALSE))
    }

    if (vline_spec == "median") {
      return(data.frame(xintercept = median_value, stringsAsFactors = FALSE))
    }

    mad_value <- stats::mad(values, na.rm = TRUE)
    if (is.na(mad_value)) {
      return(data.frame(xintercept = numeric(0), stringsAsFactors = FALSE))
    }

    if (vline_spec == "upper") {
      return(data.frame(xintercept = median_value + nmad_value * mad_value, stringsAsFactors = FALSE))
    }

    if (vline_spec == "lower") {
      return(data.frame(xintercept = median_value - nmad_value * mad_value, stringsAsFactors = FALSE))
    }

    if (vline_spec == "both") {
      return(data.frame(
        xintercept = c(median_value - nmad_value * mad_value, median_value + nmad_value * mad_value),
        stringsAsFactors = FALSE
      ))
    }

    return(data.frame(xintercept = numeric(0), stringsAsFactors = FALSE))
  }

  ComputeMedianPositions <- function(values) {
    values <- values[!is.na(values)]
    if (length(values) == 0) {
      return(data.frame(xintercept = numeric(0), stringsAsFactors = FALSE))
    }
    return(data.frame(xintercept = stats::median(values), stringsAsFactors = FALSE))
  }

  BuildSingleDensityPanel <- function(feature_df, x_label, panel_title) {
    current_plot <- ggplot2::ggplot(feature_df, ggplot2::aes(x = value)) +
      ggplot2::geom_density(fill = "lightblue", alpha = alpha, color = "black")

    if (pt.size > 0) {
      current_plot <- current_plot + ggplot2::geom_rug(sides = "b", linewidth = pt.size, alpha = 0.35)
    }

    vline_df <- ComputeVlinePositions(feature_df$value, vline, nmad)
    if (nrow(vline_df) > 0) {
      current_plot <- current_plot + ggplot2::geom_vline(
        data = vline_df,
        ggplot2::aes(xintercept = xintercept),
        color = "red",
        linetype = "dashed",
        linewidth = 0.6,
        show.legend = FALSE
      )
    }

    if (plot.median) {
      median_df <- ComputeMedianPositions(feature_df$value)
      if (nrow(median_df) > 0) {
        current_plot <- current_plot + ggplot2::geom_vline(
          data = median_df,
          ggplot2::aes(xintercept = xintercept),
          color = "black",
          linetype = "dashed",
          linewidth = 0.7,
          show.legend = FALSE
        )
      }
    }

    return(
      current_plot +
        ggplot2::theme_bw() +
        ggplot2::labs(title = panel_title, x = x_label, y = "Density") +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    )
  }

  # ─────────────────────────────────────────────────────────────────────────────
  # 5) Create a single-level parallel task index
  # This avoids nested parallel loops and lets future scheduling handle load.
  # ─────────────────────────────────────────────────────────────────────────────
  if (has_grouping && split.plot) {
    task_index <- expand.grid(
      feature_id = seq_along(features),
      group_id = seq_along(group_levels),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
  } else {
    task_index <- data.frame(
      feature_id = seq_along(features),
      group_id = NA_integer_,
      stringsAsFactors = FALSE
    )
  }

  task_results <- future.apply::future_lapply(
    X = seq_len(nrow(task_index)),
    FUN = function(task_row_id) {
      feature_id <- task_index$feature_id[[task_row_id]]
      feature_name <- features[[feature_id]]
      feature_title <- feature_titles[[feature_id]]

      feature_df <- plot_data[, c(feature_name, ".group"), drop = FALSE]
      names(feature_df)[1] <- "value"
      feature_df <- feature_df[!is.na(feature_df$value), , drop = FALSE]

      if (nrow(feature_df) == 0) {
        stop(paste0("Feature '", feature_name, "' has no non-missing numeric values to plot."))
      }

      if (!has_grouping) {
        return(list(
          feature_id = feature_id,
          group_id = NA_integer_,
          plot = BuildSingleDensityPanel(
            feature_df = feature_df,
            x_label = feature_name,
            panel_title = feature_title
          )
        ))
      }

      if (!split.plot) {
        current_plot <- ggplot2::ggplot(
          feature_df,
          ggplot2::aes(x = value, color = .group, fill = .group)
        ) +
          ggplot2::geom_density(alpha = alpha)

        if (pt.size > 0) {
          current_plot <- current_plot + ggplot2::geom_rug(
            ggplot2::aes(color = .group),
            sides = "b",
            linewidth = pt.size,
            alpha = 0.35,
            show.legend = FALSE
          )
        }

        if (!is.null(vline)) {
          vline_df <- feature_df |>
            dplyr::group_by(.group) |>
            dplyr::group_modify(~ComputeVlinePositions(.x$value, vline, nmad)) |>
            dplyr::ungroup()

          if (nrow(vline_df) > 0) {
            current_plot <- current_plot + ggplot2::geom_vline(
              data = vline_df,
              ggplot2::aes(xintercept = xintercept, group = .group),
              color = "red",
              linetype = "dashed",
              linewidth = 0.55,
              show.legend = FALSE
            )
          }
        }

        if (plot.median) {
          median_df <- feature_df |>
            dplyr::group_by(.group) |>
            dplyr::summarise(xintercept = stats::median(value, na.rm = TRUE), .groups = "drop")

          if (nrow(median_df) > 0) {
            current_plot <- current_plot + ggplot2::geom_vline(
              data = median_df,
              ggplot2::aes(xintercept = xintercept, group = .group),
              color = "black",
              linetype = "dashed",
              linewidth = 0.65,
              show.legend = FALSE
            )
          }
        }

        return(list(
          feature_id = feature_id,
          group_id = NA_integer_,
          plot = current_plot +
            ggplot2::scale_color_viridis_d(option = scale.colors, name = group_label) +
            ggplot2::scale_fill_viridis_d(option = scale.colors, name = group_label) +
            ggplot2::theme_bw() +
            ggplot2::labs(title = feature_title, x = feature_name, y = "Density") +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
        ))
      }

      group_id <- task_index$group_id[[task_row_id]]
      group_name <- group_levels[[group_id]]
      group_df <- feature_df[feature_df$.group == group_name, , drop = FALSE]

      return(list(
        feature_id = feature_id,
        group_id = group_id,
        plot = BuildSingleDensityPanel(
          feature_df = group_df,
          x_label = feature_name,
          panel_title = group_name
        )
      ))
    },
    future.seed = TRUE,
    future.scheduling = 2
  )

  # ─────────────────────────────────────────────────────────────────────────────
  # 6) Assemble per-feature outputs
  # ─────────────────────────────────────────────────────────────────────────────
  feature_plots <- vector("list", length(features))

  for (feature_id in seq_along(features)) {
    feature_title <- feature_titles[[feature_id]]
    feature_results <- task_results[vapply(task_results, function(item) item$feature_id == feature_id, logical(1))]

    if (has_grouping && split.plot) {
      ordered_results <- feature_results[order(vapply(feature_results, function(item) item$group_id, integer(1)))]
      group_plots <- lapply(ordered_results, function(item) item$plot)
      ncol_groups <- if (!is.null(ncol)) ncol else ceiling(sqrt(length(group_plots)))

      feature_plots[[feature_id]] <- patchwork::wrap_plots(group_plots, ncol = ncol_groups) +
        patchwork::plot_annotation(
          title = feature_title,
          theme = ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
          )
        )
    } else {
      feature_plots[[feature_id]] <- feature_results[[1]]$plot
    }
  }

  names(feature_plots) <- features

  if (length(feature_plots) == 1) {
    return(feature_plots[[1]])
  }

  return(feature_plots)
}
