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
#' - optional vertical reference lines (`vline`) in red,
#' - optional independent median overlays (`plot.median`) in black,
#' - custom plot titles via `plot.title`,
#' - multi-feature output as a named list (one plot per feature).
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
#' @param vline Optional reference-line specification (drawn in red). Accepted values:
#'   `NULL`, `"mean"`, `"median"`, `"upper"`, `"lower"`, `"both"`, or a
#'   numeric value. `"upper"`, `"lower"`, and `"both"` use median +/- `nmad`*MAD.
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
#' @param mc.cores Integer or `NULL`. Number of cores for parallel plotting with
#'   `parallel::mclapply()`. If `NULL`, defaults to one core per feature (or per
#'   group when `split.plot = TRUE`). On Windows, forced to `1`.
#'
#' @return If `length(features) == 1`, returns a `ggplot2`/`patchwork` plot object.
#'   If `length(features) > 1`, returns a named list of plot objects (one per feature;
#'   list names equal `features`).
#'
#' @examples
#' \dontrun{
#' # Single feature with default split-by-group mode and median line.
#' FeatureDensityPlot(
#'   SeuratObject,
#'   features = "percent.mt"
#' )
#'
#' # Multi-feature output returns a named list of plots.
#' plt_list <- FeatureDensityPlot(
#'   SeuratObject,
#'   features = c("percent.mt", "nCount_RNA"),
#'   group.by = "batch",
#'   split.plot = FALSE,
#'   vline = "upper",
#'   plot.median = TRUE,
#'   nmad = 2.5,
#'   mc.cores = 4
#' )
#'
#' # Custom title for a single feature.
#' FeatureDensityPlot(
#'   SeuratObject,
#'   features = "percent.mt",
#'   group.by = "library",
#'   plot.title = "Mitochondrial Percentage Density"
#' )
#' }
#'
#' @import ggplot2
#' @import dplyr
#' @import patchwork
#' @import Seurat
#' @importFrom parallel mclapply
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
  # Guardrails are explicit to fail early with clear error messages.
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

  # vline contract is intentionally strict because plotting semantics depend on type.
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
  # No SplitObject() is used; all grouping is resolved directly from metadata.
  # ─────────────────────────────────────────────────────────────────────────────
  md <- SeuratObject@meta.data

  missing_features <- setdiff(features, colnames(md))
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

    # `active.ident` is treated as a virtual grouping column for convenience.
    if (group.by == "active.ident") {
      group_vec <- as.character(Seurat::Idents(SeuratObject))
      group_label <- "active.ident"
    } else {
      if (!group.by %in% colnames(md)) {
        stop(paste0("'group.by' column '", group.by, "' not found in SeuratObject@meta.data."))
      }
      group_vec <- as.character(md[[group.by]])
      group_label <- group.by
    }
  } else {
    group_vec <- rep("all", nrow(md))
    group_label <- "all"
  }

  # Keep only required data in memory to reduce overhead on large objects.
  plot_df <- md[, features, drop = FALSE]
  plot_df$.group <- group_vec
  plot_df <- plot_df[!is.na(plot_df$.group), , drop = FALSE]

  # Density plots require numeric x-values; coerce metadata columns when needed.
  for (feat in features) {
    if (!is.numeric(plot_df[[feat]])) {
      suppressWarnings(plot_df[[feat]] <- as.numeric(plot_df[[feat]]))
    }
  }

  group_levels <- unique(plot_df$.group)
  if (length(group_levels) == 0) {
    stop("No groups available to plot after filtering missing grouping values.")
  }

  # ─────────────────────────────────────────────────────────────────────────────
  # 3) Parallel strategy
  # One core per feature by default, or per group when split plotting is enabled.
  # Windows falls back to one core due to mclapply limitations.
  # ─────────────────────────────────────────────────────────────────────────────
  if (is.null(mc.cores)) {
    mc.cores <- if (has_grouping && split.plot) length(group_levels) else length(features)
  }

  if (!is.numeric(mc.cores) || length(mc.cores) != 1 || is.na(mc.cores) || mc.cores < 1) {
    stop("'mc.cores' must be NULL or a positive integer.")
  }

  mc.cores <- as.integer(mc.cores)
  if (.Platform$OS.type == "windows") {
    mc.cores <- 1
  }

  # ─────────────────────────────────────────────────────────────────────────────
  # 4) Helper functions for line positions
  # These helpers isolate statistical line calculations and NA handling.
  # ─────────────────────────────────────────────────────────────────────────────
  compute_vline_positions <- function(x, vline_spec, nmad_value) {
    x <- x[!is.na(x)]
    if (length(x) == 0 || is.null(vline_spec)) {
      return(data.frame(type = character(0), xintercept = numeric(0), stringsAsFactors = FALSE))
    }

    if (is.numeric(vline_spec)) {
      return(data.frame(type = "fixed", xintercept = vline_spec, stringsAsFactors = FALSE))
    }

    med <- stats::median(x)

    if (vline_spec == "mean") {
      return(data.frame(type = "mean", xintercept = mean(x), stringsAsFactors = FALSE))
    }

    if (vline_spec == "median") {
      return(data.frame(type = "median", xintercept = med, stringsAsFactors = FALSE))
    }

    mad_val <- stats::mad(x, na.rm = TRUE)
    if (is.na(mad_val)) {
      return(data.frame(type = character(0), xintercept = numeric(0), stringsAsFactors = FALSE))
    }

    if (vline_spec == "upper") {
      return(data.frame(type = "upper", xintercept = med + nmad_value * mad_val, stringsAsFactors = FALSE))
    }

    if (vline_spec == "lower") {
      return(data.frame(type = "lower", xintercept = med - nmad_value * mad_val, stringsAsFactors = FALSE))
    }

    if (vline_spec == "both") {
      return(data.frame(
        type = c("lower", "upper"),
        xintercept = c(med - nmad_value * mad_val, med + nmad_value * mad_val),
        stringsAsFactors = FALSE
      ))
    }

    data.frame(type = character(0), xintercept = numeric(0), stringsAsFactors = FALSE)
  }

  # Median helper is intentionally separate from vline so `plot.median` remains
  # independent and can coexist with any vline setting.
  compute_median_positions <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) {
      return(data.frame(xintercept = numeric(0), stringsAsFactors = FALSE))
    }
    data.frame(xintercept = stats::median(x), stringsAsFactors = FALSE)
  }

  # Expand titles to per-feature vector after validation.
  feature_titles <- if (is.null(plot.title)) {
    features
  } else if (length(plot.title) == 1) {
    rep(plot.title, length(features))
  } else {
    plot.title
  }

  # ─────────────────────────────────────────────────────────────────────────────
  # 5) Build one plot per feature (parallelized)
  # Branches: no-grouping, grouped-overlay, grouped-split.
  # Color contract: vline is always red; median overlay is always black.
  # ─────────────────────────────────────────────────────────────────────────────
  feature_plots <- parallel::mclapply(
    X = seq_along(features),
    FUN = function(i) {
      feat <- features[[i]]
      feat_title <- feature_titles[[i]]

      feat_df <- plot_df[, c(feat, ".group"), drop = FALSE]
      names(feat_df)[1] <- "value"
      feat_df <- feat_df[!is.na(feat_df$value), , drop = FALSE]

      if (nrow(feat_df) == 0) {
        stop(paste0("Feature '", feat, "' has no non-missing numeric values to plot."))
      }

      if (!has_grouping) {
        p <- ggplot2::ggplot(feat_df, ggplot2::aes(x = value)) +
          ggplot2::geom_density(fill = "lightblue", alpha = alpha, color = "black")

        if (pt.size > 0) {
          p <- p + ggplot2::geom_rug(sides = "b", linewidth = pt.size, alpha = 0.35)
        }

        # vline rendered in red by design.
        ref_df <- compute_vline_positions(feat_df$value, vline, nmad)
        if (nrow(ref_df) > 0) {
          p <- p + ggplot2::geom_vline(
            data = ref_df,
            ggplot2::aes(xintercept = xintercept, linetype = type),
            color = "red",
            linewidth = 0.6,
            show.legend = FALSE
          )
        }

        # Median overlay rendered in black and added last for visibility.
        if (plot.median) {
          med_df <- compute_median_positions(feat_df$value)
          if (nrow(med_df) > 0) {
            p <- p + ggplot2::geom_vline(
              data = med_df,
              ggplot2::aes(xintercept = xintercept),
              color = "black",
              linewidth = 0.7,
              linetype = "dashed",
              show.legend = FALSE
            )
          }
        }

        return(
          p +
            ggplot2::theme_bw() +
            ggplot2::labs(title = feat_title, x = feat, y = "Density") +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
        )
      }

      if (!split.plot) {
        p <- ggplot2::ggplot(
          feat_df,
          ggplot2::aes(x = value, color = .group, fill = .group)
        ) +
          ggplot2::geom_density(alpha = alpha)

        if (pt.size > 0) {
          p <- p + ggplot2::geom_rug(
            ggplot2::aes(color = .group),
            sides = "b",
            linewidth = pt.size,
            alpha = 0.35,
            show.legend = FALSE
          )
        }

        # Group-level vlines are computed per group and drawn in red.
        if (!is.null(vline)) {
          ref_df <- feat_df |>
            dplyr::group_by(.group) |>
            dplyr::group_modify(~compute_vline_positions(.x$value, vline, nmad)) |>
            dplyr::ungroup()

          if (nrow(ref_df) > 0) {
            p <- p + ggplot2::geom_vline(
              data = ref_df,
              ggplot2::aes(xintercept = xintercept, group = .group),
              color = "red",
              linetype = "dashed",
              linewidth = 0.55,
              show.legend = FALSE
            )
          }
        }

        # Group-level medians are drawn in black.
        if (plot.median) {
          med_df <- feat_df |>
            dplyr::group_by(.group) |>
            dplyr::summarise(xintercept = stats::median(value, na.rm = TRUE), .groups = "drop")

          if (nrow(med_df) > 0) {
            p <- p + ggplot2::geom_vline(
              data = med_df,
              ggplot2::aes(xintercept = xintercept, group = .group),
              color = "black",
              linetype = "dashed",
              linewidth = 0.65,
              show.legend = FALSE
            )
          }
        }

        return(
          p +
            ggplot2::scale_color_viridis_d(option = scale.colors, name = group_label) +
            ggplot2::scale_fill_viridis_d(option = scale.colors, name = group_label) +
            ggplot2::theme_bw() +
            ggplot2::labs(title = feat_title, x = feat, y = "Density") +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
        )
      }

      # Split mode: one panel per group for the current feature.
      split_groups <- unique(feat_df$.group)
      group_plots <- parallel::mclapply(
        X = split_groups,
        FUN = function(grp) {
          grp_df <- feat_df[feat_df$.group == grp, , drop = FALSE]

          p <- ggplot2::ggplot(grp_df, ggplot2::aes(x = value)) +
            ggplot2::geom_density(fill = "lightblue", alpha = alpha, color = "black")

          if (pt.size > 0) {
            p <- p + ggplot2::geom_rug(sides = "b", linewidth = pt.size, alpha = 0.35)
          }

          ref_df <- compute_vline_positions(grp_df$value, vline, nmad)
          if (nrow(ref_df) > 0) {
            p <- p + ggplot2::geom_vline(
              data = ref_df,
              ggplot2::aes(xintercept = xintercept, linetype = type),
              color = "red",
              linewidth = 0.6,
              show.legend = FALSE
            )
          }

          if (plot.median) {
            med_df <- compute_median_positions(grp_df$value)
            if (nrow(med_df) > 0) {
              p <- p + ggplot2::geom_vline(
                data = med_df,
                ggplot2::aes(xintercept = xintercept),
                color = "black",
                linewidth = 0.7,
                linetype = "dashed",
                show.legend = FALSE
              )
            }
          }

          p +
            ggplot2::theme_bw() +
            ggplot2::labs(title = grp, x = feat, y = "Density") +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 10))
        },
        mc.cores = max(1L, min(mc.cores, length(split_groups)))
      )

      ncol_groups <- if (!is.null(ncol)) ncol else ceiling(sqrt(length(split_groups)))

      # Global split title is explicitly centered.
      patchwork::wrap_plots(group_plots, ncol = ncol_groups) +
        patchwork::plot_annotation(
          title = feat_title,
          theme = ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
          )
        )
    },
    mc.cores = max(1L, min(mc.cores, length(features)))
  )

  # ─────────────────────────────────────────────────────────────────────────────
  # 6) Return shape
  # Single feature -> single plot object.
  # Multiple features -> named list keyed by feature names.
  # ─────────────────────────────────────────────────────────────────────────────
  names(feature_plots) <- features

  if (length(feature_plots) == 1) {
    return(feature_plots[[1]])
  }

  return(feature_plots)
}
