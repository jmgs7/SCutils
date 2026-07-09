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
#' - optional vertical reference lines using mean/median, median +/- n*MAD, or a
#'   fixed numeric threshold,
#' - multi-feature plot composition with patchwork.
#'
#' @param SeuratObject A Seurat object.
#' @param features Character vector of metadata columns to plot on the x-axis.
#' @param group.by Character scalar. Grouping variable. Use a metadata column name,
#'   `"active.ident"` (default), or `NULL` for no grouping.
#' @param split.plot Logical. If `FALSE` (default), groups are overlaid in one panel
#'   per feature. If `TRUE`, creates one panel per group level for each feature.
#' @param scale.colors Character scalar. Viridis palette option used for grouped
#'   density colors (`"viridis"`, `"magma"`, `"plasma"`, `"inferno"`,
#'   `"cividis"`, `"rocket"`, `"mako"`, `"turbo"`). Default is `"viridis"`.
#' @param ncol Integer or `NULL`. Number of columns when combining multi-feature
#'   outputs with patchwork. If `NULL`, uses `length(features)`.
#' @param vline Optional reference-line specification. Accepted values:
#'   `NULL`, `"mean"`, `"median"`, `"upper"`, `"lower"`, `"both"`, or a
#'   numeric value. `"upper"`, `"lower"`, and `"both"` use median +/- `nmad`*MAD.
#' @param nmad Numeric. Number of MADs used when `vline` is `"upper"`, `"lower"`,
#'   or `"both"`. Default is `2`.
#' @param alpha Numeric in `[0, 1]`. Fill alpha for density geometries.
#' @param pt.size Numeric. If `0` (default), no rug is drawn. If `> 0`, adds a rug
#'   (`geom_rug()`) with this line width.
#' @param mc.cores Integer or `NULL`. Number of cores for parallel plotting with
#'   `parallel::mclapply()`. If `NULL`, defaults to one core per feature (or per
#'   group when `split.plot = TRUE`). On Windows, forced to `1`.
#'
#' @return A `ggplot2` or `patchwork` object.
#'
#' @examples
#' \dontrun{
#' FeatureDensityPlot(
#'   SeuratObject,
#'   features = c("nCount_RNA", "nFeature_RNA", "percent.mt")
#' )
#'
#' FeatureDensityPlot(
#'   SeuratObject,
#'   features = c("percent.mt", "complexity"),
#'   group.by = "batch",
#'   split.plot = FALSE,
#'   vline = "upper",
#'   nmad = 2.5,
#'   mc.cores = 4
#' )
#'
#' FeatureDensityPlot(
#'   SeuratObject,
#'   features = "percent.mt",
#'   group.by = "library",
#'   split.plot = TRUE,
#'   vline = 10
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
    split.plot = FALSE,
    scale.colors = "viridis",
    ncol = NULL,
    vline = NULL,
    nmad = 2,
    alpha = 0.3,
    pt.size = 0,
    mc.cores = NULL
) {
  # Validate basic inputs early to fail fast with clear messages.
  if (!inherits(SeuratObject, "Seurat")) {
    stop("'SeuratObject' must be a Seurat object.")
  }

  if (!is.character(features) || length(features) == 0) {
    stop("'features' must be a non-empty character vector.")
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

  # Validate vline type/specification.
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

  md <- SeuratObject@meta.data

  # Restrict to metadata features by design.
  missing_features <- setdiff(features, colnames(md))
  if (length(missing_features) > 0) {
    stop(
      paste0(
        "These features are not metadata columns: ",
        paste(missing_features, collapse = ", ")
      )
    )
  }

  # Build grouping vector directly from metadata or active identities.
  has_grouping <- !is.null(group.by)
  if (has_grouping) {
    if (!is.character(group.by) || length(group.by) != 1) {
      stop("'group.by' must be NULL or a single character value.")
    }

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

  # Keep only required columns for plotting to reduce memory footprint.
  plot_df <- md[, features, drop = FALSE]
  plot_df$.group <- group_vec
  plot_df <- plot_df[!is.na(plot_df$.group), , drop = FALSE]

  # Convert each feature to numeric as density requires continuous data.
  for (feat in features) {
    if (!is.numeric(plot_df[[feat]])) {
      suppressWarnings(plot_df[[feat]] <- as.numeric(plot_df[[feat]]))
    }
  }

  # Remove groups with no cells after filtering.
  group_levels <- unique(plot_df$.group)
  if (length(group_levels) == 0) {
    stop("No groups available to plot after filtering missing grouping values.")
  }

  # Default core strategy: per-group for split mode, per-feature otherwise.
  if (is.null(mc.cores)) {
    mc.cores <- if (has_grouping && split.plot) length(group_levels) else length(features)
  }

  if (!is.numeric(mc.cores) || length(mc.cores) != 1 || is.na(mc.cores) || mc.cores < 1) {
    stop("'mc.cores' must be NULL or a positive integer.")
  }

  mc.cores <- as.integer(mc.cores)

  # mclapply runs as single-core on Windows.
  if (.Platform$OS.type == "windows") {
    mc.cores <- 1
  }

  # Helper to compute vline positions from a numeric vector.
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

  # Build one plot per feature; this is parallelized with mclapply.
  feature_plots <- parallel::mclapply(
    X = features,
    FUN = function(feat) {
      feat_df <- plot_df[, c(feat, ".group"), drop = FALSE]
      names(feat_df)[1] <- "value"
      feat_df <- feat_df[!is.na(feat_df$value), , drop = FALSE]

      if (nrow(feat_df) == 0) {
        stop(paste0("Feature '", feat, "' has no non-missing numeric values to plot."))
      }

      if (!has_grouping) {
        # No grouping: a single density trace.
        p <- ggplot2::ggplot(feat_df, ggplot2::aes(x = value)) +
          ggplot2::geom_density(fill = "lightblue", alpha = alpha, color = "black")

        if (pt.size > 0) {
          p <- p + ggplot2::geom_rug(sides = "b", linewidth = pt.size, alpha = 0.35)
        }

        ref_df <- compute_vline_positions(feat_df$value, vline, nmad)
        if (nrow(ref_df) > 0) {
          p <- p + ggplot2::geom_vline(
            data = ref_df,
            ggplot2::aes(xintercept = xintercept, linetype = type),
            color = if (is.character(vline) && vline %in% c("upper", "lower", "both")) "red" else "black",
            linewidth = 0.6,
            show.legend = FALSE
          )
        }

        return(
          p +
            ggplot2::theme_bw() +
            ggplot2::labs(title = feat, x = feat, y = "Density") +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
        )
      }

      if (!split.plot) {
        # Grouped overlay: one panel with densities per group.
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

        # Compute reference lines per group so thresholds match each group's distribution.
        if (!is.null(vline)) {
          ref_df <- feat_df |>
            dplyr::group_by(.group) |>
            dplyr::group_modify(~compute_vline_positions(.x$value, vline, nmad)) |>
            dplyr::ungroup()

          if (nrow(ref_df) > 0) {
            p <- p + ggplot2::geom_vline(
              data = ref_df,
              ggplot2::aes(xintercept = xintercept, color = .group),
              linetype = "dashed",
              linewidth = 0.55,
              show.legend = FALSE
            )
          }
        }

        return(
          p +
            ggplot2::scale_color_viridis_d(option = scale.colors, name = group_label) +
            ggplot2::scale_fill_viridis_d(option = scale.colors, name = group_label) +
            ggplot2::theme_bw() +
            ggplot2::labs(title = feat, x = feat, y = "Density") +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
        )
      }

      # Split mode: one panel per group level for this feature.
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
              color = if (is.character(vline) && vline %in% c("upper", "lower", "both")) "red" else "black",
              linewidth = 0.6,
              show.legend = FALSE
            )
          }

          p +
            ggplot2::theme_bw() +
            ggplot2::labs(title = grp, x = feat, y = "Density") +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 10))
        },
        mc.cores = max(1L, min(mc.cores, length(split_groups)))
      )

      # Layout within each feature's split panels.
      ncol_groups <- ceiling(sqrt(length(split_groups)))

      patchwork::wrap_plots(group_plots, ncol = ncol_groups) +
        patchwork::plot_annotation(title = feat)
    },
    mc.cores = max(1L, min(mc.cores, length(features)))
  )

  if (length(feature_plots) == 1) {
    return(feature_plots[[1]])
  }

  # Layout across features.
  n_cols <- if (is.null(ncol)) length(features) else ncol

  patchwork::wrap_plots(feature_plots, ncol = n_cols) +
    patchwork::plot_layout(guides = "collect")
}
