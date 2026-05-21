#' @title VlnPlotGradient
#' @description `VlnPlotGradient()` extends `Seurat::VlnPlot()` functionality
#'   by coloring each violin plot with a gradient based on the per-identity value
#'   of a given feature. Particularly useful for visualizing QC metrics (nCount_RNA,
#'   nFeature_RNA, percent.mt) simultaneously while coloring by number of cells per
#'   identity or any other aggregated feature value.
#'
#' @param SeuratObject A Seurat object.
#' @param features Character vector of features to plot (gene names or metadata columns).
#' @param gradient Character string. The feature used to compute the per-identity
#'   color value. Use \code{"nCells"} to color by the number of cells per identity.
#'   Any metadata column or gene is also accepted (its mean per identity is used).
#' @param group.by Character string. Metadata column to group cells by. If NULL
#'   (default), uses the active identity (\code{Idents(SeuratObject)}).
#' @param scale.colors Character string. Viridis palette to use for the gradient.
#'   One of \code{"magma"} / \code{"A"}, \code{"inferno"} / \code{"B"},
#'   \code{"plasma"} / \code{"C"}, \code{"viridis"} / \code{"D"},
#'   \code{"cividis"} / \code{"E"}, \code{"rocket"} / \code{"F"},
#'   \code{"mako"} / \code{"G"}, \code{"turbo"} / \code{"H"}.
#'   Default is \code{"viridis"}.
#' @param lower.limit Numeric. Lower limit of the gradient scale.
#'   Default is \code{0}. Only applied when \code{upper.limit} is specified.
#' @param upper.limit Numeric or NULL. Upper limit of the gradient scale.
#'   When \code{NULL} (default), limits are set automatically.
#' @param pt.size Numeric. Size of the jittered points overlaid on the violins.
#'   Set to \code{0} to hide points. Default is \code{0.1}.
#' @param ncol Integer or NULL. Number of columns in the combined plot.
#'   If NULL, Seurat's default (number of features) is used.
#'
#' @return A \code{ggplot2} / \code{patchwork} object.
#'
#' @examples
#' \dontrun{
#' VlnPlotGradient(
#'   SeuratObject,
#'   features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
#'   gradient = "nCells"
#' )
#'
#' VlnPlotGradient(
#'   SeuratObject,
#'   features = c("CD3D", "CD8A"),
#'   gradient = "nCount_RNA",
#'   group.by = "seurat_clusters",
#'   scale.colors = "plasma",
#'   upper.limit = 5000
#' )
#'}
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
    group.by     = NULL,
    scale.colors = "viridis",
    lower.limit  = 0,
    upper.limit  = NULL,
    pt.size      = 0.1,
    ncol         = NULL
) {

  # ── 1. Resolve identity ──────────────────────────────────────────────────────────────────────
  original_idents <- NULL
  if (!is.null(group.by)) {
    if (!group.by %in% colnames(SeuratObject@meta.data)) {
      stop(paste0("'group.by' column '", group.by, "' not found in meta.data."))
    }
    original_idents <- Seurat::Idents(SeuratObject)
    Seurat::Idents(SeuratObject) <- SeuratObject@meta.data[[group.by]]
  }
  cell_idents <- as.character(Seurat::Idents(SeuratObject))

  # ── 2. Compute per-identity gradient value ───────────────────────────────────────────────────────
  if (gradient == "nCells") {
    ident_df <- data.frame(identity = cell_idents, stringsAsFactors = FALSE)
    gradient_values <- ident_df %>%
      dplyr::count(identity, name = "nCells") %>%
      dplyr::rename(gradient_val = nCells)
    gradient_label <- "nCells"

  } else {
    feat_data <- tryCatch(
      Seurat::FetchData(SeuratObject, vars = gradient),
      error = function(e) stop(paste0("Cannot fetch gradient feature '", gradient, "': ", e$message))
    )
    gradient_values <- data.frame(
      identity     = cell_idents,
      feature_val  = feat_data[[gradient]],
      stringsAsFactors = FALSE
    ) %>%
      dplyr::group_by(identity) %>%
      dplyr::summarise(gradient_val = mean(feature_val, na.rm = TRUE), .groups = "drop")
    gradient_label <- gradient
  }

  identity_levels <- levels(Seurat::Idents(SeuratObject))
  if (is.null(identity_levels)) identity_levels <- sort(unique(cell_idents))

  gradient_values <- gradient_values %>%
    dplyr::filter(identity %in% identity_levels) %>%
    dplyr::arrange(match(identity, identity_levels))

  color_values <- setNames(gradient_values$gradient_val, gradient_values$identity)

  # ── 3. Build per-feature violin plots ───────────────────────────────────────────────────────────
  plots <- lapply(features, function(feat) {

    p <- Seurat::VlnPlot(
      SeuratObject,
      features = feat,
      pt.size  = pt.size,
      cols     = scales::rescale(color_values[identity_levels])
    )

    if (inherits(p, "patchwork") || inherits(p, "gg")) {
      layer_data <- p$data
      if (is.null(layer_data)) {
        p_inner <- p[[1]]
      } else {
        p_inner <- p
      }
    } else {
      p_inner <- p
    }

    p_inner <- p_inner +
      aes(fill = ident) +
      scale_fill_manual(
        values = if (!is.null(upper.limit)) {
          scales::col_numeric(
            palette = viridisLite::viridis(256, option = scale.colors),
            domain  = c(lower.limit, upper.limit)
          )(color_values[identity_levels])
        } else {
          scales::col_numeric(
            palette = viridisLite::viridis(256, option = scale.colors),
            domain  = range(color_values, na.rm = TRUE)
          )(color_values[identity_levels])
        },
        name = gradient_label
      ) +
      theme(legend.position = "none")

    p_inner
  })

  # ── 4. Assemble and add shared legend ─────────────────────────────────────────────────────────
  grad_range <- if (!is.null(upper.limit)) c(lower.limit, upper.limit) else range(color_values, na.rm = TRUE)

  legend_plot <- ggplot2::ggplot(
    data.frame(x = 1, y = grad_range[1]:grad_range[2], fill = grad_range[1]:grad_range[2]),
    ggplot2::aes(x = x, y = y, fill = fill)
  ) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(name = gradient_label, option = scale.colors,
                                   limits = grad_range) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "right")

  legend_grob <- cowplot::get_legend(legend_plot)

  ncol_val <- if (is.null(ncol)) length(features) else ncol

  combined <- patchwork::wrap_plots(plots, ncol = ncol_val) +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "none")

  final_plot <- cowplot::plot_grid(combined, legend_grob, rel_widths = c(1, 0.1))

  if (!is.null(original_idents)) {
    Seurat::Idents(SeuratObject) <- original_idents
  }

  return(final_plot)
}
