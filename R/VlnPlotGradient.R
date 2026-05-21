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
  # If group.by is provided, temporarily set the active identity to that column.
  # This mirrors Seurat::VlnPlot() behaviour exactly.
  original_idents <- NULL
  if (!is.null(group.by)) {
    if (!group.by %in% colnames(SeuratObject@meta.data)) {
      stop(paste0("'group.by' column '", group.by, "' not found in meta.data."))
    }
    original_idents <- Seurat::Idents(SeuratObject)
    Seurat::Idents(SeuratObject) <- SeuratObject@meta.data[[group.by]]
  }
  # Capture the identity labels for each cell
  cell_idents <- as.character(Seurat::Idents(SeuratObject))

  # ── 2. Compute per-identity gradient value ───────────────────────────────────────────────────────
  if (gradient == "nCells") {
    # Special case: count cells per identity using dplyr
    ident_df <- data.frame(identity = cell_idents, stringsAsFactors = FALSE)
    gradient_values <- ident_df %>%
      dplyr::count(identity, name = "nCells") %>%
      dplyr::rename(gradient_val = nCells)
    gradient_label <- "nCells"

  } else {
    # General case: use the mean of the feature across cells in each identity.
    # Works for metadata columns and genes (Seurat::FetchData handles both).
    feat_data <- tryCatch(
      Seurat::FetchData(SeuratObject, vars = gradient),
      error = function(e) stop(paste0("Cannot fetch gradient feature '", gradient, "': ", e$message))
    )
    ident_df <- data.frame(
      identity     = cell_idents,
      gradient_raw = feat_data[[1]],
      stringsAsFactors = FALSE
    )
    gradient_values <- ident_df %>%
      dplyr::group_by(identity) %>%
      dplyr::summarise(gradient_val = mean(gradient_raw, na.rm = TRUE), .groups = "drop")
    gradient_label <- paste0("mean(", gradient, ")")
  }

  # ── 2b. Order identities by gradient value descending ──────────────────────────────────────────
  # Sort the gradient table so the identity with the highest gradient value comes
  # first. The resulting order vector is used as the factor level order in ggplot,
  # which controls the left-to-right position of violins on the x-axis.
  gradient_values <- gradient_values %>%
    dplyr::arrange(dplyr::desc(gradient_val))

  # Ordered identity levels (high → low gradient, left → right on x-axis)
  ordered_levels <- gradient_values$identity

  # ── 3. Restore original identity if changed ─────────────────────────────────────────────────────
  if (!is.null(original_idents)) {
    Seurat::Idents(SeuratObject) <- original_idents
  }

  # ── 4. Build per-feature violin plots ─────────────────────────────────────────────────────────────
  # Fetch feature data for all requested features at once
  feat_matrix <- tryCatch(
    Seurat::FetchData(SeuratObject, vars = features),
    error = function(e) stop(paste0("Cannot fetch one or more features: ", e$message))
  )

  # Global color scale limits (consistent across all panels)
  if (!is.null(upper.limit)) {
    scale_limits <- c(lower.limit, upper.limit)
  } else {
    scale_limits <- NULL  # let ggplot auto-scale
  }

  # Build one ggplot per feature
  plot_list <- lapply(features, function(feat) {

    # Per-cell data frame; join gradient values first as plain character,
    # then re-apply the ordered factor AFTER the join to prevent left_join()
    # from silently dropping the factor class and reverting to alphabetical order.
    cell_df <- data.frame(
      identity = cell_idents,
      value    = feat_matrix[[feat]],
      stringsAsFactors = FALSE
    ) %>%
      dplyr::left_join(gradient_values, by = "identity") %>%
      dplyr::mutate(identity = factor(identity, levels = ordered_levels))

    # ── ggplot2 violin, styled to resemble Seurat::VlnPlot() ─────────────────
    p <- ggplot2::ggplot(cell_df,
           ggplot2::aes(x = identity, y = value, fill = gradient_val)) +

      # Violin body
      ggplot2::geom_violin(
        scale     = "width",
        trim      = TRUE,
        adjust    = 1,
        linewidth = 0.3,
        color     = "black"
      ) +

      # Jittered individual points (optional)
      {
        if (pt.size > 0) {
          ggplot2::geom_jitter(
            ggplot2::aes(color = gradient_val),
            width  = 0.3,
            size   = pt.size,
            alpha  = 0.6,
            show.legend = FALSE
          )
        }
      } +

      # Viridis gradient fill for violins
      ggplot2::scale_fill_viridis_c(
        name     = gradient_label,
        option   = scale.colors,
        limits   = scale_limits,
        na.value = "grey70"
      ) +

      # Matching color scale for points
      ggplot2::scale_color_viridis_c(
        option = scale.colors,
        limits = scale_limits,
        guide  = "none"
      ) +

      # Seurat-like axis labels
      ggplot2::labs(
        title = feat,
        x     = NULL,
        y     = NULL
      ) +

      # Seurat-inspired theme
      ggplot2::theme_classic(base_size = 11) +
      ggplot2::theme(
        # Title
        plot.title         = ggplot2::element_text(
          face = "bold", size = 11, hjust = 0.5
        ),
        # X-axis identity labels — rotated as in Seurat
        axis.text.x        = ggplot2::element_text(
          angle = 45, hjust = 1, vjust = 1, size = 9
        ),
        axis.text.y        = ggplot2::element_text(size = 9),
        axis.title.y       = ggplot2::element_text(size = 9),
        # Legend
        legend.title       = ggplot2::element_text(size = 9),
        legend.text        = ggplot2::element_text(size = 8),
        legend.key.height  = ggplot2::unit(0.8, "cm"),
        legend.key.width   = ggplot2::unit(0.25, "cm"),
        # Panel
        panel.grid.major.y = ggplot2::element_line(
          color = "grey90", linewidth = 0.25
        ),
        panel.border       = ggplot2::element_blank(),
        axis.line          = ggplot2::element_line(linewidth = 0.4)
      )

    return(p)
  })

  # ── 5. Combine panels using patchwork ──────────────────────────────────────────────────────────────────
  # Share a single colour legend across all panels (collect_guides)
  n_cols <- if (!is.null(ncol)) ncol else length(features)

  combined <- patchwork::wrap_plots(plot_list, ncol = n_cols) +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "right")

  return(combined)
}
