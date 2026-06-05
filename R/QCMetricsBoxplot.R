#' @title QCMetricsBoxplot
#' @description `QCMetricsBoxplot()` creates boxplots to visualize the distribution of 
#'   QC metrics (nFeature, nCount, percent.mt) for cells within a specific entity or 
#'   across multiple entities. Individual cells are represented as points with colors 
#'   reflecting the magnitude of each QC metric, allowing for rapid visual assessment 
#'   of QC metric distributions by sample, batch, cluster, or other grouping variables.
#'
#' @param SeuratObject A Seurat object.
#' @param entity_name Character string. The name of the metadata column to group cells by 
#'   (e.g., "orig.ident", "batch", "seurat_clusters").
#' @param entity_type Character string or NULL. If specified, filters the plot to show only 
#'   cells from this specific entity. If NULL (default), displays all entities in the 
#'   metadata column.
#' @param qc_metrics Character vector. QC metrics to visualize. Default is 
#'   `c("nFeature_RNA", "nCount_RNA", "percent.mt")`. Any metadata columns can be used.
#' @param gradient_col Character string or NULL. Metadata column/feature to use for point 
#'   gradient coloring. If NULL (default), each QC metric is colored by its own values, 
#'   creating a per-metric gradient where higher values appear in warmer colors.
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
#' @param pt.size Numeric. Size of points representing individual cells. Default is \code{1}.
#' @param pt.alpha Numeric. Transparency of points (0-1). Default is \code{0.6}.
#' @param fill_color Character. Color for the boxplot fill. Default is \code{"lightblue"}.
#'   Set to NA to remove fill color.
#' @param outlier.size Numeric. Size of boxplot outliers. Default is \code{1}.
#' @param ncol Integer or NULL. Number of columns in the combined plot.
#'   If NULL (default), uses 3 columns for layout.
#'
#' @return A \code{ggplot2} / \code{patchwork} object combining boxplots for each QC metric.
#'
#' @examples
#' \dontrun{
#' # Visualize QC metrics across all samples
#' QCMetricsBoxplot(
#'   SeuratObject,
#'   entity_name = "orig.ident"
#' )
#'
#' # Visualize QC metrics for a single sample
#' QCMetricsBoxplot(
#'   SeuratObject,
#'   entity_name = "orig.ident",
#'   entity_type = "Sample_1"
#' )
#'
#' # Customize colors and aesthetics
#' QCMetricsBoxplot(
#'   SeuratObject,
#'   entity_name = "seurat_clusters",
#'   qc_metrics = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
#'   scale.colors = "plasma",
#'   pt.size = 1.5,
#'   pt.alpha = 0.7,
#'   upper.limit = 5000
#' )
#' }
#'
#' @import ggplot2
#' @import dplyr
#' @import patchwork
#' @import Seurat
#'
#' @export

QCMetricsBoxplot <- function(
    SeuratObject,
    entity_name,
    entity_type = NULL,
    qc_metrics = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    gradient_col = NULL,
    scale.colors = "viridis",
    lower.limit = 0,
    upper.limit = NULL,
    pt.size = 1,
    pt.alpha = 0.6,
    fill_color = "lightblue",
    outlier.size = 1,
    ncol = NULL
) {

  # ── 1. Input validation ──────────────────────────────────────────────────────────────────────
  
  # Check if entity_name exists in metadata
  if (!entity_name %in% colnames(SeuratObject@meta.data)) {
    stop(paste0("'entity_name' column '", entity_name, "' not found in meta.data."))
  }

  # Check if all QC metrics exist
  missing_metrics <- setdiff(qc_metrics, colnames(SeuratObject@meta.data))
  if (length(missing_metrics) > 0) {
    stop(paste0("The following QC metrics not found in meta.data: ", 
                paste(missing_metrics, collapse = ", ")))
  }

  # Check if gradient_col exists (if specified)
  if (!is.null(gradient_col)) {
    if (!gradient_col %in% colnames(SeuratObject@meta.data)) {
      stop(paste0("'gradient_col' column '", gradient_col, "' not found in meta.data."))
    }
  }

  # Validate numeric parameters
  if (!is.numeric(pt.size) || pt.size < 0) {
    stop("'pt.size' must be a positive numeric value.")
  }
  if (!is.numeric(pt.alpha) || pt.alpha < 0 || pt.alpha > 1) {
    stop("'pt.alpha' must be a numeric value between 0 and 1.")
  }
  if (!is.numeric(outlier.size) || outlier.size < 0) {
    stop("'outlier.size' must be a positive numeric value.")
  }

  # ── 2. Prepare data ──────────────────────────────────────────────────────────────────────────

  # Get metadata
  metadata <- SeuratObject@meta.data

  # Filter to entity_type if specified
  if (!is.null(entity_type)) {
    if (!entity_type %in% metadata[[entity_name]]) {
      stop(paste0("'entity_type' value '", entity_type, 
                  "' not found in '", entity_name, "' column."))
    }
    metadata <- metadata[metadata[[entity_name]] == entity_type, , drop = FALSE]
  }

  # Create list to store plots
  plot_list <- list()

  # Set default ncol if not specified
  if (is.null(ncol)) {
    ncol <- min(3, length(qc_metrics))
  }

  # ── 3. Create one boxplot per QC metric ──────────────────────────────────────────────────────

  for (metric in qc_metrics) {

    # Prepare data frame for this metric
    if (is.null(gradient_col)) {
      # Use the metric itself for coloring (per-value gradient)
      plot_data <- data.frame(
        entity = metadata[[entity_name]],
        value = metadata[[metric]],
        gradient = metadata[[metric]],
        stringsAsFactors = FALSE
      )
      gradient_label <- metric
    } else {
      # Use specified gradient column
      plot_data <- data.frame(
        entity = metadata[[entity_name]],
        value = metadata[[metric]],
        gradient = metadata[[gradient_col]],
        stringsAsFactors = FALSE
      )
      gradient_label <- gradient_col
    }

    # Remove any rows with NA values in value or entity
    plot_data <- plot_data[!is.na(plot_data$value) & !is.na(plot_data$entity), ]
    plot_data <- plot_data[!is.na(plot_data$gradient), ]

    # Convert entity to factor for proper ordering
    plot_data$entity <- factor(plot_data$entity)

    # Set gradient scale limits
    if (!is.null(upper.limit)) {
      scale_limits <- c(lower.limit, upper.limit)
    } else {
      scale_limits <- NULL
    }

    # ── Build ggplot for this metric ───────────────────────────────────────────────────────────

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = entity, y = value)) +

      # Boxplot layer (without points to avoid duplication)
      ggplot2::geom_boxplot(
        fill = fill_color,
        color = "black",
        alpha = 0.7,
        outlier.size = outlier.size,
        outlier.shape = NA  # Hide outliers here; they'll show as gradient points
      ) +

      # Points layer with gradient coloring
      ggplot2::geom_jitter(
        ggplot2::aes(color = gradient),
        width = 0.2,
        size = pt.size,
        alpha = pt.alpha,
        show.legend = (metric == qc_metrics[1])  # Show legend only for first metric
      ) +

      # Viridis gradient scale for points
      ggplot2::scale_color_viridis_c(
        name = gradient_label,
        option = scale.colors,
        limits = scale_limits,
        na.value = "grey70"
      ) +

      # Labels
      ggplot2::labs(
        title = metric,
        x = entity_name,
        y = "Value"
      ) +

      # Theme
      ggplot2::theme_classic(base_size = 11) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 11, hjust = 0.5),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
        axis.text.y = ggplot2::element_text(size = 9),
        axis.title.x = ggplot2::element_text(size = 10),
        axis.title.y = ggplot2::element_text(size = 10),
        legend.title = ggplot2::element_text(size = 9),
        legend.text = ggplot2::element_text(size = 8),
        legend.key.height = ggplot2::unit(0.8, "cm"),
        legend.key.width = ggplot2::unit(0.25, "cm"),
        panel.grid.major.y = ggplot2::element_line(color = "grey90", linewidth = 0.25),
        panel.border = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(linewidth = 0.4)
      )

    plot_list[[metric]] <- p
  }

  # ── 4. Combine panels using patchwork ────────────────────────────────────────────────────────

  if (length(plot_list) == 1) {
    combined <- plot_list[[1]]
  } else {
    combined <- patchwork::wrap_plots(plot_list, ncol = ncol) +
      patchwork::plot_layout(guides = "collect") &
      ggplot2::theme(legend.position = "right")
  }

  return(combined)
}
