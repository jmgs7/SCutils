#' @title VariableFeaturePlot2
#'
#' @description
#' This function is a modified version of Seurat's VariableFeaturePlot.
#' It allows users to specify a custom list of highly variable features (hvf) for plotting,
#' rather than relying solely on the variable features stored in the Seurat object. It also forces
#' the use of the output of SeuratObject::VariableFeatures() to highligh the high variable features
#' when the hvf parameter is not provided. This can be useful for visualizing specific sets of
#' features, for example, those which have been filtered based on certain criteria.
#'
#' @inheritParams Seurat::FeatureScatter
#' @inheritParams SeuratObject::HVFInfo
#' @param cols Colors to specify non-variable/variable status
#' @param assay Assay to pull variable features from
#' @param log Plot the x-axis in log scale
#' @param raster Convert points to raster format, default is \code{NULL}
#' which will automatically use raster if the number of points plotted is greater than
#' 100,000
#' @param hvf The list of highly variable features to use for plotting. If \code{NULL}, the function will
#' use the variable features stored in the Seurat object.
#' @return A ggplot object
#'
#' @import Seurat
#' @import SeuratObject
#' @importFrom ggplot2 labs scale_color_manual scale_x_log10
#' @export
#'
#' @examples
#' \dontrun{
#'   data("pbmc_small")
#'   pbmc_small <- data("pbmc_small") |>
#'     Seurat::NormalizeData() |>
#'     Seurat::FindVariableFeatures() |>
#'     SCutils::FilterVariableFeatures() |>
#'     Seurat::ScaleData()
#'   hvf <- SeuratObject::VariableFeatures(pbmc_small)
#'   # Plot a unique variable feature plot for all layers but highligh
#'   # the consensus variable features across all layers instead of the r
#'   # recalculated variable features for each layer.
#'   pbmc_small |>
#'     JoinLayers() |>
#'     Seurat::FindVariableFeatures() |>
#'     SCutils::VariableFeaturePlot2(hvf = hvf)
#' }
#'
VariableFeaturePlot2 <- function(
  object,
  cols = c('black', 'red'),
  pt.size = 1,
  log = NULL,
  selection.method = NULL,
  assay = NULL,
  raster = NULL,
  raster.dpi = c(512, 512),
  hvf = NULL
) {
  if (length(x = cols) != 2) {
    stop("'cols' must be of length 2")
  }
  hvf.info <- SeuratObject::HVFInfo(
    object = object,
    assay = assay,
    method = selection.method,
    status = TRUE
  )

  status.col <- colnames(hvf.info)[grepl("variable", colnames(hvf.info))][[1]]

  if (is.null(x = hvf)) {
    hvf <- SeuratObject::VariableFeatures(
      object = object,
      assay = assay,
      method = selection.method
    )
    # Ensure that the filtered variables have their status set to 'no'.
    hvf.info[!rownames(x = hvf.info) %in% hvf, status.col] <- FALSE
  } else {
    if (!all(hvf %in% rownames(x = hvf.info))) {
      stop("Some features in 'hvf' are not present in the hvf info dataframe.")
    }
    # Reset the status of all features to FALSE.
    hvf.info[rownames(x = hvf.info), status.col] <- FALSE
    # Set the status of the specified hvf features to TRUE.
    hvf.info[rownames(x = hvf.info) %in% hvf, status.col] <- TRUE
  }

  var.status <- c('no', 'yes')[unlist(hvf.info[[status.col]]) + 1]

  if (colnames(x = hvf.info)[3] == 'dispersion.scaled') {
    hvf.info <- hvf.info[, c(1, 2)]
  } else if (colnames(x = hvf.info)[3] == 'variance.expected') {
    hvf.info <- hvf.info[, c(1, 4)]
  } else {
    hvf.info <- hvf.info[, c(1, 3)]
  }
  axis.labels <- switch(
    EXPR = colnames(x = hvf.info)[2],
    'variance.standardized' = c('Average Expression', 'Standardized Variance'),
    'dispersion' = c('Average Expression', 'Dispersion'),
    'residual_variance' = c('Geometric Mean of Expression', 'Residual Variance')
  )
  log <- log %||%
    (any(
      c('variance.standardized', 'residual_variance') %in%
        colnames(x = hvf.info)
    ))

  plot <- Seurat::SingleCorPlot(
    data = hvf.info,
    col.by = var.status,
    pt.size = pt.size,
    raster = raster,
    raster.dpi = raster.dpi
  )
  if (length(x = unique(x = var.status)) == 1) {
    switch(
      EXPR = var.status[1],
      'yes' = {
        cols <- cols[2]
        labels.legend <- 'Variable'
      },
      'no' = {
        cols <- cols[1]
        labels.legend <- 'Non-variable'
      }
    )
  } else {
    labels.legend <- c('Non-variable', 'Variable')
  }
  plot <- plot +
    labs(title = NULL, x = axis.labels[1], y = axis.labels[2]) +
    scale_color_manual(
      labels = paste(labels.legend, 'count:', table(var.status)),
      values = cols
    )
  if (log) {
    plot <- plot + scale_x_log10()
  }
  return(plot)

  # TODO: Integrate labeling of top variable features as an option.
}
